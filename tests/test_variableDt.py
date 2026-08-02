#!/usr/bin/env python3
"""
test_variableDt.py -- bin-edge-straddling and variable-time-step edges
(roadmap step 4 / the numerical rule in step-plan section 4).

Two runs of the bin-constant prescribe trajectory whose steps do NOT land on
bin edges:
  V1: fixed dt = 0.05 (does not divide binDeltaT = 0.125, and is not exactly
      representable in binary -> float drift on top of straddling);
  V2: genuinely VARIABLE dt via the setTimeStep function object
      (deltaT = table of (t, dt), linearly interpolated).

Under the v2 OVERLAP contract each step is split at every bin boundary it
crosses and each bin is credited its exact overlap; the step's single
sample (classified by its midpoint) is treated as constant over the
interval. Both runs assert:

  1. implementation: on-disk == replica bit-faithfully (the replica
     implements the identical splitting);
  2. weight conservation, the v2 identities: per completed cycle
     |sum_j W_{c,j} - T| < 1e-9 (machine-level), and the total assigned
     mass equals the logged mass clipped at the phase origin (the FO's
     startup execute covers the phantom interval (-dt, 0], which the
     overlap contract correctly excludes -- the old whole-step rule
     silently counted it);
  3. discretisation: on-disk vs the continuous closed form now shows the
     HONEST first-order quadrature error of one sample per step under
     straddling (reported, sanity-bounded < 0.5) -- under the old rule the
     bin-constant field masked it by construction.

Runnable standalone (exit 1 on failure, for Alltest) and via pytest.
"""
from __future__ import annotations
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bbpa_verify_lib as rv
import bbpa_offline as bo

TOL = 1e-12

SET_DT_FO = """
    setDt
    {
        type            setTimeStep;
        libs            ("libutilityFunctionObjects.so");
        deltaT          table ((0 0.028) (2 0.0415) (4 0.019) (6.5 0.033));
    }
"""


def _err(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def max_field_err(a, b, bins, fields=("PA", "PA_UU", "PAVariance")):
    out = {}
    for fld in fields:
        e = 0.0
        for j in bins:
            va, vb = a.get(j, {}).get(fld), b.get(j, {}).get(fld)
            if va is None or vb is None:
                continue
            e = max(e, _err(va, vb))
        out[fld] = e
    return out


def _weight_report(replica, samples):
    """Total-assignment conservation + per-cycle deviation (v2: ~0)."""
    # logged mass clipped at the phase origin: the startup execute's step
    # covers (-dt, 0], which is not integrable data
    total_logged = sum(min(dt, max(0.0, t)) for (t, dt, _) in samples)
    total_assigned = float(sum(w.sum() for w in replica.cycleWeightHist)
                           + replica.binTime.sum())
    max_dt = max(dt for (_, dt, _) in samples)
    devs = [abs(float(w.sum()) - rv.T) for w in replica.cycleWeightHist]
    return dict(total_logged=total_logged, total_assigned=total_assigned,
                conservation=abs(total_assigned - total_logged),
                max_cycle_dev=max(devs) if devs else 0.0, max_dt=max_dt)


def run_one(name, dt, extra_fo="", end_time=6.3):
    os.makedirs(rv.RUNROOT, exist_ok=True)
    rundir, logpath = rv.make_case(name, "bin_constant", dt, end_time, 1.0,
                                   extra_fo=extra_fo)
    assert rv.block_mesh(rundir).returncode == 0, name + ": blockMesh failed"
    assert rv.run_serial(rundir).returncode == 0, name + ": foamRun failed"

    samples = rv.read_log(logpath)
    cs = rv.latest_cycle_start(rundir)
    assert cs is not None, name + ": no full phase-aligned bin set"
    on_disk = rv.read_bbpa_output(rundir, cs)
    allb = [j for j in range(rv.I) if j in on_disk]
    assert len(allb) == rv.I, f"{name}: expected {rv.I} bins, got {len(allb)}"

    replica = bo.BBPAReplica(rv.T, rv.I).replay(samples)
    closed = bo.closed_form_bin_constant(rv.A, rv.c_of_n, int(replica.nCycles))

    e_impl = max_field_err(on_disk, replica.write_fields(), allb)
    e_disc = max_field_err(on_disk, closed, allb)
    wrep = _weight_report(replica, samples)
    return dict(name=name, n_cycles=int(replica.nCycles), cycle_start=cs,
                impl=e_impl, disc=e_disc, weights=wrep,
                n_steps=len(samples))


def _check(r, errors):
    n = r["name"]
    if max(r["impl"].values()) >= TOL:
        errors.append(f"{n}: impl {r['impl']}")
    # the honest first-order quadrature error: reported, sanity-bounded
    if max(r["disc"].values()) >= 0.5:
        errors.append(f"{n}: disc beyond first-order sanity {r['disc']}")
    w = r["weights"]
    if w["conservation"] >= 1e-9 * max(1.0, w["total_logged"]):
        errors.append(f"{n}: total-weight conservation {w['conservation']:.3e}")
    # the v2 identity: sum_j W_{c,j} == T per completed cycle
    if w["max_cycle_dev"] >= 1e-9:
        errors.append(f"{n}: per-cycle weight deviation {w['max_cycle_dev']:.3e}"
                      f" (v2 requires ~0)")


def run_matrix():
    errors = []
    results = []
    results.append(run_one("V1_nondividing", 0.05))
    results.append(run_one("V2_variable", 0.028, extra_fo=SET_DT_FO))
    for r in results:
        _check(r, errors)
    return results, errors


def test_variableDt():
    results, errors = run_matrix()
    assert not errors, "\n".join(errors)


if __name__ == "__main__":
    results, errors = run_matrix()
    for r in results:
        w = r["weights"]
        print(f"{r['name']}: steps={r['n_steps']} cycles={r['n_cycles']} "
              f"cycleStart={r['cycle_start']}")
        print(f"    impl  (on-disk vs replica)     max = "
              f"{max(r['impl'].values()):.3e}")
        print(f"    disc  (on-disk vs closed form) max = "
              f"{max(r['disc'].values()):.3e}")
        print(f"    weights: total conservation = {w['conservation']:.3e}; "
              f"per-cycle |sum W - T| = {w['max_cycle_dev']:.3e} "
              f"(v2 overlap identity)")
    if errors:
        print("FAIL:")
        for e in errors:
            print("  ", e)
        sys.exit(1)
    print("RESULT: PASS")
    sys.exit(0)
