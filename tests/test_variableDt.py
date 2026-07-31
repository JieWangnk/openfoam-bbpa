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

The bin-constant field is STRADDLE-IMMUNE by construction: a step crossing a
bin edge carries U evaluated at its END time -- the constant of the very bin
v1's whole-step rule assigns it to -- so bin MEANS (and hence variance) stay
exact even though bin WEIGHTS are misassigned.  This lets both runs assert:

  1. implementation: on-disk == replica bit-faithfully (<= 1e-12), i.e. BBPA
     behaves deterministically per its documented whole-step rule;
  2. closed form: on-disk finalised _PA/_PA_UU/_PAVariance == the exact
     analytic ensemble (<= 1e-12) despite straddling;
  3. weight bookkeeping: every logged dt is assigned exactly once (total
     conservation), and the per-cycle deviation |sum_j W_{c,j} - T| is bounded
     by the largest step in that cycle -- the quantified v1 straddle effect
     that the v2 split rule (Delta t_{n,j} overlap integration) will drive
     to zero.  The deviation is reported so the v2 change is measurable.

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
    """Total-assignment conservation + per-cycle straddle deviation."""
    total_logged = sum(dt for (_, dt, _) in samples)
    total_assigned = float(sum(w.sum() for w in replica.cycleWeightHist)
                           + replica.binTime.sum())
    # per-cycle |sum_j W_{c,j} - T|, bounded by the max step in the run
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
    if max(r["disc"].values()) >= TOL:
        errors.append(f"{n}: disc {r['disc']}")
    w = r["weights"]
    if w["conservation"] >= TOL * max(1.0, w["total_logged"]):
        errors.append(f"{n}: total-weight conservation {w['conservation']:.3e}")
    if w["max_cycle_dev"] > w["max_dt"] + TOL:
        errors.append(f"{n}: per-cycle weight deviation {w['max_cycle_dev']:.3e}"
                      f" exceeds max step {w['max_dt']:.3e}")


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
              f"per-cycle straddle deviation = {w['max_cycle_dev']:.3e} "
              f"(v1 whole-step rule; bound = max dt {w['max_dt']:.3e})")
    if errors:
        print("FAIL:")
        for e in errors:
            print("  ", e)
        sys.exit(1)
    print("RESULT: PASS")
    sys.exit(0)
