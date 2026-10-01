#!/usr/bin/env python3
"""
test_fixedDt.py -- BBPA ground-truth arithmetic test (roadmap Tutorial-1 baseline).

Prescribe case: bin-constant analytic U(t) = A[bin] + c[cycle].  Two conditions
make the bin integral EXACT (no quadrature error), so every product is checkable
to machine precision:
  * the field is constant within each phase bin, and
  * dt is a power-of-two fraction of binDeltaT (dt = binDeltaT/16 here) so steps
    land exactly on bin edges with no floating-point drift.
    NB: a non-representable dt like 0.0125 (=1/80) drifts and mixes adjacent bins.

Asserts, to machine precision (roadmap Tutorial-1 baseline row + §6 tolerances):
  1. implementation: on-disk _PA/_PA_UU/_PAVariance == BBPAReplica.write_fields()
     over ALL bins (BBPA == its own documented algorithm, bit-faithful replica).
  2. discretisation: on-disk == closed_form_bin_constant over PRISTINE bins (the
     ones with no in-progress current-cycle data) -- the algorithm-independent
     closed form.  The single in-progress bin is validated by check (1) only,
     since the on-disk blend of a partial bin is a write-layout artefact, not a
     phase-mean.

Runnable standalone (exit 1 on failure, for Alltest) and importable by pytest.
"""
from __future__ import annotations
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bbpa_verify_lib as rv          # proven case-gen + OF field readers
import bbpa_offline as bo             # replica + closed-form references

TOL = 1e-12
DT = rv.BIN_DT / 16                   # 16 steps/bin, exact 2**-7 -> no drift
END_TIME = 6.0                        # 6 completed cycles -> variance defined
WRITE_INTERVAL = 3.0                  # proven cadence


def max_field_err(a, b, bins, fields=("PA", "PA_UU", "PAVariance")):
    out = {}
    for fld in fields:
        e = 0.0
        for j in bins:
            va, vb = a[j].get(fld), b[j].get(fld)
            if va is None or vb is None:
                continue
            e = max(e, float(np.max(np.abs(np.asarray(va, float) - np.asarray(vb, float)))))
        out[fld] = e
    return out


def run_case():
    os.makedirs(rv.RUNROOT, exist_ok=True)
    rundir, logpath = rv.make_case("fixedDt", "bin_constant", DT, END_TIME,
                                   WRITE_INTERVAL)
    assert rv.block_mesh(rundir).returncode == 0, "blockMesh failed"
    assert rv.run_serial(rundir).returncode == 0, "foamRun failed"
    return rundir, logpath


def compute(rundir, logpath):
    samples = rv.read_log(logpath)
    cs = rv.latest_cycle_start(rundir)
    assert cs is not None, "no full phase-aligned bin set on disk"
    on_disk = rv.read_bbpa_output(rundir, cs)

    replica = bo.BBPAReplica(rv.T, rv.I, start_time=0.0, start_cycle=0)
    replica.replay(samples)
    wref = replica.write_fields()
    n_cycles = int(replica.nCycles)
    closed = bo.closed_form_bin_constant(rv.A, rv.c_of_n, n_cycles)

    allb = [j for j in range(rv.I) if j in on_disk]
    prist = [j for j in allb if replica.binTime[j] == 0.0]   # no in-progress data

    e_impl = max_field_err(on_disk, wref, allb)
    e_disc = max_field_err(on_disk, closed, prist)
    return dict(n_cycles=n_cycles, cycle_start=cs, allb=allb, prist=prist,
                impl=e_impl, disc=e_disc)


def test_fixedDt():
    """Ground-truth arithmetic, all three layers, to machine precision:
      implementation: on-disk == replica.write_fields() over ALL bins
                      (including the in-progress bin's pure-finalised _PA);
      discretisation: on-disk == algorithm-independent closed form over
                      PRISTINE bins, for _PA, _PA_UU AND _PAVariance.
    Reconciliation history: the ported replica initially modelled two
    superseded behaviours -- the pre-fix 0.5*(magSqr+magSqr) pseudo-Welford
    M2 and the removed in-progress-cycle write blend. Both were stale
    reference models (the library's own comments document each fix); with
    the true Welford inner product and the no-blend write rule the errors
    fell from 1.7e-3 / 5e-2 to ~1e-16 / 0."""
    rundir, logpath = run_case()
    r = compute(rundir, logpath)
    assert max(r["impl"].values()) < TOL, ("impl", r["impl"])
    assert max(r["disc"].values()) < TOL, ("disc", r["disc"])


if __name__ == "__main__":
    rundir, logpath = run_case()
    r = compute(rundir, logpath)
    print(f"fixedDt: dt={DT} cycles={r['n_cycles']} cycleStart={r['cycle_start']} "
          f"bins={len(r['allb'])} pristine={len(r['prist'])}")
    print("  implementation (on-disk vs replica, all bins):")
    for k, v in r["impl"].items():
        print(f"      {k:<11} max|err| = {v:.3e}")
    print("  discretisation (on-disk vs closed form, pristine bins):")
    for k, v in r["disc"].items():
        print(f"      {k:<11} max|err| = {v:.3e}")
    ok = max(r["impl"].values()) < TOL and max(r["disc"].values()) < TOL
    print("  RESULT:", "PASS" if ok else "FAIL", f"(tol {TOL:.0e})")
    sys.exit(0 if ok else 1)
