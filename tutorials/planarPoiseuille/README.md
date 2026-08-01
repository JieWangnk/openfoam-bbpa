# planarPoiseuille — cross-cycle variance and restart

**One job: prove the cross-cycle variance and the restart are correct.**
The Waters–King startup is a DEVELOPING flow: successive "cycles" genuinely
differ, so the cross-cycle variance (`U_PAVariance`) is non-zero and
checkable. On a perfectly periodic flow the variance is ~0 and a variance
bug is invisible — the pre-fix pseudo-Welford (biased 25% high on a
two-cycle example) survived for exactly that reason, and this case is what
exercises the true Welford inner-product recurrence.

`period 0.9997` is deliberately not a multiple of `deltaT`, so no state
lands exactly on a cycle boundary or bin edge.

Run:

```bash
./Allrun
```

Restart check (the second half of the job):

```bash
./Allrun                      # uninterrupted reference, endTime 5
# then: edit controlDict startTime to a mid-run write (e.g. 2.5),
# rerun foamRun, and compare the final U_PA / U_PAVariance — they match
# the uninterrupted run exactly. Changing period/nBins/binAlignment (or
# any accumulator-shaping entry) across the restart ABORTS with a message
# naming the entry: the checkpoint carries a configuration signature.
```

The machine-precision versions of both checks (variance == closed form ==
bit-faithful replica; continuous == restarted == restore-replica, exactly)
run in `../../tests/` (`test_fixedDt.py`, `test_restart.py`). The
dictionary controls are documented in `../bbpaVerify/README.md`.
