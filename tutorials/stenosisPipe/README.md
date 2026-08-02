# stenosisPipe — production-like parallel operation

**One job: prove serial and parallel BBPA agree.** A pulsatile stenosed
pipe (18k cells, LES, 3 forcing periods), decomposed onto 4 ranks by the
Allrun. The phase-averaged products involve per-bin reductions across
processor boundaries; this case is where a broken reduction or a
processor-local write path would show.

`binDeltaT = period/nBins = 0.02 = writeInterval`, so every phase-aligned
directory IS a solver write directory: `U` and `U_PA` land side by side
with no extra directories, serial or parallel.
The first (transient) cycle is discarded (`discardCycles 1`). The IPA
function object runs alongside for a classical pointwise-estimator
comparison.

Run:

```bash
./Allrun            # blockMesh, decomposePar, mpirun -np 4, reconstructPar
```

The serial-equals-parallel assertion at machine precision (and the
restart-in-parallel path) lives in `../../tests/`. Dictionary reference:
`../bbpaVerify/README.md`.
