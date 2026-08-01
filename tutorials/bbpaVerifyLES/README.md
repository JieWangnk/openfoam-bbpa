# bbpaVerifyLES — adaptive time stepping and second-order statistics

**One job: stress BBPA with a genuinely VARIABLE time step on an unsteady
LES field.** `adjustTimeStep yes` (maxCo 0.6): steps grow and shrink with
the pulsatile Courant number and routinely straddle bin edges — the
situation the whole-step assignment rule and the boundary snap exist to
make deterministic. 11 520 cells, forcing period 0.5 s, ~5 cycles, an
inlet record whose ~0.07 rms perturbation differs every cycle, so the
phase mean, the second moment `<U U>` and the cross-cycle variance are
all non-trivial — for `U` and for the derived `wallShearStress`, both
with the per-cycle `_Bin` capture on (restart-safe at write granularity).

`writeControl adjustableRunTime` keeps the solver's writes exactly on the
0.02 grid, and `phaseSelection { mode solverWriteTimes; }` derives the 25
tracked bins from it (stride 4 of 100), each landing in a write folder,
centred on its label.

Run (~15–25 min):

```bash
./Allrun
```

The reference for correctness under variable dt is an offline
recomputation from the same logged trajectory, never a second turbulent
run (chaotic divergence would mix function-object accuracy with flow
realisation differences) — that machinery, at machine precision, is
`../../tests/test_variableDt.py`. Dictionary reference:
`../bbpaVerify/README.md`.
