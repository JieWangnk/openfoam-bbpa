# bbpaVerify — the dictionary reference

**One job: demonstrate every BBPA dictionary control on a case small enough
to check by hand.** 320 cells, fixed `deltaT`, a sinusoidal inlet with
period 1 s, 10 bins, 2 cycles, ~10 s wall clock. This is the only tutorial
that explains the controls; every other tutorial has one focused job and
refers here.

Run:

```bash
./Allrun          # blockMesh + foamRun (OpenFOAM 14, libPhaseAveraging built)
```

## The v2 dictionary model

* `fields ( … )` decides **what exists** — per-field products.
* `output { … }` decides **where and when** it is written.
* The accumulator has one numerical meaning throughout: nothing in
  `output{}` can change a number.

### Per-field products

```
fields
(
    U                                    // plain name: every product
    p { secondMoment off; }              // tune per field
);
```

| option | default | kind | meaning |
|---|---|---|---|
| `writePhaseMean` | on | output | write the `_PA` mean file. The mean itself **always** accumulates (the variance folds against it); off suppresses only the file. With the variance on, the mean then persists as the checkpoint-only `_PAmeanState`. |
| `secondMoment` | on | calculation | allocate + accumulate `<phi·phi>`, write `_PA_UU`. Off skips the dominant per-bin memory (a symmTensor field per bin for vectors). |
| `crossCycleVariance` | on | calculation | allocate + accumulate the Welford M2, write `_PAVariance` (defined for ≥ 2 cycles). |
| `writePerCycle` | off | output | keep and write the completed cycle's own bin mean as `_Bin` — this cycle's value, never blended into the cross-cycle mean. |

### The phase grid

`period`, `nBins` (bin width = period/nBins), `phaseOrigin` (where phase 0
sits; omit to auto-align to the cycle boundary at or before the start; the
old spelling `startTime` is deprecated — it collided with the standard
function-object scheduling entry), `discardCycles` (skip startup cycles),
`maxCycles` (stop after N complete cycles; −1 unlimited).

### Which bins

```
phaseSelection { mode all; }                      // every bin (default)
phaseSelection { mode list; bins (2 5); }         // a sparse subset
phaseSelection { mode solverWriteTimes; }         // derived from the solver
```

`solverWriteTimes` (used here) requires the solver's `writeControl
runTime|adjustableRunTime` with `writeInterval` an integer multiple of the
bin width — then every bin label coincides with a solver write directory
and the phase-averaged fields land beside `U` and `p`. `binAlignment
centred` centres each bin's window on its label (removing the half-bin
phase lag); `left` starts it there.

### Output

```
output
{
    layout                  phaseAligned;  // | companion | both
    cycleWriteInterval      1;   // 0 = every write; N = every N cycles
    meanWriteIntervalCycles 1;   // <= 0 = final write only
    writeFinal              true;
    writeProperties         true;  // off = no restart checkpoint
}
```

`phaseAligned` writes the bins to their own label-time directories;
`companion` writes the current bin beside the solver's fields at every
write; `both` does both. The per-cycle `_Bin` capture is **not** a layout —
it is the per-field `writePerCycle` product.

Deprecated spellings (`startCycle`, `cycles`, `binsOfInterest`,
`binsAtWriteTimes`, `binCentred`, `PAwriteInterval`, `writeMode`) still
work for one release, each with a warning naming its replacement; giving
old and new together aborts.

## What to look at

With this setup every write directory `0.1, 0.2, …, 2` carries `U` and
`U_PA` side by side, `U_PA` is the running cross-cycle phase mean of the
sinusoidal inlet response, and the run restarts exactly (kill it, restart
from the last write, the numbers match the uninterrupted run — the strict
restart metadata aborts if the configuration changed meanwhile).

The machine-precision ground truth for all of this lives in `../../tests/`
(`./Alltest`): a prescribed analytic field, a bit-faithful offline replica,
and closed-form references at 1e-12.
