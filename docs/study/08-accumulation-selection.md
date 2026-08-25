# 08 — `accumulationSelection`: which bins exist

```c
accumulationSelection { mode all; }                // every bin (default)
accumulationSelection { mode list; bins (2 5); }   // a sparse subset
accumulationSelection { mode solverWriteTimes; }   // derived from the solver
```

`accumulationSelection` decides which bins exist in memory. A bin that
is out of the selection costs nothing: no fields are allocated for it
and its solver-step segments are dropped after classification. This is
the memory control for large cases: per tracked bin and vector field,
BBPA holds a mean field plus (products permitting) a symmetric-tensor
second moment and a variance accumulator.

Because it shapes the statistics themselves, the selection is part of
the restart signature: a restart with a different selection is refused.
Chapter 09's `outputSelection` is the free-to-change counterpart on the
disk side.

## `mode list`

`bins ( ... )` names the bin indices to keep, each in `[0, nBins)`,
each once. Typical use: a fine grid (`nBins 1000`) for sharp phase
resolution, with only the phases you care about (peak systole, early
deceleration) actually tracked.

`study_d08_accumulation_list` keeps bins 2 and 5 of 8:

<!-- demo:study_d08_accumulation_list -->
```text
demo: study_d08_accumulation_list
chapter: 08-accumulation-selection.md

case: study_d08_accumulation_list  (endTime 5, solver writeInterval 1)
BBPA dictionary body:
    fields          (U);
    period          1.0;
    nBins           8;
    discardCycles   0;
    phaseOrigin     0;
    accumulationSelection { mode list; bins (2 5); }
    executeControl  timeStep;
    executeInterval 1;
    writeControl    writeTime;
run: completed

    file                  dirs  first    last
    U_Bin                    0  -        -
    U_Bin_UU                 0  -        -
    U_PA                    10  0.25     4.625
    U_PAVarianceMag          8  1.25     4.625
    U_PA_UU                 10  0.25     4.625
    U_PAmeanState            0  -        -
    uniform/*Properties      6
folded cycles at the newest checkpoint: U=5
deprecation warnings in the log: 0
```

`U_PA` appears in 10 directories: 2 bins × 5 cycles, at labels
`c + 0.25` and `c + 0.625` only. Nothing else was computed at all.

## `mode solverWriteTimes`

Derives the list from the solver's write schedule so that every
tracked bin's label coincides with a solver write directory, and the
phase-averaged fields land beside `U` and `p`. It requires the
solver's `writeControl` to be `runTime` or `adjustableRunTime`, with
`writeInterval` an integer multiple of the bin width `period/nBins`;
the object then tracks every `writeInterval/binWidth`-th bin and logs
the stride. Anything else is refused at startup, because the write
times cannot be predicted from the dictionary.

This mode is a co-location convenience, used by `tutorials/bbpaVerify`
and the LES tutorial. For publication runs, prefer an explicit
`mode list`: the tracked phases are then stated in the dictionary
itself rather than inherited from the solver's I/O settings.

The deprecated spellings `binsOfInterest` (= `mode list`) and
`binsAtWriteTimes` (= `mode solverWriteTimes`) still work for one
release, with warnings. `phaseSelection`, which mixed this control
with output selection, is rejected outright.
