# 05 — `fields` and the per-field products

```c
fields
(
    U                                    // plain name: every default
    p { secondMoment off; }              // a { } after a name tunes it
);
```

`fields` lists what to average. Volume scalar and vector fields are
supported: solver fields (`U`, `p`) and derived fields made by other
function objects (`wallShearStress`), as long as their producer runs
before BBPA (chapter 11). A plain name gets the default products; a
dictionary after a name switches individual products for that field.

The four per-field switches:

| Switch | Default | Kind | Effect |
|---|---|---|---|
| `writePhaseMean` | on | output | write the `_PA` mean file. The mean itself always accumulates, because the variance is computed against it. Off suppresses only the file; the mean then persists as the checkpoint-only `_PAmeanState` |
| `secondMoment` | on | calculation | allocate and accumulate `⟨φ⊗φ⟩`, write `_PA_UU`. This is the dominant per-bin memory for vectors (a symmetric tensor field per bin), so switch it off when you only need means |
| `crossCycleVariance` | on | calculation | allocate and accumulate the cycle-to-cycle variance, write `_PAVariance` / `_PAVarianceMag` |
| `writePerCycle` | off | output | keep and write each completed cycle's own bin mean as `_Bin` (and its second moment as `_Bin_UU`). These are that cycle's values; they are never blended into the cross-cycle mean |

The two "calculation" switches change what is allocated and computed;
turning one off saves memory and leaves every remaining product
unchanged. The two "output" switches only decide which files appear.

## Demo: calculation switches off

`study_d05a_second_moment_off` turns both calculation switches off:

<!-- demo:study_d05a_second_moment_off -->
```text
demo: study_d05a_second_moment_off
chapter: 05-fields.md

case: study_d05a_second_moment_off  (endTime 5, solver writeInterval 1)
BBPA dictionary body:
    fields          ( U { secondMoment off; crossCycleVariance off; } );
    period          1.0;
    nBins           8;
    discardCycles   0;
    phaseOrigin     0;
    executeControl  timeStep;
    executeInterval 1;
    writeControl    writeTime;
run: completed

    file                  dirs  first    last
    U_Bin                    0  -        -
    U_Bin_UU                 0  -        -
    U_PA                    40  0        4.875
    U_PAVarianceMag          0  -        -
    U_PA_UU                  0  -        -
    U_PAmeanState            0  -        -
    uniform/*Properties      6
folded cycles at the newest checkpoint: U=5
deprecation warnings in the log: 0
```

Only the mean family survives: `U_PA` is written exactly as in the
baseline (chapter 02) and nothing else appears.

## Demo: the mean file off

`study_d05b_mean_unpublished` turns off only the mean FILE:

<!-- demo:study_d05b_mean_unpublished -->
```text
demo: study_d05b_mean_unpublished
chapter: 05-fields.md

case: study_d05b_mean_unpublished  (endTime 5, solver writeInterval 1)
BBPA dictionary body:
    fields          ( U { writePhaseMean off; } );
    period          1.0;
    nBins           8;
    discardCycles   0;
    phaseOrigin     0;
    executeControl  timeStep;
    executeInterval 1;
    writeControl    writeTime;
run: completed

    file                  dirs  first    last
    U_Bin                    0  -        -
    U_Bin_UU                 0  -        -
    U_PA                     0  -        -
    U_PAVarianceMag         32  1        4.875
    U_PA_UU                 40  0        4.875
    U_PAmeanState           40  0        4.875
    uniform/*Properties      6
folded cycles at the newest checkpoint: U=5
deprecation warnings in the log: 0
```

`U_PA` vanishes and `U_PAmeanState` appears in the same 40
directories. The mean kept accumulating: the variance (which needs it)
is still produced, and the `meanState` spelling persists the mean so a
restart is complete. File output and calculation are separate
concerns; switching a file off never degrades the statistics.
