# 09 — `outputSelection`: which accumulated bins reach disk

```c
outputSelection { mode allAccumulated; }        // default
outputSelection { mode list; bins (2 5); }      // publish a subset
```

`outputSelection` gates serialisation only. Chapter 08's selection
decides which bins are computed; this one decides which of those get
product files. It can never widen the set: listing a bin that is
outside the accumulation selection is a startup error (the subset
rule). And because it touches no statistic, it may change freely
across a restart.

## The `allAccumulated` spelling

The default mode is `allAccumulated`: every accumulated bin is
written. Two tempting spellings are rejected with a fatal error:
`mode accumulated` and `mode all`. The mode was renamed before release
because an "all" distinct from the accumulated set would break the
subset rule. If a case aborts at startup right after you added this
block, check this spelling first.

## Demo: publish two bins, compute all eight

<!-- demo:study_d09_output_list -->
```text
demo: study_d09_output_list
chapter: 09-output-selection.md

case: study_d09_output_list  (endTime 5, solver writeInterval 1)
BBPA dictionary body:
    fields          (U);
    period          1.0;
    nBins           8;
    discardCycles   0;
    phaseOrigin     0;
    accumulationSelection { mode all; }
    outputSelection { mode list; bins (2 5); }
    executeControl  timeStep;
    executeInterval 1;
    writeControl    writeTime;
run: completed

    file                   dirs  first    last
    U_Bin                     0  -        -
    U_Bin_UU                  0  -        -
    U_PA                     10  0.25     4.625
    U_PAVarianceMag           8  1.25     4.625
    U_PAVarianceMagState     24  1        4.875
    U_PA_UU                  10  0.25     4.625
    U_PA_UUState             30  0        4.875
    U_PAmeanState            30  0        4.875
    uniform/*Properties       6
folded cycles at the newest checkpoint: U=5
deprecation warnings in the log: 0
```

The product files (`U_PA`, `U_PA_UU`, `U_PAVarianceMag`) appear only
at bins 2 and 5's labels, exactly as in chapter 08's sparse demo. The
new rows are the point:

* `U_PAmeanState`, `U_PA_UUState`, `U_PAVarianceMagState` appear at
  the other six bins' labels (30 = 6 bins × 5 cycles; the variance
  state needs two cycles, hence 24).

The ensemble state of every accumulated bin is checkpointed under the
`*State` spellings regardless of the output selection, so a restart
restores all eight bins and the selection can be widened later without
losing history. (The run also prints an older warning claiming
unwritten bins cold-start on restart; that text predates the `*State`
checkpoint, and the manifest above shows the state is on disk.)

Use this control when a fine accumulation grid would flood the disk:
compute everything, publish the phases you need now, and re-publish
others after a restart if it turns out you need them.
