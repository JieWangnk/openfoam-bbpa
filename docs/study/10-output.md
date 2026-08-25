# 10 — `output{}`: where and when files are written

```c
output
{
    layout                  phaseAligned;  // | companion | both
    cycleWriteInterval      1;   // 0 = at every write; N = every N cycles
    meanWriteIntervalCycles 1;   // <= 0 = final write only
    writeFinal              true;
    writeProperties         true;
}
```

Everything in `output{}` is presentation: it decides directories and
cadence, never a number. The same accumulator state produces the same
values under every setting here.

One default depends on the block's presence. With `output{}` given (even
empty), `cycleWriteInterval` defaults to 1: one dump per completed
cycle. With no `output{}` at all, the pre-v2 behaviour of dumping at
every solver write is kept.

## `layout`

* `phaseAligned` (default): each bin's files go to its own label-time
  directory (chapter 02). The full statistics set is written.
* `companion`: at each solver write, the current bin's phase mean is
  written into the solver's own time directory, beside `U`. Only the
  `_PA` mean product is written in this layout.
* `both`: both of the above.

`study_d10a_companion` uses `companion` with the solver writing every
bin width (0.125 s), so every write has a current bin to report:

<!-- demo:study_d10a_companion -->
```text
demo: study_d10a_companion
chapter: 10-output.md

case: study_d10a_companion  (endTime 5, solver writeInterval 0.125)
BBPA dictionary body:
    fields          (U);
    period          1.0;
    nBins           8;
    discardCycles   0;
    phaseOrigin     0;
    output { layout companion; }
    executeControl  timeStep;
    executeInterval 1;
    writeControl    writeTime;
run: completed

    file                  dirs  first    last
    U_Bin                    0  -        -
    U_Bin_UU                 0  -        -
    U_PA                    33  1        5
    U_PAVarianceMag          0  -        -
    U_PA_UU                  0  -        -
    U_PAmeanState            0  -        -
    uniform/*Properties     41
folded cycles at the newest checkpoint: U=5
deprecation warnings in the log: 0
```

`U_PA` sits in the solver's own directories from `1/` on: products are
finalised-only, so nothing appears until the first cycle completes.
There are no `U_PA_UU` or variance files anywhere; pair `companion`
with `phaseAligned` (i.e. `both`) when you need those.

## `cycleWriteInterval` and `meanWriteIntervalCycles`

`cycleWriteInterval N` fires a phase-aligned dump on the first solver
write after every Nth completed cycle. `meanWriteIntervalCycles`
throttles the cross-cycle `_PA` family within those dumps; `<= 0`
writes it at the final dump only.

`study_d10b_cadence` combines `cycleWriteInterval 2` with
`meanWriteIntervalCycles 0`:

<!-- demo:study_d10b_cadence -->
```text
demo: study_d10b_cadence
chapter: 10-output.md

case: study_d10b_cadence  (endTime 5, solver writeInterval 1)
BBPA dictionary body:
    fields          (U);
    period          1.0;
    nBins           8;
    discardCycles   0;
    phaseOrigin     0;
    output { cycleWriteInterval 2; meanWriteIntervalCycles 0; }
    executeControl  timeStep;
    executeInterval 1;
    writeControl    writeTime;
run: completed

    file                  dirs  first    last
    U_Bin                    0  -        -
    U_Bin_UU                 0  -        -
    U_PA                     8  4        4.875
    U_PAVarianceMag          8  4        4.875
    U_PA_UU                  8  4        4.875
    U_PAmeanState            0  -        -
    uniform/*Properties      6
folded cycles at the newest checkpoint: U=5
deprecation warnings in the log: 0
```

Only the final dump produced files: 8 directories, `4 … 4.875`, the
last completed cycle's labels. The cadence dumps after the second and
fourth completed cycles happened, but with the `_PA` family deferred
to the end and no per-cycle capture on, they had nothing to write. With 5 folded cycles
inside, the run wrote 8 small directories instead of 40.

## `writeFinal` and `writeProperties`

`writeFinal true` (default) forces a dump at the run's last write, on
top of any cadence. `writeProperties true` (default) writes the
restart checkpoint (`uniform/<name>Properties` plus the `*State`
files); `off` disables checkpointing and the run warns that it cannot
be restarted.

`study_d10c_no_final` switches both off, with a cadence that fires
once:

<!-- demo:study_d10c_no_final -->
```text
demo: study_d10c_no_final
chapter: 10-output.md

case: study_d10c_no_final  (endTime 5, solver writeInterval 1)
BBPA dictionary body:
    fields          (U);
    period          1.0;
    nBins           8;
    discardCycles   0;
    phaseOrigin     0;
    output { writeFinal off; writeProperties off; cycleWriteInterval 4; }
    executeControl  timeStep;
    executeInterval 1;
    writeControl    writeTime;
run: completed

    file                  dirs  first    last
    U_Bin                    0  -        -
    U_Bin_UU                 0  -        -
    U_PA                     8  3        3.875
    U_PAVarianceMag          8  3        3.875
    U_PA_UU                  8  3        3.875
    U_PAmeanState            0  -        -
    uniform/*Properties      0
folded cycles at the newest checkpoint: (no checkpoint on disk)
deprecation warnings in the log: 0
```

One dump exists, fired when the fourth cycle completed, at that
cycle's labels (`3 … 3.875`); the fifth cycle's fold reached no dump
because the final write was disabled. The checkpoint row reads 0: this run cannot
be restarted, by choice.
