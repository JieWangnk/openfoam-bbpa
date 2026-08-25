# 06 — The phase grid

Five entries define where the bins sit in time. They can be given flat
in the BBPA block, or grouped in a `phaseGrid{}` sub-dictionary (either
form, never both):

```c
phaseGrid
{
    period          1.0;
    nBins           8;
    phaseOrigin     0.125;
    binAlignment    left;
    timeQuadrature  endPoint;
}
```

## `period` (required)

The cycle length in seconds. The only entry with no default.

## `nBins` (default 100)

Bins per cycle. The bin width is `period/nBins`. Bin `j` of cycle `c`
gets the label time `phaseOrigin + c*period + j*period/nBins`, and the
output for that bin lands in the directory named after the label.

## `phaseOrigin` (default: auto-aligned)

Where phase 0 sits on the time axis. When omitted, BBPA aligns it to
the cycle boundary at or before the run's start time
(`floor(startTime/period)*period`), which is what you want for a run
starting at t = 0. Give it explicitly when your cycle has a known
offset, e.g. an inlet waveform whose systole starts at t = 0.125 s.

The old spelling `startTime` is deprecated; it collided with the
standard function-object scheduling entry of the same name.

## `binAlignment` (default `left`)

Where a bin's integration window sits relative to its label. `left`
starts the window at the label. `centred` centres it there, which
removes the half-bin phase lag of a left-edge window; the label then
names the middle of what was averaged. `tutorials/bbpaVerify` uses
`centred` and its README explains the consequences for reading the
output.

## `timeQuadrature` (default `endPoint`, the only value)

The rule for turning one field sample per step into an integral: the
step's single sample is treated as constant over the whole step
interval, first order in time. The entry exists so the numerical
meaning of a case is written down and checkpointed; a future
trapezoidal rule would be a different value here and a restart across
that change would be refused, never silently reinterpreted.

## Demo: shifted origin, grouped form

`study_d06_phase_grid` uses the grouped `phaseGrid{}` with
`phaseOrigin 0.125` and left alignment:

<!-- demo:study_d06_phase_grid -->
```text
demo: study_d06_phase_grid
chapter: 06-phase-grid.md

case: study_d06_phase_grid  (endTime 5, solver writeInterval 1)
BBPA dictionary body:
    fields          (U);
    phaseGrid
    {
        period          1.0;
        nBins           8;
        phaseOrigin     0.125;
        binAlignment    left;
        timeQuadrature  endPoint;
    }
    discardCycles   0;
    executeControl  timeStep;
    executeInterval 1;
    writeControl    writeTime;
run: completed

    file                  dirs  first    last
    U_Bin                    0  -        -
    U_Bin_UU                 0  -        -
    U_PA                    32  0.125    4
    U_PAVarianceMag         24  1.125    4
    U_PA_UU                 32  0.125    4
    U_PAmeanState            0  -        -
    uniform/*Properties      6
folded cycles at the newest checkpoint: U=4
deprecation warnings in the log: 0
```

Two things moved against the baseline (chapter 02):

* Every label shifted by 0.125 s: `U_PA` now starts at `0.125/` and
  the grid runs `0.125, 0.25, …`.
* Only 4 cycles folded instead of 5. With the origin at 0.125 s, the
  run's first 0.125 s sits before cycle 0 begins, and the cycle open
  at the end (starting 4.125 s) never completes. Whole cycles enter
  the ensemble; partial ones never do.

`deprecation warnings in the log: 0` confirms these five spellings are
the current interface, no compatibility path involved.
