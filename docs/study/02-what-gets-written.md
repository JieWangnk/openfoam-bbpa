# 02 — What gets written

BBPA cuts each cycle into `nBins` phase bins. Bin `j` gets a label time
inside every cycle: `phaseOrigin + cycle*period + j*(period/nBins)`.
Output lands in time directories named after those labels, so ParaView
and every OpenFOAM tool read it like ordinary solver output.

For a tracked field `U` (any vector field; scalar fields work the same
way), these files can appear:

| File | What it holds |
|---|---|
| `U_PA` | the cross-cycle phase mean of bin `j`: the equal-weight average of the completed cycles' bin means |
| `U_PA_UU` | the cross-cycle mean of the raw second moment `⟨U⊗U⟩` in that bin (a symmetric tensor). Reynolds stress `R = ⟨U⊗U⟩ − ⟨U⟩⊗⟨U⟩` and TKE `½·tr(R)` follow offline |
| `U_PAVarianceMag` | the cycle-to-cycle sample variance of the per-cycle bin means, summed over components (a scalar; needs 2 or more cycles). A scalar field `p` gets `p_PAVariance` instead |
| `U_Bin` | one cycle's own time-integrated bin mean. Only written with the per-field `writePerCycle` switch on (chapter 05); written once, when that cycle's window closes, and never updated |
| `U_Bin_UU` | the same cycle's own second moment |
| `U_PAmeanState`, `U_PA_UUState`, … | checkpoint spellings. State that must survive a restart but was excluded from the published products (chapters 05 and 09) |
| `uniform/<name>Properties` | the restart checkpoint: counters plus a configuration signature (period, nBins, alignment, …) that a restart validates before restoring |

Within one cycle, bin `j`'s value is the time-weighted mean of the field
over that bin's window: every solver step is split at bin edges and each
bin is credited its exact overlap in time. Across cycles, `U_PA` is the
plain average of those per-cycle bin means. Nothing partial is ever
written under a product name: a directory either holds a finished
statistic or nothing.

## The baseline demo

`study_d01_minimal` runs the analytic demo case (8 bins, period 1 s,
5 complete cycles, defaults everywhere) and lists what landed on disk:

<!-- demo:study_d01_minimal -->
```text
demo: study_d01_minimal
chapter: 02-what-gets-written.md

case: study_d01_minimal  (endTime 5, solver writeInterval 1)
BBPA dictionary body:
    fields          (U);
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
    U_PAVarianceMag         32  1        4.875
    U_PA_UU                 40  0        4.875
    U_PAmeanState            0  -        -
    uniform/*Properties      6
folded cycles at the newest checkpoint: U=5
deprecation warnings in the log: 0
```

Reading the table:

* `U_PA` in 40 directories = 8 bins × 5 cycles. Cycle 0's ensemble is
  dumped at labels `0 … 0.875`, cycle 1's at `1 … 1.875`, and so on.
  Each directory is written once and then left alone; walking the label
  grid forward shows the ensemble growing by one cycle per period.
* `U_PAVarianceMag` starts one cycle later (32 directories, first at
  `1/`): variance needs at least two completed cycles.
* `U_Bin` is absent because `writePerCycle` defaults to off.
* Six checkpoints: one `uniform/bbpaProperties` per solver write.
* `folded cycles: U=5`: all five completed cycles entered the ensemble.

## `U` vs `U_Bin` vs `U_PA` in one directory

The `tutorials/bbpaVerify` case (chapter 01) turns `writePerCycle` on,
so its directory `1.5/` (cycle 1's label for phase 0.5) holds three
different quantities under similar names:

* `U` is the instantaneous solution at t = 1.5 s. One realisation of
  one cycle, at one instant.
* `U_Bin` is cycle 1's own time-integrated mean over the window around
  phase 0.5.
* `U_PA` is the cross-cycle ensemble mean at phase 0.5 as of that
  cycle's dump. At `1.5/` it covers one retained cycle, so there
  `U_PA == U_Bin` identically. At `2.5/` it covers two.

Comparing them pointwise is the fastest way to build intuition. `U`
differs from `U_PA` at the same phase unless the flow repeats
perfectly; that difference is the cycle-to-cycle variability that
`U_PAVarianceMag` measures. The tutorial's `check_bbpa.py` turns these
relations into an automated test.
