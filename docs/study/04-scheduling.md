# 04 — Scheduling: `executeControl`, `executeInterval`, `writeControl`

```c
executeControl  timeStep;
executeInterval 1;
writeControl    writeTime;
```

These are the standard OpenFOAM function-object scheduling entries.
They exist on every function object; this chapter explains them once,
and the other chapters assume them.

Every function object has two hooks. `execute()` runs on a schedule set
by `executeControl` + `executeInterval`; `write()` runs on a schedule
set by `writeControl` (+ `writeInterval` where that applies). For BBPA:

* `execute()` accumulates: it takes the current solver step, splits it
  at bin boundaries, and adds each piece into its bin.
* `write()` serialises: it dumps finished statistics and the restart
  checkpoint to disk. It never changes a number.

## Why `executeControl timeStep; executeInterval 1;` is required

BBPA is a time integrator. Each call to `execute()` integrates exactly
the interval of the step that just finished, `(t − dt, t]`. With
`executeInterval 2`, `execute()` runs on every second step and the
skipped steps' intervals are never accumulated: the bin means then
describe half the signal. The same applies to any `executeControl`
that skips steps. Always accumulate every step.

This is cheap. The measured overhead of per-step accumulation is below
one percent of solver walltime (see the top-level README).

## `writeControl writeTime`

`writeTime` makes BBPA write whenever the solver writes, which keeps
the phase-averaged output and the restart checkpoint in step with the
solver's own snapshots. That is the right choice in production: a
restart from the last solver write then finds a matching BBPA
checkpoint next to it.

Other values (`runTime 0.5;`, `timeStep 100;`, …) work and only change
when `write()` fires. What each `write()` actually puts on disk is
governed by `output{}` on top of this schedule (chapter 10): with the
default `cycleWriteInterval 1`, a write call dumps the ensemble only
after a new cycle has completed since the last dump.
