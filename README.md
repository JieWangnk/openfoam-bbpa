# openfoam-bbpa

**Bin-Based Phase Averaging for OpenFOAM** — an in-situ function object that
produces phase-averaged mean, cross-cycle variance, and resolved turbulent
kinetic energy directly in memory, without writing terabytes of transient
snapshots.

This branch (`of14-port`) is **BBPA v2 on OpenFOAM Foundation v14**; `main`
carries the v1 implementation on Foundation v12. The v2 dictionary interface
below is **frozen** as of tag `of14-v2.0`.

For the theory, error analysis, and cardiovascular LES validation see the
companion paper:

> Wang, J., Skillen, A., Revell, A. *Bin-based in-situ phase averaging for
> pulsatile internal-flow LES: error budget, OpenFOAM implementation, and
> validation.* Computers & Fluids (2026, submitted). Paper repo:
> [`Bin-Phase-Average`](https://github.com/JieWangnk/Bin-Phase-Average)

## Frozen definitions

These definitions are the contract. No output option changes them; the test
suite (`tests/`) asserts them against algorithm-independent closed forms at
machine precision.

For field *q*, cycle *c*, and phase bin *j* with window *I<sub>c,j</sub>*:

```
q_Bin[c,j]   = ( ∫_{I(c,j)} q(t) dt ) / W[c,j],      W[c,j] = ∫_{I(c,j)} dt
q_BinPA[j,N] = (1/N) · Σ_{c=1..N} q_Bin[c,j]
```

| Product | File | Meaning (unambiguous) |
|---|---|---|
| per-cycle bin mean | `<q>_Bin` | the just-completed cycle's own time-integrated bin mean — **never** blended into the cross-cycle mean |
| per-cycle second moment | `<q>_Bin_UU` | the same cycle's bin mean of *q⊗q* |
| cross-cycle mean | `<q>_PA` | equal-weight mean of the completed cycles' `_Bin` values |
| cross-cycle second moment | `<q>_PA_UU` | equal-weight mean of the per-cycle second moments; Reynolds stress `R = ⟨q⊗q⟩ − ⟨q⟩⊗⟨q⟩` and TKE `½·tr(R)` follow offline |
| cross-cycle variance | `<q>_PAVariance` | Welford `M₂/(N−1)` of the per-cycle bin means (true inner-product recurrence `M₂ += δ·δ′`; defined for N ≥ 2) |

Each solver interval [t<sub>n</sub>, t<sub>n+1</sub>] is **split at every
phase-bin or cycle boundary it crosses**; its contribution to bin
*I<sub>c,j</sub>* is the exact overlap

```
Δt[n,c,j] = max( 0, min(t_{n+1}, b_{c,j+1}) − max(t_n, b_{c,j}) )
```

so that `Σ_{c,j} Δt[n,c,j] = t_{n+1} − t_n` for every solver step and
`Σ_j W[c,j] = T` for every completed cycle (both asserted at machine level
in `tests/test_variableDt.py`). The step's single field sample is treated
as constant over the whole interval — a first-order-in-time quadrature; a
step that straddles a boundary therefore carries its one sample into both
bins, and the resulting bin-mean quadrature error is O(Δt) (measured and
reported in the same test). A floating-point snap (1e-9 in bin units) is
used **only to classify** times lying numerically on a boundary; it never
assigns a straddling step whole to one bin. The in-progress cycle is
**never** blended into any cross-cycle product, in any output layout.

## Features

| | BBPA | `fieldAverage` (OpenFOAM-native) | Classical offline PA |
|---|---|---|---|
| Phase-centred ensemble mean | ✅ | ✗ (causal running filter) | ✅ |
| In-situ (no snapshot I/O) | ✅ | ✅ | ✗ (terabytes of writes) |
| Raw 2nd moment `⟨UᵢUⱼ⟩` → Reynolds stress / TKE | ✅ | ✗ | ✅ (offline) |
| Cross-cycle variance | ✅ (Welford) | partial (`UPrime2Mean`) | ✅ |
| Per-field product selection (memory-gated) | ✅ | ✅ | n/a |
| Restart-safe, with a strict config signature | ✅ | partial | n/a |
| Overhead (8.2M-cell LES, 200 cores) | **+0.67%** at I=100 (95% CI [0.2,1.1]%) | baseline | +100%+ (post-run) |

The accompanying `IPA` function object samples at exact phase crossings
without binning — a zero-width limit of BBPA for reference comparisons.

## Requirements and build

OpenFOAM Foundation **v14** (openfoam.org), GCC ≥ 9, OpenMPI or equivalent.

```bash
cd src/phaseAveraging
wmake libso              # -> $FOAM_USER_LIBBIN/libPhaseAveraging.so
```

## The v2 dictionary

```c
functions
{
    #includeFunc wallShearStress    // producers of derived fields BEFORE bbpa

    bbpa
    {
        type            BBPA;
        libs            ("libPhaseAveraging.so");

        // scheduling: accumulate EVERY step (BBPA is a time integrator,
        // not a write-time sampler); write on the solver's writes
        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;

        // WHAT exists: per-field products (a plain name = all defaults)
        fields
        (
            U
            p               { secondMoment off; crossCycleVariance off; }
            wallShearStress { crossCycleVariance off; }
        );

        // the phase grid
        period          1.0;      // cycle period [s]
        nBins           50;       // bin width = period/nBins
        phaseOrigin     0;        // omit to auto-align to a cycle boundary
        discardCycles   2;        // skip startup cycles
        maxCycles       -1;       // stop after N complete cycles (-1 = all)

        // which bins, and window placement
        phaseSelection  { mode all; }   // | list; bins (…) | solverWriteTimes
        binAlignment    left;           // | centred (window centred on label)

        // WHERE and WHEN (never changes the numbers)
        output
        {
            layout                  phaseAligned;  // | companion | both
            cycleWriteInterval      1;   // 0 = every write; N = every N cycles
            meanWriteIntervalCycles 1;   // <= 0 = final write only
            writeFinal              true;
            writeProperties         true;  // off = no restart checkpoint
        }
    }
}
```

### Per-field product options

| option | default | kind | effect |
|---|---|---|---|
| `writePhaseMean` | on | output | write `_PA`. The mean itself **always** accumulates; with the variance on and this off, the mean persists as the checkpoint-only `_PAmeanState` |
| `secondMoment` | on | calculation | allocate + accumulate `⟨q⊗q⟩`, write `_PA_UU`/`_Bin_UU`. Off releases the dominant per-bin memory (measured: 90 MB at 500 bins × 1600 cells, vector) |
| `crossCycleVariance` | on | calculation | allocate + accumulate Welford `M₂`, write `_PAVariance` |
| `writePerCycle` | off | output | keep + write the completed cycle's `_Bin`/`_Bin_UU`, flushed as each bin's window closes (restart-safe at write granularity) |

### Restart

Every checkpoint carries a configuration signature (format version, period,
phase origin, bins, alignment, selection, per-field types and products,
counters, weights). On restart, any change that alters accumulator shape or
meaning **aborts naming the entry**: period, `nBins`, `binAlignment`, the
phase selection, an explicit `phaseOrigin` that contradicts the checkpoint,
a removed field or changed type, a product turned on without history.
Turning a product off warns and discards. `tests/test_restart.py` proves
continuous == restarted exactly.

### Migration from the v1 dictionary

Deprecated spellings work for one release, each warning with its
replacement; old + new together aborts.

| deprecated | v2 |
|---|---|
| `startTime` | `phaseOrigin` (the old word collided with the standard function-object scheduling entry) |
| `startCycle` | `discardCycles` |
| `cycles` | `maxCycles` |
| `binsOfInterest (…)` | `phaseSelection { mode list; bins (…); }` |
| `binsAtWriteTimes` | `phaseSelection { mode solverWriteTimes; }` |
| `binCentred` | `binAlignment centred` |
| `PAwriteInterval` | `output { meanWriteIntervalCycles …; }` |
| `writeMode phaseAlignedDirs / perCycle / companion / both` | `output { layout …; cycleWriteInterval …; }` |
| `writeMode bin` | per-field `writePerCycle on` (one naming: `_Bin` beside `_PA`) |

## Tutorials — one job each

Each case has a README stating its single job and an `Allrun`/`Allclean`.
`tutorials/bbpaVerify/README.md` is the full dictionary reference.

| case | one job |
|---|---|
| `bbpaVerify` | every dictionary control, hand-checkable (320 cells, ~10 s) |
| `planarPoiseuille` | cross-cycle variance + restart (developing flow, non-zero variance) |
| `oscillatingLidCavity` | field types + a derived boundary field, per-field products |
| `stenosisPipe` | parallel production; serial == parallel; IPA comparison |
| `bbpaVerifyLES` | adaptive time stepping + second-order statistics |

## Verification

`tests/` is an automated ground-truth suite (`./Alltest`, or
`python3 -m pytest tests/`): a coded function object prescribes an analytic
field, a bit-faithful Python replica of the accumulator replays the logged
trajectory, and closed-form references check every product. 12 test files
assert, at 1e-12 or bit-exact: the arithmetic, output-layout invariance
(all layouts, identical numbers), restart == continuous, variable-dt
behaviour and weight conservation, every dictionary option, per-field
product manifests, serial == parallel, the strict restart metadata, the
memory reduction, and that 19 invalid dictionaries abort with messages
naming the entry.

## Implementation notes

- One `phaseClock` owns every which-bin/which-cycle decision (BBPA and
  IPA); `binIntegrator`/`cycleStatistics` hold the current-cycle and
  cross-cycle state; a const-access `writer` serialises and a `state`
  component inverts it on restart — the output layer cannot touch a
  statistic by construction.
- Products gate by **allocation**: an off product never allocates its
  accumulators, so the same sparse-bin guards skip it everywhere.
- Lazy field lookup: fields registered by other function objects
  (`wallShearStress`) are picked up when they appear and restored against
  the original checkpoint instant.

## License

GPL-3.0 (same as OpenFOAM). See `LICENSE`.

## Citation

```bibtex
@article{wang2026bbpa,
  author  = {Wang, Jie and Skillen, Alex and Revell, Alistair},
  title   = {Bin-based in-situ phase averaging for pulsatile internal-flow {LES}: error budget, {OpenFOAM} implementation, and validation},
  journal = {Computers \& Fluids},
  year    = {2026},
  note    = {submitted}
}
```
