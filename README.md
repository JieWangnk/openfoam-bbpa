# openfoam-bbpa

**Bin-Based Phase Averaging for OpenFOAM** — an in-situ function object that
produces phase-averaged mean, cross-cycle variance, and resolved turbulent
kinetic energy directly in memory, without writing terabytes of transient
snapshots.

For the theory, error analysis, and cardiovascular LES validation see the
companion paper:

> Wang, J. *Bin-Based Phase Averaging of turbulence-resolving pulsatile flows.*
> (2026, in preparation). Paper repo:
> [`Bin-Phase-Average`](https://github.com/JieWangnk/Bin-Phase-Average)

## Features

| | BBPA | `fieldAverage` (OpenFOAM-native) | Classical offline PA |
|---|---|---|---|
| Phase-centred ensemble mean | ✅ | ✗ (causal running filter) | ✅ |
| In-situ (no snapshot I/O) | ✅ | ✅ | ✗ (terabytes of writes) |
| Phase-resolved TKE | ✅ (`_PATKE`) | ✗ | ✅ (offline) |
| Cross-cycle variance | ✅ (Welford) | partial (`UPrime2Mean`) | ✅ |
| Wall shear stress binning | ✅ (lazy init) | ✅ | ✅ |
| Output per bin | phase-aligned time directories | single time dir | per-phase dir |
| Memory | `O(I)` per field | `O(1)` per field | n/a |
| Overhead (6M-cell LES, 160 cores) | **+0.13%** at I=50 | baseline | +100%+ (post-run) |

Accompanying `IPA` (Instantaneous Phase Averaging) function object samples at
exact phase crossings without binning — a zero-width limit of BBPA for
reference comparisons.

## What gets written

For each tracked vector/scalar field `φ` (e.g. `U`, `p`, `wallShearStress`),
at every solver writeTime, BBPA produces three fields in each
phase-aligned time directory `t_j = t_0 + j · (T/I)`:

| File | Contents | Units |
|---|---|---|
| `<φ>_PA` | Phase-average `⟨φ⟩_j` | same as `φ` |
| `<φ>_PAVariance` | Cross-cycle variance `M₂/(N−1)` of the bin mean | `dim(φ)²` |
| `<φ>_PATKE` | Resolved phase-ensemble TKE `½(⟨\|φ\|²⟩_j − \|⟨φ⟩_j\|²)` | `dim(φ)²` |

For a vector field the `_PATKE` output is the turbulent kinetic energy
entering the triple decomposition. For a scalar field it reduces to the
within-bin sub-timestep variance.

## Requirements

- OpenFOAM Foundation **v12** (openfoam.org)
- GCC ≥ 9, OpenMPI / equivalent

Tested on:
- Ubuntu 24.04 (system OpenFOAM 12)
- csf3 (`apps/gcc/openfoam/12`, GCC 13.3.0)
- csf4 (`openfoam/12-foss-2023a`, GCC 12.3.0)

## Build

```bash
git clone https://github.com/JieWangnk/openfoam-bbpa.git
cd openfoam-bbpa/src/phaseAveraging
wmake libso
# writes $FOAM_USER_LIBBIN/libPhaseAveraging.so
```

## Quickstart — oscillating-lid cavity tutorial

```bash
cd tutorials/oscillatingLidCavity
./Allrun
paraview case.foam &
# Time directories 0.0, 0.1, ..., 2.9 contain U_PA, p_PA, U_PATKE, etc.
```

This 400-cell laminar case ends in ≈1 min and exercises the full pipeline:
Strategy B phase-aligned output, Welford variance, TKE accumulator.

## Usage in your own case

In `system/controlDict`:

```c
libs ("libPhaseAveraging.so");

functions
{
    #includeFunc wallShearStress    // if tracking WSS: list BEFORE bbpa

    bbpa
    {
        type            BBPA;
        libs            ("libPhaseAveraging.so");
        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;       // fire on solver writeTimes
        fields          (U p wallShearStress);
        period          1.0;             // cardiac cycle (s)
        nBins           50;              // I
        startCycle      2;               // skip transient
    }
}
```

Optional keys:
- `startTime <s>` — explicit phase origin; defaults to `floor(t/period)*period`
- `cycles <N>` — stop averaging after N cycles; defaults to unlimited

See `src/phaseAveraging/exampleDict_BBPA` and `exampleDict_IPA` for full
dictionaries.

## Scripts

| `scripts/bench_overhead.py` | Parse matched-pair logs (with/without BBPA) and report per-step walltime + overhead % |
| `scripts/screenshot_WSS_BPM120.py` | pvbatch renderer for WSS-filter-width 3-panel figure (paper Fig. X) |

## Implementation notes

- **Welford M₂** for cross-cycle variance — numerically stable; avoids catastrophic cancellation when cycle-to-cycle differences are small.
- **Lazy field lookup**: the field list is re-checked every `execute()`, so fields registered by other function objects (e.g. `wallShearStress`) are picked up even if BBPA is constructed before them.
- **Strategy B output**: each bin writes to its own phase-aligned time directory, so `foamListTimes` / `foamToVTK` handle the output natively without special postprocessing.
- **Partial-cycle TKE**: `_PATKE` fires even on a single-cycle run (`binCounts_ == 0`) by falling through to the in-progress accumulators — eliminates endpoint-alignment gotchas.

## License

GPL-3.0 (same as OpenFOAM). See `LICENSE`.

## Citation

If you use BBPA in your work, please cite the paper:

```bibtex
@article{wang2026bbpa,
  author = {Wang, Jie},
  title  = {Bin-Based Phase Averaging of turbulence-resolving pulsatile flows},
  year   = {2026},
  note   = {in preparation}
}
```
