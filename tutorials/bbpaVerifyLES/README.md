# bbpaVerifyLES — periodic + turbulent BBPA verification

Third and most complete of the BBPA verification cases. It is the only one that
exercises the **phase mean, the raw second moment and the cross-cycle variance
at the same time**, in the regime the paper actually runs in (periodic forcing
with genuine turbulent cycle-to-cycle variation), and the only one that checks
**boundary** values rather than just the internal field.

## What it caught

`<field>_PA_UU` was written with its **boundary field never populated** —
every patch came out `value uniform (0 0 0 0 0 0)`. The accumulation was fine
(`+=` on a `GeometricField` updates the boundary too); the value existed in
memory and was dropped on write, because `paUU` was built as a Zero field and
then only `primitiveFieldRef()` was assigned.

For `wallShearStress` this made the field **entirely useless rather than merely
incomplete**: WSS has a uniform `(0 0 0)` internalField, so *all* of its data
lives on the wall patch. `wallShearStress_PA_UU` was therefore zero everywhere,
and any WSS fluctuation intensity or Reynolds stress derived from it would be
silently zero rather than obviously broken.

The same `primitiveFieldRef()` pattern was in `restoreBinItemFromDisk()`, so
**every restart discarded the boundary phase-average** and then folded new
cycles against a phantom zero history, biasing the result low — while the
Properties counters still claimed the cycles were all there.

A checker that only compares `internalField` cannot see any of this: it
"passes" the WSS test by comparing zeros against zeros. `check_bbpa.py` here
compares boundary patches and reports **INCONCLUSIVE rather than PASS** if it
only ever compared ~zero values.

## Why this case exists

The other two cannot do this job:

- **bbpaVerify** (clean sine inlet) is deterministic and periodic: once
  converged every cycle is identical, so the cross-cycle variance is ~0 by
  construction and cannot catch a variance bug.
- **planarPoiseuille** is a developing flow: it exercises the variance but is
  not periodic-with-turbulence.

Here the inlet `boundaryData` record is a periodic base profile **plus a ~0.070
rms perturbation that genuinely differs every cycle**, so all three outputs are
non-trivial at once. Verified properties of the record:

- period is exactly **1.0** — rms mismatch of the mean waveform at lag 1.0 is
  0.0034 (0.4% of amplitude) vs 0.565 at lag 0.5 (anti-phase).
- it is **not** an identical replay: same-phase fields one cycle apart differ by
  ~0.070 rms (correlation 0.99 — same structure, different fluctuation).

Because the variance is **inlet-driven**, it survives coarsening the mesh from
the parent `pulsatilePipeLES` (400 20 20)=160k to (40 6 6)=1440 cells. The
interior does not need to sustain turbulence for the test to be meaningful.

## Instrument settings — do NOT copy these into production runs

These three exist to make the check possible and would be wrong in a real run:

| setting | why |
|---|---|
| `writePrecision 12` | the default 6 writes ~6 significant digits, putting a **~5e-7 noise floor** under the comparison — 5 orders of magnitude above the 5e-12 agreement being detected. A real 1e-8 bug would hide inside the file's own rounding. |
| `writeInterval 1` (every step) | the FO folds **every** solver step into its bin. Re-deriving that independently needs every state the FO saw; saving every 10th step would make the "truth" a *different* average and disagree by O(1%) — looking like a bug that isn't one. This is why the mesh is tiny: per-step saving still costs ~1.3 GB for 3 cycles. |
| `period 0.9997` (true period is 1.0) | keeps every bin edge **and** cycle boundary off the 1e-3 timestep grid. On-grid edges make the FO and the checker disagree about which bin a step belongs to — an artifact that reads as a defect. The 0.03% offset drifts the phase ~0.1% of a period over 3 cycles. |

Note the awkward bin-dir names (`0.142814285714`, …) follow from `period/nBins`,
**not** from `timePrecision`. In production they are round: the aorta at
`period 0.5, I=50` bins at 0.01, 0.02, 0.03 — exactly the solver write times.
At `I=200` they fall at 0.0025 spacing, i.e. 4x more directories than write
folders, which is what `writeMode companion` / the centred-window design
address.

## Running

```sh
blockMesh
foamRun

# mean + second moment + variance, internal AND inlet
python3 check_bbpa.py --field U --period 0.9997 --nbins 7 \
                     --start-cycle 0 --include-t0 --tol 1e-8

# the wall patch, where WSS actually lives
python3 check_bbpa.py --field wallShearStress --patches walls \
                     --period 0.9997 --nbins 7 --start-cycle 0 --include-t0
```

`--include-t0` is required and is not a fudge: the FO's `execute()` is called at
the initial condition and folds it with a full dt. With the flag **every** bin
agrees to round-off; without it only bin 0 of cycle 0 is contaminated (~2e-3 in
the 3-cycle mean). A wrong hypothesis would not land all 7 bins at 5e-12.

## Result (fixed library)

Same-trajectory identity check, so agreement is round-off, not "close enough":

| field | region | `_PA` | `_PA_UU` | `_PAVariance` |
|---|---|---|---|---|
| U | internal | 5.3e-12 | 5.7e-12 | 1.3e-13 |
| U | inlet | 4.9e-12 | 5.7e-12 | 2.1e-14 |
| wallShearStress | walls | 5.8e-16 | 7.3e-19 | 4.0e-19 |

Before the fix, `U` inlet `_PA_UU` was **3.561e+00** and
`wallShearStress_PA_UU` was identically zero.

### Restart safety

Running 0→2.0, stopping, and restarting 2.0→3.0 reproduces the continuous
0→3.0 run to **1.0e-11** (the `writePrecision 12` file-rounding floor for
values of O(2)), including `wallShearStress_PA` on the wall patch at 1.5e-15.
The restart restores 7/7 bins for both tracked fields and recovers
`startTime=0` from the Properties checkpoint rather than the auto-aligned
restart time.
