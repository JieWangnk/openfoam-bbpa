# The BBPA study — every dictionary entry, explained and demonstrated

This folder teaches the BBPA function object (v2, OpenFOAM 14) one
dictionary entry at a time. Each chapter explains one top-level entry in
plain language. Sub-entries are covered inside their parent's chapter.
Most chapters end with a small runnable demo: a case that changes one
thing and a table (a "manifest") showing what that change did to the
files on disk. The manifests in the chapters are real program output,
committed under `expected/` and checked against reruns.

The meaning of every entry was taken from the source
(`src/phaseAveraging/BBPA.C` and `BBPAWriterTemplates.C`), and every
claim about files on disk comes from a demo run.

## Learning path

| Chapter | Entries covered |
|---|---|
| [01 — first run](01-first-run.md) | the smallest working dictionary, on a real tutorial |
| [02 — what gets written](02-what-gets-written.md) | the output files and how to read them |
| [03 — type and libs](03-type-and-libs.md) | `type`, `libs` |
| [04 — scheduling](04-scheduling.md) | `executeControl`, `executeInterval`, `writeControl` |
| [05 — fields](05-fields.md) | `fields` + `writePhaseMean`, `secondMoment`, `crossCycleVariance`, `writePerCycle` |
| [06 — the phase grid](06-phase-grid.md) | `period`, `nBins`, `phaseOrigin`, `binAlignment`, `timeQuadrature` (flat or in `phaseGrid{}`) |
| [07 — the cycle window](07-cycle-window.md) | `discardCycles`, `maxCycles` |
| [08 — accumulation selection](08-accumulation-selection.md) | `accumulationSelection` (`mode`, `bins`) |
| [09 — output selection](09-output-selection.md) | `outputSelection` (`mode`, `bins`) |
| [10 — output](10-output.md) | `output` + `layout`, `cycleWriteInterval`, `meanWriteIntervalCycles`, `writeFinal`, `writeProperties` |
| [11 — field availability](11-field-availability.md) | `fieldAvailability` |

That is all 26 entries of the v2 dictionary. Chapters 03 and 04 have no
demo of their own: those entries decide whether the object loads and
when it runs, and every demo in the other chapters already uses them.
The scheduling entries in chapter 04 are standard OpenFOAM
function-object controls; they are explained once, there.

## Running the demos

The demos need OpenFOAM 14 and the built library:

```bash
source /opt/openfoam14/etc/bashrc
cd src/phaseAveraging && wmake libso && cd -   # once
```

Then, from this folder:

```bash
./run_demos.py --list            # what exists
./run_demos.py --only d01_minimal
./run_demos.py                   # run everything (~1 minute), print manifests
./run_demos.py --check           # rerun and compare against expected/
```

Each demo is a tiny 8-cell case in which a coded function object
overwrites `U` every step with a known analytic signal (the test
suite's case machinery, `tests/bbpa_verify_lib.py`). Period 1 s,
8 bins, 5 cycles. Cases land in `tests/_runs/` (gitignored), so
nothing here gets dirty.

`./run_demos.py --check-docs` needs no OpenFOAM: it verifies that the
manifest printed in each chapter matches its file in `expected/` byte
for byte, and that demos and chapters map one to one.

## How to read a manifest

```
    file                  dirs  first    last
    U_PA                    40  0        4.875
```

means: a file named `U_PA` exists in 40 time directories; the earliest
is `0/` and the latest is `4.875/`. A count of 0 means the file exists
nowhere, which is often the point of the demo. The dictionary body that
produced the table is printed above it.
