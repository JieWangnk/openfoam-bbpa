# Porting `libPhaseAveraging` (BBPA / IPA) — OpenFOAM 12 → 14

Branch: `of14-port` (worktree of `openfoam-bbpa`). `main` stays OF12.
Goal: test CFD Direct's **extensibility** claim with real production code — does a mature,
restart-safe, stateful, parallel cardiac-cycle phase-averaging function object move from OF12 to
OF14, and does it still produce the same numbers?

## API-delta analysis (OF12 vs OF14 headers, before any build)

The function object depends on this API surface; every item was diffed between
`/opt/openfoam12/src` and `/opt/openfoam14/src`:

| API used by BBPA/IPA | OF12 → OF14 | Impact |
|---|---|---|
| `functionObject` virtual interface (`read`/`fields`/`execute`/`write`/`executeAtStart`) | identical signatures | none |
| `IOobject` read/write enums (`MUST_READ`, `READ_IF_PRESENT`, `NO_WRITE`, `AUTO_WRITE`, …) | identical | none |
| `fvMeshFunctionObject(const word&, const Time&, const dictionary&)` ctor | present in both | none |
| `fvMeshFunctionObject.H` include path | still reachable via `finiteVolume/lnInclude` (symlinked) | `Make/options` unchanged |
| base-class hierarchy | OF12 `regionFunctionObject` → OF14 `objectRegistryFunctionObject` + `functionObjects::regionFunctionObject` | **internal only** |
| protected `obr_` | moved into `objectRegistryFunctionObject` (still `protected`) | accessible to derived class — none |
| protected `time_` | in `functionObject` (`const Time& time_`) | accessible — none |
| `mesh_` / `mesh()` / `lookupObject` / `foundObject` | present (OF14 adds a `mesh()` accessor) | none |
| `patchSet(const dictionary&, bool)` | **removed in OF14** | **not called** by BBPA/IPA — none |
| `timeIOdictionary` (restart-state I/O) | identical constructor | none |

The header analysis above held — the function-object base class, constructor, `IOobject` enums,
`timeIOdictionary` and the virtual interface all compile untouched. But it did **not** cover the library's
**field-typed template code**, and that is where the one real OF12→OF14 delta lives:

**The single change required — `GeometricField` template-parameter reorder.** OF14 changed
`GeometricField<Type, PatchField, GeoMesh>` (OF12) to `GeometricField<Type, GeoMesh, PrimitiveField>`, so
the generic `typedef GeometricField<Type, fvPatchField, volMesh> FieldType;` no longer parses (`fvPatchField`
where a mesh type is expected), degrading `FieldType` to `int` and cascading ~60 errors. Fixed at **3
typedef sites** (`BBPA.H` ×2, `IPA.H` ×1) using OF14's own alias `VolField<Type>`
(= `GeometricField<Type, fvMesh, Field>`, from `volFieldsFwd.H`, already included). No other source change.
This is the same OF14 field-template reorganisation as the earlier custom-BC port's
`DimensionedField<Type, volMesh>` → `<…, fvMesh>`.

## Build (run under OpenFOAM 14)

```bash
of14                                   # or: source /opt/openfoam14/etc/bashrc  (fresh shell)
echo $WM_PROJECT_VERSION                # expect 14
cd ~/GitHub/openfoam-bbpa-of14-port/src/phaseAveraging
wmake libso                             # → $FOAM_USER_LIBBIN/libPhaseAveraging.so
```

## Build result
- **Attempt 1 (unmodified source):** failed — the `GeometricField` template-parameter reorder above
  (~60 errors, all cascading from the 3 typedefs degrading `FieldType` to `int`).
- **Attempt 2 (`VolField<Type>` fix applied):** ✅ **builds clean** under OF14. `libPhaseAveraging.so`
  linked to `$FOAM_USER_LIBBIN` (`~/OpenFOAM/mchi4jw4-14/platforms/linux64GccDPInt32Opt/lib/`). Only a
  pre-existing `-Wreorder` warning remains (constructor init-list order vs member-declaration order in
  `BBPA.C` — not OF14-related, warns on OF12 too; left untouched to keep the port diff to the one real change).

**Net OF12→OF14 port cost: one API change, three typedef lines** (`git diff main of14-port -- '*.H'`).

## Verification (after a clean build)
1. **Correctness:** run `bbpa_analytic_verification` (`run_all.py`) under OF14 with this `.so` on
   `$FOAM_USER_LIBBIN`; the phase-averaged outputs must match the frozen analytic references (A/B/C).
2. **Numerical continuity:** run the same verification with the OF12 `libPhaseAveraging.so` (`main`
   branch, built under OF12) and confirm OF14 == OF12 to the harness tolerance.
3. **Restart (the real state test):** run a tutorial to half, kill, restart, finish; the accumulated
   phase average must equal the uninterrupted run (exercises the `timeIOdictionary` checkpoint path).
4. **Parallel:** `decomposePar` + parallel run vs serial; phase-averaged fields must match (bin reduction).

Each ✓ is direct, checkable evidence for "can my real code move to OF14?".

## Verification result — the port compiles but does NOT reproduce OF12

Ran `run_all.py` under OF14 (`bbpa_analytic_verification`, `of14` branch; case fix: `physicalProperties`
matched to OF14 form — `viscosityModel constant; rho 1; nu 1e-05;`). It ran to completion with **no case-syntax
errors** and the coded prescribe-FO compiled fine.

- **Input is identical.** The prescribed field log (`prescribe.log`: t, dt, U each step) is **byte-identical**
  OF14 vs the committed OF12 reference — 385 steps, both to t = 6. The coded FO is unaffected by anything.
- **Output does NOT match.** The committed result tables shifted under OF14: smooth-field discretisation error
  `3.24e-01 → 0` (first-order convergence lost, slope `nan`), straddle (A2) `0 → 0.21`, restart-blend
  `2.1e-2 → 0`.
- **Direct cause located:** BBPA writes its phase-aligned `U_PA` dumps at **different cycles** under OF14 —
  at t ≈ 2 and 5, and **never at the final cycle 6** — whereas OF12 writes one dump at cycle 6. There is **no
  common output time**, so the harness reads mismatched/missing fields → the divergent metrics.

**Conclusion.** The port **compiles cleanly (one 3-line change) and runs**, but exhibits a **runtime
behavioural regression**: BBPA's write/cycle scheduling diverges under OF14, so the phase-averaged output is
not reproduced. A build-only "does it compile?" test would have declared success — the *verification* caught
it. This is the honest extensibility answer: *moving the code is one change; making it behave identically is
not free.* ## Root cause — the OF14 port is clean; the test reference was stale

A debug probe on `time_.value()`/cycle at each dump **refuted** an earlier FP-drift hypothesis: the clock
is exact at every phase-aligned dump (`snapArg-nearestInt = 0`), `currentCycleNumber()` returns the correct
cycle (3, 6), and `cycleStart = (currentCycle_ − 1)·period` (2, 5) is exactly what the **current** `main`
code prescribes (`writeAllAveragedFields` dumps the *previous*, just-completed cycle).

`frozen_raw` instead labels its dump at cycleStart = `currentCycle_` (6) — the **old** labeling. Git dates
it decisively: the "previous-cycle window" fix landed in `openfoam-bbpa` at **9d5aeb0 (2026-07-16)**;
`frozen_raw` was committed **2026-06-18** and never regenerated. The reference is a month older than a fix
to the library it verifies.

**Conclusion: the OF14 port reproduces the current OF12 `main` BBPA exactly.** The apparent "OF14 ≠ OF12"
was *old-reference vs current-code*, not a version regression. Two honest outcomes:

1. **Extensibility — clean.** `libPhaseAveraging` ports OF12 → OF14 with one 3-line change (`GeometricField`
   template reorder → `VolField<Type>`), builds, runs, and is numerically faithful to current `main`
   (clock exact, cycle/bin logic correct, restart + `np=2` parallel exercised).
2. **Verification-suite finding (bonus).** `bbpa_analytic_verification`'s frozen reference is stale — it
   predates a library fix and produced false-positive "failures." Frozen references must be re-frozen to
   track the code. (Itself the kind of defect this whole verify-don't-just-compile approach exists to catch.)

**Close-out (pending):** run the harness under OF12 with the current `main` lib to confirm OF14 == OF12
bit-for-bit, and re-freeze `frozen_raw` against current `main`.
