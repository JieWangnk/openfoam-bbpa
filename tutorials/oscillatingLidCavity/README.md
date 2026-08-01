# oscillatingLidCavity — field types and a derived field

**One job: prove BBPA handles every supported field type, including a
DERIVED wall field, with independent per-field products.** A sine-driven
lid (period 1 s) over 3 cycles:

* `U` — a `volVectorField`, every product (mean, `<U U>` for Reynolds
  stress, cross-cycle variance);
* `p` — a `volScalarField`, mean only (`secondMoment off;
  crossCycleVariance off;` — no `p_PA_UU`/`p_PAVariance` files, and no
  accumulator memory for them);
* `wallShearStress` — produced by its own function object *before* BBPA,
  with all of its data on the wall **boundary**: the case that catches
  internal-field-only bugs (a phase-averaged WSS whose patches read
  `uniform (0 0 0)` is broken, not converged). Mean + second moment, no
  variance.

Run:

```bash
./Allrun
```

Then look at any phase-aligned directory (`0.1 … 2.9`): `U_PA`, `p_PA`,
`wallShearStress_PA` (+ `_PA_UU` where the product is on) — and check the
`wallShearStress_PA` **boundaryField** on `movingWall`/`fixedWalls` is
populated. `output { layout both; }` also drops the companion copies
beside the solver's own fields.

Dictionary reference: `../bbpaVerify/README.md`. Machine-precision ground
truth: `../../tests/`.
