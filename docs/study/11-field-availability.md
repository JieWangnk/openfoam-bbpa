# 11 — `fieldAvailability`: fields that appear late, or never

```c
fieldAvailability lazy;      // default   | strict
```

Derived fields such as `wallShearStress` are created by their own
function objects, which must be listed before BBPA in `functions{}`
(function objects run in listed order):

```c
functions
{
    #includeFunc wallShearStress    // the producer, BEFORE bbpa
    bbpa { ... fields (U wallShearStress); ... }
}
```

`fieldAvailability` decides what happens when a requested field is
missing when accumulation should begin.

* `lazy` (default): keep going. BBPA retries the lookup each step and
  starts accumulating the field when it appears. Each field's ensemble
  is independent: its own cycle count covers only the cycles it
  actually observed, so a late field simply averages fewer cycles.
* `strict`: abort with an error naming the field, checked at the first
  retained accumulation interval (after `discardCycles`). Use this in
  production dictionaries: a misspelt field name, or a producer listed
  after BBPA, then stops the run at once instead of quietly producing
  nothing for that field.

## Demo: a field with no producer at all

`study_d11_field_availability` requests
`fields (U wallShearStress)` in a case that never creates
`wallShearStress`. Both runs share that mistake; only the mode
differs.

<!-- demo:study_d11_field_availability -->
```text
demo: study_d11_field_availability
chapter: 11-field-availability.md

case: study_d11_field_availability_lazy  (endTime 5, solver writeInterval 1)
BBPA dictionary body:
    fields          (U wallShearStress);
    period          1.0;
    nBins           8;
    discardCycles   0;
    phaseOrigin     0;
    fieldAvailability lazy;
    executeControl  timeStep;
    executeInterval 1;
    writeControl    writeTime;
run: completed

    file                            dirs  first    last
    U_Bin                              0  -        -
    U_Bin_UU                           0  -        -
    U_PA                              40  0        4.875
    U_PAVarianceMag                   32  1        4.875
    U_PA_UU                           40  0        4.875
    U_PAmeanState                      0  -        -
    wallShearStress_PA                 0  -        -
    wallShearStress_PAVarianceMag      0  -        -
    wallShearStress_PA_UU              0  -        -
    uniform/*Properties                6
folded cycles at the newest checkpoint: U=5
deprecation warnings in the log: 0

case: study_d11_field_availability_strict  (endTime 5, solver writeInterval 1)
BBPA dictionary body:
    fields          (U wallShearStress);
    period          1.0;
    nBins           8;
    discardCycles   0;
    phaseOrigin     0;
    fieldAvailability strict;
    executeControl  timeStep;
    executeInterval 1;
    writeControl    writeTime;
run: FATAL (expected). The error:
    fieldAvailability strict: field 'wallShearStress' is not registered at the first retained accumulation interval. Ensure its producing function object runs before BBPA, or use fieldAvailability lazy.
```

Under `lazy`, the run completes: `U` gets its full set of products,
every `wallShearStress` row is 0, and the checkpoint's folded-cycle
list names `U` only. Nothing in the log stops you from wasting a long
run this way, which is the argument for `strict`: the second run dies
immediately, and its error names the missing field and both fixes
(list the producer before BBPA, or accept lazy behaviour).
