# 03 — `type` and `libs`

```c
bbpa
{
    type            BBPA;
    libs            ("libPhaseAveraging.so");
    ...
}
```

These two entries connect the dictionary to the compiled code.

## `type`

`type BBPA;` names the function-object class to construct. The name is
registered by the library at load time. A companion class `IPA`
(instantaneous phase averaging, one sample per cycle per phase point)
lives in the same library and is selected the same way.

## `libs`

`libs ("libPhaseAveraging.so");` tells the solver to load that shared
library before constructing the object. The library comes from:

```bash
source /opt/openfoam14/etc/bashrc
cd src/phaseAveraging
wmake libso        # -> $FOAM_USER_LIBBIN/libPhaseAveraging.so
```

`$FOAM_USER_LIBBIN` is on the solver's library search path, so the
entry needs only the file name.

Two failure modes are worth recognising:

* Library missing or never built: the solver stops at startup with
  `Unknown function type BBPA`. Build it as above.
* Library built against a different OpenFOAM version: loading fails
  with unresolved symbols. Rebuild after sourcing the right
  `etc/bashrc`.

A top-level `libs (...)` line in `controlDict` (outside `functions{}`)
also works and loads the library once for all function objects; the
tutorial `stenosisPipe` shows that form. Keeping the entry inside each
BBPA block is self-contained and is what every other case here does.
