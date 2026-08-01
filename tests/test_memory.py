#!/usr/bin/env python3
"""
test_memory.py -- direct evidence that per-field products deliver the
intended memory reduction (review-matrix row).

Two construction-heavy runs of the prescribe case on an enlarged mesh
(1600 cells) with nBins = 500 and a vector field: secondMoment +
crossCycleVariance ON (500 symmTensor accumulator pairs + 500 M2 fields,
~2 x 77 KB + 13 KB per bin) vs OFF (none of them). Peak RSS is measured
with /usr/bin/time -v; the run is 2 timesteps -- allocation happens at
construction.

Asserts the ON run's peak RSS exceeds the OFF run's by a conservative
margin (> 20 MB; the expected delta is ~80 MB), and prints both so the
number is on record.

Runnable standalone (exit 1 on failure, for Alltest) and via pytest.
"""
from __future__ import annotations
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bbpa_verify_lib as rv

# No prescribe FO here: its codeExecute compilation forks g++, whose own
# ~300 MB peak RSS would swamp the measurement. The solver + BBPA alone.
CONTROL = """/*--------------------------------*- C++ -*----------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}}
application     foamRun;
solver          incompressibleFluid;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         0.02;
deltaT          0.01;
writeControl    runTime;
writeInterval   1.0;
purgeWrite      0;
writeFormat     ascii;
runTimeModifiable false;
adjustTimeStep  no;
functions
{{
    bbpa
    {{
        type            BBPA;
        libs            ("libPhaseAveraging.so");
{body}    }}
}}
"""

N_BINS = 500
MESH = "(40 40 1)"                     # 1600 cells: per-bin fields dominate
MARGIN_KB = 20*1024                    # > 20 MB (observed delta ~80 MB)


def _body(products):
    return f"""        fields          ( U {products} );
        period          1.0;
        nBins           {N_BINS};
        discardCycles   0;
        phaseOrigin     0;
        output {{ layout phaseAligned; cycleWriteInterval 1; }}
        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;
"""


def _peak_rss_kb(name, products):
    rundir, _ = rv.make_case(name, "bin_constant", 0.01, 0.02, 1.0)
    with open(os.path.join(rundir, "system", "controlDict"), "w") as f:
        f.write(CONTROL.format(body=_body(products)))
    # enlarge the mesh so the per-bin accumulators dominate the footprint
    bmd = os.path.join(rundir, "system", "blockMeshDict")
    m = open(bmd).read().replace("(4 4 1)", MESH)
    open(bmd, "w").write(m)
    assert rv.block_mesh(rundir).returncode == 0, name + ": blockMesh failed"
    r = rv.sh("/usr/bin/time -v foamRun > log.foamRun 2> log.time", rundir)
    assert r.returncode == 0, name + ": foamRun failed"
    txt = open(os.path.join(rundir, "log.time")).read()
    m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", txt)
    assert m, name + ": no RSS in time output"
    return int(m.group(1))


def run_check():
    os.makedirs(rv.RUNROOT, exist_ok=True)
    on_kb = _peak_rss_kb("mem_on", "")
    off_kb = _peak_rss_kb(
        "mem_off", "{ secondMoment off; crossCycleVariance off; }")
    delta = on_kb - off_kb
    report = (f"peak RSS, nBins={N_BINS}, vector U: all products "
              f"{on_kb/1024:.1f} MB; secondMoment+variance off "
              f"{off_kb/1024:.1f} MB; saved {delta/1024:.1f} MB")
    errors = []
    if delta < MARGIN_KB:
        errors.append(f"expected > {MARGIN_KB} kB saved, got {delta} kB "
                      f"({report})")
    return errors, report


def test_memory():
    errors, _ = run_check()
    assert not errors, "\n".join(errors)


if __name__ == "__main__":
    errors, report = run_check()
    print(report)
    if errors:
        print("FAIL:")
        for e in errors:
            print("  ", e)
        sys.exit(1)
    print("memory: PASS")
    sys.exit(0)
