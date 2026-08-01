#!/usr/bin/env python3
"""
test_invalidInput.py -- negative tests (roadmap section 7): BBPA must abort,
with a message naming the offending entry, for each invalid dictionary below.
A configuration error that runs silently is itself the defect.

Rows (expected message fragment in parentheses):
  missing fields               (keyword fields is undefined)
  empty fields ()              ('fields' is empty)
  missing period               (keyword period is undefined)
  period 0                     ('period' must be positive)
  nBins 0                      ('nBins' must be positive)
  binsOfInterest out of range  (out of range)
  binsOfInterest duplicate     (listed more than once)
  unknown writeMode            (unknown writeMode)
  binsAtWriteTimes with a writeInterval that period/nBins does not divide
                               (integer multiple)

Each case is a minimal controlDict with ONLY the BBPA FO (no prescribe FO --
the abort happens at construction, before any physics matters).

Runnable standalone (exit 1 on failure, for Alltest) and via pytest.
"""
from __future__ import annotations
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bbpa_verify_lib as rv

CONTROL = """/*--------------------------------*- C++ -*----------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     foamRun;
solver          incompressibleFluid;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         0.02;
deltaT          0.01;

writeControl    runTime;
writeInterval   {wi};
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
{body}
    }}
}}

// ************************************************************************* //
"""

BASE = """        fields          (U);
        period          1.0;
        nBins           8;
        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;
"""

ROWS = [
    ("missing_fields",
     "        period 1.0;\n        nBins 8;\n",
     "keyword fields is undefined", 1.0),
    ("empty_fields",
     "        fields ();\n        period 1.0;\n        nBins 8;\n",
     "'fields' is empty", 1.0),
    ("missing_period",
     "        fields (U);\n        nBins 8;\n",
     "keyword period is undefined", 1.0),
    ("zero_period",
     "        fields (U);\n        period 0;\n        nBins 8;\n",
     "'period' must be positive", 1.0),
    ("zero_nBins",
     "        fields (U);\n        period 1.0;\n        nBins 0;\n",
     "'nBins' must be positive", 1.0),
    ("bin_out_of_range",
     BASE + "        binsOfInterest (99);\n",
     "out of range", 1.0),
    ("bin_duplicate",
     BASE + "        binsOfInterest (2 2);\n",
     "listed more than once", 1.0),
    ("unknown_writeMode",
     BASE + "        writeMode comapnion;\n",
     "unknown writeMode", 1.0),
    ("incommensurate_binsAtWriteTimes",
     BASE + "        binsAtWriteTimes yes;\n",
     "integer multiple", 0.3),      # 0.3 / (1/8) is not an integer
    # --- step-6 rename rows ---
    ("old_and_new_name_conflict",
     BASE + "        discardCycles 1;\n        startCycle 1;\n",
     "deprecated spelling", 1.0),
    ("unknown_phaseSelection_mode",
     BASE + "        phaseSelection { mode everything; }\n",
     "unknown phaseSelection mode", 1.0),
    ("phaseSelection_and_old_bins",
     BASE + "        phaseSelection { mode list; bins (2); }\n"
            "        binsOfInterest (2);\n",
     "not both", 1.0),
    ("unknown_binAlignment",
     BASE + "        binAlignment middle;\n",
     "unknown binAlignment", 1.0),
    # --- step-7 per-field product rows ---
    ("unknown_field_option",
     "        fields ( U { secondMomnet off; } );\n"
     "        period 1.0;\n        nBins 8;\n",
     "unknown per-field option", 1.0),
    ("old_phaseMean_option",
     "        fields ( U { phaseMean off; } );\n"
     "        period 1.0;\n        nBins 8;\n",
     "unknown per-field option", 1.0),
]


def run_rows():
    errors = []
    os.makedirs(rv.RUNROOT, exist_ok=True)
    for name, body, fragment, wi in ROWS:
        rundir = os.path.join(rv.RUNROOT, "neg_" + name)
        if os.path.exists(rundir):
            shutil.rmtree(rundir)
        shutil.copytree(rv.TEMPLATE, rundir)
        with open(os.path.join(rundir, "system", "controlDict"), "w") as f:
            f.write(CONTROL.format(body=body, wi=wi))
        assert rv.block_mesh(rundir).returncode == 0, name + ": blockMesh failed"
        r = rv.sh("foamRun > log.foamRun 2>&1", rundir)
        log = open(os.path.join(rundir, "log.foamRun")).read()
        if r.returncode == 0:
            errors.append(f"{name}: ran to completion; expected an abort")
        elif fragment not in log:
            tail = "\n".join(log.splitlines()[-15:])
            errors.append(f"{name}: aborted, but without the message "
                          f"'{fragment}'. Log tail:\n{tail}")
    return errors


def test_invalidInput():
    errors = run_rows()
    assert not errors, "\n".join(errors)


if __name__ == "__main__":
    errors = run_rows()
    if errors:
        print("FAIL:")
        for e in errors:
            print("  ", e)
        sys.exit(1)
    print(f"invalid-input rows: PASS ({len(ROWS)} configurations abort with "
          "the expected message)")
    sys.exit(0)
