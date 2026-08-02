#!/usr/bin/env python3
"""
test_selection.py -- the three-control separation: the phase grid,
the bins that are CALCULATED (accumulationSelection), and the bins that
are WRITTEN (outputSelection).

Rows (same bin-constant prescribe trajectory as everywhere else):

  calcAllWriteSome   accumulation all + outputSelection list (2 5):
                     files exist ONLY for bins 2 and 5, their values equal
                     the closed form, and read() warns that unwritten bins
                     are not restart-recoverable. (The combination the old
                     single phaseSelection could not express.)
  grouped            the full three-block spelling -- phaseGrid{} +
                     accumulationSelection{} + outputSelection{} -- runs
                     with ZERO deprecation warnings and equals the closed
                     form.
  outputChangeOnRestart
                     restart with a WIDER outputSelection is allowed (it
                     gates serialisation only): segment 1 writes bins
                     (2 5); the restart writes everything accumulated. The
                     previously-written bins keep their history exactly
                     (== closed form over all 6 cycles); never-written bins
                     cold-start with a warning, accumulating only
                     post-restart cycles.

Runnable standalone (exit 1 on failure, for Alltest) and via pytest.
"""
from __future__ import annotations
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bbpa_verify_lib as rv
import bbpa_offline as bo

TOL = 1e-12
DT = rv.BIN_DT / 16
END_TIME = 6.0
FINAL_CS = 5.0


def _err(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def _body(selection):
    return f"""        fields          (U);
        period          1.0;
        nBins           8;
        discardCycles   0;
        phaseOrigin     0;
{selection}        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;
"""


GROUPED_BODY = """        fields          (U);
        phaseGrid
        {
            period          1.0;
            nBins           8;
            phaseOrigin     0;
            binAlignment    left;
        }
        discardCycles   0;
        accumulationSelection { mode all; }
        outputSelection       { mode accumulated; }
        output { layout phaseAligned; cycleWriteInterval 1; }
        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;
"""


def _run(name, body, end_time=END_TIME):
    os.makedirs(rv.RUNROOT, exist_ok=True)
    rundir, logpath = rv.make_case(name, "bin_constant", DT, end_time, 1.0,
                                   bbpa_body=body)
    assert rv.block_mesh(rundir).returncode == 0, name + ": blockMesh failed"
    assert rv.run_serial(rundir).returncode == 0, name + ": foamRun failed"
    return rundir, logpath


def _upa_at(rundir, cs, j):
    for tv, full in rv._time_dirs(rundir).items():
        if abs(tv - (cs + j * rv.BIN_DT)) < 1e-7:
            return os.path.isfile(os.path.join(full, "U_PA"))
    return False


def run_rows():
    errors = []
    closed = bo.closed_form_bin_constant(rv.A, rv.c_of_n, 6)

    # -- calculate all 8, write only (2 5) ------------------------------------
    sel = ("        outputSelection { mode list; bins (2 5); }\n")
    rd, _ = _run("sel_calcAllWriteSome", _body(sel))
    log = open(os.path.join(rd, "log.foamRun")).read()
    if "not restart-recoverable" not in log and "cold-start" not in log:
        errors.append("calcAllWriteSome: missing the sparse-output restart "
                      "warning at read()")
    got = rv.read_bbpa_output(rd, FINAL_CS)
    for j in range(rv.I):
        present = _upa_at(rd, FINAL_CS, j)
        if j in (2, 5):
            if not present:
                errors.append(f"calcAllWriteSome: bin {j} missing")
            elif _err(got[j]["PA"], closed[j]["PA"]) >= TOL:
                errors.append(f"calcAllWriteSome: bin {j} PA err "
                              f"{_err(got[j]['PA'], closed[j]['PA']):.3e}")
        elif present:
            errors.append(f"calcAllWriteSome: unselected bin {j} written")

    # -- the grouped three-block spelling -------------------------------------
    rd, _ = _run("sel_grouped", GROUPED_BODY)
    log = open(os.path.join(rd, "log.foamRun")).read()
    if "is deprecated" in log:
        errors.append("grouped: unexpected deprecation warning")
    got = rv.read_bbpa_output(rd, FINAL_CS)
    for j in range(rv.I):
        if j not in got:
            errors.append(f"grouped: bin {j} missing")
        elif _err(got[j]["PA"], closed[j]["PA"]) >= TOL:
            errors.append(f"grouped: bin {j} PA err "
                          f"{_err(got[j]['PA'], closed[j]['PA']):.3e}")

    # -- widening outputSelection across a restart is allowed -----------------
    selA = "        outputSelection { mode list; bins (2 5); }\n"
    c2, l2 = rv.make_case("sel_restart", "bin_constant", DT, 3.0, 3.0,
                          bbpa_body=_body(selA))
    assert rv.block_mesh(c2).returncode == 0
    assert rv.run_serial(c2).returncode == 0
    with open(os.path.join(c2, "system", "controlDict"), "w") as f:
        f.write(rv.control_dict("bin_constant", DT, END_TIME, 3.0,
                                "startTime", l2, start_time=3.0,
                                bbpa_body=_body("")))   # output: accumulated
    r = rv.sh("foamRun > log.foamRun2 2>&1", c2)
    log2 = open(os.path.join(c2, "log.foamRun2")).read()
    if r.returncode != 0:
        errors.append("outputChangeOnRestart: restart aborted -- widening "
                      "outputSelection must be allowed")
    else:
        if "cold start" not in log2 and "cold-start" not in log2:
            errors.append("outputChangeOnRestart: never-written bins should "
                          "warn about cold-starting")
        got = rv.read_bbpa_output(c2, FINAL_CS)
        # previously-written bins keep their full history
        for j in (2, 5):
            if j not in got:
                errors.append(f"outputChangeOnRestart: bin {j} missing")
            elif _err(got[j]["PA"], closed[j]["PA"]) >= TOL:
                errors.append(f"outputChangeOnRestart: bin {j} lost history "
                              f"({_err(got[j]['PA'], closed[j]['PA']):.3e})")
        # a never-written bin cold-started: cycles 3..5 only
        closed35 = bo.closed_form_bin_constant(
            rv.A, lambda n: rv.c_of_n(n + 3), 3)
        j = 0
        if j in got and _err(got[j]["PA"], closed35[j]["PA"]) >= TOL:
            errors.append(f"outputChangeOnRestart: cold-started bin {j} "
                          f"should hold the post-restart ensemble "
                          f"({_err(got[j]['PA'], closed35[j]['PA']):.3e})")
    return errors


def test_selection():
    errors = run_rows()
    assert not errors, "\n".join(errors)


if __name__ == "__main__":
    errors = run_rows()
    if errors:
        print("FAIL:")
        for e in errors:
            print("  ", e)
        sys.exit(1)
    print("selection rows: PASS -- calculate-all/write-some manifest exact,")
    print("  grouped phaseGrid/accumulationSelection/outputSelection spelling")
    print("  clean, and widening outputSelection across a restart preserves")
    print("  written history while never-written bins cold-start (warned)")
    sys.exit(0)
