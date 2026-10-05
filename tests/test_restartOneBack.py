#!/usr/bin/env python3
"""
test_restartOneBack.py -- zero-loss restart one write back (v2.2 fix).

Sparse centred bins whose labels are the solver write times (writeControl
timeStep), per-field writePerCycle, layout phaseAligned -- the aorta
production set-up. A restart from the LATEST write L straddles bin L, whose
_Bin is flushed only at the next write, so that bin is dropped for one cycle
(unbiased). A restart from the write BEFORE, L-dW, straddles bin L-dW whose
_Bin IS on disk: it must be recovered whole and the run must equal the
uninterrupted one bit for bit. v2.1 restored that bin and then dropped it
anyway (the recoverable check was gated on the legacy writeMode bin).

Runnable standalone (exit 1 on failure, for Alltest) and via pytest.
"""
from __future__ import annotations

import os, sys, shutil, glob, re
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bbpa_verify_lib as rv
DT=1/1024; END=4.0; WI=128
SP=" ".join(str(j) for j in range(0,64,8))
BODY=f"""        fields          ( U {{ writePerCycle on; }} );
        phaseGrid {{ period 1.0; nBins 64; phaseOrigin 0; binAlignment centred; }}
        discardCycles   1;
        maxCycles       2;
        accumulationSelection {{ mode list; bins ( {SP} ); }}
        outputSelection {{ mode allAccumulated; }}
        output {{ layout phaseAligned; meanWriteIntervalCycles 1; writeFinal true; writeProperties true; }}
        executeControl  timeStep; executeInterval 1; writeControl writeTime;
"""
def cd(end,start,logp,sf="startTime"):
    return rv.control_dict("smooth",DT,end,WI,sf,logp,start_time=start,bbpa_body=BODY,write_control="timeStep",time_precision=8)
def go(r,tag):
    x=rv.sh(f"foamRun > log.{tag} 2>&1", r); assert x.returncode==0, open(os.path.join(r,"log."+tag)).read()[-2000:]
def tdirs(r):
    return {round(float(os.path.basename(d)),9):d for d in glob.glob(os.path.join(r,"*")) if re.fullmatch(r"[0-9.]+",os.path.basename(d))}
def rd(p): return np.array(rv.bo.read_of_internal_first(p),float)
def compare(a,b):
    da,db=tdirs(a),tdirs(b); diff=[]; n=0
    for t,d in da.items():
        if t not in db: continue
        for f in os.listdir(d):
            if not f.startswith("U_"): continue
            pb=os.path.join(db[t],f)
            if not os.path.exists(pb): diff.append(f"{t}/{f} missing"); continue
            e=np.max(np.abs(rd(os.path.join(d,f))-rd(pb))); n+=1
            if e!=0: diff.append(f"{t}/{f} {e:.2e}")
    return n,diff


def run_rows():
    errors = []
    os.makedirs(rv.RUNROOT, exist_ok=True)
    S, lp = rv.make_case("ob_cont", "smooth", DT, END, WI, bbpa_body=BODY,
                         write_control="timeStep", time_precision=8)
    assert rv.block_mesh(S).returncode == 0; go(S, "run")
    res = {}
    for name, back in (("ob_latest", False), ("ob_oneback", True)):
        r, lp = rv.make_case(name, "smooth", DT, 2.25, WI, bbpa_body=BODY,
                             write_control="timeStep", time_precision=8)
        assert rv.block_mesh(r).returncode == 0; go(r, "leg1")
        start = 2.25
        if back:
            shutil.rmtree(tdirs(r)[2.25]); start = 2.125
        open(os.path.join(r, "system", "controlDict"), "w").write(cd(END, start, lp))
        go(r, "leg2")
        res[name] = (compare(S, r), open(os.path.join(r, "log.leg2")).read())
    (n1, d1), log1 = res["ob_oneback"]
    (n0, d0), log0 = res["ob_latest"]
    if d1:
        errors.append(f"one back: {len(d1)} of {n1} files differ from the continuous run: {d1[:4]}")
    if "recovered whole" not in log1:
        errors.append("one back: straddled bin not reported as recovered whole")
    if not d0 or "Dropping bin 16" not in log0:
        errors.append("latest: expected bin 16 dropped for one cycle (control row)")
    print(f"one back: {n1} files compared, {len(d1)} differ; latest: bin 16 dropped, {len(d0)} files differ (expected)")
    return errors


def test_restartOneBack():
    errors = run_rows()
    assert not errors, "\n".join(errors)


if __name__ == "__main__":
    errors = run_rows()
    if errors:
        print("FAIL:"); [print("  ", e) for e in errors]; sys.exit(1)
    print("restartOneBack: PASS -- a restart one write back recovers the straddled bin whole (zero loss)")
    sys.exit(0)
