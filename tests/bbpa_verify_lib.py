#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_verification.py
===================
Deterministic SOFTWARE verification of the BBPA OpenFOAM function object against
a prescribed analytic field (reviewer point 4, "strongest repair").

For each configuration this script
  1. generates a tiny OpenFOAM case in which a coded functionObject OVERWRITES
     U every solver step with a known analytic value U(t) (and logs the exact
     (t, dt, U) it set), BEFORE BBPA accumulates it;
  2. runs foamRun (serial or parallel) with the freshly-built libPhaseAveraging;
  3. parses BBPA's on-disk phase-aligned _PA / _PA_UU / _PAVariance fields;
  4. compares them against the bit-faithful offline replica (bbpa_offline.py)
     driven by the SAME logged samples  ->  isolates IMPLEMENTATION error; and
     against the continuous closed form  ->  the time-discretisation tolerance.

Outputs committed result tables under results/.

Field families
  bin-constant : U(t) = A[bin(t)] + c(cycle(t))     (Riemann-exact + straddle-free
                 when dt | binDeltaT  ->  machine precision vs closed form)
  smooth       : U(t) = amp*(cos, sin, camp*cos2) (phase) + c(cycle)
                 (the continuous bin mean is a non-trivial integral  ->  exhibits
                 the documented first-order time-discretisation error)
"""

from __future__ import annotations
import os
import sys
import math
import shutil
import subprocess
import glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bbpa_offline as bo

TEMPLATE = os.path.join(HERE, "case_template")
RUNROOT = os.path.join(HERE, "_runs")          # gitignored
RESULTS = os.path.join(HERE, "results")        # committed

# --------------------------------------------------------------------------- #
#  Shared analytic-field parameters (single source of truth: baked into the
#  coded C++ FO AND used by the Python replica/closed-form).
# --------------------------------------------------------------------------- #
T = 1.0                  # period
I = 8                    # nBins
BIN_DT = T / I           # = 0.125 = 2^-3
AMP = 1.0
CAMP = 0.3
# per-cycle offset c_n = CBASE + CSTEP*n  (gives nonzero cross-cycle variance)
CBASE = np.array([0.0, 0.0, 0.0])
CSTEP = np.array([0.10, -0.05, 0.03])


def Avec(j):
    th = 2.0 * math.pi * (j + 0.5) / I
    return np.array([math.cos(th), math.sin(th), 0.3 * math.sin(2.0 * th)])


A = [Avec(j) for j in range(I)]


def c_of_n(n):
    return CBASE + CSTEP * n


# --------------------------------------------------------------------------- #
#  Coded-FO C++ body generators (params baked as full-precision literals).
# --------------------------------------------------------------------------- #
def _f(x):
    return repr(float(x))            # 17-sig-fig round-trippable literal


def coded_fo(kind, logpath, cstep=None):
    """Return the `prescribe` coded functionObject dict text."""
    if cstep is None:
        cstep = CSTEP
    Ax = ", ".join(_f(A[j][0]) for j in range(I))
    Ay = ", ".join(_f(A[j][1]) for j in range(I))
    Az = ", ".join(_f(A[j][2]) for j in range(I))

    if kind == "bin_constant":
        ucalc = f"""
        Foam::vector Uval(Ax[binI]+cx, Ay[binI]+cy, Az[binI]+cz);
        """
    elif kind == "smooth":
        ucalc = f"""
        const Foam::scalar w = 2.0*Foam::constant::mathematical::pi/T;
        Foam::vector base
        (
            {_f(AMP)}*Foam::cos(w*phase),
            {_f(AMP)}*Foam::sin(w*phase),
            {_f(AMP)}*{_f(CAMP)}*Foam::cos(2.0*w*phase)
        );
        Foam::vector Uval(base.x()+cx, base.y()+cy, base.z()+cz);
        """
    else:
        raise ValueError(kind)

    return f"""
    prescribe
    {{
        type            coded;
        libs            ("libutilityFunctionObjects.so");
        name            prescribeAnalytic;
        executeControl  timeStep;
        executeInterval 1;
        writeControl    timeStep;
        writeInterval   1000000000;

        codeInclude
        #{{
            #include "volFields.H"
            #include "Pstream.H"
            #include "mathematicalConstants.H"
            #include <fstream>
            #include <cmath>
        #}};

        codeExecute
        #{{
            const Foam::Time& rt = mesh().time();
            const Foam::scalar t  = rt.value();
            const Foam::scalar dt = rt.deltaT().value();
            const Foam::scalar T  = {_f(T)};
            const Foam::label  I  = {I};
            const Foam::scalar t0 = 0.0;
            const Foam::scalar dtbin = T/Foam::scalar(I);
            const Foam::scalar tt = t - t0;
            const Foam::scalar phase = std::fmod(tt, T);
            Foam::label binI = Foam::min(Foam::label(phase/dtbin), I-1);
            Foam::label cyc  = Foam::label(tt/T);

            static const Foam::scalar Ax[{I}] = {{ {Ax} }};
            static const Foam::scalar Ay[{I}] = {{ {Ay} }};
            static const Foam::scalar Az[{I}] = {{ {Az} }};
            const Foam::scalar cx = {_f(CBASE[0])} + ({_f(cstep[0])})*cyc;
            const Foam::scalar cy = {_f(CBASE[1])} + ({_f(cstep[1])})*cyc;
            const Foam::scalar cz = {_f(CBASE[2])} + ({_f(cstep[2])})*cyc;
            {ucalc}
            volVectorField& U = mesh().lookupObjectRef<volVectorField>("U");
            forAll(U.primitiveFieldRef(), ci)
            {{
                U.primitiveFieldRef()[ci] = Uval;
            }}
            U.correctBoundaryConditions();

            if (Foam::Pstream::master())
            {{
                std::ofstream ofs("{logpath}", std::ios::app);
                ofs.precision(17);
                ofs << t << " " << dt << " "
                    << Uval.x() << " " << Uval.y() << " " << Uval.z() << " "
                    << binI << " " << cyc << "\\n";
            }}
        #}};
    }}
"""


def bbpa_fo(extra=""):
    """The BBPA FO dict.  `extra` appends raw dict lines (e.g.
    'writeMode bin;\\n        PAwriteInterval 2;') for option-matrix tests."""
    return f"""
    bbpa
    {{
        type            BBPA;
        libs            ("libPhaseAveraging.so");
        fields          (U);
        period          {_f(T)};
        nBins           {I};
        startCycle      0;
        startTime       0;
        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;
        {extra}
    }}
"""


def control_dict(kind, dt, end_time, write_interval, start_from, logpath,
                 start_time=0.0, cstep=None, bbpa_extra=""):
    return f"""/*--------------------------------*- C++ -*----------------------------------*/
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

startFrom       {start_from};
startTime       {_f(start_time)};
stopAt          endTime;
endTime         {_f(end_time)};
deltaT          {_f(dt)};

writeControl    runTime;
writeInterval   {_f(write_interval)};
purgeWrite      0;
writeFormat     ascii;
writePrecision  17;
writeCompression off;
timeFormat      general;
timePrecision   10;
runTimeModifiable false;
adjustTimeStep  no;

functions
{{
{coded_fo(kind, logpath, cstep)}
{bbpa_fo(bbpa_extra)}
}}

// ************************************************************************* //
"""


# --------------------------------------------------------------------------- #
#  Case generation / execution helpers.
# --------------------------------------------------------------------------- #
def sh(cmd, cwd):
    """Run a shell command (already OpenFOAM-sourced via parent env)."""
    r = subprocess.run(cmd, cwd=cwd, shell=True, capture_output=True, text=True)
    return r


def make_case(name, kind, dt, end_time, write_interval, start_from="startTime",
              cstep=None, bbpa_extra=""):
    rundir = os.path.join(RUNROOT, name)
    if os.path.exists(rundir):
        shutil.rmtree(rundir)
    shutil.copytree(TEMPLATE, rundir)
    logpath = os.path.join(rundir, "prescribe.log")
    if os.path.exists(logpath):
        os.remove(logpath)
    with open(os.path.join(rundir, "system", "controlDict"), "w") as f:
        f.write(control_dict(kind, dt, end_time, write_interval, start_from,
                             logpath, cstep=cstep, bbpa_extra=bbpa_extra))
    return rundir, logpath


def block_mesh(rundir):
    return sh("blockMesh > log.blockMesh 2>&1", rundir)


def run_serial(rundir):
    return sh("foamRun > log.foamRun 2>&1", rundir)


def run_parallel(rundir, np_):
    sh(f"decomposePar -force > log.decomposePar 2>&1", rundir)
    return sh(f"mpirun -np {np_} foamRun -parallel > log.foamRun 2>&1", rundir)


# --------------------------------------------------------------------------- #
#  Parse logged samples and BBPA on-disk output.
# --------------------------------------------------------------------------- #
def read_log(logpath):
    samples = []
    with open(logpath) as f:
        for line in f:
            p = line.split()
            if len(p) < 5:
                continue
            t, dt = float(p[0]), float(p[1])
            U = np.array([float(p[2]), float(p[3]), float(p[4])])
            samples.append((t, dt, U))
    return samples


def split_restart_log(samples, t_ckpt):
    """Split a restart prescribe.log at the checkpoint time.
    The coded FO logs t_ckpt TWICE -- once as the last step of segment 1 and
    once as segment 2's startup re-execute after restore. seg1 ends at the FIRST
    occurrence (= what the checkpoint captured); seg2 begins at the SECOND
    occurrence (= exactly the steps the restarted solver re-ran)."""
    idx = next(i for i, s in enumerate(samples) if abs(s[0] - t_ckpt) < 1e-9)
    return samples[:idx + 1], samples[idx + 1:]


def _time_dirs(base):
    out = {}
    for d in os.listdir(base):
        full = os.path.join(base, d)
        if not os.path.isdir(full):
            continue
        try:
            out[float(d)] = full
        except ValueError:
            pass
    return out


def latest_cycle_start(rundir, parallel=False, fname="U_PA"):
    """Largest k*period at which BBPA wrote a FULL set of I phase-aligned bins
    (a `fname` file in every bin dir; use fname="U_BinPA" for writeMode bin).
    Needed because a non-power-of-two dt (e.g. 0.05) accumulates float drift, so
    the final write can land at cycleStart < endTime (e.g. 5.0 rather than 6.0)."""
    base = os.path.join(rundir, "processor0") if parallel else rundir
    tdirs = _time_dirs(base)
    if not tdirs:
        return None
    kmax = int(max(tdirs) / T) + 1
    for k in range(kmax, -1, -1):
        cs = k * T
        ok = True
        for j in range(I):
            bt = cs + j * BIN_DT
            if not any(abs(tv - bt) < 1e-7 and
                       os.path.isfile(os.path.join(full, fname))
                       for tv, full in tdirs.items()):
                ok = False
                break
        if ok:
            return cs
    return None


def read_bbpa_output(rundir, cycle_start, parallel=False, suffix="_PA"):
    """
    Read BBPA's phase-aligned output at the write whose cycle starts at
    `cycle_start` (bins at cycle_start + j*BIN_DT).  Returns
    bin -> {PA: vec3, PA_UU: 3x3, PAVariance: scalar or None}.
    `suffix` selects the naming family: "_PA" (default) or "_BinPA"
    (writeMode bin).
    """
    base = os.path.join(rundir, "processor0") if parallel else rundir
    tdirs = _time_dirs(base)
    out = {}
    for j in range(I):
        bt = cycle_start + j * BIN_DT
        # match a time dir to bt
        match = None
        for tv, full in tdirs.items():
            if abs(tv - bt) < 1e-7:
                match = full
                break
        if match is None:
            continue
        rec = _read_bin_record(match, suffix)
        if rec:
            out[j] = rec
    return out


def _read_bin_record(tdir, suffix):
    """Read U<suffix>/U<suffix>_UU/U<suffix>Variance from one time dir."""
    rec = {}
    pPA = os.path.join(tdir, "U" + suffix)
    pUU = os.path.join(tdir, "U" + suffix + "_UU")
    pVar = os.path.join(tdir, "U" + suffix + "Variance")
    if os.path.isfile(pPA):
        rec["PA"] = bo.read_of_internal_first(pPA)
    if os.path.isfile(pUU):
        rec["PA_UU"] = bo.symm6_to_mat(bo.read_of_internal_first(pUU))
    if os.path.isfile(pVar):
        rec["PAVariance"] = float(bo.read_of_internal_first(pVar))
    else:
        rec["PAVariance"] = None
    return rec


def read_companion_output(rundir, cycle_start):
    """
    Read companion-mode U_PA from the SOLVER's own write dirs across the final
    cycle: dir at t = cycle_start + j*BIN_DT holds the current bin j's pure
    finalised cross-cycle mean (requires solver writeInterval == BIN_DT so each
    bin phase has a write).  Returns bin -> {PA: vec3}.
    """
    tdirs = _time_dirs(rundir)
    out = {}
    for j in range(I):
        bt = cycle_start + j * BIN_DT
        for tv, full in tdirs.items():
            if abs(tv - bt) < 1e-7:
                p = os.path.join(full, "U_PA")
                if os.path.isfile(p):
                    out[j] = {"PA": bo.read_of_internal_first(p)}
                break
    return out


def read_bin_files(rundir, cycle_start):
    """
    Read writeMode-bin per-cycle U_Bin (THIS cycle's bin mean, no cross-cycle
    averaging) from the phase-aligned dirs of the cycle starting at
    `cycle_start`.  Returns bin -> vec3.
    """
    tdirs = _time_dirs(rundir)
    out = {}
    for j in range(I):
        bt = cycle_start + j * BIN_DT
        for tv, full in tdirs.items():
            if abs(tv - bt) < 1e-7:
                p = os.path.join(full, "U_Bin")
                if os.path.isfile(p):
                    out[j] = bo.read_of_internal_first(p)
                break
    return out


# --------------------------------------------------------------------------- #
#  Comparison utilities.
# --------------------------------------------------------------------------- #
def rel(a, b, scale):
    """relative error |a-b| / max(scale, tiny)."""
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b)))) / max(scale, 1e-300)


def compare(bbpa, ref, bins, label):
    """Max relative error over `bins` for PA, PA_UU, PAVariance."""
    ePA = eUU = eVar = 0.0
    sPA = sUU = sVar = 0.0
    for j in bins:
        if j not in bbpa or j not in ref:
            continue
        sPA = max(sPA, np.max(np.abs(ref[j]["PA"])))
        sUU = max(sUU, np.max(np.abs(ref[j]["PA_UU"])))
        ePA = max(ePA, np.max(np.abs(bbpa[j]["PA"] - ref[j]["PA"])))
        eUU = max(eUU, np.max(np.abs(bbpa[j]["PA_UU"] - ref[j]["PA_UU"])))
        if ref[j]["PAVariance"] is not None and bbpa[j].get("PAVariance") is not None:
            sVar = max(sVar, abs(ref[j]["PAVariance"]))
            eVar = max(eVar, abs(bbpa[j]["PAVariance"] - ref[j]["PAVariance"]))
    return {
        "label": label, "nbins": len(bins),
        "PA_abs": ePA, "PA_rel": ePA / max(sPA, 1e-300),
        "UU_abs": eUU, "UU_rel": eUU / max(sUU, 1e-300),
        "Var_abs": eVar, "Var_rel": eVar / max(sVar, 1e-300),
    }


if __name__ == "__main__":
    os.makedirs(RUNROOT, exist_ok=True)
    os.makedirs(RESULTS, exist_ok=True)
    print("run_verification: use driver functions from run_all.py")
