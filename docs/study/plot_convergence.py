#!/usr/bin/env python3
"""plot_convergence.py -- the chapter-07 figure: what the ensemble keeps.

Runs tutorials/bbpaVerify twice with endTime 10 (ten cycles of the same
flow): once with discardCycles 0 and once with discardCycles 1. At the
liveliest cell and the phase with the largest startup imprint, it then
plots -- on a log scale, against the cycle index -- the distance from
the settled flow (cycle 9's own U_Bin) of:

  U_Bin(c)              each cycle's own mean: decays as the startup
                        transient dies out
  U_PA with discard 0   the ensemble that folded the startup cycle:
                        its error only fades as 1/N
  U_PA with discard 1   the ensemble that skipped it

Needs OpenFOAM sourced and libPhaseAveraging built; the two runs
(~1 min each, reused if already on disk) land in tests/_runs/.

Usage:  python3 plot_convergence.py [--out figures/pa-convergence.png]
                                    [--force]   # rerun the cases
"""
import argparse
import os
import re
import shutil
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
TUTORIAL = os.path.join(REPO, "tutorials", "bbpaVerify")
sys.path.insert(0, os.path.join(REPO, "tests"))
sys.path.insert(0, TUTORIAL)
import bbpa_verify_lib as rv                                    # noqa: E402
from check_bbpa import read_internal, time_dirs                 # noqa: E402

BLUE = "#2a78d6"      # U_Bin
ORANGE = "#eb6834"    # U_PA, discardCycles 0
AQUA = "#1baf7a"      # U_PA, discardCycles 1
INK = "#0b0b0b"
INK2 = "#52514e"
SURFACE = "#fcfcfb"

PERIOD, NBINS, END = 1.0, 10, 10
BIN_DT = PERIOD / NBINS


def make_case(name, discard):
    """Copy the tutorial (inputs only), set endTime 10 and discardCycles."""
    rundir = os.path.join(rv.RUNROOT, name)
    if os.path.exists(rundir):
        shutil.rmtree(rundir)
    os.makedirs(os.path.join(rundir, "0"))
    for f in ("U", "p"):
        shutil.copy(os.path.join(TUTORIAL, "0", f),
                    os.path.join(rundir, "0", f))
    shutil.copytree(os.path.join(TUTORIAL, "constant"),
                    os.path.join(rundir, "constant"),
                    ignore=shutil.ignore_patterns("polyMesh"))
    shutil.copytree(os.path.join(TUTORIAL, "system"),
                    os.path.join(rundir, "system"))
    cd = os.path.join(rundir, "system", "controlDict")
    txt = open(cd).read()
    txt, n1 = re.subn(r"(?m)^endTime\s+3;", f"endTime         {END};", txt)
    txt, n2 = re.subn(r"discardCycles\s+1;",
                      f"discardCycles   {discard};", txt)
    assert n1 == 1 and n2 == 1, f"{name}: controlDict patch failed"
    open(cd, "w").write(txt)
    return rundir


def run_case(name, discard, force):
    rundir = os.path.join(rv.RUNROOT, name)
    done = os.path.join(rundir, "%g" % (END - 0.1), "U_PA")
    if not force and os.path.isfile(done):
        print(f"{name}: reusing the existing run")
        return rundir
    rundir = make_case(name, discard)
    print(f"{name}: running (endTime {END}, discardCycles {discard}) ...")
    assert rv.block_mesh(rundir).returncode == 0, f"{name}: blockMesh failed"
    assert rv.run_serial(rundir).returncode == 0, f"{name}: foamRun failed"
    return rundir


def field_at(case, tds, t, fname, cell):
    for tv, d in tds:
        if abs(tv - t) < 1e-7:
            p = os.path.join(case, d, fname)
            if os.path.isfile(p):
                v = read_internal(p)
                return v[0, 0] if v.shape[0] == 1 else v[cell, 0]
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "figures",
                                                  "pa-convergence.png"))
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    if shutil.which("foamRun") is None:
        sys.exit("foamRun not in PATH -- source /opt/openfoam14/etc/bashrc")
    os.makedirs(rv.RUNROOT, exist_ok=True)
    d0 = run_case("study_fig_conv_d0", 0, a.force)
    d1 = run_case("study_fig_conv_d1", 1, a.force)
    tds0, tds1 = time_dirs(d0), time_dirs(d1)

    # Probe cell: where Ux moves the most over the run.
    inst = [read_internal(os.path.join(d0, d, "U"))
            for t, d in tds0 if os.path.isfile(os.path.join(d0, d, "U"))
            and os.path.isfile(os.path.join(d0, d, "p"))]
    ncell = max(U.shape[0] for U in inst)
    ux = np.array([np.repeat(U[:, 0], ncell) if U.shape[0] == 1
                   else U[:, 0] for U in inst])
    cell = int(np.argmax(np.std(ux, axis=0)))

    # Phase with the largest startup imprint; reference = cycle 9's U_Bin.
    # (Cycle 0's bin 0 is absent: its centred window straddles the run
    # start, and a partial window is never written.)
    imprint = []
    for jc in range(NBINS):
        v0 = field_at(d0, tds0, 0 + jc*BIN_DT, "U_Bin", cell)
        v9 = field_at(d0, tds0, 9 + jc*BIN_DT, "U_Bin", cell)
        if v0 is not None and v9 is not None:
            imprint.append((abs(v0 - v9), jc))
    j = max(imprint)[1]
    ref = field_at(d0, tds0, 9 + j*BIN_DT, "U_Bin", cell)

    def curve(case, tds, fname, cycles):
        # keep positive errors only: a zero cannot sit on a log axis
        out = []
        for c in cycles:
            v = field_at(case, tds, c + j*BIN_DT, fname, cell)
            if v is not None and abs(v - ref) > 0:
                out.append((c, abs(v - ref)))
        return [c for c, _ in out], [e for _, e in out]

    cb, eb = curve(d0, tds0, "U_Bin", range(0, 9))
    c0, e0 = curve(d0, tds0, "U_PA", range(0, 9))
    c1, e1 = curve(d1, tds1, "U_PA", range(1, 9))

    fig, ax = plt.subplots(figsize=(9.2, 4.6), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.semilogy(cb, eb, "o-", color=BLUE, lw=1.8, ms=5,
                label="U_Bin(c): each cycle's own mean")
    ax.semilogy(c0, e0, "D-", color=ORANGE, lw=1.8, ms=5,
                label="U_PA, discardCycles 0")
    ax.semilogy(c1, e1, "s-", color=AQUA, lw=1.8, ms=5,
                label="U_PA, discardCycles 1")

    # A visibly offset 1/N guide (the discard-0 curve itself IS 1/N).
    n = np.arange(1, 10)
    ax.semilogy(n - 1, 0.3*e0[0]/n, "--", color=INK2, lw=1.0, alpha=0.6)
    ax.text(7.6, 0.3*e0[0]/9*0.5, "∝ 1/N", fontsize=8.5, color=INK2)

    ax.annotate("the startup cycle stays in this ensemble;\n"
                "its error only fades as 1/N",
                xy=(5, e0[5]), xytext=(3.6, 1.5e-6), ha="center",
                fontsize=8.5, color=INK,
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8,
                                shrinkB=6))
    ax.annotate("one cycle after startup the flow already repeats to the\n"
                "file's write precision; the discardCycles 1 ensemble\n"
                "sits there from its first fold",
                xy=(2, e1[1]), xytext=(2.6, 1e-9), ha="center", va="top",
                fontsize=8.5, color=INK,
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8,
                                shrinkB=6))

    ax.set_xlabel("cycle index c (values read at that cycle's labels)",
                  color=INK)
    ax.set_ylabel(f"|Ux − settled Ux|  at cell {cell}, bin {j}",
                  color=INK)
    ax.set_title("tutorials/bbpaVerify over ten cycles: the ensemble mean "
                 "keeps every cycle it folds", color=INK, fontsize=11,
                 pad=34)
    ax.grid(axis="y", color=INK2, alpha=0.15, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8.5)
    leg = ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0),
                    ncol=3, fontsize=8.5, frameon=False,
                    columnspacing=1.4, handletextpad=0.6)
    for t in leg.get_texts():
        t.set_color(INK)

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.tight_layout()
    fig.savefig(a.out, facecolor=SURFACE, bbox_inches="tight")
    print(f"wrote {a.out}  (probe cell {cell}, bin {j})")


if __name__ == "__main__":
    main()
