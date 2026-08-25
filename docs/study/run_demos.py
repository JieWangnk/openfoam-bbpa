#!/usr/bin/env python3
"""run_demos.py -- runnable demos for the BBPA study (docs/study/).

Each demo is one small OpenFOAM case that changes ONE thing in the BBPA
dictionary and shows the effect on disk. The case machinery comes from
the test suite (tests/bbpa_verify_lib.py): an 8-cell box whose U is overwritten
every step by a coded function object with a bin-constant analytic signal
(period 1 s, 8 bins), so runs are fast (~25-35 s each, dominated by the
coded FO's one-off compile) and byte-reproducible.

A demo's result is a MANIFEST: the dictionary that ran, then for every
product and state file the number of time directories holding it and the
first/last such time. The manifests are committed under docs/study/expected/
and embedded verbatim in the study chapters between <!-- demo:NAME -->
markers.

Usage:
  ./run_demos.py                run every demo, print the manifests
  ./run_demos.py --list         list the demos and what each shows
  ./run_demos.py --only NAME    restrict to one demo (with or without the
                                study_ prefix)
  ./run_demos.py --update       run, then (re)write docs/study/expected/
  ./run_demos.py --check        run, then diff against docs/study/expected/
                                (exit 1 on any difference)
  ./run_demos.py --check-docs   no OpenFOAM needed: verify that every
                                chapter's embedded manifest matches
                                docs/study/expected/ byte for byte and that
                                demos and chapters map one to one

Run modes need OpenFOAM sourced and libPhaseAveraging built
(source /opt/openfoam14/etc/bashrc; wmake libso in src/phaseAveraging).
Cases land in tests/_runs/ (gitignored).
"""
from __future__ import annotations

import argparse
import difflib
import glob
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "tests"))
import bbpa_verify_lib as rv                                    # noqa: E402

EXPECTED = os.path.join(HERE, "expected")
END_TIME = 5.0
DT = rv.BIN_DT / 16          # 0.0078125: divides the bin width exactly

# Rows shown in every manifest even when absent (absence is the point of
# several demos). The scan below adds any further product/state name it
# finds, so nothing new can appear silently.
BASE_NAMES = ["U_PA", "U_PA_UU", "U_PAVarianceMag", "U_PAmeanState",
              "U_Bin", "U_Bin_UU"]


def _tail(extra=""):
    return (extra + """        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;
""")


def _grid(fields="fields          (U);", discard="discardCycles   0;"):
    return f"""        {fields}
        period          1.0;
        nBins           8;
        {discard}
        phaseOrigin     0;
"""


# --------------------------------------------------------------------------- #
#  The demos. body = the BBPA dictionary body (everything between type/libs
#  and the closing brace); chapter = the study chapter that embeds the
#  manifest; runs = list of (suffix, body, write_interval, expect) --
#  most demos have one run, d11 has two.
# --------------------------------------------------------------------------- #
DEMOS = {}


def demo(name, chapter, blurb, runs, names=()):
    DEMOS[name] = {"chapter": chapter, "blurb": blurb, "runs": runs,
                   "names": list(names)}


demo("study_d01_minimal", "02-what-gets-written.md",
     "the baseline: all default products of one field",
     [("", _grid() + _tail(), 1.0, "ok")])

demo("study_d05a_second_moment_off", "05-fields.md",
     "secondMoment off + crossCycleVariance off: mean only",
     [("", _grid("fields          ( U { secondMoment off; "
                 "crossCycleVariance off; } );") + _tail(), 1.0, "ok")])

demo("study_d05b_mean_unpublished", "05-fields.md",
     "writePhaseMean off: no U_PA file, U_PAmeanState appears",
     [("", _grid("fields          ( U { writePhaseMean off; } );")
       + _tail(), 1.0, "ok")])

demo("study_d06_phase_grid", "06-phase-grid.md",
     "grouped phaseGrid{} with phaseOrigin 0.125: labels shift",
     [("", """        fields          (U);
        phaseGrid
        {
            period          1.0;
            nBins           8;
            phaseOrigin     0.125;
            binAlignment    left;
            timeQuadrature  endPoint;
        }
        discardCycles   0;
""" + _tail(), 1.0, "ok")])

demo("study_d07a_discard_cycles", "07-cycle-window.md",
     "discardCycles 1: the startup cycle leaves no products",
     [("", _grid(discard="discardCycles   1;") + _tail(), 1.0, "ok")])

demo("study_d07b_max_cycles", "07-cycle-window.md",
     "maxCycles 2: the ensemble freezes, later dumps repeat it",
     [("", _grid() + "        maxCycles       2;\n" + _tail(), 1.0, "ok")])

demo("study_d08_accumulation_list", "08-accumulation-selection.md",
     "accumulationSelection list (2 5): only two bins exist",
     [("", _grid()
       + "        accumulationSelection { mode list; bins (2 5); }\n"
       + _tail(), 1.0, "ok")])

demo("study_d09_output_list", "09-output-selection.md",
     "outputSelection list (2 5): products for two bins, state for all",
     [("", _grid()
       + "        accumulationSelection { mode all; }\n"
       + "        outputSelection { mode list; bins (2 5); }\n"
       + _tail(), 1.0, "ok")])

demo("study_d10a_companion", "10-output.md",
     "layout companion: products beside the solver's own fields",
     [("", _grid() + "        output { layout companion; }\n" + _tail(),
       rv.BIN_DT, "ok")])

demo("study_d10b_cadence", "10-output.md",
     "cycleWriteInterval 2 + meanWriteIntervalCycles 0",
     [("", _grid() + "        output { cycleWriteInterval 2; "
       "meanWriteIntervalCycles 0; }\n" + _tail(), 1.0, "ok")])

demo("study_d10c_no_final", "10-output.md",
     "writeFinal off + writeProperties off + cycleWriteInterval 4",
     [("", _grid() + "        output { writeFinal off; "
       "writeProperties off; cycleWriteInterval 4; }\n" + _tail(),
       1.0, "ok")])

demo("study_d11_field_availability", "11-field-availability.md",
     "a field with no producer: lazy completes, strict aborts",
     [("_lazy", _grid("fields          (U wallShearStress);")
       + "        fieldAvailability lazy;\n" + _tail(), 1.0, "ok"),
      ("_strict", _grid("fields          (U wallShearStress);")
       + "        fieldAvailability strict;\n" + _tail(), 1.0, "fatal")],
     names=["wallShearStress_PA", "wallShearStress_PA_UU",
            "wallShearStress_PAVarianceMag"])


# --------------------------------------------------------------------------- #
#  Manifest generation.
# --------------------------------------------------------------------------- #
def _num(d):
    try:
        return float(d)
    except ValueError:
        return None


def scan_case(rundir):
    """All product/state files across the case's time dirs:
    name -> (count, first time dir, last time dir). Solver-owned files
    (U, p, phi) are excluded; everything else is reported."""
    tdirs = sorted(((t, d) for d in os.listdir(rundir)
                    if os.path.isdir(os.path.join(rundir, d))
                    and (t := _num(d)) is not None))
    found = {}
    nprops = 0
    props_path = None
    for t, d in tdirs:
        full = os.path.join(rundir, d)
        for f in sorted(os.listdir(full)):
            p = os.path.join(full, f)
            if f == "uniform":
                for pf in sorted(glob.glob(os.path.join(p, "*Properties"))):
                    nprops += 1
                    props_path = pf
                continue
            if not os.path.isfile(p) or f in ("U", "p", "phi"):
                continue
            c, first, _ = found.get(f, (0, d, d))
            found[f] = (c + 1, first, d)
    return found, nprops, props_path


def folded_cycles(props_path):
    """Per-field folded-cycle counts from the newest Properties file."""
    if props_path is None:
        return None
    txt = open(props_path).read()
    out = []
    for m in re.finditer(r"(\w+)\s*\{[^{}]*?nCycles\s+(\d+);", txt):
        out.append(f"{m.group(1)}={m.group(2)}")
    return " ".join(out) if out else None


def fatal_message(rundir):
    """The FatalError message from log.foamRun, whitespace-normalised."""
    log = os.path.join(rundir, "log.foamRun")
    if not os.path.isfile(log):
        return "(no log.foamRun)"
    txt = open(log).read()
    m = re.search(r"FOAM FATAL (?:IO )?ERROR.*?\n(.*?)\n\s*\n", txt,
                  re.DOTALL)
    if not m:
        return "(no FatalError found in log.foamRun)"
    return " ".join(m.group(1).split())


def run_one(name, spec):
    """Run every case of one demo; return its manifest text."""
    lines = [f"demo: {name}", f"chapter: {spec['chapter']}"]
    for suffix, body, write_interval, expect in spec["runs"]:
        case = name + suffix
        rundir, _ = rv.make_case(case, "bin_constant", DT, END_TIME,
                                 write_interval, bbpa_body=body)
        r = rv.block_mesh(rundir)
        if r.returncode != 0:
            raise RuntimeError(f"{case}: blockMesh failed:\n{r.stderr}")
        r = rv.run_serial(rundir)

        lines.append("")
        lines.append(f"case: {case}  (endTime {END_TIME:g}, solver "
                     f"writeInterval {write_interval:g})")
        lines.append("BBPA dictionary body:")
        for bl in body.rstrip("\n").split("\n"):
            bl = bl[8:] if bl.startswith(" " * 8) else bl.lstrip()
            lines.append(("    " + bl).rstrip())

        if expect == "fatal":
            if r.returncode == 0:
                raise RuntimeError(f"{case}: expected a fatal error, "
                                   "but foamRun completed")
            lines.append("run: FATAL (expected). The error:")
            lines.append("    " + fatal_message(rundir))
            continue
        if r.returncode != 0:
            raise RuntimeError(f"{case}: foamRun failed:\n"
                               + fatal_message(rundir))
        lines.append("run: completed")

        found, nprops, props_path = scan_case(rundir)
        names = sorted(set(BASE_NAMES) | set(spec.get("names", []))
                       | set(found))
        w = max(len(n) for n in names + ["uniform/*Properties"]) + 2
        lines.append("")
        lines.append(f"    {'file':<{w}}{'dirs':>5}  first    last")
        for n in names:
            if n in found:
                c, first, last = found[n]
                lines.append(f"    {n:<{w}}{c:>5}  {first:<8} {last}")
            else:
                lines.append(f"    {n:<{w}}{0:>5}  -        -")
        lines.append(f"    {'uniform/*Properties':<{w}}{nprops:>5}")
        fc = folded_cycles(props_path)
        lines.append("folded cycles at the newest checkpoint: "
                     + (fc if fc else "(no checkpoint on disk)"))
        ndep = open(os.path.join(rundir, "log.foamRun")).read() \
            .count("deprecated")
        lines.append(f"deprecation warnings in the log: {ndep}")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
#  CLI modes.
# --------------------------------------------------------------------------- #
def select(only):
    if only is None:
        return dict(DEMOS)
    for cand in (only, "study_" + only):
        if cand in DEMOS:
            return {cand: DEMOS[cand]}
    sys.exit(f"unknown demo '{only}'; --list shows the names")


def mode_list():
    for name, spec in DEMOS.items():
        print(f"{name:<32} {spec['blurb']}   [{spec['chapter']}]")


def run_all(demos):
    if shutil.which("foamRun") is None:
        sys.exit("foamRun not in PATH -- source /opt/openfoam14/etc/bashrc "
                 "(and build src/phaseAveraging with wmake libso)")
    os.makedirs(rv.RUNROOT, exist_ok=True)
    return {name: run_one(name, spec) for name, spec in demos.items()}


def mode_update(demos):
    manifests = run_all(demos)
    os.makedirs(EXPECTED, exist_ok=True)
    for name, text in manifests.items():
        with open(os.path.join(EXPECTED, name + ".txt"), "w") as f:
            f.write(text)
        print(f"wrote expected/{name}.txt")


def mode_check(demos, full):
    manifests = run_all(demos)
    bad = []
    for name, text in manifests.items():
        path = os.path.join(EXPECTED, name + ".txt")
        if not os.path.isfile(path):
            bad.append(f"{name}: expected/{name}.txt missing (--update?)")
            continue
        want = open(path).read()
        if text != want:
            diff = "".join(difflib.unified_diff(
                want.splitlines(True), text.splitlines(True),
                f"expected/{name}.txt", "rerun"))
            bad.append(f"{name}: manifest changed\n{diff}")
    if full:
        stale = {os.path.basename(p)[:-4]
                 for p in glob.glob(os.path.join(EXPECTED, "*.txt"))} \
            - set(DEMOS)
        for s in sorted(stale):
            bad.append(f"expected/{s}.txt matches no demo")
    if bad:
        print("\n".join(bad))
        print(f"\nCHECK FAILED ({len(bad)} problem(s))")
        return 1
    print(f"CHECK PASSED ({len(manifests)} demo(s) reproduce their "
          "committed manifests)")
    return 0


DOC_RE = re.compile(r"<!-- demo:(\S+) -->\n```text\n(.*?)```", re.DOTALL)


def mode_check_docs():
    bad = []
    seen = {}
    for path in sorted(glob.glob(os.path.join(HERE, "*.md"))):
        chapter = os.path.basename(path)
        for m in DOC_RE.finditer(open(path).read()):
            name, text = m.group(1), m.group(2)
            if name in seen:
                bad.append(f"{chapter}: demo {name} embedded more than once "
                           f"(also in {seen[name]})")
                continue
            seen[name] = chapter
            if name not in DEMOS:
                bad.append(f"{chapter}: unknown demo marker {name}")
                continue
            if DEMOS[name]["chapter"] != chapter:
                bad.append(f"{chapter}: demo {name} belongs in "
                           f"{DEMOS[name]['chapter']}")
            path_exp = os.path.join(EXPECTED, name + ".txt")
            if not os.path.isfile(path_exp):
                bad.append(f"{chapter}: expected/{name}.txt missing")
            elif text != open(path_exp).read():
                bad.append(f"{chapter}: embedded manifest for {name} "
                           f"differs from expected/{name}.txt")
    for name in sorted(set(DEMOS) - set(seen)):
        bad.append(f"demo {name} is embedded in no chapter "
                   f"(should be in {DEMOS[name]['chapter']})")
    if bad:
        print("\n".join(bad))
        print(f"\nDOCS CHECK FAILED ({len(bad)} problem(s))")
        return 1
    print(f"DOCS CHECK PASSED ({len(seen)} embedded manifest(s) match "
          "expected/ byte for byte; demos and chapters map one to one)")
    return 0


def main():
    ap = argparse.ArgumentParser(
        description="runnable demos for the BBPA study")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--only", metavar="NAME")
    ap.add_argument("--update", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--check-docs", action="store_true")
    a = ap.parse_args()

    if a.list:
        mode_list()
        return 0
    if a.check_docs:
        return mode_check_docs()
    demos = select(a.only)
    if a.update:
        mode_update(demos)
        return 0
    if a.check:
        return mode_check(demos, full=(a.only is None))
    for name, text in run_all(demos).items():
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
