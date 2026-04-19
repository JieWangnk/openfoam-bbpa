#!/usr/bin/env python3
"""
Parse matched-pair OpenFOAM logs (with/without BBPA) and report per-step
wall-time, per-step overhead, memory, and output size.

Usage:
    python3 bench_overhead.py \
        --nobbpa /path/to/bench_nobbpa_<jobid>.out \
        --bbpa   /path/to/bench_bbpa_<jobid>.out \
        [--timefile-nobbpa bench_nobbpa_time.txt] \
        [--timefile-bbpa   bench_bbpa_time.txt] \
        [--case-dir /scratch/.../bbpa_v2_test]

Reports:
    * n timesteps completed in each run
    * Total ExecutionTime and per-step mean / median / std
    * Overhead %  (bbpa / nobbpa - 1)
    * Output size delta
    * LaTeX-ready one-liner for Table II

Intended for the csf4 bench jobs that exercise a matched pair on the same
aorta mesh: bench_nobbpa_<id>.out vs bench_bbpa_<id>.out.
"""

from __future__ import annotations

import argparse
import os
import re
import statistics
import sys
from pathlib import Path


EXEC_RE = re.compile(r"^ExecutionTime\s*=\s*([\d.eE+-]+)\s*s\s*ClockTime\s*=\s*([\d.eE+-]+)\s*s")
TIME_RE = re.compile(r"^Time\s*=\s*([\d.eE+-]+)s")
ELAPSED_RE = re.compile(r"^ELAPSED_s=([\d.eE+-]+)")


def parse_foam_log(path: Path) -> dict:
    """Return dict with per-step Execution/Clock time lists and sim-time progression."""
    exec_times: list[float] = []
    clock_times: list[float] = []
    sim_times: list[float] = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = EXEC_RE.match(line)
            if m:
                exec_times.append(float(m.group(1)))
                clock_times.append(float(m.group(2)))
                continue
            m = TIME_RE.match(line)
            if m:
                sim_times.append(float(m.group(1)))
    return {
        "exec_cum_s": exec_times,
        "clock_cum_s": clock_times,
        "sim_time_s": sim_times,
    }


def parse_elapsed(path: Path) -> float | None:
    if not path.exists():
        return None
    with open(path) as f:
        for line in f:
            m = ELAPSED_RE.match(line.strip())
            if m:
                return float(m.group(1))
    return None


def per_step_stats(exec_cum: list[float]) -> dict:
    if len(exec_cum) < 2:
        return {"n_steps": 0, "total_s": 0.0, "mean_s": 0.0, "median_s": 0.0, "std_s": 0.0}
    # Drop the first entry (covers startup/IO), keep step-to-step deltas.
    deltas = [b - a for a, b in zip(exec_cum[:-1], exec_cum[1:])]
    return {
        "n_steps": len(deltas),
        "total_s": exec_cum[-1] - exec_cum[0],
        "mean_s": statistics.mean(deltas),
        "median_s": statistics.median(deltas),
        "std_s": statistics.pstdev(deltas) if len(deltas) > 1 else 0.0,
    }


def dir_size_kb(path: Path) -> int:
    if not path.exists():
        return 0
    total = 0
    for root, _dirs, files in os.walk(path):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return total // 1024


def fmt_row(label: str, nobbpa: float | int | str, bbpa: float | int | str, unit: str = "") -> str:
    return f"  {label:<30} {str(nobbpa):>15} {str(bbpa):>15} {unit}"


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--nobbpa", type=Path, required=True, help="OpenFOAM log from no-BBPA run")
    p.add_argument("--bbpa", type=Path, required=True, help="OpenFOAM log from BBPA run")
    p.add_argument("--timefile-nobbpa", type=Path, default=None, help="bench_nobbpa_time.txt")
    p.add_argument("--timefile-bbpa", type=Path, default=None, help="bench_bbpa_time.txt")
    p.add_argument("--case-dir", type=Path, default=None, help="case root for output-size diff")
    p.add_argument("--I", type=int, default=None, help="BBPA nBins used (for table label)")
    args = p.parse_args()

    for f in (args.nobbpa, args.bbpa):
        if not f.exists():
            print(f"ERROR: log missing: {f}", file=sys.stderr)
            return 2

    no = parse_foam_log(args.nobbpa)
    bb = parse_foam_log(args.bbpa)
    no_stats = per_step_stats(no["exec_cum_s"])
    bb_stats = per_step_stats(bb["exec_cum_s"])

    no_wall = parse_elapsed(args.timefile_nobbpa) if args.timefile_nobbpa else None
    bb_wall = parse_elapsed(args.timefile_bbpa) if args.timefile_bbpa else None

    if no_stats["mean_s"] > 0:
        overhead_per_step = (bb_stats["mean_s"] - no_stats["mean_s"]) / no_stats["mean_s"] * 100
    else:
        overhead_per_step = float("nan")

    print("=" * 72)
    print(" BBPA wall-time overhead benchmark")
    print("=" * 72)
    print(f"  Log (no BBPA): {args.nobbpa}")
    print(f"  Log (BBPA):    {args.bbpa}")
    print()
    print(fmt_row("Metric", "No BBPA", "BBPA", ""))
    print("  " + "-" * 62)
    print(fmt_row("Timesteps completed", no_stats["n_steps"], bb_stats["n_steps"]))
    print(fmt_row("Total ExecutionTime",
                  f"{no_stats['total_s']:.1f}",
                  f"{bb_stats['total_s']:.1f}", "s"))
    print(fmt_row("Per-step mean",
                  f"{no_stats['mean_s']:.4f}",
                  f"{bb_stats['mean_s']:.4f}", "s"))
    print(fmt_row("Per-step median",
                  f"{no_stats['median_s']:.4f}",
                  f"{bb_stats['median_s']:.4f}", "s"))
    print(fmt_row("Per-step stdev",
                  f"{no_stats['std_s']:.4f}",
                  f"{bb_stats['std_s']:.4f}", "s"))
    if no_wall is not None and bb_wall is not None:
        print(fmt_row("SLURM wall time",
                      f"{no_wall:.1f}",
                      f"{bb_wall:.1f}", "s"))
    if no["sim_time_s"] and bb["sim_time_s"]:
        print(fmt_row("Sim time reached",
                      f"{no['sim_time_s'][-1]:.4f}",
                      f"{bb['sim_time_s'][-1]:.4f}", "s (sim)"))
    if args.case_dir:
        out_no = 0
        out_bb = 0
        for d in sorted(args.case_dir.glob("[0-9]*")):
            if d.is_dir():
                out_bb += dir_size_kb(d)
        print(fmt_row("Case output dir", "(baseline)", f"{out_bb/1024:.1f}", "MB"))
    print("  " + "-" * 62)
    print(fmt_row("Overhead (mean step)",
                  "—",
                  f"+{overhead_per_step:.2f}", "%"))
    print()
    label = f"I={args.I}" if args.I else "I=?"
    print("LaTeX row for Table II:")
    print(f"  Solver wall-time overhead        & —  & "
          f"\\bbpaOverhead{{placeholder}} ({label}) & ... \\\\")
    print()
    print("Suggested macro fill-in (update \\bbpaOverheadHundred or "
          "\\bbpaOverheadFiveHundred in interacttfqsample.tex):")
    print(f"  \\newcommand{{\\bbpaOverhead{label.replace('=','')}}}{{{overhead_per_step:+.1f}\\%}}")
    print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
