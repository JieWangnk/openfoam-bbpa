#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bbpa_offline.py
===============
Exact offline replica of the BBPA OpenFOAM function-object accumulator, plus
closed-form analytic references, for the deterministic software verification
(reviewer point 4, "strongest repair").

The replica reproduces -- operation for operation, in the same order -- the
arithmetic in
    src/functionObjects/field/phaseAveraging/{BBPA.C, BBPATemplates.C}
so that, fed the *exact* per-step (t, dt, U) samples that BBPA accumulated,
its outputs equal BBPA's on-disk _PA / _PA_UU / _PAVariance fields to machine
precision.  This isolates *implementation* error (replica == compiled object)
from the *discretisation* tolerance (replica == continuous closed form, only
when dt | binDeltaT and the field is bin-constant).

Authoritative arithmetic being replicated
------------------------------------------
phaseBin(t)        : phase = fmod(t - startTime, period);
                     binI  = min(int(phase/binDeltaT), nBins-1)
cycleNumber(t)     : int((t - startTime)/period)        [trunc toward zero]
accumulate(j,U,dt) : binTime[j]   += dt
                     binSum[j]    += dt * U
                     binSqrSum[j] += dt * (U outer U)         (symmTensor)
finalizeCycle()    : for each bin with binTime>0:
                       cycleMean    = binSum/binTime
                       cycleMeanSqr = binSqrSum/binTime
                       count        = ++binCounts[j]
                       if count == 1: meanFields=cycleMean; meanSqr=cycleMeanSqr; M2=0
                       else:
                         delta  = cycleMean - meanFields
                         meanFields = (meanFields*(count-1) + cycleMean)/count
                         delta2 = cycleMean - meanFields
                         M2 += delta & delta2          # true Welford inner product (scalar)
                         meanSqr = (meanSqr*(count-1) + cycleMeanSqr)/count
write() per bin (what lands on disk) -- see BBPATemplates.C write():
   phi_PA          = meanFields (PURE finalised; running partial ONLY when
                     no complete cycle exists yet -- never blended)
   phi_PAVariance  = M2/(binCounts-1)            (PURE finalised; never blended)
   phi_PA_UU       = meanSqr    (same rule as phi_PA)
Execute() (BBPA.C, v2 overlap contract): each solver step (t-dt, t] is
split at every bin boundary it crosses; each bin is credited its exact
overlap, with the step's single sample treated as constant over the whole
interval. The previous cycle is finalised as its last segment closes
(including a step ENDING exactly on the boundary); pre-checkpoint segments
of a restart's replayed step are skipped.

Vectors are length-3 numpy arrays; the symmetric second moment U_i U_j is a
full 3x3 symmetric matrix.  magSqr(vector) = sum of squares (a scalar), exactly
as Foam::magSqr(vector).
"""

from __future__ import annotations
import math
import numpy as np


# --------------------------------------------------------------------------- #
#  Bin / cycle index maps -- replicate BBPA.C's snapped floor EXACTLY.
#  The library (BBPA.C bbpaSnapFloor + phaseBin) floors with a snap-up when the
#  argument sits within 1e-9 (bin units) below an integer, and wraps an fmod
#  result of period-1ulp back to phase 0: an accumulated solver clock like
#  fmod(2.9, 1.0) = 0.8999999999999999 would otherwise truncate one bin low
#  (documented in BBPA.C:55-72).  Plain int() truncation here diverges from the
#  compiled object on any non-representable dt (e.g. 0.05).
# --------------------------------------------------------------------------- #
def snap_floor(x):
    """Foam bbpaSnapFloor: floor, snapping up within 1e-9 below an integer."""
    r = math.floor(x)
    return int(r + 1.0) if (x - r) > (1.0 - 1e-9) else int(r)


def phase_bin(t, start_time, period, n_bins):
    bin_dt = period / n_bins
    phase = math.fmod(t - start_time, period)
    # fmod can return period-1ulp instead of 0 for an on-grid time (BBPA.C:83-85)
    if phase > period * (1.0 - 1e-12):
        phase = 0.0
    return min(snap_floor(phase / bin_dt), n_bins - 1)


def cycle_number(t, start_time, period):
    return snap_floor((t - start_time) / period)


def magsqr_vec(v):
    """Foam::magSqr(vector) -> scalar sum of squared components."""
    return float(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def outer_sym(u):
    """U outer U as a 3x3 symmetric matrix (Foam::sqr(vector) -> symmTensor)."""
    return np.outer(u, u)


# --------------------------------------------------------------------------- #
#  Exact accumulator replica for a single (uniform-field) cell.
# --------------------------------------------------------------------------- #
class BBPAReplica:
    """
    Bit-faithful replica of binItem<vector> for one cell, driven by the same
    execute()/finalizeCycle()/write() control flow as BBPA.
    """

    def __init__(self, period, n_bins, start_time=0.0, start_cycle=0):
        self.period = float(period)
        self.n_bins = int(n_bins)
        self.bin_dt = self.period / self.n_bins
        self.start_time = float(start_time)
        self.start_cycle = int(start_cycle)

        I = self.n_bins
        # cross-cycle finalised state
        self.binCounts = np.zeros(I)                  # scalarField
        self.meanFields = [np.zeros(3) for _ in range(I)]
        self.meanSqr = [np.zeros((3, 3)) for _ in range(I)]
        self.M2 = np.zeros(I)                         # scalar per bin
        self.nCycles = 0
        # current-cycle accumulators
        self.binTime = np.zeros(I)
        self.binSum = [np.zeros(3) for _ in range(I)]
        self.binSqrSum = [np.zeros((3, 3)) for _ in range(I)]
        # execute() state
        self.previousCycle = -1
        self.currentCycle = 0
        self.currentBin = -1
        # per-finalised-cycle record of the bin weights W_{c,j} (a copy of
        # binTime at each finalizeCycle), for the weight-conservation checks:
        # sum_j W_{c,j} == T holds exactly only when steps never straddle a
        # cycle boundary (v1 assigns a whole step to the bin/cycle containing
        # its END time; the v2 split rule will restore exact conservation).
        self.cycleWeightHist = []

    # -- per-step accumulate (binItem::accumulate) ------------------------- #
    def _accumulate(self, binI, U, dt):
        self.binTime[binI] += dt
        self.binSum[binI] = self.binSum[binI] + dt * U
        self.binSqrSum[binI] = self.binSqrSum[binI] + dt * outer_sym(U)

    # -- reset current cycle (resetCurrentCycle) --------------------------- #
    def _reset_cycle(self):
        I = self.n_bins
        self.binTime = np.zeros(I)
        self.binSum = [np.zeros(3) for _ in range(I)]
        self.binSqrSum = [np.zeros((3, 3)) for _ in range(I)]

    # -- finalize previous cycle (binItem::finalizeCycle) ------------------ #
    def _finalize_cycle(self):
        self.cycleWeightHist.append(self.binTime.copy())
        self.nCycles += 1
        for j in range(self.n_bins):
            if self.binTime[j] > 0:
                cycleMean = self.binSum[j] / self.binTime[j]
                cycleMeanSqr = self.binSqrSum[j] / self.binTime[j]
                self.binCounts[j] += 1.0
                count = self.binCounts[j]
                if count < 1.5:
                    self.meanFields[j] = cycleMean.copy()
                    self.meanSqr[j] = cycleMeanSqr.copy()
                    self.M2[j] = 0.0
                else:
                    delta = cycleMean - self.meanFields[j]
                    self.meanFields[j] = (self.meanFields[j] * (count - 1.0) + cycleMean) / count
                    delta2 = cycleMean - self.meanFields[j]
                    # True Welford: M2 += delta . delta2 (inner product), matching
                    # BBPATemplates.C finalizeCycle(). The pre-fix library used
                    # 0.5*(magSqr(delta)+magSqr(delta2)), which is not Welford.
                    self.M2[j] += float(np.dot(delta, delta2))
                    self.meanSqr[j] = (self.meanSqr[j] * (count - 1.0) + cycleMeanSqr) / count

    # -- overlap splitting (phaseClock::split) ----------------------------- #
    def _split(self, t0, t1):
        """Split (t0, t1] at every bin boundary it crosses; yields
        (global_bin_index, dt) with exact overlap weights. Mirrors
        phaseClock::split including the endpoint snap semantics (a segment
        ending ON a boundary belongs to the bin BELOW it)."""
        tol = 1e-9 * self.bin_dt
        m = snap_floor((t0 - self.start_time) / self.bin_dt)
        cur = t0
        b = self.start_time + (m + 1) * self.bin_dt
        segs = []
        while b < t1 - tol:
            segs.append((m, b - cur))
            cur = b
            m += 1
            b = self.start_time + (m + 1) * self.bin_dt
        segs.append((m, t1 - cur))
        return segs

    def _advance_to(self, cyc):
        if cyc <= self.previousCycle:
            return
        if self.previousCycle >= self.start_cycle:
            self._finalize_cycle()
        self._reset_cycle()
        self.previousCycle = cyc

    # -- one solver-step execute() (BBPA::execute, v2 overlap contract) ---- #
    def execute(self, t, dt, U):
        U = np.asarray(U, dtype=float)

        for g, dts in self._split(t - dt, t):
            cyc = g // self.n_bins           # floor division (negatives ok)
            binI = g - cyc * self.n_bins
            # pre-checkpoint segments are already folded into the ensemble
            if cyc < self.previousCycle:
                continue
            self._advance_to(cyc)
            if cyc < self.start_cycle:
                continue
            self.currentBin = binI
            self._accumulate(binI, U, dts)

        # a step ending exactly on a cycle boundary completes that cycle
        end_cycle = cycle_number(t, self.start_time, self.period)
        if end_cycle > self.previousCycle:
            self._advance_to(end_cycle)
        self.currentCycle = end_cycle

    # -- replay a list of (t, dt, U) samples ------------------------------- #
    def replay(self, samples):
        for (t, dt, U) in samples:
            self.execute(t, dt, U)
        return self

    # -- write() outputs that land on disk (binItem::write) ---------------- #
    def write_fields(self):
        """
        Returns dict bin -> {PA: vec3, PA_UU: 3x3, PAVariance: scalar or None}
        reproducing exactly what BBPA writes to the phase-aligned dirs
        (BBPATemplates.C binItem::write): _PA and _PA_UU are the PURE finalised
        means over the complete cycles whenever at least one complete cycle
        exists; the running partial is written ONLY when no complete cycle
        exists yet. The in-progress cycle is never blended in -- the old
        (mean*nOld + partial)/(nOld+1) blend weighted a one-timestep partial
        bin as a whole cycle and corrupted the restart checkpoint.
        """
        out = {}
        for j in range(self.n_bins):
            hasCross = self.binCounts[j] > 0
            hasCurrent = self.binTime[j] > 0
            if not (hasCross or hasCurrent):
                continue
            # --- mean (phi_PA): pure finalised; partial only if no cross data ---
            if hasCross:
                PA = self.meanFields[j].copy()
            else:
                inv = 1.0 / self.binTime[j]
                PA = self.binSum[j] * inv
            # --- variance (phi_PAVariance): PURE finalised, only if count>1 ---
            if self.binCounts[j] > 1:
                PAVar = self.M2[j] / (self.binCounts[j] - 1.0)
            else:
                PAVar = None
            # --- second moment (phi_PA_UU): same rule as the mean ---
            if hasCross:
                PAUU = self.meanSqr[j].copy()
            else:
                invT = 1.0 / self.binTime[j]
                PAUU = self.binSqrSum[j] * invT
            out[j] = {"PA": PA, "PA_UU": PAUU, "PAVariance": PAVar}
        return out

    # -- restart: emulate readCheckpoint() restore from written fields ----- #
    @classmethod
    def restore(cls, written, binCounts, nCycles, currentCycle, previousCycle,
                period, n_bins, start_time=0.0, start_cycle=0):
        """
        Reconstruct an accumulator from the on-disk fields exactly as
        BBPA::readCheckpoint() / restoreBinItemFromDisk() does:
            meanFields[j] = _PA[j]
            M2[j]         = _PAVariance[j]*(nCyc-1)   (only if nCyc>1, file present)
            meanSqr[j]    = _PA_UU[j]
        Counters (binCounts, nCycles, currentCycle, previousCycle) come from the
        Properties dict.  Current-cycle accumulators start at zero.
        """
        r = cls(period, n_bins, start_time, start_cycle)
        r.binCounts = np.array(binCounts, dtype=float)
        r.nCycles = int(nCycles)
        r.currentCycle = int(currentCycle)
        r.previousCycle = int(previousCycle)
        for j in range(n_bins):
            if r.binCounts[j] <= 0:
                continue
            w = written.get(j)
            if w is None:
                continue
            r.meanFields[j] = np.array(w["PA"], dtype=float)
            r.meanSqr[j] = np.array(w["PA_UU"], dtype=float)
            if r.binCounts[j] > 1 and w["PAVariance"] is not None:
                r.M2[j] = float(w["PAVariance"]) * (r.binCounts[j] - 1.0)
            else:
                r.M2[j] = 0.0
        return r


# --------------------------------------------------------------------------- #
#  Closed-form analytic references (independent cross-check of the replica).
# --------------------------------------------------------------------------- #
def closed_form_bin_constant(A, c_of_n, n_cycles):
    """
    Bin-constant field U(t)=A[bin]+c[cycle].  Returns, per bin j, the EXACT
    finalised cross-cycle ensemble over cycles 0..n_cycles-1:
       meanFields[j] = A[j] + mean_n c_n
       meanSqr[j]    = mean_n (A[j]+c_n) outer (A[j]+c_n)
       M2[j]         = true-Welford (inner-product) accumulation over x_n = A[j]+c_n
       variance[j]   = M2[j]/(n_cycles-1)
    The Welford M2 here is the SAME inner-product recurrence BBPA uses
    (BBPATemplates.C finalizeCycle: M2 += delta & delta2), so it is the exact
    analytic value of what BBPA computes. For exact per-cycle bin means this
    equals the textbook sum of squared deviations.
    """
    I = len(A)
    cs = [np.asarray(c_of_n(n), dtype=float) for n in range(n_cycles)]
    cbar = np.mean(cs, axis=0)
    res = {}
    for j in range(I):
        Aj = np.asarray(A[j], dtype=float)
        xs = [Aj + cn for cn in cs]                      # per-cycle bin mean
        meanF = Aj + cbar
        meanSqr = np.mean([outer_sym(x) for x in xs], axis=0)
        # true Welford inner-product recurrence, identical to BBPA
        mean_run = None
        M2 = 0.0
        for k, x in enumerate(xs):
            count = k + 1
            if count == 1:
                mean_run = x.copy()
                M2 = 0.0
            else:
                delta = x - mean_run
                mean_run = (mean_run * (count - 1) + x) / count
                delta2 = x - mean_run
                # True Welford inner-product form, matching BBPATemplates.C.
                M2 += float(np.dot(delta, delta2))
        var = M2 / (n_cycles - 1) if n_cycles > 1 else None
        res[j] = {"PA": meanF, "PA_UU": meanSqr, "PAVariance": var, "M2": M2}
    return res


def smooth_field(t, period, amp, camp, c_of_n):
    """Smooth analytic field U(t) used for the discretisation (first-order) study."""
    ph = math.fmod(t, period)
    w = 2.0 * math.pi / period
    n = int(t / period)
    base = np.array([amp * math.cos(w * ph),
                     amp * math.sin(w * ph),
                     amp * camp * math.cos(2.0 * w * ph)])
    return base + np.asarray(c_of_n(n), dtype=float)


def closed_form_smooth_binmean(period, n_bins, amp, camp, c_of_n, n_cycles):
    r"""
    Exact continuous bin mean of the smooth field, averaged equally across
    cycles 0..n_cycles-1 (matches BBPA's equal-weight cross-cycle mean):
        meanFields[j] = (1/Delta) \int_{jD}^{(j+1)D} base(ph) dph + mean_n c_n
    The phase integral is done analytically.
    """
    I = n_bins
    D = period / I
    w = 2.0 * math.pi / period
    cs = [np.asarray(c_of_n(n), dtype=float) for n in range(n_cycles)]
    cbar = np.mean(cs, axis=0)
    res = {}
    for j in range(I):
        a, b = j * D, (j + 1) * D
        # (1/D) * integral of [cos(w ph), sin(w ph), camp cos(2w ph)]
        mx = amp * (math.sin(w * b) - math.sin(w * a)) / (w * D)
        my = amp * (-(math.cos(w * b) - math.cos(w * a))) / (w * D)
        mz = amp * camp * (math.sin(2 * w * b) - math.sin(2 * w * a)) / (2 * w * D)
        res[j] = np.array([mx, my, mz]) + cbar
    return res


# --------------------------------------------------------------------------- #
#  Minimal OpenFOAM field reader: first internal-cell value of a uniform or
#  nonuniform internalField (scalar / vector / symmTensor).
# --------------------------------------------------------------------------- #
import re


def read_of_internal_first(path):
    """Return the first internal-cell value as a float (scalar) or np.array."""
    with open(path, "r") as f:
        txt = f.read()
    m = re.search(r"internalField\s+(.*?);", txt, re.DOTALL)
    if not m:
        raise ValueError(f"no internalField in {path}")
    body = m.group(1).strip()
    if body.startswith("uniform"):
        rest = body[len("uniform"):].strip()
        return _parse_of_value(rest)
    # nonuniform List<...> N ( v0 v1 ... )
    paren = re.search(r"\((.*)\)", body, re.DOTALL)
    inner = paren.group(1).strip()
    # first token group: either a number (scalar list) or a (...) tuple
    if inner.startswith("("):
        first = re.search(r"\(([^()]*)\)", inner).group(1)
        return np.array([float(x) for x in first.split()])
    else:
        return float(inner.split()[0])


def _parse_of_value(s):
    s = s.strip()
    if s.startswith("("):
        nums = [float(x) for x in re.search(r"\((.*)\)", s).group(1).split()]
        return np.array(nums)
    return float(s)


def symm6_to_mat(v6):
    """OF symmTensor order (xx xy xz yy yz zz) -> 3x3 symmetric matrix."""
    xx, xy, xz, yy, yz, zz = v6
    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])


def mat_to_symm6(m):
    return np.array([m[0, 0], m[0, 1], m[0, 2], m[1, 1], m[1, 2], m[2, 2]])
