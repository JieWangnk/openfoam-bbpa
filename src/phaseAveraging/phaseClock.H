/*---------------------------------------------------------------------------*\
Class
    Foam::phaseClock

Description
    The single home for phase/cycle/bin boundary arithmetic shared by the
    phase-averaging function objects (BBPA, IPA).

    Every quantity that decides WHICH bin or cycle an instant belongs to is
    computed here and nowhere else. The motivation is drift: this arithmetic
    previously existed in four independent copies (BBPA, IPA, the Python
    verification replica, and the test harness's prescribe function object),
    and a boundary-handling fix applied to one copy silently left the others
    wrong. One class, one behaviour.

    The boundary handling itself (see snapFloor and wrapPhase):

    The solver clock accumulates dt in floating point, so an on-grid time
    arrives as e.g. fmod(2.9, 1.0) = 0.8999999999999999, and a plain
    truncation files the step one bin (or one cycle) LOW. The error grows
    with t (fewer mantissa bits left for the fraction), so a run drifts from
    51 steps per bin in the first cycle to 49 by the third instead of a
    constant 50. At I = 200 a bin can be ~25 timesteps wide, so one misplaced
    step is ~4% of the bin.

    snapFloor therefore floors with a snap-up when the argument sits within
    1e-9 (in bin units) below an integer -- ~7 orders of magnitude smaller
    than any sane timestep-to-bin ratio, so it can only ever absorb
    floating-point noise, never a step that genuinely belongs to the lower
    bin. wrapPhase additionally maps an fmod result of period-1ulp back to
    phase 0, which would otherwise put the first step of a cycle in the LAST
    bin of that cycle.

    Two origins are kept:
      - labelOrigin: bin LABELS (output times) sit at
        labelOrigin + cycle*period + bin*binDeltaT.
      - phaseOrigin: the integration WINDOW origin used to classify an
        instant; equals labelOrigin, or labelOrigin - binDeltaT/2 when the
        bins are centred on their labels (binCentred), so that a bin centred
        on phase 0 lies wholly inside one cycle instead of straddling the
        boundary and losing its first half to the cycle reset.

\*---------------------------------------------------------------------------*/

#ifndef phaseClock_H
#define phaseClock_H

#include "scalar.H"
#include "label.H"
#include "DynamicList.H"

#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class phaseClock Declaration
\*---------------------------------------------------------------------------*/

class phaseClock
{
    // Private Data

        //- Cycle period
        scalar period_;

        //- Number of phase bins per period
        label nBins_;

        //- Bin width: period/nBins
        scalar binDeltaT_;

        //- Origin of the bin LABELS (output times)
        scalar labelOrigin_;

        //- Origin used to classify an instant into a phase/cycle;
        //  labelOrigin_ - binDeltaT_/2 when centred, else labelOrigin_
        scalar phaseOrigin_;

        //- Are bins centred on their labels?
        bool centred_;


public:

    // Constructors

        //- Construct null; reset() before use
        phaseClock()
        :
            period_(0),
            nBins_(0),
            binDeltaT_(0),
            labelOrigin_(0),
            phaseOrigin_(0),
            centred_(false)
        {}


    // Member Functions

        //- Configure. Called whenever period/nBins/origin/centring change
        //  (dictionary read and checkpoint restore).
        void reset
        (
            const scalar period,
            const label nBins,
            const scalar labelOrigin,
            const bool centred
        )
        {
            period_ = period;
            nBins_ = nBins;
            binDeltaT_ = period/nBins;
            labelOrigin_ = labelOrigin;
            centred_ = centred;
            phaseOrigin_ = centred ? labelOrigin - 0.5*binDeltaT_ : labelOrigin;
        }


        // Boundary-safe primitives

        //- Floor with a snap-up within 1e-9 below an integer
        static label snapFloor(const scalar x)
        {
            const scalar r = std::floor(x);
            return label((x - r) > (1.0 - 1e-9) ? r + 1.0 : r);
        }

        //- Phase of instant t in [0, period), with the period-1ulp fmod
        //  wrap-around guard
        scalar wrapPhase(const scalar t) const
        {
            scalar phase = std::fmod(t - phaseOrigin_, period_);
            if (phase > period_*(1.0 - 1e-12)) phase = 0;
            return phase;
        }


        // Classification

        //- Phase-bin index of instant t, in [0, nBins)
        label phaseBin(const scalar t) const
        {
            return
                std::min
                (
                    snapFloor(wrapPhase(t)/binDeltaT_),
                    nBins_ - 1
                );
        }

        //- Cycle number of instant t
        label cycleNumber(const scalar t) const
        {
            return snapFloor((t - phaseOrigin_)/period_);
        }


        // Label times

        //- Start time (label) of a cycle
        scalar cycleStart(const label cycle) const
        {
            return labelOrigin_ + cycle*period_;
        }

        //- Label time of bin binI of a cycle
        scalar binLabelTime(const label cycle, const label binI) const
        {
            return cycleStart(cycle) + binI*binDeltaT_;
        }

        //- Half-bin offset between a centred bin's label and its window start
        scalar centreOffset() const
        {
            return centred_ ? 0.5*binDeltaT_ : 0.0;
        }


        // Overlap integration (the v2 numerical contract)

        //- One integration segment of a solver step: the part of the
        //  interval lying inside a single bin of a single cycle
        struct segment
        {
            label cycle;
            label bin;
            scalar dt;
        };

        //- Split the solver interval (t0, t1] at every bin boundary it
        //  crosses, so each returned segment is the exact overlap
        //  dt(n,c,j) = max(0, min(t1,b(j+1)) - max(t0,b(j))) with one
        //  (cycle, bin). Consequently sum(segments.dt) == t1 - t0 for
        //  every step and sum_j W(c,j) == period for every completed
        //  cycle. The snap classifies an ENDPOINT lying numerically on a
        //  boundary (a segment ending on a boundary belongs to the bin
        //  BELOW it; one starting there to the bin above) -- it never
        //  assigns a genuinely straddling interval whole to one bin.
        void split
        (
            const scalar t0,
            const scalar t1,
            DynamicList<segment>& segs
        ) const
        {
            segs.clear();

            const scalar tol = 1e-9*binDeltaT_;

            // global bin index just above t0 (snap: t0 on a boundary
            // starts the interval in the bin ABOVE it, consistent with
            // phaseBin of points immediately after t0)
            label m = snapFloor((t0 - phaseOrigin_)/binDeltaT_);
            scalar cur = t0;
            scalar b = phaseOrigin_ + (m + 1)*binDeltaT_;

            while (b < t1 - tol)
            {
                segs.append(globalSegment(m, b - cur));
                cur = b;
                ++m;
                b = phaseOrigin_ + (m + 1)*binDeltaT_;
            }
            segs.append(globalSegment(m, t1 - cur));
        }


        // Access

        scalar period() const { return period_; }
        label nBins() const { return nBins_; }
        scalar binDeltaT() const { return binDeltaT_; }
        scalar labelOrigin() const { return labelOrigin_; }
        scalar phaseOrigin() const { return phaseOrigin_; }
        bool centred() const { return centred_; }


private:

        //- Map a global bin index (bins counted from the phase origin
        //  across all cycles) to a (cycle, bin, dt) segment
        segment globalSegment(const label g, const scalar dt) const
        {
            // floor division, correct for negative g
            label cyc = g/nBins_;
            label bin = g - cyc*nBins_;
            if (bin < 0)
            {
                bin += nBins_;
                cyc -= 1;
            }
            return segment{cyc, bin, dt};
        }
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
