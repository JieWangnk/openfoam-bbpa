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


        // Access

        scalar period() const { return period_; }
        label nBins() const { return nBins_; }
        scalar binDeltaT() const { return binDeltaT_; }
        scalar labelOrigin() const { return labelOrigin_; }
        scalar phaseOrigin() const { return phaseOrigin_; }
        bool centred() const { return centred_; }
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
