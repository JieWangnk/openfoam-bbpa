/*---------------------------------------------------------------------------*\
Class
    Foam::functionObjects::binIntegrator
    Foam::functionObjects::cycleStatistics

Description
    The two data components of bin-based phase averaging, split by
    responsibility (v2 internal separation):

    binIntegrator<Type> -- the CURRENT cycle's time-weighted accumulation:
        W_j   = sum dt                over steps in bin j this cycle
        S_j   = sum dt * phi          (first moment)
        Q_j   = sum dt * sqr(phi)     (second moment; symmTensor for vector)
    Reset at every cycle boundary. Knows nothing about cycles, ensembles,
    files, or restarts.

    cycleStatistics<Type> -- the CROSS-CYCLE ensemble over completed cycles:
        binCounts_[j]  completed cycles folded into bin j
        mean_[j]       running equal-weight mean of the per-cycle bin means
        meanSqr_[j]    running mean of the per-cycle bin second moments
        M2_[j]         Welford sum of squared deviations (inner-product form)
    Updated only by the fold at a cycle boundary (binItem::finalizeCycle,
    which reads the integrator and writes here). Knows nothing about the
    solver clock or output layout.

    Both are plain data holders with their smallest possible behaviour;
    allocation (field sizing, IOobject naming, sparse-bin activation) stays
    with their owner, BBPA::binItem. Separating the data this way makes the
    boundary explicit: execute() may touch the integrator, the cycle fold is
    the only writer of the statistics, and write() reads both without
    modifying either.

\*---------------------------------------------------------------------------*/

#ifndef binAccumulators_H
#define binAccumulators_H

#include "volFieldsFwd.H"
#include "scalarField.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

//- The squared type of an averaged field: sqr(vector) is a symmTensor
//  (<U_i U_j>); scalars square to scalars.
template<class T> struct bbpaSqrTypeOf { typedef T type; };
template<> struct bbpaSqrTypeOf<vector> { typedef symmTensor type; };


/*---------------------------------------------------------------------------*\
                        Struct binIntegrator Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
struct binIntegrator
{
    typedef VolField<Type> FieldType;
    typedef typename bbpaSqrTypeOf<Type>::type SqrType;
    typedef VolField<SqrType> SqrFieldType;

    //- Time accumulated in each bin this cycle (the bin weights W_j)
    scalarField binTime_;

    //- Time-weighted first-moment sum per bin
    PtrList<FieldType> binSum_;

    //- Time-weighted second-moment sum per bin
    PtrList<SqrFieldType> binSqrSum_;

    //- Size the containers; the owner allocates the per-bin fields
    binIntegrator(const label nBins)
    :
        binTime_(nBins, 0.0),
        binSum_(nBins),
        binSqrSum_(nBins)
    {}

    //- Is this bin allocated? (sparse-bin mode leaves inactive bins null)
    bool active(const label binI) const
    {
        return binSum_.set(binI);
    }

    //- One solver step's contribution to bin binI
    void integrate(const label binI, const FieldType& field, const scalar dt)
    {
        binTime_[binI] += dt;
        binSum_[binI] += dt*field;
        binSqrSum_[binI] += dt*sqr(field);
    }

    //- Zero everything for a new cycle
    void reset()
    {
        binTime_ = 0.0;

        forAll(binSum_, binI)
        {
            if (!binSum_.set(binI)) continue;

            binSum_[binI] =
                dimensioned<Type>(binSum_[binI].dimensions(), Zero);
            binSqrSum_[binI] =
                dimensioned<SqrType>(binSqrSum_[binI].dimensions(), Zero);
        }
    }
};


/*---------------------------------------------------------------------------*\
                       Struct cycleStatistics Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
struct cycleStatistics
{
    typedef VolField<Type> FieldType;
    typedef typename bbpaSqrTypeOf<Type>::type SqrType;
    typedef VolField<SqrType> SqrFieldType;

    //- Completed cycles folded in (all bins; the Info-line counter)
    label nCycles_;

    //- Completed cycles folded into each bin
    scalarField binCounts_;

    //- Running equal-weight mean of the per-cycle bin means
    PtrList<FieldType> mean_;

    //- Welford M2: sum of squared deviations of the per-cycle bin means
    //  from their running mean (inner-product form; variance = M2/(n-1))
    PtrList<volScalarField> M2_;

    //- Running mean of the per-cycle bin second moments <phi*phi>
    PtrList<SqrFieldType> meanSqr_;

    //- Size the containers; the owner allocates the per-bin fields
    cycleStatistics(const label nBins)
    :
        nCycles_(0),
        binCounts_(nBins, 0.0),
        mean_(nBins),
        M2_(nBins),
        meanSqr_(nBins)
    {}
};

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
