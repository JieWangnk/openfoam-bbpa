/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "BBPA.H"
#include "volFields.H"
#include "OSspecific.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Inner product of two fields of the tracked type, returning the scalar field
//  used by the Welford cycle-variance accumulator  M2 = sum |x - <x>|^2.
//  Overloaded rather than templated: for a vector field the update is the dot
//  product, while for a scalar field operator& is not defined and it is just
//  the product. Defined here (not in BBPA.H) because BBPA.H carries only
//  volFieldsFwd.H, so the field operators are not visible there.
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
namespace functionObjects
{

inline tmp<volScalarField> bbpaInner
(
    const volScalarField& a,
    const volScalarField& b
)
{
    return a*b;
}

inline tmp<volScalarField> bbpaInner
(
    const volVectorField& a,
    const volVectorField& b
)
{
    return a & b;
}

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::functionObjects::BBPA::binItem<Type>::binItem
(
    const word& name,
    const fvMesh& mesh,
    const label nBins,
    const labelList& activeBins,
    const bool keepLastCycle
)
:
    fieldName_(name),
    integ_(nBins),
    stats_(nBins),
    lastCycleFields_(keepLastCycle ? nBins : 0),
    lastCycleBinTime_(keepLastCycle ? nBins : 0, 0.0),
    binFlushed_(keepLastCycle ? nBins : 0, false)
{
    // Look up field to get dimensions
    const FieldType& field = mesh.lookupObject<FieldType>(name);

    // Determine which bins to allocate. An EMPTY activeBins list means
    // "all bins" (original behaviour); a non-empty list activates only
    // the listed bins, leaving the rest as null PtrList entries that the
    // accumulate / finalize / write paths skip.
    boolList isActive(nBins, activeBins.empty());
    forAll(activeBins, k)
    {
        isActive[activeBins[k]] = true;
    }

    forAll(stats_.mean_, binI)
    {
        if (!isActive[binI])
        {
            // Leave this bin's accumulators unallocated (null). Skipped
            // everywhere downstream by the stats_.binCounts_/integ_.binTime_
            // and PtrList::set() guards.
            continue;
        }
        stats_.mean_.set
        (
            binI,
            new FieldType
            (
                IOobject
                (
                    name + "_mean_bin" + Foam::name(binI),
                    mesh.time().name(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensioned<Type>(field.dimensions(), Zero)
            )
        );

        stats_.M2_.set
        (
            binI,
            new volScalarField
            (
                IOobject
                (
                    name + "_M2_bin" + Foam::name(binI),
                    mesh.time().name(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar(sqr(field.dimensions()), Zero)
            )
        );

        integ_.binSum_.set
        (
            binI,
            new FieldType
            (
                IOobject
                (
                    name + "_currentBinSum_bin" + Foam::name(binI),
                    mesh.time().name(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensioned<Type>(field.dimensions(), Zero)
            )
        );

        stats_.meanSqr_.set
        (
            binI,
            new SqrFieldType
            (
                IOobject
                (
                    name + "_meanSqr_bin" + Foam::name(binI),
                    mesh.time().name(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensioned<SqrType>(sqr(field.dimensions()), Zero)
            )
        );

        integ_.binSqrSum_.set
        (
            binI,
            new SqrFieldType
            (
                IOobject
                (
                    name + "_currentBinSqr_bin" + Foam::name(binI),
                    mesh.time().name(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensioned<SqrType>(sqr(field.dimensions()), Zero)
            )
        );

        // Note: internal fields may appear in output dirs due to OF12
        // auto-write. This is cosmetic — they can be cleaned up with
        // foamListTimes or a post-processing script.
    }
}


template<class Type>
Foam::label
Foam::functionObjects::BBPA::binItem<Type>::restoreCurrentCycleBins
(
    const fvMesh& mesh,
    const scalar cycleStartTime,
    const scalar binDeltaT,
    const scalar now,
    const scalar centreOffset
)
{
    if (binFlushed_.empty()) return 0;

    label recovered = 0;

    forAll(integ_.binSum_, binI)
    {
        if (!integ_.binSum_.set(binI)) continue;

        // Recover any bin whose window STARTED before the restart instant and
        // whose complete <field>_Bin is on disk. This deliberately INCLUDES the
        // bin that CONTAINS the restart instant (the straddled bin): the flush
        // only ever writes a CLOSED window, so an on-disk _Bin is the WHOLE
        // window average -- strictly better than re-accumulating only its
        // post-restart part. This is what lets a restart from an EARLIER write
        // (e.g. one write before the kill) lose nothing even with binCentred:
        // the bin centred on that write has its complete _Bin on disk from
        // before the kill, so it is recovered rather than straddle-dropped.
        //
        // A bin whose window starts at or after the restart has no pre-restart
        // data (its dir was never written, or was trimmed) and is re-accumulated
        // by the resumed solver instead.
        const scalar binStart =
            cycleStartTime + binI*binDeltaT - centreOffset;
        if (binStart >= now - 1e-9*binDeltaT) continue;

        const scalar binTime = cycleStartTime + binI*binDeltaT;
        const word binTimeName = Time::timeName(binTime);

        const fileName binPath =
            mesh.time().path()/binTimeName/(fieldName_ + "_Bin");
        const fileName uuPath =
            mesh.time().path()/binTimeName/(fieldName_ + "_Bin_UU");

        // Both files or neither: a bin recovered with a mean but no second
        // moment would fold a zero <UiUj> in at full weight.
        if (!Foam::isFile(binPath) || !Foam::isFile(uuPath)) continue;

        FieldType loaded
        (
            IOobject
            (
                fieldName_ + "_Bin",
                binTimeName,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh
        );

        SqrFieldType loadedUU
        (
            IOobject
            (
                fieldName_ + "_Bin_UU",
                binTimeName,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh
        );

        // The files hold the bin MEANS; the accumulators hold time-weighted
        // SUMS. The bin is complete, so its weight is exactly binDeltaT.
        integ_.binSum_[binI] = loaded*binDeltaT;
        integ_.binSqrSum_[binI] = loadedUU*binDeltaT;
        integ_.binTime_[binI] = binDeltaT;
        binFlushed_[binI] = true;      // already on disk; do not rewrite
        recovered++;
    }

    return recovered;
}


template<class Type>
void Foam::functionObjects::BBPA::binItem<Type>::resetCurrentCycle()
{
    binFlushed_ = false;
    integ_.reset();
}


template<class Type>
void Foam::functionObjects::BBPA::binItem<Type>::accumulate
(
    const label binI,
    const FieldType& field,
    const scalar dt
)
{
    // Sparse-bin mode: if the current phase lands in a bin that was not
    // allocated, do nothing -- and crucially do NOT touch the bin weight,
    // so finalizeCycle()/write() skip it too.
    if (!integ_.active(binI)) return;

    // A bin already recovered from disk this cycle (binFlushed_) holds its
    // COMPLETE window average. When a restart resumes from an EARLIER write the
    // solver re-runs that bin's phase; adding those steps again would
    // double-count the recovered window, so skip it. binFlushed_ is empty
    // outside "bin" mode (no-op there), and in normal forward running phaseBin
    // never revisits a flushed bin within a cycle (no-op there too) -- this
    // only bites on a restart-from-earlier re-run.
    if (!binFlushed_.empty() && binFlushed_[binI]) return;

    // Time-weighted accumulation: phi(t)*dt (first moment) and
    // sqr(phi(t))*dt (second moment; symmTensor phi_i phi_j for vector phi,
    // enabling full Reynolds-stress recovery offline).
    integ_.integrate(binI, field, dt);
}


template<class Type>
void Foam::functionObjects::BBPA::binItem<Type>::finalizeCycle
(
    const label skipBin
)
{
    stats_.nCycles_++;

    // Update overall statistics for each bin that has data
    forAll(integ_.binSum_, binI)
    {
        // Skip bins not allocated in sparse-bin mode (also have time==0).
        if (!integ_.binSum_.set(binI)) continue;

        // Restart-straddled bin: a restart resumed at an instant strictly
        // inside this bin's phase window, so the cycle covered only the part
        // of the window lying after that instant. Its mean is a partial-window
        // mean, biased toward whichever phases happen to follow the restart,
        // and folding it at full weight would carry that bias into the phase
        // average -- while looking entirely healthy on disk.
        //
        // Drop it exactly as a bin that received no data at all is dropped:
        // stats_.binCounts_ is not incremented, so no cycle is ever credited without
        // the data to back it. One cycle noisier, but UNBIASED -- the same
        // trade the missing-checkpoint recovery makes.
        //
        // lastCycleBinTime_ is zeroed too, so "bin" mode does not emit a
        // <field>_Bin file built from the partial window either.
        if (binI == skipBin)
        {
            if (lastCycleFields_.size())
            {
                lastCycleBinTime_[binI] = 0;
            }
            continue;
        }

        // "bin" mode: capture THIS cycle's bin before it is folded into the
        // running mean and the accumulators are reset. write() runs after
        // resetAllBinsForNewCycle(), so this is the only point at which the
        // completed cycle's bin still exists.
        if (lastCycleFields_.size())
        {
            lastCycleBinTime_[binI] = integ_.binTime_[binI];
            if (integ_.binTime_[binI] > 0)
            {
                if (!lastCycleFields_.set(binI))
                {
                    lastCycleFields_.set
                    (
                        binI,
                        new FieldType
                        (
                            IOobject
                            (
                                fieldName_ + "_lastCycle_bin"
                              + Foam::name(binI),
                                integ_.binSum_[binI].mesh().time().name(),
                                integ_.binSum_[binI].mesh(),
                                IOobject::NO_READ,
                                IOobject::NO_WRITE,
                                false
                            ),
                            integ_.binSum_[binI]/integ_.binTime_[binI]
                        )
                    );
                }
                else
                {
                    lastCycleFields_[binI] =
                        integ_.binSum_[binI]/integ_.binTime_[binI];
                }
            }
        }

        if (integ_.binTime_[binI] > 0)
        {
            const FieldType cycleMean =
                integ_.binSum_[binI] / integ_.binTime_[binI];

            const SqrFieldType cycleMeanSqr =
                integ_.binSqrSum_[binI] / integ_.binTime_[binI];

            stats_.binCounts_[binI] += 1.0;
            const scalar count = stats_.binCounts_[binI];

            if (count < 1.5)
            {
                stats_.mean_[binI] = cycleMean;
                stats_.meanSqr_[binI] = cycleMeanSqr;
                stats_.M2_[binI] = dimensionedScalar(stats_.M2_[binI].dimensions(), 0);
            }
            else
            {
                // Welford's online algorithm for cycle-mean variance:
                // M2 accumulates sum of squared deviations |x - <x>|^2
                // of the BIN-AVERAGED value across cycles.
                const FieldType delta(cycleMean - stats_.mean_[binI]);

                stats_.mean_[binI] =
                    (stats_.mean_[binI] * (count - 1.0) + cycleMean) / count;

                const FieldType delta2(cycleMean - stats_.mean_[binI]);

                // Welford: M2 += (x - mean_old) . (x - mean_new). For a vector
                // field, M2 is the scalar sum |x - <x>|^2, so the update is the
                // INNER PRODUCT of the two deltas.
                //
                // This was previously 0.5*(magSqr(delta) + magSqr(delta2)), the
                // mean of the squared magnitudes, which is not Welford and is
                // biased HIGH: by AM-GM, 0.5*(|a|^2 + |b|^2) >= a.b for all a,b,
                // with equality only when a == b. Worked example for two cycles
                // v1, v2: the true sum of squared deviations is |v2-v1|^2/2,
                // which delta.delta2 reproduces exactly, whereas the old form
                // gave 0.625*|v2-v1|^2 -- 25% too large. The error is invisible
                // on a periodic flow (every cycle identical => variance ~ 0),
                // which is why it survived; tutorials/planarPoiseuille exercises
                // it with a developing flow whose cycles genuinely differ.
                stats_.M2_[binI] += bbpaInner(delta, delta2);

                // Running cross-cycle mean of <phi*phi>_bin. Combined
                // with stats_.mean_ this recovers Reynolds stress
                // R_ij = <U_i U_j> - <U_i><U_j> and TKE = 0.5 trace(R)
                // offline from the written fields.
                stats_.meanSqr_[binI] =
                    (stats_.meanSqr_[binI] * (count - 1.0) + cycleMeanSqr) / count;
            }
        }
    }
}


// * * * * * * * * * * * * Restore from checkpoint * * * * * * * * * * * * * //
//
// At restart the existing _PA / _PAVariance / _PA_UU files written
// at the most recent solver writeTime carry the cycle-N ensemble
// state.  Per binItem<Type>::write() at the post-finalisation cycle
// boundary (stats_.binCounts_[binI] > 0 && integ_.binTime_[binI] == 0):
//
//     paField     = stats_.mean_[binI]                  (no blending)
//     paVariance  = stats_.M2_[binI] / (n - 1)          (only if n > 1)
//     paUU        = stats_.meanSqr_[binI]               (no blending)
//
// Reverse-engineer the three internal arrays from the public outputs.
template<class Type>
void Foam::functionObjects::BBPA::restoreBinItemFromDisk
(
    binItem<Type>& bi,
    const word& fieldName,
    const scalar cycleStartTime,
    const word& paSuffix
)
{
    typedef typename binItem<Type>::FieldType FieldType;
    typedef typename binItem<Type>::SqrFieldType SqrFieldType;

    // Mesh handle from the first allocated bin (sparse-bin mode may
    // leave bin 0 unallocated).
    const fvMesh* meshPtr = nullptr;
    forAll(bi.stats_.mean_, b)
    {
        if (bi.stats_.mean_.set(b)) { meshPtr = &bi.stats_.mean_[b].mesh(); break; }
    }
    if (meshPtr == nullptr) return;
    const fvMesh& mesh = *meshPtr;
    label restored = 0;

    forAll(bi.stats_.mean_, binI)
    {
        // Skip bins not allocated in sparse-bin mode.
        if (!bi.stats_.mean_.set(binI)) continue;

        const scalar nCyc = bi.stats_.binCounts_[binI];
        if (nCyc <= 0) continue;  // bin never received data

        const scalar binTime = cycleStartTime + binI * binDeltaT_;
        const word binTimeName = Time::timeName(binTime);

        // 1. mean: read fieldName_PA into stats_.mean_[binI].
        {
            IOobject paIO
            (
                fieldName + paSuffix,
                binTimeName,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            const fileName paPath =
                mesh.time().path() / binTimeName / paIO.name();
            if (!Foam::isFile(paPath))
            {
                // The mean for this bin is gone but its cycle count survived in
                // the Properties dict.  Keeping the count would make Welford
                // average the remaining cycles against a phantom zero -- with
                // count n0 and k further cycles the bin converges to
                // k*x/(n0+k), i.e. biased LOW by n0/(n0+k) and indistinguishable
                // from a healthy value (non-zero, full-size file).  Reset the
                // bin's history instead: it loses the lost cycles and is
                // noisier, but it is UNBIASED.
                WarningInFunction
                    << "Checkpoint mean field missing: " << binTimeName
                    << "/" << fieldName << paSuffix << " -- bin " << binI
                    << " reset to a clean cold start (count " << nCyc
                    << " -> 0) to avoid a low-biased average." << endl;
                bi.stats_.binCounts_[binI] = 0;
                bi.stats_.mean_[binI] =
                    dimensioned<Type>(bi.stats_.mean_[binI].dimensions(), Zero);
                bi.stats_.M2_[binI] =
                    dimensionedScalar(bi.stats_.M2_[binI].dimensions(), 0);
                bi.stats_.meanSqr_[binI] =
                    dimensioned<typename binItem<Type>::SqrType>
                    (bi.stats_.meanSqr_[binI].dimensions(), Zero);
                continue;
            }

            // Whole-field assignment, not primitiveFieldRef(): restoring only
            // the internal field leaves the accumulator's BOUNDARY at its
            // constructed zero, so every restart silently discarded the
            // boundary phase-average and then folded new cycles against a
            // phantom zero history. For wallShearStress -- whose data is
            // entirely on the wall patch -- that meant the restart threw away
            // the whole accumulated field while the Properties counters still
            // claimed the cycles were there.
            FieldType loaded(paIO, mesh);
            bi.stats_.mean_[binI] = loaded;
        }

        // 2. M2: read fieldName_PAVariance and rescale by (n-1).
        if (nCyc > 1)
        {
            IOobject varIO
            (
                fieldName + paSuffix + "Variance",
                binTimeName,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            const fileName varPath =
                mesh.time().path() / binTimeName / varIO.name();
            if (Foam::isFile(varPath))
            {
                volScalarField loadedVar(varIO, mesh);
                bi.stats_.M2_[binI] = loadedVar * (nCyc - 1.0);
            }
        }

        // 3. meanSqr: read fieldName_PA_UU into stats_.meanSqr_[binI].
        {
            IOobject uuIO
            (
                fieldName + paSuffix + "_UU",
                binTimeName,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            const fileName uuPath =
                mesh.time().path() / binTimeName / uuIO.name();
            if (Foam::isFile(uuPath))
            {
                SqrFieldType loadedUU(uuIO, mesh);
                bi.stats_.meanSqr_[binI] = loadedUU;
            }
        }

        restored++;
    }

    Info<< "    " << fieldName << ": restored "
        << restored << "/" << bi.stats_.mean_.size()
        << " bins from phase-aligned dirs (cycleStart="
        << cycleStartTime << ", nCycles=" << bi.stats_.nCycles_ << ")"
        << endl;
}


// Explicit instantiation for common types
template struct Foam::functionObjects::BBPA::binItem<Foam::scalar>;
template struct Foam::functionObjects::BBPA::binItem<Foam::vector>;
// Tensor/symmTensor instantiations removed: OpenFOAM doesn't provide a
// generic sqr(tensor) primitive, and BBPA is only used for volScalarField
// (e.g. pressure) and volVectorField (e.g. velocity, wallShearStress).

// Explicit instantiation of the restore helper.
template void
Foam::functionObjects::BBPA::restoreBinItemFromDisk<Foam::scalar>
(
    binItem<Foam::scalar>&,
    const word&,
    const scalar,
    const word&
);
template void
Foam::functionObjects::BBPA::restoreBinItemFromDisk<Foam::vector>
(
    binItem<Foam::vector>&,
    const word&,
    const scalar,
    const word&
);


// ************************************************************************* //
