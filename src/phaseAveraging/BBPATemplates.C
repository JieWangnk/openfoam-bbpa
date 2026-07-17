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
//  writePAField: write a volField to a phase-aligned time directory using an
//  explicit OFstream. regIOobject::write() in OpenFOAM is hard-wired to the
//  current solver time: writing with instance != time.timeName() silently
//  fails for time dirs that don't match a solver writeTime. This helper
//  bypasses that restriction by constructing the path explicitly and writing
//  the full FoamFile format (header + data) via an OFstream. It handles the
//  parallel case (prepending processorN to the path) and the serial case.
//
//  This is a function-object-local write; the file does not register with
//  the objectRegistry, so repeated writes are idempotent.
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
template<class FieldT>
static void writePAField
(
    const FieldT& paField,
    const Foam::fvMesh& mesh,
    const Foam::word& binTimeName
)
{
    // mesh.time().path() == rootPath()/caseName(); in a parallel run
    // caseName() is already the processor-local case (".../processorN"),
    // so the processorN segment must NOT be appended again. The previous
    // code did rootPath/caseName + processorN, producing the
    // processorN/processorN/<time> nesting -- and the write location then
    // disagreed with restoreBinItemFromDisk()'s read path
    // (mesh.time().path()/binTime), silently breaking parallel restart.
    const Foam::fileName timeDir = mesh.time().path() / binTimeName;
    Foam::mkDir(timeDir);

    Foam::OFstream os(timeDir / paField.name());
    paField.regIOobject::writeHeader(os);
    os << paField;
    Foam::IOobject::writeEndDivider(os);
}

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
    nCycles_(0),
    binCounts_(nBins, 0.0),
    meanFields_(nBins),
    M2Fields_(nBins),
    currentCycleBinTime_(nBins, 0.0),
    currentCycleBinSum_(nBins),
    meanSqrFields_(nBins),
    currentCycleBinSqrSum_(nBins),
    lastCycleFields_(keepLastCycle ? nBins : 0),
    lastCycleBinTime_(keepLastCycle ? nBins : 0, 0.0)
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

    forAll(meanFields_, binI)
    {
        if (!isActive[binI])
        {
            // Leave this bin's accumulators unallocated (null). Skipped
            // everywhere downstream by the binCounts_/currentCycleBinTime_
            // and PtrList::set() guards.
            continue;
        }
        meanFields_.set
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

        M2Fields_.set
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

        currentCycleBinSum_.set
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

        meanSqrFields_.set
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

        currentCycleBinSqrSum_.set
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
void Foam::functionObjects::BBPA::binItem<Type>::resetCurrentCycle()
{
    currentCycleBinTime_ = 0.0;

    forAll(currentCycleBinSum_, binI)
    {
        // Skip bins not allocated in sparse-bin mode.
        if (!currentCycleBinSum_.set(binI)) continue;

        currentCycleBinSum_[binI] =
            dimensioned<Type>(currentCycleBinSum_[binI].dimensions(), Zero);
        currentCycleBinSqrSum_[binI] =
            dimensioned<SqrType>(currentCycleBinSqrSum_[binI].dimensions(), Zero);
    }
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
    // allocated, do nothing -- and crucially do NOT touch
    // currentCycleBinTime_, so finalizeCycle()/write() skip it too.
    if (!currentCycleBinSum_.set(binI)) return;

    // Time-weighted accumulation: phi(t)*dt (first moment) and
    // sqr(phi(t))*dt (second moment). For vector phi, sqr() returns
    // a symmTensor (phi_i phi_j) enabling full Reynolds-stress recovery
    // offline. For scalar phi, sqr() is the scalar square.
    currentCycleBinTime_[binI] += dt;
    currentCycleBinSum_[binI] += dt * field;
    currentCycleBinSqrSum_[binI] += dt * sqr(field);
}


template<class Type>
void Foam::functionObjects::BBPA::binItem<Type>::finalizeCycle
(
    const label skipBin
)
{
    nCycles_++;

    // Update overall statistics for each bin that has data
    forAll(currentCycleBinSum_, binI)
    {
        // Skip bins not allocated in sparse-bin mode (also have time==0).
        if (!currentCycleBinSum_.set(binI)) continue;

        // Restart-straddled bin: a restart resumed at an instant strictly
        // inside this bin's phase window, so the cycle covered only the part
        // of the window lying after that instant. Its mean is a partial-window
        // mean, biased toward whichever phases happen to follow the restart,
        // and folding it at full weight would carry that bias into the phase
        // average -- while looking entirely healthy on disk.
        //
        // Drop it exactly as a bin that received no data at all is dropped:
        // binCounts_ is not incremented, so no cycle is ever credited without
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
            lastCycleBinTime_[binI] = currentCycleBinTime_[binI];
            if (currentCycleBinTime_[binI] > 0)
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
                                currentCycleBinSum_[binI].mesh().time().name(),
                                currentCycleBinSum_[binI].mesh(),
                                IOobject::NO_READ,
                                IOobject::NO_WRITE,
                                false
                            ),
                            currentCycleBinSum_[binI]/currentCycleBinTime_[binI]
                        )
                    );
                }
                else
                {
                    lastCycleFields_[binI] =
                        currentCycleBinSum_[binI]/currentCycleBinTime_[binI];
                }
            }
        }

        if (currentCycleBinTime_[binI] > 0)
        {
            const FieldType cycleMean =
                currentCycleBinSum_[binI] / currentCycleBinTime_[binI];

            const SqrFieldType cycleMeanSqr =
                currentCycleBinSqrSum_[binI] / currentCycleBinTime_[binI];

            binCounts_[binI] += 1.0;
            const scalar count = binCounts_[binI];

            if (count < 1.5)
            {
                meanFields_[binI] = cycleMean;
                meanSqrFields_[binI] = cycleMeanSqr;
                M2Fields_[binI] = dimensionedScalar(M2Fields_[binI].dimensions(), 0);
            }
            else
            {
                // Welford's online algorithm for cycle-mean variance:
                // M2 accumulates sum of squared deviations |x - <x>|^2
                // of the BIN-AVERAGED value across cycles.
                const FieldType delta(cycleMean - meanFields_[binI]);

                meanFields_[binI] =
                    (meanFields_[binI] * (count - 1.0) + cycleMean) / count;

                const FieldType delta2(cycleMean - meanFields_[binI]);

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
                M2Fields_[binI] += bbpaInner(delta, delta2);

                // Running cross-cycle mean of <phi*phi>_bin. Combined
                // with meanFields_ this recovers Reynolds stress
                // R_ij = <U_i U_j> - <U_i><U_j> and TKE = 0.5 trace(R)
                // offline from the written fields.
                meanSqrFields_[binI] =
                    (meanSqrFields_[binI] * (count - 1.0) + cycleMeanSqr) / count;
            }
        }
    }
}


template<class Type>
void Foam::functionObjects::BBPA::binItem<Type>::write
(
    const Time& time,
    const label nBins,
    const label currentCycle,
    const label currentBin,
    const scalar cycleStartTime,
    const scalar binDeltaT,
    const word& paSuffix,
    const bool emitBin,
    const bool emitPA
) const
{
    // Find a mesh handle from the first ALLOCATED bin (in sparse-bin mode
    // bin 0 may be unallocated).
    const fvMesh* meshPtr = nullptr;
    forAll(meanFields_, b)
    {
        if (meanFields_.set(b)) { meshPtr = &meanFields_[b].mesh(); break; }
    }
    if (meshPtr == nullptr) return;   // nothing allocated -> nothing to write
    const fvMesh& mesh = *meshPtr;
    label binsWritten = 0;

    forAll(meanFields_, binI)
    {
        // Skip bins not allocated in sparse-bin mode.
        if (!meanFields_.set(binI)) continue;

        const bool hasCrossData = binCounts_[binI] > 0;
        const bool hasCurrentData = currentCycleBinTime_[binI] > 0;

        if (hasCrossData || hasCurrentData)
        {
            // Phase-aligned time for this bin. When binDeltaT equals the
            // solver's writeInterval this IS one of the solver's own write
            // directories, so the bin lands next to the instantaneous field
            // rather than creating a directory of its own.
            const scalar binTime = cycleStartTime + binI * binDeltaT;
            const word binTimeName = Time::timeName(binTime);

            // <field>_Bin -- THIS cycle's bin value, a distinct quantity from
            // the cross-cycle mean below. Writing the two separately is what
            // removes any reason to blend the in-progress cycle into the mean.
            // It is also what makes the output restart-safe at write
            // granularity: every write directory carries a complete bin.
            if
            (
                emitBin
             && lastCycleFields_.size()
             && lastCycleFields_.set(binI)
             && lastCycleBinTime_[binI] > 0
            )
            {
                FieldType binField
                (
                    IOobject
                    (
                        fieldName_ + "_Bin",
                        binTimeName,
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    lastCycleFields_[binI]
                );

                writePAField(binField, mesh, binTimeName);
            }

            if (!emitPA)
            {
                binsWritten++;
                continue;
            }

            // Build the PA field values without triggering tmp<> expressions
            // on registered fields. We compute in a copy.
            FieldType paField
            (
                IOobject
                (
                    fieldName_ + paSuffix,
                    binTimeName,
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                meanFields_[binI]
            );

            // paField was constructed FROM meanFields_[binI], i.e. the mean over
            // the binCounts_[binI] COMPLETE cycles folded so far. That is what
            // <name>_PA must contain.
            //
            // It must NOT be blended with the in-progress cycle. The old code did
            //     (meanFields_*nOld + currentCycleBinSum_/currentCycleBinTime_)
            //         / (nOld + 1)
            // which is wrong twice over:
            //  (1) the in-progress bin is PARTIAL -- with writeMode perCycle the
            //      write fires on the cycle transition, when the new cycle has
            //      accumulated exactly ONE timestep into bin 0 -- yet it is
            //      weighted as a whole cycle. Measured on tutorials/bbpaVerify:
            //      bin 0 came out 0.632 against a true 0.701, i.e. (mean*2 + U_inst)/3,
            //      while every other bin (just reset, so no current data) was exact.
            //  (2) these files ARE the restart checkpoint: restoreBinItemFromDisk()
            //      reads meanFields_ back out of <name>_PA. A blended value is
            //      therefore restored AS the folded mean and compounds at every
            //      restart.
            //
            // Only when no complete cycle exists yet is there nothing to report but
            // the running partial. binCounts_[binI] is 0 in that case, and the
            // restore skips bins with count 0, so it can never be read back as a mean.
            //
            // Assign the WHOLE field, not primitiveFieldRef(): an internal-only
            // assignment leaves the boundary at whatever the accumulator was
            // constructed with (zero), silently dropping the boundary values
            // of the partial bin. See the _PA_UU note below.
            if (!hasCrossData && hasCurrentData)
            {
                const scalar inv = 1.0 / currentCycleBinTime_[binI];
                paField = currentCycleBinSum_[binI] * inv;
            }

            writePAField(paField, mesh, binTimeName);

            binsWritten++;

            // Welford cycle-variance (cycle-to-cycle variability of bin mean)
            if (binCounts_[binI] > 1)
            {
                volScalarField paVariance
                (
                    IOobject
                    (
                        fieldName_ + paSuffix + "Variance",
                        binTimeName,
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    M2Fields_[binI] / (binCounts_[binI] - 1.0)
                );

                writePAField(paVariance, mesh, binTimeName);
            }

            // Second-moment output: <phi*phi>_bin as a phase-aligned field.
            // For vector U this is the 6-component symmetric tensor
            // <U_i U_j> from which the Reynolds stress R_ij =
            // <U_i U_j> - <U_i><U_j> and TKE = 0.5*trace(R) are computed
            // offline. For scalar p this reduces to <p*p>.
            // Fires whenever we have data (finalised or in-progress bin).
            if (hasCrossData || hasCurrentData)
            {
                // Construct FROM meanSqrFields_[binI], exactly as paField is
                // constructed from meanFields_[binI] above. This is not a
                // stylistic change: the previous form built a Zero field and
                // then assigned only primitiveFieldRef(), so the BOUNDARY of
                // <name>_PA_UU was never written and every patch came out
                // "value uniform (0 0 0 0 0 0)".
                //
                // The accumulation itself is fine -- operator+= on a
                // GeometricField updates the boundary too -- so the correct
                // boundary existed in memory and was simply dropped on write.
                //
                // For wallShearStress this made the field entirely useless
                // rather than merely incomplete: WSS has a uniform (0 0 0)
                // internalField, so ALL of its data is on the wall patch.
                // wallShearStress_PA_UU was therefore zero everywhere, and any
                // WSS fluctuation intensity or Reynolds stress derived from it
                // would be silently zero rather than obviously broken.
                // Verified by tutorials/bbpaVerifyLES, which compares boundary
                // patches and not just the internalField.
                SqrFieldType paUU
                (
                    IOobject
                    (
                        fieldName_ + paSuffix + "_UU",
                        binTimeName,
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    meanSqrFields_[binI]
                );

                // The raw second moment over the COMPLETE cycles only -- as for
                // _PA above, and for the same two reasons: blending the
                // in-progress cycle weights a partial bin (one timestep, on a
                // perCycle write) as a whole cycle, and this file is the
                // checkpoint that restoreBinItemFromDisk() reads meanSqrFields_
                // back out of, so a blended value would be restored AS the
                // second moment and compound every restart.
                //
                // This matters more than _PA: R_ij^bin = <U_iU_j> - <U_i><U_j>
                // and k_bin = 0.5*tr(R_ij^bin) are reconstructed FROM this
                // field, so a contaminated bin propagates straight into the
                // reported Reynolds stress and resolved energy.
                //
                // Only when no complete cycle exists yet is there nothing to
                // report but the running partial.
                if (!hasCrossData && hasCurrentData)
                {
                    const scalar invT = 1.0 / currentCycleBinTime_[binI];
                    paUU = currentCycleBinSqrSum_[binI] * invT;
                }

                writePAField(paUU, mesh, binTimeName);
            }
        }
    }

    Info<< "    " << fieldName_ << " phase-averaged output:" << nl
        << "        Cycles accumulated: " << nCycles_ << nl
        << "        Bins written: " << binsWritten << "/" << nBins << nl
        << "        Time range: " << Time::timeName(cycleStartTime)
        << " to " << Time::timeName(cycleStartTime + (nBins-1)*binDeltaT)
        << endl;
}


template<class Type>
void Foam::functionObjects::BBPA::binItem<Type>::writeCompanion
(
    const Time& time,
    const label currentBinIndex
) const
{
    // Companion write: one phase (the current bin) into the solver's
    // current time directory. Intended to be paired with the solver's
    // own instantaneous field write at the same time, giving the
    // reader both U and U_PA in the same directory.
    const label binI = currentBinIndex;
    if (binI < 0 || binI >= meanFields_.size()) return;

    // Sparse-bin mode: nothing to write if the current bin is not tracked.
    if (!meanFields_.set(binI)) return;

    const bool hasCrossData = binCounts_[binI] > 0;
    const bool hasCurrentData = currentCycleBinTime_[binI] > 0;
    if (!hasCrossData && !hasCurrentData) return;

    const fvMesh& mesh = meanFields_[binI].mesh();
    const word solverTimeName = time.timeName(time.value());

    // Build PA field for this bin, labelled with solver's current
    // time (not the phase-aligned time).
    FieldType paField
    (
        IOobject
        (
            fieldName_ + "_PA",
            solverTimeName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        meanFields_[binI]
    );

    if (hasCrossData && hasCurrentData)
    {
        const scalar inv = 1.0 / currentCycleBinTime_[binI];
        const scalar nOld = binCounts_[binI];
        const scalar total = nOld + 1.0;
        paField.primitiveFieldRef() =
            (
                meanFields_[binI].primitiveField() * nOld
              + currentCycleBinSum_[binI].primitiveField() * inv
            ) / total;
    }
    else if (!hasCrossData)
    {
        const scalar inv = 1.0 / currentCycleBinTime_[binI];
        paField.primitiveFieldRef() =
            currentCycleBinSum_[binI].primitiveField() * inv;
    }

    mkDir(paField.objectPath().path());
    paField.regIOobject::write();
}


// * * * * * * * * * * * * Restore from checkpoint * * * * * * * * * * * * * //
//
// At restart the existing _PA / _PAVariance / _PA_UU files written
// at the most recent solver writeTime carry the cycle-N ensemble
// state.  Per binItem<Type>::write() at the post-finalisation cycle
// boundary (binCounts_[binI] > 0 && currentCycleBinTime_[binI] == 0):
//
//     paField     = meanFields_[binI]                  (no blending)
//     paVariance  = M2Fields_[binI] / (n - 1)          (only if n > 1)
//     paUU        = meanSqrFields_[binI]               (no blending)
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
    forAll(bi.meanFields_, b)
    {
        if (bi.meanFields_.set(b)) { meshPtr = &bi.meanFields_[b].mesh(); break; }
    }
    if (meshPtr == nullptr) return;
    const fvMesh& mesh = *meshPtr;
    label restored = 0;

    forAll(bi.meanFields_, binI)
    {
        // Skip bins not allocated in sparse-bin mode.
        if (!bi.meanFields_.set(binI)) continue;

        const scalar nCyc = bi.binCounts_[binI];
        if (nCyc <= 0) continue;  // bin never received data

        const scalar binTime = cycleStartTime + binI * binDeltaT_;
        const word binTimeName = Time::timeName(binTime);

        // 1. mean: read fieldName_PA into meanFields_[binI].
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
                bi.binCounts_[binI] = 0;
                bi.meanFields_[binI] =
                    dimensioned<Type>(bi.meanFields_[binI].dimensions(), Zero);
                bi.M2Fields_[binI] =
                    dimensionedScalar(bi.M2Fields_[binI].dimensions(), 0);
                bi.meanSqrFields_[binI] =
                    dimensioned<typename binItem<Type>::SqrType>
                    (bi.meanSqrFields_[binI].dimensions(), Zero);
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
            bi.meanFields_[binI] = loaded;
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
                bi.M2Fields_[binI] = loadedVar * (nCyc - 1.0);
            }
        }

        // 3. meanSqr: read fieldName_PA_UU into meanSqrFields_[binI].
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
                bi.meanSqrFields_[binI] = loadedUU;
            }
        }

        restored++;
    }

    Info<< "    " << fieldName << ": restored "
        << restored << "/" << bi.meanFields_.size()
        << " bins from phase-aligned dirs (cycleStart="
        << cycleStartTime << ", nCycles=" << bi.nCycles_ << ")"
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
