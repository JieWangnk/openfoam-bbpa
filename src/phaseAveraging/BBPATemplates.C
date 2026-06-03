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
    Foam::fileName caseDir =
        mesh.time().rootPath() / mesh.time().caseName();

    Foam::fileName timeDir;
    if (Foam::Pstream::parRun())
    {
        timeDir = caseDir
            / ("processor" + Foam::name(Foam::Pstream::myProcNo()))
            / binTimeName;
    }
    else
    {
        timeDir = caseDir / binTimeName;
    }
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
    const labelList& activeBins
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
    currentCycleBinSqrSum_(nBins)
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
void Foam::functionObjects::BBPA::binItem<Type>::finalizeCycle()
{
    nCycles_++;

    // Update overall statistics for each bin that has data
    forAll(currentCycleBinSum_, binI)
    {
        // Skip bins not allocated in sparse-bin mode (also have time==0).
        if (!currentCycleBinSum_.set(binI)) continue;

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

                M2Fields_[binI] += 0.5*(magSqr(delta) + magSqr(delta2));

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
    const scalar binDeltaT
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
            // Phase-aligned time for this bin
            const scalar binTime = cycleStartTime + binI * binDeltaT;
            const word binTimeName = Time::timeName(binTime);

            // Build the PA field values without triggering tmp<> expressions
            // on registered fields. We compute in a copy.
            FieldType paField
            (
                IOobject
                (
                    fieldName_ + "_PA",
                    binTimeName,
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

            writePAField(paField, mesh, binTimeName);

            binsWritten++;

            // Welford cycle-variance (cycle-to-cycle variability of bin mean)
            if (binCounts_[binI] > 1)
            {
                volScalarField paVariance
                (
                    IOobject
                    (
                        fieldName_ + "_PAVariance",
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
                SqrFieldType paUU
                (
                    IOobject
                    (
                        fieldName_ + "_PA_UU",
                        binTimeName,
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh,
                    dimensioned<SqrType>(meanSqrFields_[binI].dimensions(), Zero)
                );

                if (hasCrossData && hasCurrentData)
                {
                    // Blend finalised cross-cycle stats with current cycle.
                    const scalar invT = 1.0 / currentCycleBinTime_[binI];
                    const scalar nOld = binCounts_[binI];
                    const scalar total = nOld + 1.0;
                    paUU.primitiveFieldRef() =
                        (
                            meanSqrFields_[binI].primitiveField() * nOld
                          + currentCycleBinSqrSum_[binI].primitiveField() * invT
                        ) / total;
                }
                else if (hasCrossData)
                {
                    paUU.primitiveFieldRef() =
                        meanSqrFields_[binI].primitiveField();
                }
                else
                {
                    // Partial-cycle only: single-realisation second moment.
                    const scalar invT = 1.0 / currentCycleBinTime_[binI];
                    paUU.primitiveFieldRef() =
                        currentCycleBinSqrSum_[binI].primitiveField() * invT;
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
    const scalar cycleStartTime
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
                fieldName + "_PA",
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
                WarningInFunction
                    << "Checkpoint mean field missing: " << binTimeName
                    << "/" << fieldName << "_PA -- bin " << binI
                    << " left at cold-start zero." << endl;
                continue;
            }

            FieldType loaded(paIO, mesh);
            bi.meanFields_[binI].primitiveFieldRef() =
                loaded.primitiveField();
        }

        // 2. M2: read fieldName_PAVariance and rescale by (n-1).
        if (nCyc > 1)
        {
            IOobject varIO
            (
                fieldName + "_PAVariance",
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
                bi.M2Fields_[binI].primitiveFieldRef() =
                    loadedVar.primitiveField() * (nCyc - 1.0);
            }
        }

        // 3. meanSqr: read fieldName_PA_UU into meanSqrFields_[binI].
        {
            IOobject uuIO
            (
                fieldName + "_PA_UU",
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
                bi.meanSqrFields_[binI].primitiveFieldRef() =
                    loadedUU.primitiveField();
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
    const scalar
);
template void
Foam::functionObjects::BBPA::restoreBinItemFromDisk<Foam::vector>
(
    binItem<Foam::vector>&,
    const word&,
    const scalar
);


// ************************************************************************* //
