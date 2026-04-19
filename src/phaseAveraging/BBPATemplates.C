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

template<class Type>
Foam::functionObjects::BBPA::binItem<Type>::binItem
(
    const word& name,
    const fvMesh& mesh,
    const label nBins
)
:
    fieldName_(name),
    nCycles_(0),
    binCounts_(nBins, 0.0),
    meanFields_(nBins),
    M2Fields_(nBins),
    currentCycleBinTime_(nBins, 0.0),
    currentCycleBinSum_(nBins),
    meanUUFields_(nBins),
    currentCycleBinUUSum_(nBins)
{
    // Look up field to get dimensions
    const FieldType& field = mesh.lookupObject<FieldType>(name);

    forAll(meanFields_, binI)
    {
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

        meanUUFields_.set
        (
            binI,
            new volScalarField
            (
                IOobject
                (
                    name + "_meanUU_bin" + Foam::name(binI),
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

        currentCycleBinUUSum_.set
        (
            binI,
            new volScalarField
            (
                IOobject
                (
                    name + "_currentBinUU_bin" + Foam::name(binI),
                    mesh.time().name(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar(sqr(field.dimensions()), Zero)
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
        currentCycleBinSum_[binI] =
            dimensioned<Type>(currentCycleBinSum_[binI].dimensions(), Zero);
        currentCycleBinUUSum_[binI] =
            dimensionedScalar(currentCycleBinUUSum_[binI].dimensions(), Zero);
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
    // Time-weighted accumulation: phi(t) * dt, and |u|^2 * dt for TKE
    currentCycleBinTime_[binI] += dt;
    currentCycleBinSum_[binI] += dt * field;
    currentCycleBinUUSum_[binI] += dt * magSqr(field);
}


template<class Type>
void Foam::functionObjects::BBPA::binItem<Type>::finalizeCycle()
{
    nCycles_++;

    // Update overall statistics for each bin that has data
    forAll(currentCycleBinSum_, binI)
    {
        if (currentCycleBinTime_[binI] > 0)
        {
            const FieldType cycleMean =
                currentCycleBinSum_[binI] / currentCycleBinTime_[binI];

            const volScalarField cycleMeanUU =
                currentCycleBinUUSum_[binI] / currentCycleBinTime_[binI];

            binCounts_[binI] += 1.0;
            const scalar count = binCounts_[binI];

            if (count < 1.5)
            {
                meanFields_[binI] = cycleMean;
                meanUUFields_[binI] = cycleMeanUU;
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

                // Running cross-cycle mean of <|u|^2>_bin: time-weighted
                // intra-bin second moment, then ensemble-averaged over
                // cycles. Combined with meanFields_ this yields the true
                // LES resolved TKE per phase bin.
                meanUUFields_[binI] =
                    (meanUUFields_[binI] * (count - 1.0) + cycleMeanUU) / count;
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
    const fvMesh& mesh = meanFields_[0].mesh();
    label binsWritten = 0;

    forAll(meanFields_, binI)
    {
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

            mkDir(paField.objectPath().path());
            paField.regIOobject::write();

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

                mkDir(paVariance.objectPath().path());
                paVariance.regIOobject::write();
            }

            // True LES resolved TKE per phase bin:
            //   TKE = 0.5 * ( <|u|^2>_bin - |<u>_bin|^2 )
            // Fires whenever we have data, either from a finalised cycle
            // (hasCrossData) or an in-progress one (hasCurrentData). The
            // partial-cycle fall-through avoids a silent miss when the
            // simulation endTime lands exactly on a cycle boundary (the
            // cycle-transition finalisation never runs in that case).
            if (hasCrossData || hasCurrentData)
            {
                volScalarField paTKE
                (
                    IOobject
                    (
                        fieldName_ + "_PATKE",
                        binTimeName,
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh,
                    dimensionedScalar(meanUUFields_[binI].dimensions(), Zero)
                );

                if (hasCrossData && hasCurrentData)
                {
                    // Blend finalised cross-cycle stats with current cycle,
                    // matching the weighting used for _PA above.
                    const scalar invT = 1.0 / currentCycleBinTime_[binI];
                    const scalar nOld = binCounts_[binI];
                    const scalar total = nOld + 1.0;
                    const scalarField meanUU =
                        (
                            meanUUFields_[binI].primitiveField() * nOld
                          + currentCycleBinUUSum_[binI].primitiveField() * invT
                        ) / total;
                    paTKE.primitiveFieldRef() =
                        0.5 * (meanUU - magSqr(paField.primitiveField()));
                }
                else if (hasCrossData)
                {
                    paTKE.primitiveFieldRef() =
                        0.5 *
                        (
                            meanUUFields_[binI].primitiveField()
                          - magSqr(meanFields_[binI].primitiveField())
                        );
                }
                else
                {
                    // Partial-cycle only: single-realisation TKE.
                    // Statistically noisier than N>=2 but mathematically
                    // valid (identity holds at any N).
                    const scalar invT = 1.0 / currentCycleBinTime_[binI];
                    paTKE.primitiveFieldRef() =
                        0.5 *
                        (
                            currentCycleBinUUSum_[binI].primitiveField() * invT
                          - magSqr(paField.primitiveField())
                        );
                }

                mkDir(paTKE.objectPath().path());
                paTKE.regIOobject::write();
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


// Explicit instantiation for common types
template struct Foam::functionObjects::BBPA::binItem<Foam::scalar>;
template struct Foam::functionObjects::BBPA::binItem<Foam::vector>;
template struct Foam::functionObjects::BBPA::binItem<Foam::symmTensor>;
template struct Foam::functionObjects::BBPA::binItem<Foam::tensor>;


// ************************************************************************* //
