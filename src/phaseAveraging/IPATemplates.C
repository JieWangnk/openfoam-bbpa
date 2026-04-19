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

#include "IPA.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::functionObjects::IPA::phaseItem<Type>::phaseItem
(
    const word& name,
    const fvMesh& mesh,
    const label nPhases
)
:
    fieldName_(name),
    nCycles_(0),
    phaseAveragedFields_(nPhases),
    phaseSampleCounts_(nPhases, 0),
    currentCycleSamples_(nPhases),
    currentCycleSampled_(nPhases, false)
{
    // Look up field to get dimensions
    const FieldType& field = mesh.lookupObject<FieldType>(name);
    
    // Initialize phase-averaged fields
    forAll(phaseAveragedFields_, phaseI)
    {
        phaseAveragedFields_.set
        (
            phaseI,
            new FieldType
            (
                IOobject
                (
                    name + "_phaseAvg_" + Foam::name(phaseI),
                    mesh.time().name(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensioned<Type>(field.dimensions(), Zero)
            )
        );
        
        currentCycleSamples_.set
        (
            phaseI,
            new FieldType
            (
                IOobject
                (
                    name + "_currentSample_" + Foam::name(phaseI),
                    mesh.time().name(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensioned<Type>(field.dimensions(), Zero)
            )
        );
    }
}


template<class Type>
void Foam::functionObjects::IPA::phaseItem<Type>::resetCurrentCycle()
{
    // Reset current cycle tracking
    currentCycleSampled_ = false;
    
    // Reset current cycle samples to zero
    forAll(currentCycleSamples_, phaseI)
    {
        currentCycleSamples_[phaseI] = dimensioned<Type>(currentCycleSamples_[phaseI].dimensions(), Zero);
    }
}


template<class Type>
void Foam::functionObjects::IPA::phaseItem<Type>::sampleAtPhase
(
    const label phaseI,
    const FieldType& field
)
{
    // Only sample if we haven't already sampled this phase in this cycle
    if (!currentCycleSampled_[phaseI])
    {
        currentCycleSamples_[phaseI] = field;
        currentCycleSampled_[phaseI] = true;
    }
}


template<class Type>
void Foam::functionObjects::IPA::phaseItem<Type>::finalizeCycle()
{
    nCycles_++;
    
    // Update phase averages for all sampled phases
    forAll(currentCycleSamples_, phaseI)
    {
        if (currentCycleSampled_[phaseI])
        {
            phaseSampleCounts_[phaseI]++;
            
            const FieldType& sampleValue = currentCycleSamples_[phaseI];
            
            // Update phase average: PA_φ,i^(N) = (1/N) Σ φ_inst,n,i
            if (phaseSampleCounts_[phaseI] == 1)
            {
                phaseAveragedFields_[phaseI] = sampleValue;
            }
            else
            {
                // Incremental average update
                const scalar count = scalar(phaseSampleCounts_[phaseI]);
                phaseAveragedFields_[phaseI] = 
                    (phaseAveragedFields_[phaseI] * (count - 1.0) + sampleValue) / count;
            }
        }
    }
}


template<class Type>
void Foam::functionObjects::IPA::phaseItem<Type>::write(const Time& time) const
{
    label phasesWritten = 0;

    // Write phase-averaged fields, including current cycle samples
    // in the output for best-estimate (non-destructive to accumulators)
    forAll(phaseAveragedFields_, phaseI)
    {
        const bool hasAvgData = phaseSampleCounts_[phaseI] > 0;
        const bool hasCurrentSample = currentCycleSampled_[phaseI];

        if (hasAvgData || hasCurrentSample)
        {
            FieldType phaseField
            (
                IOobject
                (
                    fieldName_ + "_IPA_phase" + Foam::name(phaseI),
                    time.name(),
                    phaseAveragedFields_[phaseI].mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                phaseAveragedFields_[phaseI]
            );

            if (hasAvgData && hasCurrentSample)
            {
                // Blend finalized average with current cycle sample
                const scalar count =
                    scalar(phaseSampleCounts_[phaseI]) + 1.0;
                phaseField =
                (
                    phaseAveragedFields_[phaseI]
                  * scalar(phaseSampleCounts_[phaseI])
                  + currentCycleSamples_[phaseI]
                ) / count;
            }
            else if (!hasAvgData)
            {
                // Only current cycle sample available
                phaseField = currentCycleSamples_[phaseI];
            }
            // else: hasAvgData only, already set from phaseAveragedFields_

            phaseField.write();
            phasesWritten++;
        }
    }
    
    // Print summary
    Info<< "        " << fieldName_ << " (instantaneous): wrote " << phasesWritten 
        << "/" << phaseAveragedFields_.size() << " phase points" << nl
        << "        Total cycles processed: " << nCycles_ << nl;
    
    // Show sampling statistics
    label totalSamples = 0;
    forAll(phaseSampleCounts_, phaseI)
    {
        totalSamples += phaseSampleCounts_[phaseI];
    }
    
    Info<< "        Total samples collected: " << totalSamples << nl
        << "        Average samples per phase: " << (totalSamples / max(1, label(phaseAveragedFields_.size()))) << endl;
}


// Explicit instantiation for common types
template class Foam::functionObjects::IPA::phaseItem<Foam::scalar>;
template class Foam::functionObjects::IPA::phaseItem<Foam::vector>;
template class Foam::functionObjects::IPA::phaseItem<Foam::symmTensor>;
template class Foam::functionObjects::IPA::phaseItem<Foam::tensor>;


// ************************************************************************* //