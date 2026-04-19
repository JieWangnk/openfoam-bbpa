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
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(IPA, 0);
    addToRunTimeSelectionTable(functionObject, IPA, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::functionObjects::IPA::currentCycleNumber() const
{
    const scalar t = time_.value() - startTime_;
    return label(t / period_);
}


Foam::label Foam::functionObjects::IPA::currentPhaseIndex() const
{
    const scalar t = time_.value() - startTime_;
    const scalar phase = std::fmod(t, period_);
    return label(phase / phaseInterval_);
}


bool Foam::functionObjects::IPA::isAtPhasePoint(label& phaseIndex) const
{
    const scalar t = time_.value() - startTime_;
    const scalar phase = std::fmod(t, period_);
    
    // Check each phase point to see if we're within tolerance
    for (label i = 0; i < nPhases_; i++)
    {
        const scalar targetPhase = scalar(i) * phaseInterval_;
        const scalar phaseDiff = mag(phase - targetPhase);
        
        // Also check wrapped difference (near period boundary)
        const scalar wrappedDiff = min(phaseDiff, period_ - phaseDiff);
        
        if (wrappedDiff <= phaseTolerance_)
        {
            phaseIndex = i;
            return true;
        }
    }
    
    return false;
}


void Foam::functionObjects::IPA::resetPhaseItemsForNewCycle()
{
    forAllIter(HashTable<autoPtr<phaseItem<scalar>>>, scalarPhaseItems_, iter)
    {
        iter()->resetCurrentCycle();
    }
    forAllIter(HashTable<autoPtr<phaseItem<vector>>>, vectorPhaseItems_, iter)
    {
        iter()->resetCurrentCycle();
    }
    forAllIter(HashTable<autoPtr<phaseItem<symmTensor>>>, symmTensorPhaseItems_, iter)
    {
        iter()->resetCurrentCycle();
    }
    forAllIter(HashTable<autoPtr<phaseItem<tensor>>>, tensorPhaseItems_, iter)
    {
        iter()->resetCurrentCycle();
    }
}


void Foam::functionObjects::IPA::finalizePhaseItemsForCycle()
{
    forAllIter(HashTable<autoPtr<phaseItem<scalar>>>, scalarPhaseItems_, iter)
    {
        iter()->finalizeCycle();
    }
    forAllIter(HashTable<autoPtr<phaseItem<vector>>>, vectorPhaseItems_, iter)
    {
        iter()->finalizeCycle();
    }
    forAllIter(HashTable<autoPtr<phaseItem<symmTensor>>>, symmTensorPhaseItems_, iter)
    {
        iter()->finalizeCycle();
    }
    forAllIter(HashTable<autoPtr<phaseItem<tensor>>>, tensorPhaseItems_, iter)
    {
        iter()->finalizeCycle();
    }
}


void Foam::functionObjects::IPA::sampleFieldsAtPhase(const label phaseI)
{
    // Sample scalar fields
    forAllIter(HashTable<autoPtr<phaseItem<scalar>>>, scalarPhaseItems_, iter)
    {
        const word& fieldName = iter.key();
        const volScalarField& field = mesh_.lookupObject<volScalarField>(fieldName);
        iter()->sampleAtPhase(phaseI, field);
    }
    
    // Sample vector fields
    forAllIter(HashTable<autoPtr<phaseItem<vector>>>, vectorPhaseItems_, iter)
    {
        const word& fieldName = iter.key();
        const volVectorField& field = mesh_.lookupObject<volVectorField>(fieldName);
        iter()->sampleAtPhase(phaseI, field);
    }
    
    // Sample symmTensor fields
    forAllIter(HashTable<autoPtr<phaseItem<symmTensor>>>, symmTensorPhaseItems_, iter)
    {
        const word& fieldName = iter.key();
        const volSymmTensorField& field = mesh_.lookupObject<volSymmTensorField>(fieldName);
        iter()->sampleAtPhase(phaseI, field);
    }
    
    // Sample tensor fields
    forAllIter(HashTable<autoPtr<phaseItem<tensor>>>, tensorPhaseItems_, iter)
    {
        const word& fieldName = iter.key();
        const volTensorField& field = mesh_.lookupObject<volTensorField>(fieldName);
        iter()->sampleAtPhase(phaseI, field);
    }
}


void Foam::functionObjects::IPA::initializePhaseItems()
{
    forAll(fieldNames_, i)
    {
        const word& fieldName = fieldNames_[i];
        
        if (mesh_.foundObject<volScalarField>(fieldName))
        {
            scalarPhaseItems_.insert
            (
                fieldName,
                autoPtr<phaseItem<scalar>>(new phaseItem<scalar>(fieldName, mesh_, nPhases_))
            );
        }
        else if (mesh_.foundObject<volVectorField>(fieldName))
        {
            vectorPhaseItems_.insert
            (
                fieldName,
                autoPtr<phaseItem<vector>>(new phaseItem<vector>(fieldName, mesh_, nPhases_))
            );
        }
        else if (mesh_.foundObject<volSymmTensorField>(fieldName))
        {
            symmTensorPhaseItems_.insert
            (
                fieldName,
                autoPtr<phaseItem<symmTensor>>(new phaseItem<symmTensor>(fieldName, mesh_, nPhases_))
            );
        }
        else if (mesh_.foundObject<volTensorField>(fieldName))
        {
            tensorPhaseItems_.insert
            (
                fieldName,
                autoPtr<phaseItem<tensor>>(new phaseItem<tensor>(fieldName, mesh_, nPhases_))
            );
        }
        else
        {
            FatalErrorInFunction
                << "Field " << fieldName << " not found or unsupported type"
                << exit(FatalError);
        }
    }
}


void Foam::functionObjects::IPA::writePhaseAveragedFields()
{
    // Write scalar fields
    forAllConstIter(HashTable<autoPtr<phaseItem<scalar>>>, scalarPhaseItems_, iter)
    {
        iter()->write(time_);
    }
    
    // Write vector fields
    forAllConstIter(HashTable<autoPtr<phaseItem<vector>>>, vectorPhaseItems_, iter)
    {
        iter()->write(time_);
    }
    
    // Write symmTensor fields
    forAllConstIter(HashTable<autoPtr<phaseItem<symmTensor>>>, symmTensorPhaseItems_, iter)
    {
        iter()->write(time_);
    }
    
    // Write tensor fields
    forAllConstIter(HashTable<autoPtr<phaseItem<tensor>>>, tensorPhaseItems_, iter)
    {
        iter()->write(time_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::IPA::IPA
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldNames_(),
    period_(0),
    nPhases_(50),
    startCycle_(0),
    startTime_(runTime.value()),
    phaseInterval_(0),
    currentCycle_(0),
    previousCycle_(-1),
    phaseTolerance_(0.01),
    scalarPhaseItems_(),
    vectorPhaseItems_(),
    symmTensorPhaseItems_(),
    tensorPhaseItems_()
{
    read(dict);
    initializePhaseItems();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::IPA::~IPA()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::IPA::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    
    dict.lookup("fields") >> fieldNames_;
    dict.lookup("period") >> period_;
    nPhases_ = dict.lookupOrDefault("nPhases", 50);
    startCycle_ = dict.lookupOrDefault("startCycle", 0);
    phaseTolerance_ = dict.lookupOrDefault("phaseTolerance", 0.01);
    
    phaseInterval_ = period_ / scalar(nPhases_);
    
    Info<< type() << " " << name() << ":" << nl
        << "    Period: " << period_ << " s" << nl
        << "    Number of phase points: " << nPhases_ << nl
        << "    Phase interval: " << phaseInterval_ << " s" << nl
        << "    Phase tolerance: " << phaseTolerance_ << " s" << nl
        << "    Start cycle: " << startCycle_ << nl
        << "    Fields: " << fieldNames_ << endl;
    
    return true;
}


bool Foam::functionObjects::IPA::execute()
{
    currentCycle_ = currentCycleNumber();
    
    // Check if we should start averaging yet
    if (currentCycle_ < startCycle_)
    {
        return true;
    }
    
    // Detect cycle transition
    if (currentCycle_ != previousCycle_)
    {
        if (previousCycle_ >= startCycle_)
        {
            // Finalize the previous cycle (only if it was being averaged)
            finalizePhaseItemsForCycle();
        }
        
        // Reset phase items for new cycle
        resetPhaseItemsForNewCycle();
        
        Info<< "Starting cycle " << currentCycle_;
        if (currentCycle_ == startCycle_)
        {
            Info<< " (instantaneous phase averaging begins)";
        }
        Info<< endl;
        
        previousCycle_ = currentCycle_;
    }
    
    // Check if we're at a phase point
    label phaseIndex;
    if (isAtPhasePoint(phaseIndex))
    {
        Info<< "    Sampling phase " << phaseIndex << "/" << (nPhases_-1) 
            << " at time " << time_.value() << " (cycle " << currentCycle_ << ")" << endl;
        
        sampleFieldsAtPhase(phaseIndex);
    }
    
    return true;
}


bool Foam::functionObjects::IPA::write()
{
    Info<< type() << " " << name() << " write:" << nl;
    
    // Only write if we're in the averaging phase
    if (currentCycle_ >= startCycle_)
    {
        writePhaseAveragedFields();
    }
    else
    {
        Info<< "    Pre-averaging phase - no output written" << endl;
    }
    
    return true;
}


// ************************************************************************* //