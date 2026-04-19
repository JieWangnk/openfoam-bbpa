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
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "timeIOdictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(BBPA, 0);
    addToRunTimeSelectionTable(functionObject, BBPA, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::functionObjects::BBPA::phaseBin() const
{
    const scalar t = time_.value() - startTime_;
    const scalar phase = std::fmod(t, period_);
    return min(label(phase/binDeltaT_), nBins_ - 1);
}


Foam::label Foam::functionObjects::BBPA::currentCycleNumber() const
{
    const scalar t = time_.value() - startTime_;
    return label(t/period_);
}


void Foam::functionObjects::BBPA::initializeBinItems()
{
    // Idempotent: creates binItems for fields available now that don't
    // already have one. Fields registered later (e.g. wallShearStress,
    // which is populated by its own function object) are picked up on a
    // later call from execute().
    forAll(fieldNames_, i)
    {
        const word& fieldName = fieldNames_[i];

        if
        (
            scalarBinItems_.found(fieldName)
         || vectorBinItems_.found(fieldName)
         || symmTensorBinItems_.found(fieldName)
         || tensorBinItems_.found(fieldName)
        )
        {
            continue;
        }

        if (mesh_.foundObject<volScalarField>(fieldName))
        {
            scalarBinItems_.insert
            (
                fieldName,
                autoPtr<binItem<scalar>>
                (
                    new binItem<scalar>(fieldName, mesh_, nBins_)
                )
            );
        }
        else if (mesh_.foundObject<volVectorField>(fieldName))
        {
            vectorBinItems_.insert
            (
                fieldName,
                autoPtr<binItem<vector>>
                (
                    new binItem<vector>(fieldName, mesh_, nBins_)
                )
            );
        }
        else if (mesh_.foundObject<volSymmTensorField>(fieldName))
        {
            symmTensorBinItems_.insert
            (
                fieldName,
                autoPtr<binItem<symmTensor>>
                (
                    new binItem<symmTensor>(fieldName, mesh_, nBins_)
                )
            );
        }
        else if (mesh_.foundObject<volTensorField>(fieldName))
        {
            tensorBinItems_.insert
            (
                fieldName,
                autoPtr<binItem<tensor>>
                (
                    new binItem<tensor>(fieldName, mesh_, nBins_)
                )
            );
        }
        // Silently skip fields not yet registered; they'll be retried
        // on the next execute() call.
    }
}


void Foam::functionObjects::BBPA::resetAllBinsForNewCycle()
{
    forAllIter(HashTable<autoPtr<binItem<scalar>>>, scalarBinItems_, iter)
    {
        iter()->resetCurrentCycle();
    }
    forAllIter(HashTable<autoPtr<binItem<vector>>>, vectorBinItems_, iter)
    {
        iter()->resetCurrentCycle();
    }
    forAllIter
    (
        HashTable<autoPtr<binItem<symmTensor>>>,
        symmTensorBinItems_,
        iter
    )
    {
        iter()->resetCurrentCycle();
    }
    forAllIter(HashTable<autoPtr<binItem<tensor>>>, tensorBinItems_, iter)
    {
        iter()->resetCurrentCycle();
    }
}


void Foam::functionObjects::BBPA::finalizeAllForCycle()
{
    forAllIter(HashTable<autoPtr<binItem<scalar>>>, scalarBinItems_, iter)
    {
        iter()->finalizeCycle();
    }
    forAllIter(HashTable<autoPtr<binItem<vector>>>, vectorBinItems_, iter)
    {
        iter()->finalizeCycle();
    }
    forAllIter
    (
        HashTable<autoPtr<binItem<symmTensor>>>,
        symmTensorBinItems_,
        iter
    )
    {
        iter()->finalizeCycle();
    }
    forAllIter(HashTable<autoPtr<binItem<tensor>>>, tensorBinItems_, iter)
    {
        iter()->finalizeCycle();
    }
}


void Foam::functionObjects::BBPA::accumulateAllFields
(
    const label binI,
    const scalar dt
)
{
    forAllIter(HashTable<autoPtr<binItem<scalar>>>, scalarBinItems_, iter)
    {
        const word& fieldName = iter.key();
        const volScalarField& field =
            mesh_.lookupObject<volScalarField>(fieldName);
        iter()->accumulate(binI, field, dt);
    }

    forAllIter(HashTable<autoPtr<binItem<vector>>>, vectorBinItems_, iter)
    {
        const word& fieldName = iter.key();
        const volVectorField& field =
            mesh_.lookupObject<volVectorField>(fieldName);
        iter()->accumulate(binI, field, dt);
    }

    forAllIter
    (
        HashTable<autoPtr<binItem<symmTensor>>>,
        symmTensorBinItems_,
        iter
    )
    {
        const word& fieldName = iter.key();
        const volSymmTensorField& field =
            mesh_.lookupObject<volSymmTensorField>(fieldName);
        iter()->accumulate(binI, field, dt);
    }

    forAllIter(HashTable<autoPtr<binItem<tensor>>>, tensorBinItems_, iter)
    {
        const word& fieldName = iter.key();
        const volTensorField& field =
            mesh_.lookupObject<volTensorField>(fieldName);
        iter()->accumulate(binI, field, dt);
    }
}


void Foam::functionObjects::BBPA::writeAllAveragedFields()
{
    // Compute the start time of the current (or most recent) cycle
    // so bins are written to phase-aligned time directories.
    const scalar cycleStart =
        startTime_ + currentCycle_ * period_;

    forAllConstIter
    (
        HashTable<autoPtr<binItem<scalar>>>,
        scalarBinItems_,
        iter
    )
    {
        iter()->write
        (
            time_, nBins_, currentCycle_, currentBin_,
            cycleStart, binDeltaT_
        );
    }

    forAllConstIter
    (
        HashTable<autoPtr<binItem<vector>>>,
        vectorBinItems_,
        iter
    )
    {
        iter()->write
        (
            time_, nBins_, currentCycle_, currentBin_,
            cycleStart, binDeltaT_
        );
    }

    forAllConstIter
    (
        HashTable<autoPtr<binItem<symmTensor>>>,
        symmTensorBinItems_,
        iter
    )
    {
        iter()->write
        (
            time_, nBins_, currentCycle_, currentBin_,
            cycleStart, binDeltaT_
        );
    }

    forAllConstIter
    (
        HashTable<autoPtr<binItem<tensor>>>,
        tensorBinItems_,
        iter
    )
    {
        iter()->write
        (
            time_, nBins_, currentCycle_, currentBin_,
            cycleStart, binDeltaT_
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::BBPA::BBPA
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldNames_(),
    period_(0),
    nBins_(100),
    cycles_(-1),
    startCycle_(0),
    startTime_(0),
    binDeltaT_(0),
    currentCycle_(0),
    previousCycle_(-1),
    currentBin_(-1),
    previousBin_(-1),
    scalarBinItems_(),
    vectorBinItems_(),
    symmTensorBinItems_(),
    tensorBinItems_()
{
    read(dict);
    initializeBinItems();

    // Checkpoint restore: meanFields_ and M2Fields_ are restored via
    // READ_IF_PRESENT in their IOobject constructors. Scalar counters
    // (nCycles, binCounts) are restored from the Properties file if present.
    // TODO: implement Properties file read for scalar counters on restart.
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::BBPA::~BBPA()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::BBPA::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.lookup("fields") >> fieldNames_;
    dict.lookup("period") >> period_;
    nBins_ = dict.lookupOrDefault("nBins", 100);
    cycles_ = dict.lookupOrDefault("cycles", -1);
    startCycle_ = dict.lookupOrDefault("startCycle", 0);

    // Phase origin: explicit override or auto-align to nearest
    // cycle boundary at or before the current time.
    if (dict.found("startTime"))
    {
        startTime_ = dict.lookup<scalar>("startTime");
    }
    else
    {
        startTime_ = std::floor(time_.value() / period_) * period_;
    }

    binDeltaT_ = period_/nBins_;

    Info<< type() << " " << name() << ":" << nl
        << "    Period: " << period_ << " s" << nl
        << "    Number of bins: " << nBins_ << nl
        << "    Bin time step: " << binDeltaT_ << " s" << nl
        << "    Start cycle: " << startCycle_ << nl
        << "    Phase origin: " << startTime_ << " s" << nl
        << "    Number of cycles: "
        << (cycles_ == -1 ? "unlimited" : Foam::name(cycles_)) << nl
        << "    Fields: " << fieldNames_ << endl;

    return true;
}


bool Foam::functionObjects::BBPA::execute()
{
    // Retry initialisation for fields that weren't registered at
    // construction time (e.g. wallShearStress, populated by its own
    // function object that runs before BBPA).
    initializeBinItems();

    const label binI = phaseBin();
    currentCycle_ = currentCycleNumber();

    // Check if we should start averaging yet
    if (currentCycle_ < startCycle_)
    {
        return true;
    }

    // Check if within averaging cycle limit
    if (cycles_ > 0 && currentCycle_ >= (startCycle_ + cycles_))
    {
        return true;
    }

    // Detect cycle transition
    if (currentCycle_ != previousCycle_)
    {
        if (previousCycle_ >= startCycle_)
        {
            // Finalize the previous cycle (only if it was being averaged)
            finalizeAllForCycle();
        }

        // Reset bins for new cycle
        resetAllBinsForNewCycle();

        Info<< "Starting cycle " << currentCycle_;
        if (currentCycle_ == startCycle_)
        {
            Info<< " (phase averaging begins)";
        }
        Info<< endl;

        previousCycle_ = currentCycle_;
    }

    // Update current bin tracking
    currentBin_ = binI;

    // Accumulate all fields for current time step
    const scalar dt = time_.deltaT().value();
    accumulateAllFields(binI, dt);

    previousBin_ = currentBin_;

    return true;
}


bool Foam::functionObjects::BBPA::write()
{
    Info<< type() << " " << name() << " write:" << nl;

    writeAllAveragedFields();

    // Checkpoint: save scalar counters so the run can be restarted.
    // The meanFields_ and M2Fields_ are written by AUTO_WRITE; this
    // dict saves the counters that cannot be reconstructed from the
    // fields alone.
    {
        timeIOdictionary propsDict
        (
            IOobject
            (
                name() + "Properties",
                time_.name(),
                "uniform",
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        propsDict.add("startTime", startTime_);
        propsDict.add("currentCycle", currentCycle_);
        propsDict.add("previousCycle", previousCycle_);

        // Per-field counters
        forAllConstIter
        (
            HashTable<autoPtr<binItem<scalar>>>,
            scalarBinItems_,
            iter
        )
        {
            dictionary fieldDict;
            fieldDict.add("nCycles", iter()->nCycles_);
            fieldDict.add("binCounts", iter()->binCounts_);
            propsDict.add(iter.key(), fieldDict);
        }
        forAllConstIter
        (
            HashTable<autoPtr<binItem<vector>>>,
            vectorBinItems_,
            iter
        )
        {
            dictionary fieldDict;
            fieldDict.add("nCycles", iter()->nCycles_);
            fieldDict.add("binCounts", iter()->binCounts_);
            propsDict.add(iter.key(), fieldDict);
        }
        forAllConstIter
        (
            HashTable<autoPtr<binItem<symmTensor>>>,
            symmTensorBinItems_,
            iter
        )
        {
            dictionary fieldDict;
            fieldDict.add("nCycles", iter()->nCycles_);
            fieldDict.add("binCounts", iter()->binCounts_);
            propsDict.add(iter.key(), fieldDict);
        }
        forAllConstIter
        (
            HashTable<autoPtr<binItem<tensor>>>,
            tensorBinItems_,
            iter
        )
        {
            dictionary fieldDict;
            fieldDict.add("nCycles", iter()->nCycles_);
            fieldDict.add("binCounts", iter()->binCounts_);
            propsDict.add(iter.key(), fieldDict);
        }

        propsDict.regIOobject::write();
    }

    return true;
}


// ************************************************************************* //
