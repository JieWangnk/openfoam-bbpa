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
#include "IFstream.H"

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


bool Foam::functionObjects::BBPA::initializeBinItems()
{
    // Idempotent: creates binItems for fields available now that don't
    // already have one. Fields registered later (e.g. wallShearStress,
    // which is populated by its own function object) are picked up on a
    // later call from execute().  Returns true if anything was created,
    // which signals execute() to re-apply the checkpoint to the newcomer
    // (otherwise it would silently cold-start on every restart).
    bool created = false;

    forAll(fieldNames_, i)
    {
        const word& fieldName = fieldNames_[i];

        if
        (
            scalarBinItems_.found(fieldName)
         || vectorBinItems_.found(fieldName)
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
                    new binItem<scalar>
                    (
                        fieldName, mesh_, nBins_, binsOfInterest_
                    )
                )
            );
            created = true;
        }
        else if (mesh_.foundObject<volVectorField>(fieldName))
        {
            vectorBinItems_.insert
            (
                fieldName,
                autoPtr<binItem<vector>>
                (
                    new binItem<vector>
                    (
                        fieldName, mesh_, nBins_, binsOfInterest_
                    )
                )
            );
            created = true;
        }
        // Other field types (symmTensor, tensor) are not supported: the
        // second-moment accumulator phi*phi isn't defined in OpenFOAM for
        // generic tensor fields. Silently skip unknown or unsupported
        // fields; they'll be retried on the next execute() call in case
        // a vector/scalar version appears.
    }

    return created;
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
}


void Foam::functionObjects::BBPA::writeCompanionAllFields()
{
    // Companion mode: write only the current-phase bin for every
    // tracked field, into the solver's current time directory.
    forAllConstIter
    (
        HashTable<autoPtr<binItem<scalar>>>, scalarBinItems_, iter
    )
    {
        iter()->writeCompanion(time_, currentBin_);
    }
    forAllConstIter
    (
        HashTable<autoPtr<binItem<vector>>>, vectorBinItems_, iter
    )
    {
        iter()->writeCompanion(time_, currentBin_);
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
    lastFullWriteCycle_(-1),
    checkpointSearched_(false),
    checkpointPath_(fileName::null),
    countersRestored_(false),
    checkpointCycleStart_(0),
    restoredFields_(),
    scalarBinItems_(),
    vectorBinItems_()
{
    read(dict);
    initializeBinItems();

    // Checkpoint restore.  Reads <name>Properties and the
    // phase-aligned _PA / _PAVariance / _PA_UU fields from the most
    // recent solver writeTime to reconstruct the internal accumulator
    // state.  No-op on cold start (no Properties file found).
    readCheckpoint();
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

    // Optional sparse-bin mode: only track the listed bin indices.
    // Default empty list = track all nBins_ bins (original behaviour).
    binsOfInterest_ =
        dict.lookupOrDefault<labelList>("binsOfInterest", labelList());

    // Validate any supplied bin indices against [0, nBins_).
    forAll(binsOfInterest_, k)
    {
        const label b = binsOfInterest_[k];
        if (b < 0 || b >= nBins_)
        {
            FatalErrorInFunction
                << "binsOfInterest entry " << b
                << " is out of range [0," << nBins_ << ")"
                << exit(FatalError);
        }
    }

    // Output mode: "phaseAlignedDirs" (default, Strategy B),
    // "companion" (one bin into solver's time dir), or "both".
    writeMode_ = dict.lookupOrDefault<word>
    (
        "writeMode", word("phaseAlignedDirs")
    );

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
        << "    Fields: " << fieldNames_ << nl
        << "    Bin mode: "
        << (binsOfInterest_.empty()
                ? word("all " + Foam::name(nBins_) + " bins")
                : word("sparse, "
                     + Foam::name(binsOfInterest_.size()) + " bins"))
        << endl;
    if (!binsOfInterest_.empty())
    {
        Info<< "    binsOfInterest: " << binsOfInterest_ << endl;
    }

    return true;
}


bool Foam::functionObjects::BBPA::execute()
{
    // Retry initialisation for fields that weren't registered at
    // construction time (e.g. wallShearStress, populated by its own
    // function object that runs before BBPA).  If one appears now, apply
    // the checkpoint to it immediately -- otherwise it would cold-start on
    // every restart while U (a solver field, present at construction)
    // restored correctly.  Done BEFORE currentCycle_ is advanced below, so
    // the restore sees the checkpoint's cycle state.
    if (initializeBinItems())
    {
        readCheckpoint();
    }

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
    Info<< type() << " " << name() << " write (mode="
        << writeMode_ << "):" << nl;

    const bool isLastWrite =
        (time_.value() + 0.5*time_.deltaT().value() >= time_.endTime().value());

    if (writeMode_ == "companion")
    {
        writeCompanionAllFields();
    }
    else if (writeMode_ == "both")
    {
        writeCompanionAllFields();
        if (isLastWrite)
        {
            writeAllAveragedFields();
        }
    }
    else if (writeMode_ == "perCycle")
    {
        // Full phase-aligned dump, but only once per completed cycle (the
        // first write after a cycle boundary) instead of every solver write.
        // Cuts the bin-output I/O by ~period/writeInterval while staying a
        // full, restart-safe checkpoint at each cycle (restart granularity =
        // one cycle). lastFullWriteCycle_ starts at -1 and is not persisted,
        // so the first write after any (re)start re-establishes the checkpoint.
        if (currentCycle_ > lastFullWriteCycle_)
        {
            writeAllAveragedFields();
            lastFullWriteCycle_ = currentCycle_;
        }
    }
    else
    {
        // default: "phaseAlignedDirs" (Strategy B, back-compat)
        writeAllAveragedFields();
    }

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
        propsDict.regIOobject::write();
    }

    return true;
}


void Foam::functionObjects::BBPA::readCheckpoint()
{
    // Look for the Properties file in any uniform/ subdir, latest first.
    // findInstance returns the time directory containing the file, or
    // throws if missing -- so we ask it to be quiet about a missing
    // file via NO_READ_IF_PRESENT semantics by checking explicitly.
    const word propsName = name() + "Properties";

    // Locate the checkpoint ONCE and remember it.  Re-locating on a later
    // call would find a Properties file written by THIS run, so a
    // late-registering field would be restored against the wrong state.
    if (!checkpointSearched_)
    {
        checkpointSearched_ = true;

        // Search all time dirs in the registry, take the latest
        // that contains "uniform/<propsName>".  Using direct
        // filesystem checks rather than IOobject::typeHeaderOk so we
        // never require constructing the dictionary at the wrong time
        // instance.
        const fileName casePath = obr_.time().path();
        const instantList times = obr_.time().times();
        forAllReverse(times, i)
        {
            const fileName candidate =
                casePath / times[i].name() / "uniform" / propsName;
            if (Foam::isFile(candidate))
            {
                // Construct the full processor-local path directly from the
                // already-verified casePath (== obr_.time().path(), which is
                // processor-aware in parallel runs).  Going through
                // IOdictionary(IOobject(..., obr_, MUST_READ)) here was
                // resolving to the case-level path and crashing in parallel.
                checkpointPath_ = candidate;
                break;
            }
        }

        if (checkpointPath_.empty())
        {
            Info<< type() << " " << name()
                << ": no checkpoint Properties found; cold start."
                << endl;
        }
    }

    // Genuine cold start (no checkpoint on disk): nothing to restore, and
    // every later call is a cheap no-op.
    if (checkpointPath_.empty())
    {
        return;
    }

    const fileName& propsFullPath = checkpointPath_;

    Info<< type() << " " << name()
        << ": restoring checkpoint state from "
        << propsFullPath << endl;

    // Read directly via IFstream + dictionary -- bypasses the IOobject
    // path-resolution machinery that mishandles uniform/ in parallel.
    IFstream propsFile(propsFullPath);
    if (!propsFile.good())
    {
        WarningInFunction
            << "Could not open " << propsFullPath
            << " for read; falling back to cold start." << endl;
        return;
    }
    dictionary propsDict(propsFile);

    // Scalar counters.  Phase origin always taken from the checkpoint:
    // a manual dict override after a restart would desynchronise the
    // accumulator state (which was indexed against the saved
    // startTime_) and is therefore not honoured here.  Restored ONCE:
    // on a later call (a field registering mid-run) currentCycle_ has
    // already advanced and must not be rewound to the checkpoint value.
    if (!countersRestored_)
    {
        startTime_ = propsDict.lookup<scalar>("startTime");
        currentCycle_ = propsDict.lookup<label>("currentCycle");
        previousCycle_ = propsDict.lookup<label>("previousCycle");

        // Capture the checkpoint's cycle start so a late-registering field
        // is restored from the phase-aligned dirs the checkpoint actually
        // wrote, not from an already-advanced cycle.
        checkpointCycleStart_ = startTime_ + currentCycle_ * period_;
        countersRestored_ = true;
    }

    // Reconstruct meanFields_, M2Fields_ and meanSqrFields_ from the
    // phase-aligned _PA / _PAVariance / _PA_UU files.  At the most
    // recent writeTime (just after cycle finalisation), the public
    // bin fields satisfy:
    //
    //     <name>_PA[binI]         at t = cycleStart + binI*binDeltaT
    //                              = meanFields_[binI]
    //     <name>_PAVariance[binI] at the same t
    //                              = M2Fields_[binI] / (nCycles - 1)
    //     <name>_PA_UU[binI]      at the same t
    //                              = meanSqrFields_[binI]
    //
    // Restore by reading from the appropriate phase-aligned dirs.  Uses the
    // checkpoint's own cycle start (captured above), so this is correct both
    // at construction and for a field that registers later.
    const scalar cycleStart = checkpointCycleStart_;

    // Per-field restore: counters + bin fields, skipping any field already
    // restored.  This is what lets a derived field (wallShearStress) that
    // only appears once its own function object has run pick up the very
    // same checkpoint U was restored from, instead of cold-starting.
    forAllIter
    (
        HashTable<autoPtr<binItem<scalar>>>,
        scalarBinItems_,
        iter
    )
    {
        const word& fname = iter.key();
        if (restoredFields_.found(fname)) continue;

        if (propsDict.found(fname))
        {
            const dictionary& fd = propsDict.subDict(fname);
            iter()->nCycles_ = fd.lookup<label>("nCycles");
            iter()->binCounts_ = fd.lookup<scalarField>("binCounts");
            restoreBinItemFromDisk(*(iter()), fname, cycleStart);

            Info<< "    Restored field " << fname
                << ": nCycles=" << iter()->nCycles_ << endl;
        }
        // Field absent from the checkpoint (not tracked in the previous
        // run): a cold start is correct for it. Mark it done so we do not
        // retry the lookup on every step.
        restoredFields_.insert(fname);
    }
    forAllIter
    (
        HashTable<autoPtr<binItem<vector>>>,
        vectorBinItems_,
        iter
    )
    {
        const word& fname = iter.key();
        if (restoredFields_.found(fname)) continue;

        if (propsDict.found(fname))
        {
            const dictionary& fd = propsDict.subDict(fname);
            iter()->nCycles_ = fd.lookup<label>("nCycles");
            iter()->binCounts_ = fd.lookup<scalarField>("binCounts");
            restoreBinItemFromDisk(*(iter()), fname, cycleStart);

            Info<< "    Restored field " << fname
                << ": nCycles=" << iter()->nCycles_ << endl;
        }
        restoredFields_.insert(fname);
    }

    Info<< "    Restored: startTime=" << startTime_
        << "  currentCycle=" << currentCycle_
        << "  previousCycle=" << previousCycle_ << endl;
}


// ************************************************************************* //
