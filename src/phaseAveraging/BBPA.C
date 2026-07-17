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

// Snap a value that is within a whisker of a whole number onto it, then
// truncate. Needed because the natural expression label(x) is exact only when
// x is not sitting on an integer boundary -- and in production it almost always
// IS: a case with period 0.5 and I=50 has a bin width of 0.01, equal to its
// writeInterval, so EVERY bin edge lands on a solver time.
//
// The failure is not hypothetical. fmod(2.9, 1.0) returns 0.8999999999999999,
// not 0.9, so 0.8999999999999999/0.05 = 17.999999999999996 truncates to bin 17
// when the step belongs in bin 18. The error grows with t (the mantissa has
// fewer bits left for the fraction), which is why a run drifts from 51 steps
// per bin in the first cycle to 49 by the third instead of a constant 50.
//
// Whether a given step lands right therefore depended on the last bit of
// OpenFOAM's accumulated clock. At I=200 a bin can be only ~25 timesteps wide,
// so one misplaced step is ~4% of the bin.
//
// The tolerance is in BIN UNITS and is ~7 orders of magnitude smaller than the
// timestep-to-bin ratio of any sane configuration, so it can only ever absorb
// floating-point noise, never a step that genuinely belongs to the lower bin.
static inline Foam::label bbpaSnapFloor(const Foam::scalar x)
{
    const Foam::scalar r = std::floor(x);
    return Foam::label((x - r) > (1.0 - 1e-9) ? r + 1.0 : r);
}


Foam::label Foam::functionObjects::BBPA::phaseBin() const
{
    const scalar t = time_.value() - startTime_;
    scalar phase = std::fmod(t, period_);

    // fmod can return period-1ulp instead of 0 for an on-grid time, which
    // would put the first step of a cycle in the LAST bin of that cycle.
    if (phase > period_*(1.0 - 1e-12)) phase = 0;

    return min(bbpaSnapFloor(phase/binDeltaT_), nBins_ - 1);
}


Foam::label Foam::functionObjects::BBPA::currentCycleNumber() const
{
    const scalar t = time_.value() - startTime_;
    return bbpaSnapFloor(t/period_);
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
                        fieldName, mesh_, nBins_, binsOfInterest_, binNaming_
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
                        fieldName, mesh_, nBins_, binsOfInterest_, binNaming_
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


void Foam::functionObjects::BBPA::finalizeAllForCycle(const label skipBin)
{
    forAllIter(HashTable<autoPtr<binItem<scalar>>>, scalarBinItems_, iter)
    {
        iter()->finalizeCycle(skipBin);
    }
    forAllIter(HashTable<autoPtr<binItem<vector>>>, vectorBinItems_, iter)
    {
        iter()->finalizeCycle(skipBin);
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


void Foam::functionObjects::BBPA::flushCurrentCycleBins(const label cycle)
{
    // The in-progress cycle's window. Only bins whose upper edge the solver has
    // already passed are written, so every directory touched is at or behind
    // the solver clock -- never ahead of it, which is what made restart cleanup
    // delete live bin data in an earlier version.
    const scalar cycleStart = startTime_ + cycle*period_;
    const scalar now = time_.value();

    // Never persist the bin a restart landed inside: it covers only part of its
    // window this cycle, and a file written from it would be recovered later as
    // a complete bin.
    const label skip = (cycle == straddledCycle_) ? straddledBin_ : -1;

    forAllIter(HashTable<autoPtr<binItem<scalar>>>, scalarBinItems_, iter)
    {
        iter()->writeCurrentCycleBins(time_, cycleStart, binDeltaT_, now, skip);
    }
    forAllIter(HashTable<autoPtr<binItem<vector>>>, vectorBinItems_, iter)
    {
        iter()->writeCurrentCycleBins(time_, cycleStart, binDeltaT_, now, skip);
    }
}


Foam::label Foam::functionObjects::BBPA::restoreCurrentCycle()
{
    const scalar cycleStart = startTime_ + currentCycle_*period_;
    const scalar now = time_.value();
    label n = 0;

    forAllIter(HashTable<autoPtr<binItem<scalar>>>, scalarBinItems_, iter)
    {
        n += iter()->restoreCurrentCycleBins
        (
            mesh_, cycleStart, binDeltaT_, now
        );
    }
    forAllIter(HashTable<autoPtr<binItem<vector>>>, vectorBinItems_, iter)
    {
        n += iter()->restoreCurrentCycleBins
        (
            mesh_, cycleStart, binDeltaT_, now
        );
    }

    return n;
}


void Foam::functionObjects::BBPA::writeAllAveragedFields(const bool emitPA)
{
    // Phase-aligned time directories for the cycle the accumulator actually
    // DESCRIBES, which is the one that just finished (finalizeAllForCycle()
    // folds cycle currentCycle_-1 in at the transition, and the perCycle dump
    // fires on the first write after that boundary).
    //
    // This must be the PREVIOUS cycle's window, not the current one.  Writing
    // to the current cycle's window put the bin dirs AHEAD of the solver clock,
    // where they were indistinguishable from a half-written aborted run: the
    // restart cleanup deleted them (it must, or `latestTime` lands on a
    // bin-only dir with no p and the case will not start), and the restore then
    // found the mean missing while the counters had already been restored --
    // silently averaging real cycles against a phantom zero.  Stamping the
    // window the data belongs to keeps every bin dir at or behind the solver
    // clock, so the cleanup never touches it and no cleanup change is needed.
    const scalar cycleStart =
        startTime_ + max(currentCycle_ - 1, 0) * period_;

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
            cycleStart, binDeltaT_,
            binNaming_ ? word("_BinPA") : word("_PA"),
            binNaming_, emitPA
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
            cycleStart, binDeltaT_,
            binNaming_ ? word("_BinPA") : word("_PA"),
            binNaming_, emitPA
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
    binNaming_(false),
    paWriteInterval_(1),
    checkpointSearched_(false),
    checkpointPath_(fileName::null),
    countersRestored_(false),
    checkpointCycleStart_(0),
    restoredFields_(),
    legacyStart_probed_(false),
    straddledBin_(-1),
    straddledCycle_(-1),
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
    // "companion" (one bin into solver's time dir), "both", "perCycle",
    // or "bin".
    writeMode_ = dict.lookupOrDefault<word>
    (
        "writeMode", word("phaseAlignedDirs")
    );

    // "bin" mode: <field>_Bin / _BinPA / _BinPA_UU / _BinPAVariance naming,
    // plus the extra per-cycle _Bin file. Everything else behaves as
    // "perCycle" (one full dump per completed cycle).
    binNaming_ = (writeMode_ == "bin");

    // How often the cross-cycle _BinPA* files are written, in completed
    // cycles. The per-cycle _Bin file is always written.
    paWriteInterval_ = dict.lookupOrDefault<label>("PAwriteInterval", 1);

    if (binNaming_)
    {
        // Warn when the bins will NOT land on the solver's write times, since
        // the whole point of this mode is that <time>/U and <time>/U_Bin sit
        // side by side. This is a configuration smell, not an error: the mode
        // still works, it just creates directories of its own.
        //
        // Bin times are k*binWidth; solver write times are m*writeInterval.
        // The bins land in directories the solver already wrote iff every
        // k*binWidth is some m*writeInterval, i.e. iff binWidth is an integer
        // multiple of writeInterval -- NOT the other way round.
        //
        // Getting this backwards inverts the test exactly where it matters:
        // period 0.5 with I=200 gives binWidth 0.0025 against a writeInterval
        // of 0.01, so writeInterval/binWidth = 4 is a clean integer while the
        // bins fall at 0.0025, 0.005, 0.0075 -- times the solver never writes.
        // That is the very configuration this warning exists to catch.
        //
        // NB: binDeltaT_ is not assigned until below, so use period_/nBins_.
        // Time has no public writeInterval() accessor, so read controlDict.
        const scalar binWidth = period_/nBins_;
        const scalar wi =
            time_.controlDict().lookupOrDefault<scalar>("writeInterval", -1);
        const word wc =
            time_.controlDict().lookupOrDefault<word>("writeControl", word::null);

        if (wi > 0 && wc == "runTime")
        {
            const scalar ratio = binWidth/wi;
            if (mag(ratio - std::round(ratio)) > 1e-6)
            {
                WarningInFunction
                    << "writeMode 'bin': the bin width period/nBins ("
                    << binWidth << ") is not an integer multiple of the solver "
                    << "writeInterval (" << wi << ")." << nl
                    << "    The bins fall at times the solver does not write, "
                    << "so the _Bin files will create directories of their own "
                    << "rather than landing beside the instantaneous fields."
                    << nl
                    << "    For co-located output choose nBins so that "
                    << "period/nBins is a multiple of " << wi << "."
                    << endl;
            }
        }
    }

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
            // The completing cycle's LAST bin closes exactly at this boundary,
            // so it has never been flushed by a write() call inside the cycle.
            // Persist it before resetAllBinsForNewCycle() wipes the
            // accumulator, or its _Bin_UU would never be written at all.
            if (binNaming_)
            {
                flushCurrentCycleBins(previousCycle_);
            }

            // Finalize the previous cycle (only if it was being averaged).
            // If a restart landed mid-bin during THIS cycle, that one bin saw
            // only part of its phase window and must not be folded in as a
            // whole cycle -- drop it, exactly as a bin with no data is dropped.
            const label skip =
                (previousCycle_ == straddledCycle_) ? straddledBin_ : -1;

            if (skip >= 0)
            {
                Info<< "    dropping restart-straddled bin " << skip
                    << " from cycle " << previousCycle_
                    << " (partial window)" << endl;
            }

            finalizeAllForCycle(skip);

            // consumed: later cycles cover this bin fully again
            if (skip >= 0)
            {
                straddledBin_ = -1;
                straddledCycle_ = -1;
            }
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

    // "bin" mode: persist the CURRENT cycle's newly-closed bins before
    // anything else. This is what makes the output restart-safe at write
    // granularity: a kill then loses only the bin it lands in, instead of the
    // entire in-progress cycle (measured: 15 of 25 bins lost to a SIGKILL
    // before this, 1 after).
    if (binNaming_)
    {
        flushCurrentCycleBins(currentCycle_);
    }

    if (writeMode_ == "companion")
    {
        writeCompanionAllFields();
    }
    else if (writeMode_ == "both")
    {
        writeCompanionAllFields();
        if (isLastWrite)
        {
            writeAllAveragedFields(true);
        }
    }
    else if (writeMode_ == "bin" || writeMode_ == "perCycle")
    {
        // Full phase-aligned dump, but only once per completed cycle (the
        // first write after a cycle boundary) instead of every solver write.
        // Cuts the bin-output I/O by ~period/writeInterval while staying a
        // full, restart-safe checkpoint at each cycle (restart granularity =
        // one cycle). lastFullWriteCycle_ starts at -1 and is not persisted,
        // so the first write after any (re)start re-establishes the checkpoint.
        if (currentCycle_ > lastFullWriteCycle_)
        {
            // "bin" mode: the per-cycle _Bin file is written every cycle,
            // but the cross-cycle _BinPA* only every PAwriteInterval cycles
            // (or, when PAwriteInterval <= 0, on the final write only).
            // The completed-cycle count is currentCycle_ (cycle currentCycle_-1
            // has just been folded in).
            const bool emitPA =
                !binNaming_
             || (paWriteInterval_ <= 0
                    ? isLastWrite
                    : (currentCycle_ % paWriteInterval_) == 0)
             || isLastWrite;

            writeAllAveragedFields(emitPA);
            lastFullWriteCycle_ = currentCycle_;
        }
    }
    else
    {
        // default: "phaseAlignedDirs" (Strategy B, back-compat)
        writeAllAveragedFields(true);
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

        // Straddled-bin guard. The run resumes at the last WRITTEN time, which
        // need not be a bin boundary. The bin containing that instant will now
        // accumulate only the part of its phase window lying after the restart,
        // and finalizeCycle() would fold that partial-window mean in as a
        // complete cycle at full weight -- biased toward whichever phases
        // happen to follow the restart, and indistinguishable from a healthy
        // value on disk.
        //
        // Detect it exactly rather than by a coverage heuristic: the bin is
        // straddled iff the restart phase lies STRICTLY INSIDE it. Landing
        // exactly on a bin edge (the aorta's binWidth == writeInterval case)
        // covers the bin fully and must not be dropped.
        {
            const scalar tRes = time_.value() - startTime_;
            scalar phase = std::fmod(tRes, period_);
            if (phase > period_*(1.0 - 1e-12)) phase = 0;

            const label binI = min(bbpaSnapFloor(phase/binDeltaT_), nBins_ - 1);
            const scalar binStart = binI*binDeltaT_;

            // strictly inside => the earlier part of this bin's window is gone
            if (phase - binStart > 1e-9*binDeltaT_)
            {
                straddledBin_ = binI;
                straddledCycle_ = currentCycle_;

                WarningInFunction
                    << "Restart at t=" << time_.value() << " (phase "
                    << phase << ") lands INSIDE bin " << binI
                    << " [" << binStart << ", " << binStart + binDeltaT_
                    << "): that bin's window is only partially covered in cycle "
                    << currentCycle_ << "." << nl
                    << "    Dropping bin " << binI << " for cycle "
                    << currentCycle_ << " rather than folding a partial-window "
                    << "mean in at full weight (unbiased, one cycle noisier)."
                    << nl
                    << "    To avoid this entirely, make the solver "
                    << "writeInterval an integer multiple of period/nBins ("
                    << binDeltaT_ << ") so restarts always land on a bin edge."
                    << endl;
            }
            else
            {
                straddledBin_ = -1;
                straddledCycle_ = -1;
            }
        }

        // Capture the checkpoint's cycle start so a late-registering field
        // is restored from the phase-aligned dirs the checkpoint actually
        // wrote, not from an already-advanced cycle.  Must mirror
        // writeAllAveragedFields(): the bins live in the window of the cycle
        // they describe (currentCycle_-1), which is at/behind the solver clock.
        checkpointCycleStart_ = startTime_ + max(currentCycle_ - 1, 0) * period_;

        // Backward compatibility: checkpoints written before this fix stamped
        // the bins one cycle LATER (currentCycle_*period_). If the current
        // convention finds nothing on disk but the legacy window does, fall
        // back to it so an in-flight run restarts against its own history
        // instead of silently cold-starting.
        {
            const scalar legacyStart = startTime_ + currentCycle_ * period_;
            const word probe = fieldNames_.size() ? fieldNames_[0] : word::null;
            if (probe != word::null && !legacyStart_probed_)
            {
                legacyStart_probed_ = true;
                const fileName base = obr_.time().path();
                const bool newOk = Foam::isFile
                (
                    base/Time::timeName(checkpointCycleStart_)/(probe + "_PA")
                );
                const bool oldOk = Foam::isFile
                (
                    base/Time::timeName(legacyStart)/(probe + "_PA")
                );
                if (!newOk && oldOk)
                {
                    Info<< "    legacy checkpoint layout detected: reading bins "
                        << "from the pre-fix window (cycleStart=" << legacyStart
                        << ")" << endl;
                    checkpointCycleStart_ = legacyStart;
                }
            }
        }
        countersRestored_ = true;

        // Recover the in-progress cycle's completed bins, written as each bin
        // closed. Without this the whole partial cycle is discarded even though
        // every completed bin of it is sitting on disk. Runs AFTER the counters
        // are restored (it needs currentCycle_ to locate the window) and only
        // covers bins whose window closed strictly before the restart instant;
        // the bin containing the restart is partial and is dropped by the
        // straddle guard instead.
        if (binNaming_)
        {
            const label nrec = restoreCurrentCycle();
            if (nrec > 0)
            {
                Info<< "    recovered " << nrec
                    << " completed bin(s) of the in-progress cycle "
                    << currentCycle_ << " from <field>_Bin on disk"
                    << " (partial cycle preserved across the restart)" << endl;
            }
        }
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
            restoreBinItemFromDisk
            (
                *(iter()), fname, cycleStart,
                binNaming_ ? word("_BinPA") : word("_PA")
            );

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
            restoreBinItemFromDisk
            (
                *(iter()), fname, cycleStart,
                binNaming_ ? word("_BinPA") : word("_PA")
            );

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
