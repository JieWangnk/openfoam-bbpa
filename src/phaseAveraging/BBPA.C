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
#include "DynamicList.H"
#include "Switch.H"

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

// All phase/cycle/bin boundary arithmetic -- the snapped floor, the fmod
// wrap-around guard, the label/phase origins -- lives in phaseClock (see
// phaseClock.H for the full rationale). These wrappers exist so every caller
// classifies the CURRENT solver instant the same way.

Foam::label Foam::functionObjects::BBPA::phaseBin() const
{
    return clock_.phaseBin(time_.value());
}


Foam::label Foam::functionObjects::BBPA::currentCycleNumber() const
{
    return clock_.cycleNumber(time_.value());
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
    const scalar cycleStart = clock_.cycleStart(cycle);
    const scalar now = time_.value();
    const scalar centreOffset = clock_.centreOffset();

    // Never persist the bin a restart landed inside: it covers only part of its
    // window this cycle, and a file written from it would be recovered later as
    // a complete bin.
    const label skip = (cycle == straddledCycle_) ? straddledBin_ : -1;

    forAllIter(HashTable<autoPtr<binItem<scalar>>>, scalarBinItems_, iter)
    {
        writer::flushClosedBins
        (
            *iter(), time_, cycleStart, binDeltaT_, now, skip, centreOffset
        );
    }
    forAllIter(HashTable<autoPtr<binItem<vector>>>, vectorBinItems_, iter)
    {
        writer::flushClosedBins
        (
            *iter(), time_, cycleStart, binDeltaT_, now, skip, centreOffset
        );
    }
}


Foam::label Foam::functionObjects::BBPA::restoreCurrentCycle()
{
    const scalar cycleStart = clock_.cycleStart(currentCycle_);
    const scalar now = time_.value();
    const scalar centreOffset = clock_.centreOffset();
    label n = 0;

    forAllIter(HashTable<autoPtr<binItem<scalar>>>, scalarBinItems_, iter)
    {
        n += state::restoreClosedBins
        (
            *iter(), mesh_, cycleStart, binDeltaT_, now, centreOffset
        );
    }
    forAllIter(HashTable<autoPtr<binItem<vector>>>, vectorBinItems_, iter)
    {
        n += state::restoreClosedBins
        (
            *iter(), mesh_, cycleStart, binDeltaT_, now, centreOffset
        );
    }

    return n;
}


void Foam::functionObjects::BBPA::flagStraddledBin()
{
    // The bin containing the instant this run BEGAN covers only the part of its
    // phase window lying after that instant, yet finalizeCycle() would fold it
    // in as a complete cycle at full weight -- a partial-window mean, biased
    // toward whichever phases follow the start, and indistinguishable from a
    // healthy value on disk.
    //
    // This is a property of the START INSTANT, so it applies to a cold start
    // exactly as much as to a restart. With binCentred, bin 0 of the first
    // cycle wants [-dT/2, +dT/2) but the run begins at startTime: half that
    // window never existed. The check below catches both.
    //
    // Detect it exactly, not by a coverage heuristic: the bin is straddled iff
    // the start phase lies STRICTLY INSIDE it. Beginning exactly on a bin edge
    // (the ordinary uncentred cold start) covers the bin fully and must not be
    // dropped.
    const scalar phase = clock_.wrapPhase(runStartTime_);
    const label binI = clock_.phaseBin(runStartTime_);
    const scalar binStart = binI*binDeltaT_;
    const label cyc = clock_.cycleNumber(runStartTime_);

    if (phase - binStart > 1e-9*binDeltaT_)
    {
        // If this straddled bin's COMPLETE _Bin is already on disk, the run
        // resumed from a write EARLIER than the one that closed this bin, so
        // state::restoreClosedBins recovered the whole window -- do NOT drop it.
        // Only drop when it is genuinely unrecoverable: the cold-start bin 0
        // (nothing on disk yet), or a restart from the very write that straddles
        // a bin whose window had not closed (hence not flushed) before the kill.
        bool recoverable = false;
        if (binNaming_ && !fieldNames_.empty())
        {
            const scalar labelTime = clock_.binLabelTime(cyc, binI);
            const fileName binFile =
                mesh_.time().path()/Time::timeName(labelTime)
               /(fieldNames_[0] + "_Bin");
            recoverable = Foam::isFile(binFile);
        }

        if (recoverable)
        {
            straddledBin_ = -1;
            straddledCycle_ = -1;
            Info<< "    bin " << binI << " straddled by the restart, but its "
                << "complete _Bin is on disk -> recovered whole, not dropped "
                << "(zero loss)." << endl;
        }
        else
        {
            straddledBin_ = binI;
            straddledCycle_ = cyc;

            WarningInFunction
                << "run begins at t=" << runStartTime_ << " (phase " << phase
                << "), strictly inside bin " << binI << " [" << binStart << ", "
                << binStart + binDeltaT_ << "): that bin's window is only "
                << "partially covered in cycle " << cyc << " and no complete "
                << "_Bin is on disk." << nl
                << "    Dropping bin " << binI << " for cycle " << cyc
                << " rather than folding a partial-window mean in at full weight"
                << " (unbiased, one cycle noisier). To avoid it, resume from one"
                << " write EARLIER so this bin's window is recovered whole."
                << endl;
        }
    }
    else
    {
        straddledBin_ = -1;
        straddledCycle_ = -1;
    }
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
    const scalar cycleStart = clock_.cycleStart(max(currentCycle_ - 1, 0));

    forAllConstIter
    (
        HashTable<autoPtr<binItem<scalar>>>,
        scalarBinItems_,
        iter
    )
    {
        writer::phaseAligned
        (
            *iter(), 
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
        writer::phaseAligned
        (
            *iter(), 
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
        writer::companion(*iter(), time_, currentBin_);
    }
    forAllConstIter
    (
        HashTable<autoPtr<binItem<vector>>>, vectorBinItems_, iter
    )
    {
        writer::companion(*iter(), time_, currentBin_);
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
    binsAtWriteTimes_(false),
    binCentred_(false),
    phaseOrigin_(0),
    paWriteInterval_(1),
    checkpointSearched_(false),
    checkpointPath_(fileName::null),
    countersRestored_(false),
    checkpointCycleStart_(0),
    restoredFields_(),
    inProgressCycleStart_(0),
    inProgressNow_(0),
    inProgressCaptured_(false),
    legacyStart_probed_(false),
    straddledBin_(-1),
    straddledCycle_(-1),
    straddleChecked_(false),
    runStartTime_(runTime.value()),
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

    if (fieldNames_.empty())
    {
        FatalErrorInFunction
            << "'fields' is empty; give at least one field to phase-average, "
            << "e.g. fields (U)."
            << exit(FatalError);
    }

    if (period_ <= 0)
    {
        FatalErrorInFunction
            << "'period' must be positive; got " << period_
            << ". Set it to the cycle period in seconds."
            << exit(FatalError);
    }

    if (nBins_ <= 0)
    {
        FatalErrorInFunction
            << "'nBins' must be positive; got " << nBins_
            << ". It is the number of phase bins per period."
            << exit(FatalError);
    }

    // Optional sparse-bin mode: only track the listed bin indices.
    // Default empty list = track all nBins_ bins (original behaviour).
    binsOfInterest_ =
        dict.lookupOrDefault<labelList>("binsOfInterest", labelList());

    // Validate any supplied bin indices against [0, nBins_), rejecting
    // duplicates (a duplicate would double-account that bin's statistics).
    {
        labelHashSet seen;
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
            if (!seen.insert(b))
            {
                FatalErrorInFunction
                    << "binsOfInterest entry " << b
                    << " is listed more than once; give each bin index once."
                    << exit(FatalError);
            }
        }
    }

    // Output mode: "phaseAlignedDirs" (default, Strategy B),
    // "companion" (one bin into solver's time dir), "both", "perCycle",
    // or "bin".
    writeMode_ = dict.lookupOrDefault<word>
    (
        "writeMode", word("phaseAlignedDirs")
    );

    // A misspelt writeMode must not silently fall through to the default
    // layout: the write() dispatch treats every unknown word as
    // "phaseAlignedDirs".
    if
    (
        writeMode_ != "phaseAlignedDirs"
     && writeMode_ != "perCycle"
     && writeMode_ != "companion"
     && writeMode_ != "both"
     && writeMode_ != "bin"
    )
    {
        FatalErrorInFunction
            << "unknown writeMode '" << writeMode_ << "'; valid values: "
            << "phaseAlignedDirs, perCycle, companion, both, bin."
            << exit(FatalError);
    }

    // "bin" mode: <field>_Bin / _BinPA / _BinPA_UU / _BinPAVariance naming,
    // plus the extra per-cycle _Bin file. Everything else behaves as
    // "perCycle" (one full dump per completed cycle).
    binNaming_ = (writeMode_ == "bin");

    // How often the cross-cycle _BinPA* files are written, in completed
    // cycles. The per-cycle _Bin file is always written.
    paWriteInterval_ = dict.lookupOrDefault<label>("PAwriteInterval", 1);


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

    // "bin" mode conveniences: derive the tracked bins from the solver's write
    // times, and centre each bin's window on its label.
    binsAtWriteTimes_ = dict.lookupOrDefault<Switch>("binsAtWriteTimes", false);
    binCentred_ = dict.lookupOrDefault<Switch>("binCentred", false);

    if (binsAtWriteTimes_)
    {
        if (!binsOfInterest_.empty())
        {
            WarningInFunction
                << "both binsAtWriteTimes and an explicit binsOfInterest were "
                << "given; keeping the explicit list." << endl;
        }
        else
        {
            const scalar wi =
                time_.controlDict().lookupOrDefault<scalar>("writeInterval", -1);
            const word wc =
                time_.controlDict().lookupOrDefault<word>
                (
                    "writeControl", word::null
                );

            // Both runTime and adjustableRunTime write at fixed writeInterval
            // spacing (adjustableRunTime just nudges the timestep to land on
            // the write times exactly) -- either lets us know where the writes
            // are. The aorta production runs use adjustableRunTime. clockTime /
            // cpuTime / timeStep do NOT give a fixed time spacing and are
            // rejected.
            const bool fixedTimeWrites =
                (wc == "runTime" || wc == "adjustableRunTime");
            if (!fixedTimeWrites || wi <= 0)
            {
                FatalErrorInFunction
                    << "binsAtWriteTimes needs writeControl runTime or "
                    << "adjustableRunTime with a positive writeInterval to know "
                    << "where the write times are (found writeControl " << wc
                    << ", writeInterval " << wi << ")." << exit(FatalError);
            }

            // A bin is labelled j*binDeltaT_ and coincides with a write time
            // k*wi exactly when j is a multiple of wi/binDeltaT_. That ratio
            // must be a whole number or no bin label ever lands on a write.
            const scalar ratio = wi/binDeltaT_;
            const label stride = label(std::round(ratio));

            if (mag(ratio - stride) > 1e-6 || stride < 1)
            {
                FatalErrorInFunction
                    << "binsAtWriteTimes requires the solver writeInterval ("
                    << wi << ") to be an integer multiple of the bin width "
                    << "period/nBins (" << binDeltaT_ << "); the ratio is "
                    << ratio << "." << nl
                    << "    No bin label can land on a write time otherwise. "
                    << "Choose nBins so that period/nBins divides " << wi << "."
                    << exit(FatalError);
            }

            DynamicList<label> bins;
            for (label j = 0; j < nBins_; j += stride)
            {
                bins.append(j);
            }
            binsOfInterest_ = bins;

            Info<< "    binsAtWriteTimes: writeInterval " << wi
                << " / binWidth " << binDeltaT_ << " = stride " << stride
                << " -> tracking " << binsOfInterest_.size() << " of "
                << nBins_ << " bins, all landing on solver write times"
                << endl;
        }
    }

    // Centring is a shift of the phase AND cycle origin by half a bin. The bin
    // LABELS keep using startTime_, so they stay on the write times; only the
    // window they describe moves to be centred on them.
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
        // Runs AFTER binsAtWriteTimes derived binsOfInterest_: checking before
        // it would warn about bins that are about to be aligned.
        // Time has no public writeInterval() accessor, so read controlDict.
        const scalar binWidth = binDeltaT_;
        const scalar wi =
            time_.controlDict().lookupOrDefault<scalar>("writeInterval", -1);
        const word wc =
            time_.controlDict().lookupOrDefault<word>("writeControl", word::null);

        // In sparse-bin mode only the LISTED bins are written, so the test is
        // whether each of THOSE lands on a write time -- not whether every bin
        // would. The recommended configuration is exactly this: a fine bin width
        // for sharpness, with binsOfInterest picking the subset that coincides
        // with the solver's writes (e.g. period 0.5, nBins 100 -> bin 0.005,
        // writeInterval 0.02, bins 0 4 8 ... -> 25 co-located bins at I=100
        // sharpness). Testing all 100 bins there would warn about directories
        // that are never created.
        bool misaligned = false;
        if (wi > 0 && (wc == "runTime" || wc == "adjustableRunTime"))
        {
            if (binsOfInterest_.empty())
            {
                const scalar ratio = binWidth/wi;
                misaligned = mag(ratio - std::round(ratio)) > 1e-6;
            }
            else
            {
                forAll(binsOfInterest_, k)
                {
                    const scalar r = binsOfInterest_[k]*binWidth/wi;
                    if (mag(r - std::round(r)) > 1e-6)
                    {
                        misaligned = true;
                        break;
                    }
                }
            }
        }

        if (misaligned)
        {
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

    phaseOrigin_ = binCentred_ ? startTime_ - 0.5*binDeltaT_ : startTime_;
    clock_.reset(period_, nBins_, startTime_, binCentred_);

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

    // Once, after startTime_/phaseOrigin_ are final (readCheckpoint may have
    // changed them). Covers the cold start too, which never reaches
    // readCheckpoint's guard because there is no checkpoint to read.
    if (!straddleChecked_)
    {
        straddleChecked_ = true;
        flagStraddledBin();
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
        // The final averaged cycle completes exactly at this boundary, and
        // the transition block below (which normally finalises the previous
        // cycle on the first step of the next one) is no longer reached.
        // Finalise it once here, or it stays stranded in the current-cycle
        // accumulators and "cycles N" silently delivers an (N-1)-cycle
        // ensemble. One-shot: previousCycle_ is moved past the limit after.
        if
        (
            previousCycle_ >= startCycle_
         && previousCycle_ < (startCycle_ + cycles_)
        )
        {
            if (binNaming_)
            {
                flushCurrentCycleBins(previousCycle_);
            }

            const label skip =
                (previousCycle_ == straddledCycle_) ? straddledBin_ : -1;

            finalizeAllForCycle(skip);

            if (skip >= 0)
            {
                straddledBin_ = -1;
                straddledCycle_ = -1;
            }

            resetAllBinsForNewCycle();
            previousCycle_ = currentCycle_;

            Info<< "Cycle limit reached: " << cycles_
                << " cycles averaged; accumulation stopped" << endl;
        }

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
    // The stats_.mean_ and stats_.M2_ are written by AUTO_WRITE; this
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
            fieldDict.add("nCycles", iter()->stats_.nCycles_);
            fieldDict.add("binCounts", iter()->stats_.binCounts_);
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
            fieldDict.add("nCycles", iter()->stats_.nCycles_);
            fieldDict.add("binCounts", iter()->stats_.binCounts_);
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

        // startTime_ has just changed, so the derived phase origin must follow
        // it. Missing this would leave phaseOrigin_ at the value read() derived
        // from the RESTART time, silently un-centring every bin after a restart
        // while the labels stayed put.
        phaseOrigin_ = binCentred_ ? startTime_ - 0.5*binDeltaT_ : startTime_;
        clock_.reset(period_, nBins_, startTime_, binCentred_);

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
        // (straddle is flagged by flagStraddledBin(), called from execute())

        // Capture the checkpoint's cycle start so a late-registering field
        // is restored from the phase-aligned dirs the checkpoint actually
        // wrote, not from an already-advanced cycle.  Must mirror
        // writeAllAveragedFields(): the bins live in the window of the cycle
        // they describe (currentCycle_-1), which is at/behind the solver clock.
        checkpointCycleStart_ = clock_.cycleStart(max(currentCycle_ - 1, 0));

        // Backward compatibility: checkpoints written before this fix stamped
        // the bins one cycle LATER (currentCycle_*period_). If the current
        // convention finds nothing on disk but the legacy window does, fall
        // back to it so an in-flight run restarts against its own history
        // instead of silently cold-starting.
        {
            const scalar legacyStart = clock_.cycleStart(currentCycle_);
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

        // Capture the restart instant. The in-progress cycle recovery itself
        // is done PER FIELD in the loops below (guarded by restoredFields_),
        // so a late-registering derived field (wallShearStress) recovers its
        // partial cycle against this SAME instant -- not against the later time
        // it happens to register, and not skipped entirely (the bug: the
        // recovery used to run once here, for U only, so WSS lost a whole cycle
        // on every bin whose window had closed before the restart).
        inProgressCycleStart_ = clock_.cycleStart(currentCycle_);
        inProgressNow_ = time_.value();
        inProgressCaptured_ = true;
    }

    // Reconstruct stats_.mean_, stats_.M2_ and stats_.meanSqr_ from the
    // phase-aligned _PA / _PAVariance / _PA_UU files.  At the most
    // recent writeTime (just after cycle finalisation), the public
    // bin fields satisfy:
    //
    //     <name>_PA[binI]         at t = cycleStart + binI*binDeltaT
    //                              = stats_.mean_[binI]
    //     <name>_PAVariance[binI] at the same t
    //                              = stats_.M2_[binI] / (nCycles - 1)
    //     <name>_PA_UU[binI]      at the same t
    //                              = stats_.meanSqr_[binI]
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
            // Guard against a checkpoint from a run with a DIFFERENT nBins.
            // stats_.binCounts_ is sized nBins_; scalarField::operator= would RESIZE it
            // to the checkpoint's size, after which every downstream loop over
            // nBins_ bins indexes it OUT OF BOUNDS -> garbage counts and a
            // free(): invalid pointer crash. Reject a mismatched checkpoint and
            // cold-start this field instead (the usual cause is a stale
            // <name>Properties left in the start-time dir after nBins changed;
            // clean it, or this is exactly the right thing to do).
            const scalarField loadedCounts(fd.lookup<scalarField>("binCounts"));
            if (loadedCounts.size() != nBins_)
            {
                WarningInFunction
                    << "checkpoint for field " << fname << " has "
                    << loadedCounts.size() << " bins but nBins is " << nBins_
                    << " -- ignoring this checkpoint and cold-starting the field"
                    << " (a stale " << name() << "Properties from a run with a"
                    << " different nBins? remove it from the start-time dir)."
                    << endl;
                restoredFields_.insert(fname);
                continue;
            }
            iter()->stats_.nCycles_ = fd.lookup<label>("nCycles");
            iter()->stats_.binCounts_ = loadedCounts;
            state::restoreStatistics
            (
                *(iter()), fname, cycleStart, binDeltaT_,
                binNaming_ ? word("_BinPA") : word("_PA")
            );

            // Recover THIS field's in-progress cycle (the completed bins of the
            // partial cycle, written as each bin closed). Per-field and guarded
            // by restoredFields_, so it runs exactly once per field -- including
            // a late-registering derived field -- and never re-clobbers a field
            // that has already accumulated new steps since the restart.
            if (binNaming_ && inProgressCaptured_)
            {
                const scalar cOff = binCentred_ ? 0.5*binDeltaT_ : 0.0;
                const label nrec = state::restoreClosedBins
        (
                    *iter(), mesh_, inProgressCycleStart_, binDeltaT_, inProgressNow_, cOff
                );
                if (nrec > 0)
                {
                    Info<< "    " << fname << ": recovered " << nrec
                        << " completed bin(s) of the in-progress cycle "
                        << currentCycle_ << " from " << fname
                        << "_Bin on disk" << endl;
                }
            }

            Info<< "    Restored field " << fname
                << ": nCycles=" << iter()->stats_.nCycles_ << endl;
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
            // Guard against a checkpoint from a run with a DIFFERENT nBins.
            // stats_.binCounts_ is sized nBins_; scalarField::operator= would RESIZE it
            // to the checkpoint's size, after which every downstream loop over
            // nBins_ bins indexes it OUT OF BOUNDS -> garbage counts and a
            // free(): invalid pointer crash. Reject a mismatched checkpoint and
            // cold-start this field instead (the usual cause is a stale
            // <name>Properties left in the start-time dir after nBins changed;
            // clean it, or this is exactly the right thing to do).
            const scalarField loadedCounts(fd.lookup<scalarField>("binCounts"));
            if (loadedCounts.size() != nBins_)
            {
                WarningInFunction
                    << "checkpoint for field " << fname << " has "
                    << loadedCounts.size() << " bins but nBins is " << nBins_
                    << " -- ignoring this checkpoint and cold-starting the field"
                    << " (a stale " << name() << "Properties from a run with a"
                    << " different nBins? remove it from the start-time dir)."
                    << endl;
                restoredFields_.insert(fname);
                continue;
            }
            iter()->stats_.nCycles_ = fd.lookup<label>("nCycles");
            iter()->stats_.binCounts_ = loadedCounts;
            state::restoreStatistics
            (
                *(iter()), fname, cycleStart, binDeltaT_,
                binNaming_ ? word("_BinPA") : word("_PA")
            );

            // Recover THIS field's in-progress cycle (the completed bins of the
            // partial cycle, written as each bin closed). Per-field and guarded
            // by restoredFields_, so it runs exactly once per field -- including
            // a late-registering derived field -- and never re-clobbers a field
            // that has already accumulated new steps since the restart.
            if (binNaming_ && inProgressCaptured_)
            {
                const scalar cOff = binCentred_ ? 0.5*binDeltaT_ : 0.0;
                const label nrec = state::restoreClosedBins
        (
                    *iter(), mesh_, inProgressCycleStart_, binDeltaT_, inProgressNow_, cOff
                );
                if (nrec > 0)
                {
                    Info<< "    " << fname << ": recovered " << nrec
                        << " completed bin(s) of the in-progress cycle "
                        << currentCycle_ << " from " << fname
                        << "_Bin on disk" << endl;
                }
            }

            Info<< "    Restored field " << fname
                << ": nCycles=" << iter()->stats_.nCycles_ << endl;
        }
        restoredFields_.insert(fname);
    }

    Info<< "    Restored: startTime=" << startTime_
        << "  currentCycle=" << currentCycle_
        << "  previousCycle=" << previousCycle_ << endl;
}


// ************************************************************************* //
