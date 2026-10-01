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

Description
    BBPA::writer -- the output-layout concern (v2 step 5c).

    Everything here decides file NAMES, DIRECTORIES and what gets written
    WHEN; nothing here computes or mutates a statistic (flushClosedBins
    only marks bins as flushed). The calculation lives in binIntegrator /
    cycleStatistics (binAccumulators.H) and binItem::finalizeCycle
    (BBPATemplates.C); the writer READS those components through a
    const binItem and serialises them. Keeping the two in separate files
    is what makes "the output mode never changes the numbers" a property
    of the structure rather than only of the tests.

\*---------------------------------------------------------------------------*/

#include "BBPA.H"
#include "volFields.H"
#include "OFstream.H"
#include "timeDirName.H"

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


template<class Type>
void Foam::functionObjects::BBPA::writer::flushClosedBins
(
    binItem<Type>& bi,
    const boolList& writeBins,
    const Time& time,
    const scalar cycleStartTime,
    const scalar binDeltaT,
    const scalar now,
    const label skipBin,
    const scalar centreOffset
)
{
    typedef typename binItem<Type>::FieldType FieldType;
    typedef typename binItem<Type>::SqrFieldType SqrFieldType;

    if (bi.binFlushed_.empty()) return;

    forAll(bi.integ_.binSum_, binI)
    {
        if (!bi.integ_.binSum_.set(binI)) continue;   // sparse-bin mode
        if (binI == skipBin) continue;                  // restart-straddled
        if (bi.binFlushed_[binI]) continue;                // already on disk
        if (bi.integ_.binTime_[binI] <= 0) continue;  // no data yet

        // The bin is complete only once the solver has passed its UPPER edge.
        // Writing earlier would persist a partial window that a restart would
        // then restore as if it were whole -- the very bug the straddle guard
        // exists to prevent.
        const scalar binEnd =
            cycleStartTime + (binI + 1)*binDeltaT - centreOffset;
        if (now < binEnd - 1e-9*binDeltaT) continue;

        const scalar binTime = cycleStartTime + binI*binDeltaT;
        const fvMesh& mesh = bi.integ_.binSum_[binI].mesh();
        const word binTimeName = timeDirName
        (
            mesh.time(), binTime,
            timeDirTol(binDeltaT, mesh.time().deltaTValue())
        );
        const scalar inv = 1.0/bi.integ_.binTime_[binI];

        FieldType binField
        (
            IOobject
            (
                bi.fieldName_ + "_Bin"
              + (writeBins[binI] ? word::null : word("State")),
                binTimeName,
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            bi.integ_.binSum_[binI]*inv
        );
        writePAField(binField, mesh, binTimeName);

        // The per-cycle second moment must be persisted too: the restore
        // rebuilds bi.integ_.binSqrSum_ from it. Without it a recovered bin
        // would carry a full bin time but a zero second moment, and
        // finalizeCycle() would fold that zero into <UiUj> as real data.
        // (Skipped when the secondMoment product is off -- nothing
        // accumulates, so there is nothing to persist or restore.)
        if (bi.integ_.binSqrSum_.set(binI))
        {
            SqrFieldType binUU
            (
                IOobject
                (
                    bi.fieldName_ + "_Bin_UU"
                  + (writeBins[binI] ? word::null : word("State")),
                    binTimeName,
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                bi.integ_.binSqrSum_[binI]*inv
            );
            writePAField(binUU, mesh, binTimeName);
        }

        bi.binFlushed_[binI] = true;
    }
}


template<class Type>
void Foam::functionObjects::BBPA::writer::phaseAligned
(
    const binItem<Type>& bi,
    const boolList& writeBins,
    const Time& time,
    const label nBins,
    const label currentCycle,
    const label currentBin,
    const scalar cycleStartTime,
    const scalar binDeltaT,
    const word& paSuffix,
    const bool emitBin,
    const bool emitPA
)
{
    typedef typename binItem<Type>::FieldType FieldType;
    typedef typename binItem<Type>::SqrFieldType SqrFieldType;

    // Find a mesh handle from the first ALLOCATED bin (in sparse-bin mode
    // bin 0 may be unallocated).
    const fvMesh* meshPtr = nullptr;
    forAll(bi.stats_.mean_, b)
    {
        if (bi.stats_.mean_.set(b)) { meshPtr = &bi.stats_.mean_[b].mesh(); break; }
    }
    if (meshPtr == nullptr) return;   // nothing allocated -> nothing to write
    const fvMesh& mesh = *meshPtr;
    label binsWritten = 0;

    forAll(bi.stats_.mean_, binI)
    {
        // Skip bins not allocated in sparse-bin mode.
        if (!bi.stats_.mean_.set(binI)) continue;



        const bool hasCrossData = bi.stats_.binCounts_[binI] > 0;
        const bool hasCurrentData = bi.integ_.binTime_[binI] > 0;

        if (hasCrossData || hasCurrentData)
        {
            // Phase-aligned time for this bin. When binDeltaT equals the
            // solver's writeInterval this IS one of the solver's own write
            // directories, so the bin lands next to the instantaneous field
            // rather than creating a directory of its own.
            const scalar binTime = cycleStartTime + binI * binDeltaT;
            const word binTimeName = timeDirName
            (
                time, binTime, timeDirTol(binDeltaT, time.deltaTValue())
            );

            // <field>_Bin -- THIS cycle's bin value, a distinct quantity from
            // the cross-cycle mean below. Writing the two separately is what
            // removes any reason to blend the in-progress cycle into the mean.
            // It is also what makes the output restart-safe at write
            // granularity: every write directory carries a complete bin.
            if
            (
                emitBin
             && writeBins[binI]
             && bi.lastCycleFields_.size()
             && bi.lastCycleFields_.set(binI)
             && bi.lastCycleBinTime_[binI] > 0
            )
            {
                FieldType binField
                (
                    IOobject
                    (
                        bi.fieldName_ + "_Bin",
                        binTimeName,
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    bi.lastCycleFields_[binI]
                );

                writePAField(binField, mesh, binTimeName);
            }

            if (!emitPA)
            {
                binsWritten++;
                continue;
            }

            // Products are FINALISED-ONLY: nothing is published (or
            // checkpointed) for a bin before its first complete cycle.
            // A partial integral written under a product name would be
            // indistinguishable on disk from a converged value.
            if (!hasCrossData)
            {
                continue;
            }

            // Product files obey outputSelection and the per-field output
            // flags. The ensemble STATE of every accumulated bin is
            // checkpointed regardless, under the *State spellings, so
            // outputSelection is serialisation-only: it may change across
            // a restart without touching restart continuity.
            const bool meanAsProduct =
                writeBins[binI] && bi.products_.writeMean;

            {
                FieldType paField
                (
                    IOobject
                    (
                        bi.fieldName_ + paSuffix
                      + (meanAsProduct ? word::null : word("meanState")),
                        binTimeName,
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    bi.stats_.mean_[binI]
                );

                writePAField(paField, mesh, binTimeName);
            }

            binsWritten++;

            // Welford cycle-variance (cycle-to-cycle variability of the
            // bin means; defined for two or more cycles)
            if (bi.stats_.M2_.set(binI) && bi.stats_.binCounts_[binI] > 1)
            {
                volScalarField paVariance
                (
                    IOobject
                    (
                        bi.fieldName_ + paSuffix
                      + bbpaVarianceSuffix<Type>()
                      + (writeBins[binI] ? word::null : word("State")),
                        binTimeName,
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    bi.stats_.M2_[binI] / (bi.stats_.binCounts_[binI] - 1.0)
                );

                writePAField(paVariance, mesh, binTimeName);
            }

            // Second moment <phi*phi> over the COMPLETE cycles: the raw
            // tensor from which Reynolds stress and TKE are reconstructed
            // offline. Whole-field construction so the boundary is
            // written (wallShearStress carries ALL its data there).
            if (bi.stats_.meanSqr_.set(binI))
            {
                SqrFieldType paUU
                (
                    IOobject
                    (
                        bi.fieldName_ + paSuffix + "_UU"
                      + (writeBins[binI] ? word::null : word("State")),
                        binTimeName,
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    bi.stats_.meanSqr_[binI]
                );

                writePAField(paUU, mesh, binTimeName);
            }
        }
    }

    Info<< "    " << bi.fieldName_ << " phase-averaged output:" << nl
        << "        Cycles accumulated: " << bi.stats_.nCycles_ << nl
        << "        Bins written: " << binsWritten << "/" << nBins << nl
        << "        Time range: " << Time::timeName(cycleStartTime)
        << " to " << Time::timeName(cycleStartTime + (nBins-1)*binDeltaT)
        << endl;
}


template<class Type>
void Foam::functionObjects::BBPA::writer::companion
(
    const binItem<Type>& bi,
    const boolList& writeBins,
    const Time& time,
    const label currentBinIndex
)
{
    typedef typename binItem<Type>::FieldType FieldType;
    typedef typename binItem<Type>::SqrFieldType SqrFieldType;

    // Companion write: one phase (the current bin) into the solver's
    // current time directory. Intended to be paired with the solver's
    // own instantaneous field write at the same time, giving the
    // reader both U and U_PA in the same directory.
    const label binI = currentBinIndex;
    if (binI < 0 || binI >= bi.stats_.mean_.size()) return;

    // Sparse-bin mode: nothing to write if the current bin is not tracked.
    if (!bi.stats_.mean_.set(binI)) return;

    // The companion file IS the phase-mean product.
    if (!bi.products_.writeMean) return;
    if (!writeBins[currentBinIndex]) return;    // outputSelection

    // Finalised-only: the companion file is a product and is never
    // built from a partial window.
    const bool hasCrossData = bi.stats_.binCounts_[binI] > 0;
    if (!hasCrossData) return;

    const fvMesh& mesh = bi.stats_.mean_[binI].mesh();
    const word solverTimeName = time.timeName(time.value());

    // Build PA field for this bin, labelled with solver's current
    // time (not the phase-aligned time).
    FieldType paField
    (
        IOobject
        (
            bi.fieldName_ + "_PA",
            solverTimeName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        bi.stats_.mean_[binI]
    );

    // The mean over the COMPLETE cycles only -- the same no-blend rule as the
    // phase-aligned write() above, for the same reason: the in-progress bin is
    // PARTIAL (often a single timestep) yet the old
    //     (bi.stats_.mean_*nOld + bi.integ_.binSum_/bi.integ_.binTime_)
    //         / (nOld + 1)
    // blend weighted it as a whole cycle. Removing it here makes the companion
    // output numerically identical to every other writeMode -- the output
    // layout must never change the numbers. Only when no complete cycle exists
    // yet is there nothing to report but the running partial. Whole-field
    // assignment (not primitiveFieldRef) so the boundary is written too.
    if (!hasCrossData)
    {
        const scalar inv = 1.0 / bi.integ_.binTime_[binI];
        paField = bi.integ_.binSum_[binI] * inv;
    }

    mkDir(paField.objectPath().path());
    paField.regIOobject::write();
}


// ************************************************************************* //
