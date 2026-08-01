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
        const word binTimeName = Time::timeName(binTime);
        const fvMesh& mesh = bi.integ_.binSum_[binI].mesh();
        const scalar inv = 1.0/bi.integ_.binTime_[binI];

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
                    bi.fieldName_ + "_Bin_UU",
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
            const word binTimeName = Time::timeName(binTime);

            // <field>_Bin -- THIS cycle's bin value, a distinct quantity from
            // the cross-cycle mean below. Writing the two separately is what
            // removes any reason to blend the in-progress cycle into the mean.
            // It is also what makes the output restart-safe at write
            // granularity: every write directory carries a complete bin.
            if
            (
                emitBin
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

            // Build the PA field values without triggering tmp<> expressions
            // on registered fields. We compute in a copy.
            // phaseMean off suppresses only this FILE; the mean chain is
            // still accumulated (the Welford variance needs it).
            if (bi.products_.mean)
            {
            FieldType paField
            (
                IOobject
                (
                    bi.fieldName_ + paSuffix,
                    binTimeName,
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                bi.stats_.mean_[binI]
            );

            // paField was constructed FROM bi.stats_.mean_[binI], i.e. the mean over
            // the bi.stats_.binCounts_[binI] COMPLETE cycles folded so far. That is what
            // <name>_PA must contain.
            //
            // It must NOT be blended with the in-progress cycle. The old code did
            //     (bi.stats_.mean_*nOld + bi.integ_.binSum_/bi.integ_.binTime_)
            //         / (nOld + 1)
            // which is wrong twice over:
            //  (1) the in-progress bin is PARTIAL -- with writeMode perCycle the
            //      write fires on the cycle transition, when the new cycle has
            //      accumulated exactly ONE timestep into bin 0 -- yet it is
            //      weighted as a whole cycle. Measured on tutorials/bbpaVerify:
            //      bin 0 came out 0.632 against a true 0.701, i.e. (mean*2 + U_inst)/3,
            //      while every other bin (just reset, so no current data) was exact.
            //  (2) these files ARE the restart checkpoint: restoreBinItemFromDisk()
            //      reads bi.stats_.mean_ back out of <name>_PA. A blended value is
            //      therefore restored AS the folded mean and compounds at every
            //      restart.
            //
            // Only when no complete cycle exists yet is there nothing to report but
            // the running partial. bi.stats_.binCounts_[binI] is 0 in that case, and the
            // restore skips bins with count 0, so it can never be read back as a mean.
            //
            // Assign the WHOLE field, not primitiveFieldRef(): an internal-only
            // assignment leaves the boundary at whatever the accumulator was
            // constructed with (zero), silently dropping the boundary values
            // of the partial bin. See the _PA_UU note below.
            if (!hasCrossData && hasCurrentData)
            {
                const scalar inv = 1.0 / bi.integ_.binTime_[binI];
                paField = bi.integ_.binSum_[binI] * inv;
            }

            writePAField(paField, mesh, binTimeName);
            }

            binsWritten++;

            // Welford cycle-variance (cycle-to-cycle variability of bin mean)
            if (bi.stats_.M2_.set(binI) && bi.stats_.binCounts_[binI] > 1)
            {
                volScalarField paVariance
                (
                    IOobject
                    (
                        bi.fieldName_ + paSuffix + "Variance",
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

            // Second-moment output: <phi*phi>_bin as a phase-aligned field.
            // For vector U this is the 6-component symmetric tensor
            // <U_i U_j> from which the Reynolds stress R_ij =
            // <U_i U_j> - <U_i><U_j> and TKE = 0.5*trace(R) are computed
            // offline. For scalar p this reduces to <p*p>.
            // Fires whenever we have data (finalised or in-progress bin)
            // and the product is on (allocation is the gate).
            if (bi.stats_.meanSqr_.set(binI) && (hasCrossData || hasCurrentData))
            {
                // Construct FROM bi.stats_.meanSqr_[binI], exactly as paField is
                // constructed from bi.stats_.mean_[binI] above. This is not a
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
                        bi.fieldName_ + paSuffix + "_UU",
                        binTimeName,
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    bi.stats_.meanSqr_[binI]
                );

                // The raw second moment over the COMPLETE cycles only -- as for
                // _PA above, and for the same two reasons: blending the
                // in-progress cycle weights a partial bin (one timestep, on a
                // perCycle write) as a whole cycle, and this file is the
                // checkpoint that restoreBinItemFromDisk() reads bi.stats_.meanSqr_
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
                    const scalar invT = 1.0 / bi.integ_.binTime_[binI];
                    paUU = bi.integ_.binSqrSum_[binI] * invT;
                }

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
    if (!bi.products_.mean) return;

    const bool hasCrossData = bi.stats_.binCounts_[binI] > 0;
    const bool hasCurrentData = bi.integ_.binTime_[binI] > 0;
    if (!hasCrossData && !hasCurrentData) return;

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
