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
    BBPA::state -- the restart-reconstruction concern (v2 step 5d).

    The mirror of BBPA::writer: where the writer serialises a binItem's
    components to disk, state rebuilds them FROM disk on restart --
    restoreStatistics reverse-engineers the cross-cycle ensemble
    (stats_) from the public _PA / _PAVariance / _PA_UU files the writer
    emitted, and restoreClosedBins recovers the in-progress cycle's
    completed bins (integ_) from the bin-mode _Bin / _Bin_UU flush.
    Neither computes a statistic; both only invert what the writer wrote.
    The checkpoint COUNTERS (cycle numbers, per-bin counts) stay with
    BBPA::readCheckpoint()/write(), which orchestrate these calls.

\*---------------------------------------------------------------------------*/

#include "BBPA.H"
#include "timeDirName.H"
#include "volFields.H"
#include "IFstream.H"

template<class Type>
Foam::label
Foam::functionObjects::BBPA::state::restoreClosedBins
(
    binItem<Type>& bi,
    const fvMesh& mesh,
    const scalar cycleStartTime,
    const scalar binDeltaT,
    const scalar now,
    const scalar centreOffset
)
{
    typedef typename binItem<Type>::FieldType FieldType;
    typedef typename binItem<Type>::SqrFieldType SqrFieldType;

    if (bi.binFlushed_.empty()) return 0;

    label recovered = 0;

    forAll(bi.integ_.binSum_, binI)
    {
        if (!bi.integ_.binSum_.set(binI)) continue;

        // Recover any bin whose window STARTED before the restart instant and
        // whose complete <field>_Bin is on disk. This deliberately INCLUDES the
        // bin that CONTAINS the restart instant (the straddled bin): the flush
        // only ever writes a CLOSED window, so an on-disk _Bin is the WHOLE
        // window average -- strictly better than re-accumulating only its
        // post-restart part. This is what lets a restart from an EARLIER write
        // (e.g. one write before the kill) lose nothing even with binCentred:
        // the bin centred on that write has its complete _Bin on disk from
        // before the kill, so it is recovered rather than straddle-dropped.
        //
        // A bin whose window starts at or after the restart has no pre-restart
        // data (its dir was never written, or was trimmed) and is re-accumulated
        // by the resumed solver instead.
        const scalar binStart =
            cycleStartTime + binI*binDeltaT - centreOffset;
        if (binStart >= now - 1e-9*binDeltaT) continue;

        const scalar binTime = cycleStartTime + binI*binDeltaT;
        const word binTimeName = timeDirName
        (
            mesh.time(), binTime,
            timeDirTol(binDeltaT, mesh.time().deltaTValue())
        );

        // Restore only the products that accumulate (allocation is the
        // gate): with the second moment off there is no _Bin_UU on disk and
        // no accumulator to fill. When the product is on, both files or
        // neither -- a bin recovered with a mean but no second moment would
        // fold a zero <UiUj> in at full weight.
        const bool withSqr = bi.integ_.binSqrSum_.set(binI);

        // outputSelection may have changed since the flush: closed bins
        // persist under _Bin (published) or _BinState (state-only) -- try
        // the product name first, then the state spelling.
        word binName = bi.fieldName_ + "_Bin";
        if (!Foam::isFile(mesh.time().path()/binTimeName/binName))
        {
            binName = bi.fieldName_ + "_BinState";
        }
        word uuName = bi.fieldName_ + "_Bin_UU";
        if (!Foam::isFile(mesh.time().path()/binTimeName/uuName))
        {
            uuName = bi.fieldName_ + "_Bin_UUState";
        }

        const fileName binPath = mesh.time().path()/binTimeName/binName;
        const fileName uuPath = mesh.time().path()/binTimeName/uuName;

        if (!Foam::isFile(binPath)) continue;
        if (withSqr && !Foam::isFile(uuPath)) continue;

        FieldType loaded
        (
            IOobject
            (
                binName,
                binTimeName,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh
        );

        // The files hold the bin MEANS; the accumulators hold time-weighted
        // SUMS. The bin is complete, so its weight is exactly binDeltaT.
        bi.integ_.binSum_[binI] = loaded*binDeltaT;

        if (withSqr)
        {
            SqrFieldType loadedUU
            (
                IOobject
                (
                    uuName,
                    binTimeName,
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh
            );

            bi.integ_.binSqrSum_[binI] = loadedUU*binDeltaT;
        }

        bi.integ_.binTime_[binI] = binDeltaT;
        bi.binFlushed_[binI] = true;      // already on disk; do not rewrite
        recovered++;
    }

    return recovered;
}


// * * * * * * * * * * * * Restore from checkpoint * * * * * * * * * * * * * //
//
// At restart the existing _PA / _PAVariance / _PA_UU files written
// at the most recent solver writeTime carry the cycle-N ensemble
// state.  Per binItem<Type>::write() at the post-finalisation cycle
// boundary (stats_.binCounts_[binI] > 0 && integ_.binTime_[binI] == 0):
//
//     paField     = stats_.mean_[binI]                  (no blending)
//     paVariance  = stats_.M2_[binI] / (n - 1)          (only if n > 1)
//     paUU        = stats_.meanSqr_[binI]               (no blending)
//
// Reverse-engineer the three internal arrays from the public outputs.
template<class Type>
void Foam::functionObjects::BBPA::state::restoreStatistics
(
    binItem<Type>& bi,
    const word& fieldName,
    const scalar cycleStartTime,
    const scalar binDeltaT,
    const word& paSuffix
)
{
    typedef typename binItem<Type>::FieldType FieldType;
    typedef typename binItem<Type>::SqrFieldType SqrFieldType;

    // Mesh handle from the first allocated bin (sparse-bin mode may
    // leave bin 0 unallocated).
    const fvMesh* meshPtr = nullptr;
    forAll(bi.stats_.mean_, b)
    {
        if (bi.stats_.mean_.set(b)) { meshPtr = &bi.stats_.mean_[b].mesh(); break; }
    }
    if (meshPtr == nullptr) return;
    const fvMesh& mesh = *meshPtr;
    label restored = 0;

    forAll(bi.stats_.mean_, binI)
    {
        // Skip bins not allocated in sparse-bin mode.
        if (!bi.stats_.mean_.set(binI)) continue;

        const scalar nCyc = bi.stats_.binCounts_[binI];
        if (nCyc <= 0) continue;  // bin never received data

        const scalar binTime = cycleStartTime + binI * binDeltaT;
        const word binTimeName = timeDirName
        (
            mesh.time(), binTime,
            timeDirTol(binDeltaT, mesh.time().deltaTValue())
        );

        // 1. mean: the writer emits it for every finalised accumulated
        //    bin -- as the _PA product when published, as the meanState
        //    spelling otherwise. outputSelection may have changed since
        //    the checkpoint, so try the product name first and fall back
        //    to the state spelling.
        {
            word meanName = fieldName + paSuffix;
            if
            (
                !Foam::isFile(mesh.time().path()/binTimeName/meanName)
             && Foam::isFile
                (
                    mesh.time().path()/binTimeName
                   /(fieldName + paSuffix + "meanState")
                )
            )
            {
                meanName = fieldName + paSuffix + "meanState";
            }

            IOobject paIO
            (
                meanName,
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
                // The mean for this bin is gone but its cycle count survived in
                // the Properties dict.  Keeping the count would make Welford
                // average the remaining cycles against a phantom zero -- with
                // count n0 and k further cycles the bin converges to
                // k*x/(n0+k), i.e. biased LOW by n0/(n0+k) and indistinguishable
                // from a healthy value (non-zero, full-size file).  Reset the
                // bin's history instead: it loses the lost cycles and is
                // noisier, but it is UNBIASED.
                WarningInFunction
                    << "Checkpoint mean field missing: " << binTimeName
                    << "/" << fieldName << paSuffix << " -- bin " << binI
                    << " reset to a clean cold start (count " << nCyc
                    << " -> 0) to avoid a low-biased average." << endl;
                bi.stats_.binCounts_[binI] = 0;
                bi.stats_.mean_[binI] =
                    dimensioned<Type>(bi.stats_.mean_[binI].dimensions(), Zero);
                if (bi.stats_.M2_.set(binI))
                {
                    bi.stats_.M2_[binI] =
                        dimensionedScalar(bi.stats_.M2_[binI].dimensions(), 0);
                }
                if (bi.stats_.meanSqr_.set(binI))
                {
                    bi.stats_.meanSqr_[binI] =
                        dimensioned<typename binItem<Type>::SqrType>
                        (bi.stats_.meanSqr_[binI].dimensions(), Zero);
                }
                continue;
            }

            // Whole-field assignment, not primitiveFieldRef(): restoring only
            // the internal field leaves the accumulator's BOUNDARY at its
            // constructed zero, so every restart silently discarded the
            // boundary phase-average and then folded new cycles against a
            // phantom zero history. For wallShearStress -- whose data is
            // entirely on the wall patch -- that meant the restart threw away
            // the whole accumulated field while the Properties counters still
            // claimed the cycles were there.
            FieldType loaded(paIO, mesh);
            bi.stats_.mean_[binI] = loaded;
        }

        // 2. M2: read the variance file and rescale by (n-1)
        //    (only when the variance product is allocated; product name
        //    first, then the state spelling).
        if (bi.stats_.M2_.set(binI) && nCyc > 1)
        {
            word varName = fieldName + paSuffix + bbpaVarianceSuffix<Type>();
            if (!Foam::isFile(mesh.time().path()/binTimeName/varName))
            {
                varName =
                    fieldName + paSuffix + bbpaVarianceSuffix<Type>()
                  + "State";
            }

            IOobject varIO
            (
                varName,
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
                bi.stats_.M2_[binI] = loadedVar * (nCyc - 1.0);
            }
        }

        // 3. meanSqr: read the second-moment file into stats_.meanSqr_
        //    (only when allocated; product name first, then the state
        //    spelling).
        if (bi.stats_.meanSqr_.set(binI))
        {
            word uuName = fieldName + paSuffix + "_UU";
            if (!Foam::isFile(mesh.time().path()/binTimeName/uuName))
            {
                uuName = fieldName + paSuffix + "_UUState";
            }

            IOobject uuIO
            (
                uuName,
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
                bi.stats_.meanSqr_[binI] = loadedUU;
            }
        }

        restored++;
    }

    Info<< "    " << fieldName << ": restored "
        << restored << "/" << bi.stats_.mean_.size()
        << " bins from phase-aligned dirs (cycleStart="
        << cycleStartTime << ", nCycles=" << bi.stats_.nCycles_ << ")"
        << endl;
}


// ************************************************************************* //
