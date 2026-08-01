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

Class
    Foam::functionObjects::BBPA

Description
    Bin-Based Phase Averaging (BBPA) functionObject

    Performs time integration of field values within phase bins and averages
    across multiple cycles. Preserves full vector/tensor field structure.

    Mathematical formulation:
    - Temporal quantization: Divides period T into I equal bins of width
      dt_i = T/I
    - Bin integration: phi_Bin,j = (1/dt_j) integral phi(t) dt for each bin j
    - Phase averaging: PA_phi,j^(N) = (1/N) sum[n=1 to N] phi_Bin,n,j

\*---------------------------------------------------------------------------*/

#ifndef BBPA_H
#define BBPA_H

#include "fvMeshFunctionObject.H"
#include "phaseClock.H"
#include "binAccumulators.H"
#include "volFieldsFwd.H"
#include "HashPtrTable.H"
#include "HashSet.H"
#include "fileName.H"
#include "scalarField.H"
#include "boolList.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

// bbpaSqrTypeOf (the phi*phi trait) and the two data components of a binItem
// -- binIntegrator (current cycle) and cycleStatistics (cross-cycle ensemble)
// -- live in binAccumulators.H (v2 separation of calculation concerns).


/*---------------------------------------------------------------------------*\
                          Class BBPA Declaration
\*---------------------------------------------------------------------------*/

class BBPA
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Fields to average
        wordList fieldNames_;

        //- Per-field product options (v2): the optional dictionary after
        //  each entry in
        //      fields ( U { secondMoment off; } p );
        //  keyed by field name; a missing/empty dictionary keeps every
        //  product at its default. Resolved to a binProducts when the
        //  field's binItem is created.
        HashTable<dictionary> fieldOptions_;

        //- Period of oscillation
        scalar period_;

        //- Number of bins per cycle
        label nBins_;

        //- Optional list of bin indices to track ("sparse-bin" mode).
        //   When empty (default) ALL nBins_ bins are allocated and
        //   accumulated exactly as before -- the original full-bin
        //   behaviour is unchanged. When non-empty, only the listed
        //   bins allocate field-sized accumulators and are
        //   accumulated/written; every other bin is skipped. This
        //   reduces the in-memory cost from O(nBins) to
        //   O(binsOfInterest) field-sized accumulators per tracked
        //   variable -- useful when only a few phases (e.g. peak
        //   systole, post-peak deceleration) are needed at large I.
        labelList binsOfInterest_;

        //- Number of cycles to average (-1 for all)
        label cycles_;

        //- Starting cycle number (phase averaging begins from this cycle)
        label startCycle_;

        //- Start time for phase averaging
        scalar startTime_;

        //- Time step size for bins
        scalar binDeltaT_;

        //- Output mode for phase-averaged fields.
        //   "phaseAlignedDirs" (default, back-compat): at each write()
        //     call, dump all I bins to phase-aligned time directories
        //     t_j = cycleStart + j*binDeltaT. ParaView-friendly but
        //     produces I dirs per write.
        //   "companion": at each write() call, write only the current
        //     phase's bin into the solver's current time directory,
        //     alongside the instantaneous fields. One directory per
        //     write, no extra dirs created.
        //   "both": do both -- companion every write + a full
        //     phase-aligned dump on the final call.
        //   "perCycle": like "phaseAlignedDirs" but the full dump is
        //     written only once per completed cycle (the first write
        //     after a cycle boundary), not every write -- cuts the
        //     bin-output I/O by ~period/writeInterval while staying a
        //     full, restart-safe checkpoint at each cycle boundary.
        //   "bin": like "perCycle", but names the outputs <field>_Bin*
        //     and additionally emits <field>_Bin -- THIS cycle's bin
        //     value -- alongside the cross-cycle <field>_BinPA. See
        //     binNaming_ below for why the extra file matters.
        word writeMode_;

        //- "bin" mode: use the <field>_Bin / _BinPA / _BinPA_UU /
        //  _BinPAVariance naming instead of <field>_PA / _PA_UU /
        //  _PAVariance, and emit the extra _Bin file.
        //
        //  The point of _Bin is not cosmetic. The old output tried to serve
        //  two masters with one file: <field>_PA was the mean over complete
        //  cycles, but a reader also wanted to know what THIS cycle did, so
        //  the companion path blended the in-progress cycle into it --
        //  weighting a partial bin (often a single timestep) as a whole
        //  cycle, and corrupting the very file the restart reads its mean
        //  back out of. Emitting the two separately -- _Bin for this cycle,
        //  _BinPA for the average across cycles -- removes the reason to
        //  blend at all.
        //
        //  When binDeltaT_ equals the solver's writeInterval (the usual
        //  production case: period 0.5 with I=50 gives 0.01), the
        //  phase-aligned bin directories ARE the solver's own write
        //  directories, so <time>/U, <time>/U_Bin and <time>/U_BinPA land
        //  side by side and no extra directories are created.
        bool binNaming_;

        //- "bin" mode: derive binsOfInterest automatically from the solver's
        //  writeInterval instead of requiring it by hand.
        //
        //  A bin is labelled at j*binDeltaT_, so it coincides with a solver
        //  write time k*writeInterval exactly when j is a multiple of
        //  writeInterval/binDeltaT_. Listing those by hand
        //  ("binsOfInterest (0 4 8 ... 96)") is busywork that silently goes
        //  wrong the moment writeInterval or nBins changes -- the FO would keep
        //  binning phases nobody writes, and create directories of its own.
        //  Deriving it removes that whole class of mistake.
        //
        //  Requires writeInterval to be an integer multiple of binDeltaT_,
        //  which is exactly the condition for a bin label to land on a write
        //  time at all.
        bool binsAtWriteTimes_;

        //- "bin" mode: centre each bin's window on its label time instead of
        //  starting it there.
        //
        //  A rectangular average over [t_j, t_j+dT) represents the signal at
        //  the window's CENTRE, not its left edge, so labelling it t_j imposes
        //  a half-bin phase lag -- a first-order error that is pure
        //  bookkeeping. main.tex already removes it in post-processing by
        //  referencing each bin to t_j + dT/2; this does it at source, so
        //  <time>/U_Bin is the average of the window CENTRED on <time>.
        //
        //  Implemented as a phase-origin shift of -dT/2 (see phaseOrigin_).
        //  The cycle origin must shift with it: a bin centred on phase 0 spans
        //  [-dT/2, +dT/2), straddling the cycle boundary, and its first half
        //  would be discarded by the cycle reset. Shifting both keeps every bin
        //  wholly inside one cycle. Bin LABELS are unaffected -- they stay on
        //  the solver's write times.
        bool binCentred_;

        //- Phase/cycle origin actually used by phaseBin() and
        //  currentCycleNumber(): startTime_ - dT/2 when binCentred_, else
        //  startTime_. Kept separate so bin LABELS (which use startTime_)
        //  remain on the write times.
        scalar phaseOrigin_;

        //- The phase/cycle/bin boundary arithmetic, configured from the
        //  members above wherever phaseOrigin_ is (re)computed. Every
        //  WHICH-bin/WHICH-cycle decision goes through this one object --
        //  see phaseClock.H for why (four copies of this logic once drifted
        //  apart).
        phaseClock clock_;

        //- "bin" mode: write the cross-cycle _BinPA* files only every Nth
        //  completed cycle. The per-cycle _Bin file is always written (it is
        //  what makes the output restart-safe at write granularity).
        //   1  = every cycle (default)
        //   N  = every N cycles
        //  <=0 = only on the final write
        label paWriteInterval_;

        //- Current cycle number
        label currentCycle_;

        //- Previous cycle number
        label previousCycle_;

        //- Current bin number
        label currentBin_;

        //- Previous bin number
        label previousBin_;

        //- Last cycle whose full phase-aligned set was written
        //  (writeMode "perCycle" only). -1 = none yet; resets to -1 on
        //  (re)start, so the first write after any restart re-checkpoints.
        label lastFullWriteCycle_;

        //- Lazy checkpoint restore state.
        //   Derived fields (e.g. wallShearStress, produced by their own
        //   function object) are NOT yet in the objectRegistry when this
        //   FO is constructed, so they get no binItem at construction and
        //   were therefore missed by the one-shot checkpoint restore --
        //   they cold-started on every restart while U (a genuine solver
        //   field) restored correctly. readCheckpoint() is now idempotent
        //   and per-field: initializeBinItems() reports when a field
        //   registers late, and the checkpoint is re-applied to just that
        //   field, against the SAME checkpoint the run started from.

            //- Has the checkpoint file been located yet?
            bool checkpointSearched_;

            //- Path of the Properties checkpoint restored from
            //  (empty => genuine cold start; no checkpoint on disk).
            fileName checkpointPath_;

            //- Have the run-level counters been restored?
            bool countersRestored_;

            //- Cycle start time of the restored checkpoint. Captured once so
            //  a late-registering field is restored against the checkpoint's
            //  phase-aligned dirs, not against an already-advanced cycle.
            scalar checkpointCycleStart_;

            //- Fields whose checkpoint state has been restored.
            HashSet<word> restoredFields_;

            //- Restart instant, captured at the FIRST readCheckpoint so a
            //  late-registering field (e.g. wallShearStress) recovers its
            //  in-progress cycle against the SAME instant U did, not against
            //  the later time at which it happens to register.
            scalar inProgressCycleStart_;
            scalar inProgressNow_;
            bool   inProgressCaptured_;

            //- Has the legacy (pre-fix) bin-window layout been probed for?
            bool legacyStart_probed_;

        //- Straddled-bin guard.
        //   A restart resumes at the last written time, which need not be a bin
        //   boundary. The bin CONTAINING that instant then accumulates only the
        //   part of its phase window that lies after the restart, yet
        //   finalizeCycle() would fold it in as a COMPLETE cycle at full
        //   weight: a partial-window mean, biased toward the phases that
        //   happen to fall after the restart, averaged in as if it were the
        //   whole bin. Nothing about the value looks wrong on disk.
        //
        //   Measured on tutorials/bbpaVerifyLES with a SIGKILL at t=3.54: with
        //   nBins=25 (bin width 0.04) phase 0.54 falls at bin 13.5 -- mid-bin --
        //   and the restarted run's matched bins differed from an identical
        //   continuous run by 5e-4. With nBins=50 (bin width 0.02) the same
        //   restart lands exactly on bin 27's edge, no bin is straddled, and the
        //   matched bins agree to 1e-12. Only nBins changed.
        //
        //   This is EXACTLY detectable rather than a coverage heuristic: the bin
        //   is straddled if and only if the restart phase lies strictly inside
        //   it. Treat that one bin, for that one cycle, exactly as a bin that
        //   received no data at all -- drop it. It loses one cycle and is
        //   noisier, but it is UNBIASED, which is the same trade the
        //   missing-checkpoint recovery already makes.
        //
        //   Note the two constraints this interacts with, which conflict:
        //     co-location (bins land in solver write dirs): binWidth = k*writeInterval
        //     restart safety (no straddle possible):        writeInterval = m*binWidth
        //   Both hold only when binWidth == writeInterval (the aorta's I=50
        //   config). This guard removes the need to rely on that alignment.

            //- Bin straddled by the restart instant (-1 = none).
            label straddledBin_;

            //- Cycle in which straddledBin_ is contaminated (-1 = none).
            label straddledCycle_;

            //- Has the straddle check been done for this run's start instant?
            bool straddleChecked_;

            //- The instant this run actually began (restart time, or the cold
            //  start time). The straddle is a property of THAT instant, not of
            //  whether a checkpoint happened to exist.
            scalar runStartTime_;

        //- Bin item structure for each field (templated to preserve type)
        template<class Type>
        struct binItem
        {
            // OF14: GeometricField template params reordered
            // <Type, PatchField, GeoMesh> -> <Type, GeoMesh, PrimitiveField>.
            // VolField<Type> (volFieldsFwd.H) == GeometricField<Type, fvMesh, Field>.
            typedef VolField<Type> FieldType;
            typedef typename bbpaSqrTypeOf<Type>::type SqrType;
            typedef VolField<SqrType> SqrFieldType;

            word fieldName_;

            //- CURRENT cycle's time-weighted accumulation (weights, first
            //  and second moments per bin). execute() writes here; reset at
            //  every cycle boundary. See binAccumulators.H.
            binIntegrator<Type> integ_;

            //- CROSS-CYCLE ensemble over completed cycles (counts, running
            //  mean, Welford M2, mean second moment <phi*phi> from which
            //  Reynolds stress and TKE are recovered offline). Written only
            //  by the finalizeCycle() fold. See binAccumulators.H.
            cycleStatistics<Type> stats_;

            //- "bin" mode only (empty otherwise): the bin mean of the cycle
            //  that JUST completed, captured by finalizeCycle() before the
            //  accumulators are reset.
            //
            //  It has to be captured rather than read from
            //  integ_.binSum_ at write time, because execute() runs
            //  finalizeAllForCycle() and then resetAllBinsForNewCycle() on the
            //  cycle transition, and only afterwards does the framework call
            //  write(). By then integ_.binTime_ is zero and the completed
            //  cycle's bin is gone -- folded into stats_.mean_ and unrecoverable
            //  from it.
            PtrList<FieldType> lastCycleFields_;

            //- Bin time accumulated by the cycle captured in lastCycleFields_.
            //  Zero means that bin received no data in that cycle.
            scalarField lastCycleBinTime_;

            //- "bin" mode: bins of the CURRENT cycle already flushed to disk as
            //  <field>_Bin / <field>_Bin_UU. Reset every cycle.
            //
            //  Writing each bin the moment its phase window closes -- rather
            //  than batching all of them at the cycle boundary -- is what makes
            //  the output restart-safe at WRITE granularity instead of CYCLE
            //  granularity. Measured on tutorials/bbpaVerifyLES: a SIGKILL at
            //  t=3.602 threw away cycle 3's bins 0-13 even though every one of
            //  them was complete and sitting in memory, purely because nothing
            //  had persisted them yet -- 15 of 25 bins lost a cycle. With the
            //  incremental flush the same kill costs only the straddled bin.
            //
            //  It is close to free: the same one file per bin per cycle, just
            //  written as the bin closes instead of all at once afterwards.
            boolList binFlushed_;

            //- Which products this field carries (v2 per-field controls).
            //  Off calculation products never allocate their accumulators,
            //  so the PtrList::set() guards gate them everywhere downstream;
            //  products_.writeMean gates only the _PA file in the writer
            //  (the mean itself always accumulates).
            binProducts products_;

            //- Construct. activeBins lists the bin indices that should
            //  allocate field-sized accumulators; an EMPTY list means
            //  "all nBins bins" (original behaviour). Non-active bins
            //  are left as null entries in the PtrLists and are skipped
            //  by accumulate()/finalizeCycle()/write().
            //  prods selects the products: off ones are never allocated
            //  (writePerCycle allocates lastCycleFields_, one extra field
            //  per active bin).
            binItem
            (
                const word& name,
                const fvMesh& mesh,
                const label nBins,
                const labelList& activeBins,
                const binProducts& prods
            );

            void resetCurrentCycle();
            void accumulate
            (
                const label binI,
                const FieldType& field,
                const scalar dt
            );
            //- Fold the completed cycle into the running mean.
            //  skipBin (-1 = none) drops one bin for this cycle: used for the
            //  bin straddled by a restart instant, whose window is only
            //  partially covered and must not be weighted as a whole cycle.
            void finalizeCycle(const label skipBin = -1);

        };


        //- The output-layout concern (v2 separation): file names,
        //  directories, and what is written when. Reads a binItem's
        //  components; never computes or mutates a statistic
        //  (flushClosedBins only marks bins as flushed). Definitions in
        //  BBPAWriterTemplates.C.
        struct writer
        {
            //- Write the phase-aligned bin outputs.
            //  paSuffix selects the naming: "_PA" (default/back-compat) or
            //  "_BinPA" ("bin" mode). emitBin additionally writes
            //  <field>_Bin, the CURRENT cycle's bin value, which is a
            //  different quantity from the cross-cycle mean and must never
            //  be blended into it. emitPA gates the cross-cycle files so
            //  they can be written less often than the per-cycle _Bin.
            template<class Type>
            static void phaseAligned
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
            );

            //- Companion write: dump only the current-phase bin's PA
            //  field into the solver's current time directory (no
            //  phase-aligned sub-directories). Paired with an
            //  instantaneous snapshot in the same time dir.
            template<class Type>
            static void companion
            (
                const binItem<Type>& bi,
                const Time& time,
                const label currentBinIndex
            );

            //- "bin" mode: flush every bin of the CURRENT cycle whose phase
            //  window has closed at or before 'now' and that has not been
            //  written yet, as <field>_Bin and <field>_Bin_UU in the bin's own
            //  phase-aligned directory.
            //
            //  Both files are needed, not just the mean: the restore
            //  reconstructs integ_.binSum_ AND integ_.binSqrSum_ from
            //  them. Restoring the mean alone would leave the second-moment
            //  accumulator at zero while its bin time said "full", and
            //  finalizeCycle() would then fold a genuine zero into <UiUj> --
            //  silently destroying the Reynolds stress for that bin.
            //
            //  Only bins strictly BEHIND the solver clock are written, so no
            //  directory is ever stamped ahead of the run (the failure that
            //  made restart cleanup delete live bin data).
            //  skipBin (-1 = none) suppresses the restart-straddled bin: it
            //  covers only part of its window this cycle, so persisting it
            //  would put a partial-window mean on disk under the name of a
            //  complete bin -- which a later restart would then recover AS
            //  complete. The straddled bin is dropped, not written.
            //  centreOffset (0 or binDeltaT/2) shifts the window-closed test
            //  to the CENTRED window end when binCentred is active, so a bin
            //  whose centred window closes at 'now' is treated as complete.
            template<class Type>
            static void flushClosedBins
            (
                binItem<Type>& bi,
                const Time& time,
                const scalar cycleStartTime,
                const scalar binDeltaT,
                const scalar now,
                const label skipBin,
                const scalar centreOffset
            );
        };


        //- The restart-reconstruction concern (v2 separation): the mirror
        //  of writer -- rebuilds a binItem's components FROM the files the
        //  writer emitted. Neither function computes a statistic; both only
        //  invert the writer's serialisation. The checkpoint COUNTERS stay
        //  with readCheckpoint()/write(), which orchestrate these calls.
        //  Definitions in BBPAStateTemplates.C.
        struct state
        {
            //- "bin" mode: restore the in-progress cycle's completed bins
            //  (integ_) from the <field>_Bin / <field>_Bin_UU written before
            //  a kill, so a restart resumes mid-cycle instead of discarding
            //  it. Returns the number of bins recovered.
            template<class Type>
            static label restoreClosedBins
            (
                binItem<Type>& bi,
                const fvMesh& mesh,
                const scalar cycleStartTime,
                const scalar binDeltaT,
                const scalar now,
                const scalar centreOffset
            );

            //- Rebuild the cross-cycle ensemble (stats_) from the public
            //  phase-aligned _PA / _PAVariance / _PA_UU files.
            //  paSuffix must match what the writer emitted, or the restore
            //  silently finds nothing and every bin cold-starts.
            template<class Type>
            static void restoreStatistics
            (
                binItem<Type>& bi,
                const word& fieldName,
                const scalar cycleStartTime,
                const scalar binDeltaT,
                const word& paSuffix
            );
        };

        //- Hash tables of bin items for different field types.
        //  BBPA tracks volScalarFields (e.g. pressure) and volVectorFields
        //  (e.g. velocity, wall shear stress). Tensor-valued fields are
        //  not supported because the second-moment accumulator phi*phi is
        //  not defined for generic tensor fields in OpenFOAM.
        HashTable<autoPtr<binItem<scalar>>> scalarBinItems_;
        HashTable<autoPtr<binItem<vector>>> vectorBinItems_;


    // Private Member Functions

        //- Calculate bin index for current time
        label phaseBin() const;

        //- Calculate current cycle number
        label currentCycleNumber() const;

        //- Initialize bin items for all fields.
        //  Returns true if any binItem was newly created by this call, so
        //  the caller can re-apply the checkpoint to a field that only
        //  registered now (see readCheckpoint()).
        bool initializeBinItems();

        //- Reset all bins for new cycle
        void resetAllBinsForNewCycle();

        //- Finalize all bin items for completed cycle.
        //  skipBin (-1 = none) drops the restart-straddled bin.
        void finalizeAllForCycle(const label skipBin = -1);

        //- Accumulate all fields for current time step
        void accumulateAllFields(const label binI, const scalar dt);

        //- "bin" mode: flush the CURRENT cycle's newly-closed bins to disk,
        //  so a kill loses at most the bin it lands in rather than the whole
        //  in-progress cycle.
        void flushCurrentCycleBins(const label cycle);

        //- "bin" mode: recover the in-progress cycle's completed bins written
        //  before a kill. Returns bins recovered.
        label restoreCurrentCycle();

        //- Write all phase-averaged fields (Strategy B: I dirs).
        //  emitPA=false ("bin" mode with PAwriteInterval>1) writes only the
        //  per-cycle <field>_Bin and skips the cross-cycle _BinPA* files.
        void writeAllAveragedFields(const bool emitPA);

        //- Flag the bin straddled by this run's start instant, if any.
        //  Runs for BOTH a restart and a cold start: with binCentred, bin 0 of
        //  the first cycle wants [-dT/2, +dT/2) but the run begins at
        //  startTime, so half its window never existed -- exactly the same
        //  partial-window-folded-as-a-whole-cycle defect as a mid-bin restart.
        //  Must be called only after startTime_/phaseOrigin_ are final.
        void flagStraddledBin();

        //- Restore checkpoint state on restart from disk.
        //  Reads <name>Properties (counters: startTime, currentCycle,
        //  previousCycle, per-field nCycles + binCounts) and the
        //  phase-aligned _PA / _PAVariance / _PA_UU fields written
        //  at the most recent solver writeTime, reconstructing the
        //  internal stats_.mean_, stats_.M2_ and stats_.meanSqr_ arrays.
        //  No-op when no Properties file is found (= cold start).
        void readCheckpoint();

        //- Companion write: one bin (current phase) into solver's
        //  current time directory, for all tracked fields
        void writeCompanionAllFields();


public:

    //- Runtime type information
    TypeName("BBPA");


    // Constructors

        //- Construct from Time and dictionary
        BBPA
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        BBPA(const BBPA&) = delete;


    //- Destructor
    virtual ~BBPA();


    // Member Functions

        //- Read the BBPA data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return fieldNames_;
        }

        //- Execute, performs the phase averaging
        virtual bool execute();

        //- Write the phase-averaged fields
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const BBPA&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BBPATemplates.C"
    #include "BBPAWriterTemplates.C"
    #include "BBPAStateTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
