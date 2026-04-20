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
#include "volFieldsFwd.H"
#include "HashPtrTable.H"
#include "scalarField.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

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

        //- Period of oscillation
        scalar period_;

        //- Number of bins per cycle
        label nBins_;

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
        word writeMode_;

        //- Current cycle number
        label currentCycle_;

        //- Previous cycle number
        label previousCycle_;

        //- Current bin number
        label currentBin_;

        //- Previous bin number
        label previousBin_;

        //- Bin item structure for each field (templated to preserve type)
        template<class Type>
        struct binItem
        {
            typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

            word fieldName_;
            label nCycles_;
            scalarField binCounts_;

            //- Phase-averaged mean fields (one per bin, preserves full type)
            PtrList<FieldType> meanFields_;

            //- Welford M2 accumulator: sum of squared deviations from
            //  running mean (numerically stable variance)
            PtrList<volScalarField> M2Fields_;

            //- Current cycle: total time accumulated in each bin
            scalarField currentCycleBinTime_;

            //- Current cycle: time-weighted sum for each bin
            PtrList<FieldType> currentCycleBinSum_;

            //- Cross-cycle running mean of <|u|^2> per bin. Combined with
            //  meanFields_ gives true resolved turbulent kinetic energy:
            //    TKE(bin) = 0.5 * ( meanUUFields_[bin] - magSqr(meanFields_[bin]) )
            PtrList<volScalarField> meanUUFields_;

            //- Current cycle: time-weighted sum of |u|^2 for each bin
            PtrList<volScalarField> currentCycleBinUUSum_;

            binItem(const word& name, const fvMesh& mesh, const label nBins);

            void resetCurrentCycle();
            void accumulate
            (
                const label binI,
                const FieldType& field,
                const scalar dt
            );
            void finalizeCycle();
            void write
            (
                const Time& time,
                const label nBins,
                const label currentCycle,
                const label currentBin,
                const scalar cycleStartTime,
                const scalar binDeltaT
            ) const;

            //- Companion write: dump only the current-phase bin's PA
            //  field into the solver's current time directory (no
            //  phase-aligned sub-directories). Paired with an
            //  instantaneous snapshot in the same time dir.
            void writeCompanion
            (
                const Time& time,
                const label currentBinIndex
            ) const;
        };

        //- Hash tables of bin items for different field types
        HashTable<autoPtr<binItem<scalar>>> scalarBinItems_;
        HashTable<autoPtr<binItem<vector>>> vectorBinItems_;
        HashTable<autoPtr<binItem<symmTensor>>> symmTensorBinItems_;
        HashTable<autoPtr<binItem<tensor>>> tensorBinItems_;


    // Private Member Functions

        //- Calculate bin index for current time
        label phaseBin() const;

        //- Calculate current cycle number
        label currentCycleNumber() const;

        //- Initialize bin items for all fields
        void initializeBinItems();

        //- Reset all bins for new cycle
        void resetAllBinsForNewCycle();

        //- Finalize all bin items for completed cycle
        void finalizeAllForCycle();

        //- Accumulate all fields for current time step
        void accumulateAllFields(const label binI, const scalar dt);

        //- Write all phase-averaged fields (Strategy B: I dirs)
        void writeAllAveragedFields();

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
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
