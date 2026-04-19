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
    Foam::functionObjects::IPA

Description
    Instantaneous Phase Averaging (IPA) functionObject
    
    Samples instantaneous field values at specific phase points within a period
    and averages them across multiple cycles.
    
    Mathematical formulation:
    - Phase points: φ_i = i * (period / nPhases) for i = 0, 1, ..., nPhases-1
    - Instantaneous sampling: φ_inst,n,i = φ(t_n,i) at time t_n,i
    - Phase averaging: PA_φ,i^(N) = (1/N) Σ[n=1 to N] φ_inst,n,i
    
    where N is the number of completed cycles.

\*---------------------------------------------------------------------------*/

#ifndef IPA_H
#define IPA_H

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
                          Class IPA Declaration
\*---------------------------------------------------------------------------*/

class IPA
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Fields to average
        wordList fieldNames_;

        //- Period of oscillation
        scalar period_;

        //- Number of phase points per period
        label nPhases_;

        //- Starting cycle number (phase averaging begins from this cycle)
        label startCycle_;

        //- Start time for phase averaging
        scalar startTime_;

        //- Phase interval (period / nPhases)
        scalar phaseInterval_;

        //- Current cycle number
        label currentCycle_;

        //- Previous cycle number
        label previousCycle_;

        //- Tolerance for phase matching
        scalar phaseTolerance_;

        //- Phase item structure for each field
        template<class Type>
        struct phaseItem
        {
            typedef GeometricField<Type, fvPatchField, volMesh> FieldType;
            
            word fieldName_;
            label nCycles_;
            
            // Phase-averaged fields (one for each phase point)
            PtrList<FieldType> phaseAveragedFields_;
            
            // Sample counts for each phase point
            labelList phaseSampleCounts_;
            
            // Current cycle samples (for current cycle tracking)
            PtrList<FieldType> currentCycleSamples_;
            List<bool> currentCycleSampled_;
            
            phaseItem(const word& name, const fvMesh& mesh, const label nPhases);
            
            void resetCurrentCycle();
            void sampleAtPhase(const label phaseI, const FieldType& field);
            void finalizeCycle();
            void write(const Time& time) const;
        };

        //- Hash table of phase items for different field types
        HashTable<autoPtr<phaseItem<scalar>>> scalarPhaseItems_;
        HashTable<autoPtr<phaseItem<vector>>> vectorPhaseItems_;
        HashTable<autoPtr<phaseItem<symmTensor>>> symmTensorPhaseItems_;
        HashTable<autoPtr<phaseItem<tensor>>> tensorPhaseItems_;


    // Private Member Functions

        //- Calculate current cycle number
        label currentCycleNumber() const;

        //- Calculate phase index for current time
        label currentPhaseIndex() const;

        //- Check if current time matches a phase point
        bool isAtPhasePoint(label& phaseIndex) const;

        //- Reset all phase items for new cycle
        void resetPhaseItemsForNewCycle();

        //- Finalize all phase items for completed cycle
        void finalizePhaseItemsForCycle();

        //- Sample all fields at current phase point
        void sampleFieldsAtPhase(const label phaseI);

        //- Initialize phase items for all fields
        void initializePhaseItems();

        //- Write all phase-averaged fields
        void writePhaseAveragedFields();


public:

    //- Runtime type information
    TypeName("IPA");


    // Constructors

        //- Construct from Time and dictionary
        IPA
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        IPA(const IPA&) = delete;


    //- Destructor
    virtual ~IPA();


    // Member Functions

        //- Read the IPA data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return fieldNames_;
        }

        //- Execute, performs the instantaneous phase averaging
        virtual bool execute();

        //- Write the phase-averaged fields
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const IPA&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IPATemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //