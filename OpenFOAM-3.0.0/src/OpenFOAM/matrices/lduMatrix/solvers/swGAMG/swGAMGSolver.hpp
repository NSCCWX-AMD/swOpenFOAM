/*---------------------------------------------------------------------------*\
Copyright (C) 2011 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of CAELUS.

    CAELUS is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CAELUS is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CAELUS.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::SW_GAMGSolver

Description
    Geometric agglomerated algebraic multigrid solver.

  Characteristics:
      - Requires positive definite, diagonally dominant matrix.
      - Agglomeration algorithm: selectable and optionally cached.
      - Restriction operator: summation.
      - Prolongation operator: injection.
      - Smoother: Gauss-Seidel.
      - Coarse matrix creation: central coefficient: summation of fine grid
        central coefficients with the removal of intra-cluster face;
        off-diagonal coefficient: summation of off-diagonal faces.
      - Coarse matrix scaling: performed by correction scaling, using steepest
        descent optimisation.
      - Type of cycle: V-cycle with optional pre-smoothing.
      - Coarsest-level matrix solved using ICCG or BICCG.

SourceFiles
    GAMGSolver.cpp
    GAMGSolverCalcAgglomeration.cpp
    GAMGSolverMakeCoarseMatrix.cpp
    GAMGSolverOperations.cpp
    GAMGSolverSolve.cpp

\*---------------------------------------------------------------------------*/

#ifndef SW_GAMGSolver_H
#define SW_GAMGSolver_H

#include "GAMGAgglomeration.H"
#include "lduMatrix.H"
#include "labelField.H"
#include "primitiveFields.H"
#include "LUscalarMatrix.H"
#include "GAMGSolver.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "amulMacros.h"
#include "swAmul.h"
#define MAXCELLS (5000)

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class GAMGSolver Declaration
\*---------------------------------------------------------------------------*/



class SW_GAMGSolver : public GAMGSolver {
    private:
        int MAX_SW_USING_CORASE_LEVELS;
        int init_MAX_SW_USING_CORASE_LEVELS( );
    protected:
        static int if_first ;

//  static int **  _origin_location;
  static refilltion *_refill ;
  public:
    static amul_translate_array **_matrix_traslate;
    static amul_para *_amul_parameter ;
    static int _coarseLevel;

    static int* upperAddr_int32;
    static int* lowerAddr_int32;
    static int** upperAddrLevels_int32;
    static int** lowerAddrLevels_int32;

    //- determine using smoother or slover in the coarest level
    bool usingSmootherCoarest_;
    TypeName("swGAMG");

    SW_GAMGSolver (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
    );

    virtual ~SW_GAMGSolver();
    virtual solverPerformance solve
        (
            scalarField& psi,
            const scalarField& source,
            const direction cmpt=0
        ) const;
    void Vcycle
    (
        const PtrList<lduMatrix::smoother>& smoothers,
        scalarField& psi,
        const scalarField& source,
        scalarField& Apsi,
        scalarField& finestCorrection,
        scalarField& finestResidual,
        PtrList<scalarField>& coarseCorrFields,
        PtrList<scalarField>& coarseSources,
        const direction cmpt
    ) const;
    void initVcycle
    (
        PtrList<scalarField>& coarseCorrFields,
        PtrList<scalarField>& coarseSources,
        PtrList<lduMatrix::smoother>& smoothers
    ) const;

    scalar scalingFactor
    (
        scalarField& Acf,
        const lduMatrix& A,
        scalarField& field,
        const FieldField<Field, scalar>& interfaceLevelBouCoeffs,
        const lduInterfaceFieldPtrsList& interfaceLevel,
        const scalarField& source,
        const direction cmpt,
        int coarseLevel
    ) const;
    private:
    void SW_solveCoarsestLevel
    (
        scalarField& coarsestCorrField,
        const scalarField& coarsestSource
    ) const;

};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
