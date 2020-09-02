#include "printMatrixVector.hpp"
#include <fstream>
#include <iomanip>
#include <vector>
#include "processorFvPatchField.H"
#include "processorGAMGInterfaceField.H"
#include <sstream>

using namespace std;

void Foam::printLDUMatrix(const lduMatrix& A, const char* fileName)
{
	const label nFaces = A.upper().size();
	const label nCells = A.diag().size();
	const scalar* upperPtr = A.upper().begin();
    const scalar* lowerPtr = A.lower().begin();
    const label*  uPtr     = A.lduAddr().upperAddr().begin();
    const label*  lPtr     = A.lduAddr().lowerAddr().begin();
    const scalar* diagPtr  = A.diag().begin();

    label nnz = nCells + nFaces;
    label symm = 0;
    if(A.symmetric())
    {
    	symm = 1;
    }
    else
    {
    	symm = 0;
    	nnz += nFaces;
    }

	std::ofstream fout;
    FILEOPEN(fout, fileName);
    fout << "nCells: " << nCells << " " << "nnz: " << nnz << " " << "symm: " << symm << std::endl;

    //- diagonal
    for(int i=0; i<nCells; i++)
    {
        fout << i << " " << i << " " << setiosflags(ios::scientific) << setprecision(15) << diagPtr[i] << std::endl;
    }

    //- upper
    for(int i=0; i<nFaces; i++)
    {
        fout << lPtr[i] << " " << uPtr[i] << " "
             << setiosflags(ios::scientific) << setprecision(15)
             << upperPtr[i] << " ";

        //- if there is lower
        if(!symm)
        {
            fout << setiosflags(ios::scientific) << setprecision(15) << lowerPtr[i];
        }

        fout << std::endl;
    }

    FILECLOSE(fout, fileName);
}


void Foam::printInterfaces
(
    const lduInterfaceFieldPtrsList& interfaces,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const char* fileName
)
{
    label numInterfaces = 0;
    std::vector<label> interfacesUse;

    forAll(interfaces, interfaceI)
    {
        if (interfaces.set(interfaceI))
        {
            numInterfaces++;
            interfacesUse.push_back(interfaceI);
        }
    }

    label* locPosition = new label[numInterfaces+1];
    label* destRank = new label[numInterfaces];

    locPosition[0] = 0;

    for(int i=0; i<numInterfaces; ++i)
    {
        label interfaceI = interfacesUse[i];
        const lduInterfaceField* bp = &(interfaces[interfaceI]);
        const processorFvPatchField<scalar>* dp = static_cast<const processorFvPatchField<scalar> *>(bp);

        destRank[i] = dp->neighbProcNo();
        locPosition[i+1] = dp->size() + locPosition[i];
    }

    label nFacesOffDiag = locPosition[numInterfaces];
    label* faceCells = new label[nFacesOffDiag];
    scalar* faceCoffOffDiag = new scalar[nFacesOffDiag];

    //- interface coefficients
    for(int i=0; i<numInterfaces; ++i)
    {
        label interfaceI = interfacesUse[i];
        label locSize = locPosition[i+1] - locPosition[i];
        const lduInterfaceField* bp = &(interfaces[interfaceI]);
        const processorFvPatchField<scalar>* dp = static_cast<const processorFvPatchField<scalar> *>(bp);
        const labelUList& faceCellsLocal = dp->patch().faceCells();

        for(label facei=0; facei<locSize; facei++)
        {
            faceCoffOffDiag[locPosition[i] + facei] = interfaceBouCoeffs[interfaceI][facei];
            faceCells[locPosition[i] + facei] = faceCellsLocal[facei];
        }
    }

    std::ofstream fout;
    FILEOPEN(fout, fileName);
    fout << "neighbor processors: " << numInterfaces << " "
         << "number faces: " << nFacesOffDiag << std::endl;

    for(int i=0; i<numInterfaces; ++i)
    {
        label locSize = locPosition[i+1] - locPosition[i];
        fout << "processor: " << destRank[i] << " "
             << "size: " << locSize << std::endl;
        for(label facei=0; facei<locSize; facei++)
        {
            fout << faceCells[locPosition[i] + facei] << " "
                 << setiosflags(ios::scientific) << setprecision(15)
                 << faceCoffOffDiag[locPosition[i] + facei] << std::endl;
        }
    }

    delete [] locPosition;
    delete [] destRank;
    delete [] faceCells;
    delete [] faceCoffOffDiag;
    FILECLOSE(fout, fileName);
}


void Foam::printGAMGInterfaces
(
    const lduInterfaceFieldPtrsList& interfaces,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const char* fileName
)
{
    label numInterfaces = 0;
    std::vector<label> interfacesUse;

    forAll(interfaces, interfaceI)
    {
        if (interfaces.set(interfaceI))
        {
            numInterfaces++;
            interfacesUse.push_back(interfaceI);
        }
    }

    label* locPosition = new label[numInterfaces+1];
    label* destRank = new label[numInterfaces];

    locPosition[0] = 0;

    for(int i=0; i<numInterfaces; ++i)
    {
        label interfaceI = interfacesUse[i];
        const lduInterfaceField* bp = &(interfaces[interfaceI]);
        const processorGAMGInterfaceField *dp = static_cast<const processorGAMGInterfaceField *>(bp);

        destRank[i] = dp->neighbProcNo();
        locPosition[i+1] = dp->size() + locPosition[i];
    }

    label nFacesOffDiag = locPosition[numInterfaces];
    label* faceCells = new label[nFacesOffDiag];
    scalar* faceCoffOffDiag = new scalar[nFacesOffDiag];

    //- interface coefficients
    for(int i=0; i<numInterfaces; ++i)
    {
        label interfaceI = interfacesUse[i];
        label locSize = locPosition[i+1] - locPosition[i];
        const lduInterfaceField* bp = &(interfaces[interfaceI]);
        const processorGAMGInterfaceField* dp = static_cast<const processorGAMGInterfaceField *>(bp);
        const labelUList& faceCellsLocal = dp->faceCells();

        for(label facei=0; facei<locSize; facei++)
        {
            faceCoffOffDiag[locPosition[i] + facei] = interfaceBouCoeffs[interfaceI][facei];
            faceCells[locPosition[i] + facei] = faceCellsLocal[facei];
        }
    }

    std::ofstream fout;
    FILEOPEN(fout, fileName);
    fout << "neighbor processors: " << numInterfaces << " "
         << "number faces: " << nFacesOffDiag << std::endl;

    for(int i=0; i<numInterfaces; ++i)
    {
        label locSize = locPosition[i+1] - locPosition[i];
        fout << "processor: " << destRank[i] << " "
             << "size: " << locSize << std::endl;
        for(label facei=0; facei<locSize; facei++)
        {
            fout << faceCells[locPosition[i] + facei] << " "
                 << setiosflags(ios::scientific) << setprecision(15)
                 << faceCoffOffDiag[locPosition[i] + facei] << std::endl;
        }
    }

    delete [] locPosition;
    delete [] destRank;
    delete [] faceCells;
    delete [] faceCoffOffDiag;
    FILECLOSE(fout, fileName);
}
