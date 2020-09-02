#ifndef PRINTMATRIXVECTOR_HPP
#define PRINTMATRIXVECTOR_HPP

#include "lduMatrix.H"

#include <fstream>
#include <iomanip>
#include <sstream>

#define FILEOPEN(fout, fileName) \
std::ostringstream os; \
os << fileName << "_" << UPstream::myProcNo() << ".txt"; \
fout.open(os.str().c_str())

#define FILECLOSE(fout, fileName) fout.close()

namespace Foam
{

void printLDUMatrix(const lduMatrix& A, const char* name);

template<typename T>
void printVector(const T& b, const char* name);

void printInterfaces
(
	const lduInterfaceFieldPtrsList& interfaces,
	const FieldField<Field, scalar>& interfaceBouCoeffs,
	const char* name
);

void printGAMGInterfaces
(
	const lduInterfaceFieldPtrsList& interfaces,
	const FieldField<Field, scalar>& interfaceBouCoeffs,
	const char* name
);


} //- end namespace

template<typename T>
void Foam::printVector(const T& b, const char* fileName)
{
	const label size = b.size();
    std::ofstream fout;
	FILEOPEN(fout, fileName);
    fout << "nCells: " << size << std::endl;

    for(int i=0; i<size; i++)
    {
    	fout << std::setiosflags(std::ios::scientific) << std::setprecision(15) << b[i] << std::endl;
    }
    FILECLOSE(fout, fileName);
}



#endif //- PRINTMATRIXVECTOR_HPP
