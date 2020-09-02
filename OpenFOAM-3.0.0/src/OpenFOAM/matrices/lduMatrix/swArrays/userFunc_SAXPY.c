#include <math.h>
#include "slave.h"
#include "simd.h"
#include "userFunc_SAXPY.h"

void userFunc_aEbPk1Mua(MVM_Arrays *MVM_ArraysPtr)
{
	// a = b + k1*a
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	SCALAR k1 = MVM_ArraysPtr->k1;
	LABEL i, n = MVM_ArraysPtr->nCells;
	for(i=0; i<n; i++)
		A1Ptr[i] = A2Ptr[i] + k1*A1Ptr[i];
}


void userFunc_aEaMik1Mub(MVM_Arrays *MVM_ArraysPtr)
{
	// a -= k1*b
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	SCALAR k1 = MVM_ArraysPtr->k1;
	LABEL i, n = MVM_ArraysPtr->nCells;

	for(i=0; i<n; i++)
		A1Ptr[i] -= k1*A2Ptr[i];
}

void userFunc_aEaPk1Mub(MVM_Arrays *MVM_ArraysPtr)
{
	// a += k1*b
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	SCALAR k1 = MVM_ArraysPtr->k1;
	LABEL i, n = MVM_ArraysPtr->nCells;

	for(i=0; i<n; i++)
		A1Ptr[i] += k1*A2Ptr[i];
}

void userFunc_aEbMic(MVM_Arrays *MVM_ArraysPtr)
{
	// a = b - c
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	SCALAR *A3Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	A3Ptr = MVM_ArraysPtr->A3Ptr;
	LABEL i, n = MVM_ArraysPtr->nCells;

	for(i=0; i<n; i++)
		A1Ptr[i] = A2Ptr[i] - A3Ptr[i];
}

void userFunc_aEk1MuaPk2MubMuScMidS(MVM_Arrays *MVM_ArraysPtr)
{
	// a = k1*a + k2*b*(c-d)
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	SCALAR *A3Ptr;
	SCALAR *A4Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	A3Ptr = MVM_ArraysPtr->A3Ptr;
	A4Ptr = MVM_ArraysPtr->A4Ptr;
	SCALAR k1 = MVM_ArraysPtr->k1;
	SCALAR k2 = MVM_ArraysPtr->k2;
	LABEL i, n = MVM_ArraysPtr->nCells;

	for(i=0; i<n; i++)
		A1Ptr[i] = k1* A1Ptr[i] + k2 * A2Ptr[i] * (A3Ptr[i] - A4Ptr[i]);
}

void userFunc_aE1Db(MVM_Arrays *MVM_ArraysPtr)
{
	// a = 1/b
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	LABEL i, n = MVM_ArraysPtr->nCells;
	for(i=0; i<n; i++)
		A1Ptr[i] = 1 / A2Ptr[i];
}

void userFunc_aEbMuc(MVM_Arrays *MVM_ArraysPtr)
{
	// a = b * c
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	SCALAR *A3Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	A3Ptr = MVM_ArraysPtr->A3Ptr;
	LABEL i, n = MVM_ArraysPtr->nCells;

	for(i=0; i<n; i++)
		A1Ptr[i] = A2Ptr[i] * A3Ptr[i];
}

void userFunc_aEbDc(MVM_Arrays *MVM_ArraysPtr)
{
	// a = b / c
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	SCALAR *A3Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	A3Ptr = MVM_ArraysPtr->A3Ptr;
	LABEL i, n = MVM_ArraysPtr->nCells;

	for(i=0; i<n; i++)
		A1Ptr[i] = A2Ptr[i] / A3Ptr[i];
}

void userFunc_aEbMuScMiaSMuk1(MVM_Arrays *MVM_ArraysPtr)
{
	// a = b * (c - a) * k1
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	SCALAR *A3Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	A3Ptr = MVM_ArraysPtr->A3Ptr;
	LABEL i, n = MVM_ArraysPtr->nCells;
	SCALAR k1 = MVM_ArraysPtr->k1;

	for(i=0; i<n; i++)
		A1Ptr[i] = A2Ptr[i] * (A3Ptr[i] - A1Ptr[i]) * k1;
}

void userFunc_aEaPb(MVM_Arrays *MVM_ArraysPtr)
{
	// a += b
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	LABEL i, n = MVM_ArraysPtr->nCells;
	for(i=0; i<n; i++)
		A1Ptr[i] += A2Ptr[i];
}

void userFunc_jacobi(MVM_Arrays *MVM_ArraysPtr)
{
	// a += (b-c)/d
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	SCALAR *A3Ptr;
	SCALAR *A4Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	A3Ptr = MVM_ArraysPtr->A3Ptr;
	A4Ptr = MVM_ArraysPtr->A4Ptr;
	LABEL i, n = MVM_ArraysPtr->nCells;

    LABEL nsimd = (n / 4) * 4;

    doublev4 a1, a2, a3, a4;
    for(i=0; i<nsimd; i+=4)
	{
		simd_load(a1, &(A1Ptr[i]));
		simd_load(a2, &(A2Ptr[i]));
		simd_load(a3, &(A3Ptr[i]));
		simd_load(a4, &(A4Ptr[i]));

		a1 += (a2 - a3) / a4;

		simd_store(a1, &(A1Ptr[i]));
	}


	for(i=nsimd; i<n; i++)
	{
		A1Ptr[i] += (A2Ptr[i] - A3Ptr[i]) / A4Ptr[i];
	}

	// for(i=0; i<n; i++)
	// 	A1Ptr[i] += (A2Ptr[i] - A3Ptr[i]) / A4Ptr[i];
}

void userFunc_aEk1Mua(MVM_Arrays *MVM_ArraysPtr)
{
	// a = k1*b
	SCALAR *A1Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	SCALAR k1 = MVM_ArraysPtr->k1;
	LABEL i, n = MVM_ArraysPtr->nCells;

	for(i=0; i<n; i++)
		A1Ptr[i] = k1*A1Ptr[i];
}

void userFunc_aEbPk1MuSaMik2MucS(MVM_Arrays *MVM_ArraysPtr)
{
	// a = b + k1*(a - k2*c)
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	SCALAR *A3Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	A3Ptr = MVM_ArraysPtr->A3Ptr;
	LABEL i, n = MVM_ArraysPtr->nCells;
	SCALAR k1 = MVM_ArraysPtr->k1;
	SCALAR k2 = MVM_ArraysPtr->k2;

	for(i=0; i<n; i++)
		A1Ptr[i] = A2Ptr[i] + k1 * (A1Ptr[i] - k2 * A3Ptr[i]);
}


void userFunc_aEbMik1Muc(MVM_Arrays *MVM_ArraysPtr)
{
	// a = b - k1*c
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	SCALAR *A3Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	A3Ptr = MVM_ArraysPtr->A3Ptr;
	SCALAR k1 = MVM_ArraysPtr->k1;
	LABEL i, n = MVM_ArraysPtr->nCells;
	for(i=0; i<n; i++)
		A1Ptr[i] = A2Ptr[i] - k1*A3Ptr[i];
}


void userFunc_aEaPk1MubPk2Muc(MVM_Arrays *MVM_ArraysPtr)
{
	// a += k1*b + k2*c
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	SCALAR *A3Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	A3Ptr = MVM_ArraysPtr->A3Ptr;
	SCALAR k1 = MVM_ArraysPtr->k1;
	SCALAR k2 = MVM_ArraysPtr->k2;
	LABEL i, n = MVM_ArraysPtr->nCells;
	for(i=0; i<n; i++)
		A1Ptr[i] += k1 * A2Ptr[i] + k2 * A3Ptr[i];
}


void userFunc_residualSum(MVM_Arrays *MVM_ArraysPtr)
{
	// a = b - c
	// k = Sum(fabs(a))
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	SCALAR *A3Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	A3Ptr = MVM_ArraysPtr->A3Ptr;
	SCALAR *k1Ptr = MVM_ArraysPtr->k1Ptr;
	LABEL i, n = MVM_ArraysPtr->nCells;

	for(i=0; i<n; i++)
	{
		A1Ptr[i]  = A2Ptr[i] - A3Ptr[i];
		//k1Ptr[0] += fabs(A1Ptr[i]);
		k1Ptr[0] += fabs(A2Ptr[i] - A3Ptr[i]);
	}
}


void userFunc_residualSumK(MVM_Arrays *MVM_ArraysPtr)
{
	// a = b - k1 * c
	// k = Sum(fabs(a))
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	SCALAR *A3Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	A3Ptr = MVM_ArraysPtr->A3Ptr;
	SCALAR k1 = MVM_ArraysPtr->k1;
	SCALAR *k1Ptr = MVM_ArraysPtr->k1Ptr;
	LABEL i, n = MVM_ArraysPtr->nCells;

	for(i=0; i<n; i++)
	{
		A1Ptr[i]  = A2Ptr[i] - k1 * A3Ptr[i];
		//k1Ptr[0] += fabs(A1Ptr[i]);
		k1Ptr[0] += fabs(A2Ptr[i] - k1 * A3Ptr[i]);
	}
}


void userFunc_digPrecondSum(MVM_Arrays *MVM_ArraysPtr)
{
	// a = b * c
	// k = Sum(a * c)
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	SCALAR *A3Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	A3Ptr = MVM_ArraysPtr->A3Ptr;
	SCALAR *k1Ptr = MVM_ArraysPtr->k1Ptr;
	LABEL i, n = MVM_ArraysPtr->nCells;

	for(i=0; i<n; i++)
	{
		A1Ptr[i]  = A2Ptr[i] * A3Ptr[i];
		//k1Ptr[0] += A1Ptr[i] * A3Ptr[i];
		k1Ptr[0] += A2Ptr[i]*A3Ptr[i]*A3Ptr[i];
	}
}


void userFunc_sumProd(MVM_Arrays *MVM_ArraysPtr)
{
	// k = Sum(b * c)
	SCALAR *A2Ptr;
	SCALAR *A3Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	A3Ptr = MVM_ArraysPtr->A3Ptr;
	SCALAR *k1Ptr = MVM_ArraysPtr->k1Ptr;
	LABEL i, n = MVM_ArraysPtr->nCells;

	for(i=0; i<n; i++)
	{
		k1Ptr[0] += A2Ptr[i] * A3Ptr[i];
	}
}


void userFunc_sumSqr(MVM_Arrays *MVM_ArraysPtr)
{
	// k = Sum(b * b)
	SCALAR *A2Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	SCALAR *k1Ptr = MVM_ArraysPtr->k1Ptr;
	LABEL i, n = MVM_ArraysPtr->nCells;

	for(i=0; i<n; i++)
	{
		k1Ptr[0] += A2Ptr[i] * A2Ptr[i];
	}
}


void userFunc_aEcPk1Mua_bEdPk1Mub(MVM_Arrays *MVM_ArraysPtr)
{
	// a = c + k1*a
	// b = d + k1*b
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	SCALAR *A3Ptr;
	SCALAR *A4Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	A3Ptr = MVM_ArraysPtr->A3Ptr;
	A4Ptr = MVM_ArraysPtr->A4Ptr;
	SCALAR k1 = MVM_ArraysPtr->k1;
	LABEL i, n = MVM_ArraysPtr->nCells;
	for(i=0; i<n; i++)
	{
		A1Ptr[i] = A3Ptr[i] + k1*A1Ptr[i];
		A2Ptr[i] = A4Ptr[i] + k1*A2Ptr[i];
	}
}


void userFunc_aEaPk1Muc_bEbMik1Mud(MVM_Arrays *MVM_ArraysPtr)
{
	// a = a + k1*c
	// b = b - k1*d
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	SCALAR *A3Ptr;
	SCALAR *A4Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	A3Ptr = MVM_ArraysPtr->A3Ptr;
	A4Ptr = MVM_ArraysPtr->A4Ptr;
	SCALAR k1 = MVM_ArraysPtr->k1;
	LABEL i, n = MVM_ArraysPtr->nCells;
	for(i=0; i<n; i++)
	{
		A1Ptr[i] += k1*A3Ptr[i];
		A2Ptr[i] -= k1*A4Ptr[i];
	}
}


void userFunc_aEcMuSdMiaSMuk1_bEbPa(MVM_Arrays *MVM_ArraysPtr)
{
	// a = c * (d - a) * k1
	// b = b + a
	SCALAR *A1Ptr;
	SCALAR *A2Ptr;
	SCALAR *A3Ptr;
	SCALAR *A4Ptr;
	A1Ptr = MVM_ArraysPtr->A1Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	A3Ptr = MVM_ArraysPtr->A3Ptr;
	A4Ptr = MVM_ArraysPtr->A4Ptr;
	SCALAR k1 = MVM_ArraysPtr->k1;
	LABEL i, n = MVM_ArraysPtr->nCells;
	for(i=0; i<n; i++)
	{
		A2Ptr[i] += A3Ptr[i] * (A4Ptr[i] - A1Ptr[i]) * k1;
		A1Ptr[i]  = A3Ptr[i] * (A4Ptr[i] - A1Ptr[i]) * k1;
		//A2Ptr[i] += A1Ptr[i];
	}
}



void userFunc_sumFabs(MVM_Arrays *MVM_ArraysPtr)
{
	// k = SumMag(b)
	SCALAR *A2Ptr;
	A2Ptr = MVM_ArraysPtr->A2Ptr;
	SCALAR *k1Ptr = MVM_ArraysPtr->k1Ptr;
	LABEL i, n = MVM_ArraysPtr->nCells;

	for(i=0; i<n; i++)
	{
		k1Ptr[0] += fabs(A2Ptr[i]);
	}
}
