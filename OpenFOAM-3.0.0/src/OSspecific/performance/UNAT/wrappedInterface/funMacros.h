/*
* @Author: Gao Fei
* @Date:   2018-12-07 15:22:32
* @Last Modified by:   Gao Fei
* @Last Modified time: 2018-12-07 16:44:10
*/

#ifndef funMacros_H
#define funMacros_H

#include "iterator.H"

// fvcSnGrad
extern "C"
{
	//define 2 function pointers
	define_e2v_hostFunPtr(fvcSnGrad_host)
	{
		swFloat* ssf = accessArray(frontEdgeData, 0);
		swFloat* deltaCoeffs = accessArray(frontEdgeData, 1);
		swFloat* vf  = accessArray(vertexData, 0);
		swInt* own = startVertices;
		swInt* nei = endVertices;

		swInt faceNum = getArraySize( frontEdgeData ); ??????????????????????????
		swInt* dimsPtr = getArrayDims( frontEdgeData );??????????????????????????
		swInt nCpnt = *(dimsPtr+0);??????????????????????????

		for (swInt i = 0; i < faceNum; ++i)
		{
			for (swInt j = 0; j < nCpnt; ++j)
				ssf[ i*nCpnt + j ]
					= deltaCoeffs[i]
					 * ( vf[ nei[i]*nCpnt + j ] - vf[ own[i]*nCpnt + j ] );
		}
	}

	define_e2v_slaveFunPtr(fvcSnGrad_slave)
	{

	}
}

#endif

