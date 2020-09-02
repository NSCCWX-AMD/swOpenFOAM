#include "merge4_interpolation.h"
#include "swMacro.h"
#include <assert.h>

define_e2v_FunPtr(swMerge4_lam)
{
	swFloat *weights    = accessArray(frontEdgeData,0);
	swFloat *cellFields = accessArray(data->vertexData,0);
	swFloat *faceFields = accessArray(data->frontEdgeData,0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			faceFields[idx]
				= weights[idx]*cellFields[startVertices[iedge]*dims+iDim]
				+ (1-weights[idx])
				* cellFields[endVertices[iedge]*dims+iDim];
		}
	}
}

define_e2v_FunPtr(swMerge4_dPhi)
{
	swFloat *dPhidXac = accessArray(data->frontEdgeData,0);
	swFloat *dPhi     = accessArray(frontEdgeData,0);
	swFloat *cell_x   = accessArray(data->vertexData,0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			dPhi[idx] += (cell_x[startVertices[iedge]*dims+iDim]
				- cell_x[endVertices[iedge]*dims+iDim])*dPhidXac[idx];
		}
	}
}

define_e2v_FunPtr(swMerge4_mf)
{
	swFloat *directions = accessArray(frontEdgeData,0);
	swFloat *cellFields = accessArray(data->vertexData,0);
	swFloat *faceFields = accessArray(data->frontEdgeData,0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			if(directions[idx] >= 0.0)
				faceFields[idx]
					= cellFields[endVertices[iedge]*dims+iDim];
			else
				faceFields[idx]
					= cellFields[startVertices[iedge]*dims+iDim];
		}
	}
}


