#include "separate_interpolation.h"
#include "swMacro.h"
#include <assert.h>

define_e2v_FunPtr(swInterpolation_weighted)
{
	swFloat *weights    = accessArray(frontEdgeData,0);
	swFloat *cellField = accessArray(data->vertexData,0);
	swFloat *faceField = accessArray(data->frontEdgeData,0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			faceField[idx]
				= weights[idx]*cellField[startVertices[iedge]*dims+iDim]
				+ (1-weights[idx])
				* cellField[endVertices[iedge]*dims+iDim];
		}
	}
}

define_e2v_FunPtr(swInterpolation_weightedSwap)
{
	swFloat *weights    = accessArray(frontEdgeData,0);
	swFloat *cellField = accessArray(data->vertexData,0);
	swFloat *faceField = accessArray(data->frontEdgeData,0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			faceField[idx]
				= (1-weights[idx])*cellField[startVertices[iedge]*dims+iDim]
				+ weights[idx]
				* cellField[endVertices[iedge]*dims+iDim];
		}
	}
}

define_e2v_FunPtr(swInterpolation_directed)
{
	swFloat *direction = accessArray(frontEdgeData,0);
	swFloat *cellField = accessArray(data->vertexData,0);
	swFloat *faceField = accessArray(data->frontEdgeData,0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			if(direction[idx]>=0.0)
				faceField[idx] = cellField[startVertices[iedge]*dims+iDim];
			else
				faceField[idx] = cellField[endVertices[iedge]*dims+iDim];
		}
	}
}

define_e2v_FunPtr(swInterpolation_constant_1)
{
	swFloat *cellField = accessArray(data->vertexData,0);
	swFloat *faceField = accessArray(data->frontEdgeData,0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			faceField[idx] = 1.0 * cellField[endVertices[iedge]*dims+iDim]
				- 1.0 * cellField[startVertices[iedge]*dims+iDim];
		}
	}
}

define_e2v_FunPtr(swInterpolation_constant_0_5)
{
	swFloat *cellField = accessArray(data->vertexData,0);
	swFloat *faceField = accessArray(data->frontEdgeData,0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			faceField[idx] = 0.5 * cellField[endVertices[iedge]*dims+iDim]
				+ 0.5 * cellField[startVertices[iedge]*dims+iDim];
		}
	}
}


