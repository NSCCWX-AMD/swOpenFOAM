#include "interpolateViscosity.h"
#include "swMacro.h"
#include <assert.h>

define_e2v_FunPtr(swInterpolateViscosity)
{
	swFloat *facn   = accessArray(frontEdgeData,0);
	swFloat *viseff = accessArray(vertexData,0);
	swFloat *visOut = accessArray(data->frontEdgeData,1);
	swFloat *faceIn = accessArray(data->frontEdgeData,0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim;
	swFloat visac;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			visac = facn[iedge*dims+iDim]
				* viseff[startVertices[iedge]*dims+iDim] 
				+ (1-facn[iedge*dims+iDim])
				* viseff[endVertices[iedge]*dims+iDim];
			visOut[iedge*dims+iDim] = visac*faceIn[iedge*dims+iDim];
		}
	}

}

define_e2v_FunPtr(swInterpolateViscosity_1)
{
	swFloat *facn   = accessArray(frontEdgeData,0);
	swFloat *viseff = accessArray(vertexData,0);
	swFloat *visOut = accessArray(data->frontEdgeData,0);
	swFloat *faceIn = accessArray(frontEdgeData,1);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim;
	swFloat visac;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			visac = facn[iedge*dims+iDim]
				* viseff[startVertices[iedge]*dims+iDim] 
				+ (1-facn[iedge*dims+iDim])
				* viseff[endVertices[iedge]*dims+iDim];
			visOut[iedge*dims+iDim] -= visac*faceIn[iedge*dims+iDim];
		}
	}

}
