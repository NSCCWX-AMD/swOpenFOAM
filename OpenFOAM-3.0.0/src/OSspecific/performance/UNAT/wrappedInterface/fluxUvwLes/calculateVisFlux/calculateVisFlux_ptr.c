#include "calculateVisFlux.h"
#include "swMacro.h"
#include <assert.h>

define_e2v_FunPtr(swCalculateVisFlux)
{
	swFloat *facn  = accessArray(frontEdgeData,0);
	swFloat *visSy = accessArray(frontEdgeData,1);
	swFloat *phi   = accessArray(data->vertexData,0);
	swFloat *S     = accessArray(data->vertexData,1);

	swInt edgeNumber = getArraySize(frontEdgeData);
	swInt dims = getArrayDims(frontEdgeData, 0);
	int iedge,iDim;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			faceVal = (facn[iedge*dims+iDim]
				* phi[startVertices[iedge]*dims+iDim]
				+ (1-facn[iedge*dims+iDim])
				* phi[endVertices[iedge]*dims+iDim])*visSy[iedge*dims+iDim];
			S[startVertices[iedge]*dims+iDim] += faceVal;
			S[endVertices[iedge]*dims+iDim]   -= faceVal;
		}
	}

}

define_e2v_FunPtr(swCalculateVisFlux_multiply2)
{
	swFloat *facn  = accessArray(frontEdgeData,0);
	swFloat *visSy = accessArray(frontEdgeData,1);
	swFloat *phi   = accessArray(data->vertexData,0);
	swFloat *S     = accessArray(data->vertexData,1);

	swInt edgeNumber = getArraySize(frontEdgeData);
	swInt dims = getArrayDims(frontEdgeData, 0);
	int iedge,iDim;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			faceVal = 2*(facn[iedge*dims+iDim]
				* phi[startVertices[iedge]*dims+iDim]
				+ (1-facn[iedge*dims+iDim])
				* phi[endVertices[iedge]*dims+iDim])*visSy[iedge*dims+iDim];
			S[startVertices[iedge]*dims+iDim] += faceVal;
			S[endVertices[iedge]*dims+iDim]   -= faceVal;
		}
	}

}
