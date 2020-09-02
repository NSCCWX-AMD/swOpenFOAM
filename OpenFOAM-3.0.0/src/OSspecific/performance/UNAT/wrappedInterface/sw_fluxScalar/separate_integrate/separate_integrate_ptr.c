#include "separate_integrate.h"
#include "swMacro.h"
#include <assert.h>

define_e2v_FunPtr(swSeparateIntegrate_su)
{
	swFloat *su = accessArray(data->vertexData,0);
	swFloat *fu = accessArray(data->frontEdgeData,0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			su[startVertices[iedge]*dims+iDim] += fu[idx];
			su[endVertices[iedge]*dims+iDim] -= fu[idx];
		}
	}
}


