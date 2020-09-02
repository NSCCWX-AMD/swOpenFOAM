#include "calculateFccFlux.h"
#include "swMacro.h"
#include <assert.h>

define_e2v_FunPtr(swCalculateFccFlux)
{
	swFloat *massFlux = accessArray(frontEdgeData,0);
	swFloat *fccx     = accessArray(frontEdgeData,1);
	swFloat *phi      = accessArray(data->vertexData,0);
	swFloat *S        = accessArray(data->vertexData,1);

	swInt edgeNumber = getArraySize(frontEdgeData);
	swInt dims = getArrayDims(frontEdgeData, 0);
	int iedge,iDim;
	swFloat facn,facp,fce;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			facp = massFlux[iedge*dims+iDim] >= 0.0
				? 1.0 : 0.0;
			facn = massFlux[iedge*dims+iDim] < 0.0
				? 1.0 : 0.0;
			fce = fccx[iedge*dims+iDim]
				* (facn*phi[startVertices[iedge]*dims+iDim]
				+ facp*phi[endVertices[iedge]*dims+iDim]);
			S[startVertices[iedge]*dims+iDim] += fce;
			S[endVertices[iedge]*dims+iDim]   -= fce;
		}
	}
}
