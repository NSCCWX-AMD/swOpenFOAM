#include "calculateUvwFlux.h"
#include "swMacro.h"
#include <assert.h>

define_e2v_FunPtr(swCalculateUvwFlux)
{
	swFloat *massFlux = accessArray(frontEdgeData,0);
	swFloat *gamblend = accessArray(paraData, 0);
	swFloat *u  = accessArray(data->vertexData,0);
	swFloat *Su = accessArray(data->vertexData,1);

	swInt edgeNumber = getArraySize(frontEdgeData);
	swInt dims = getArrayDims(frontEdgeData,0);
	int iedge,iDim;
	swFloat facp,facn,fuci,fuce;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			facp = massFlux[iedge*dims+iDim]>=0.0 ? 
				massFlux[iedge*dims+iDim]:0.0;
			facn = massFlux[iedge*dims+iDim]< 0.0 ? 
				massFlux[iedge*dims+iDim]:0.0;
			fuci = facn*u[startVertices[iedge]*dims+iDim]
				+ facp*u[endVertices[iedge]*dims+iDim];
			fuce = facp*u[startVertices[iedge]*dims+iDim]
				+ facn*u[endVertices[iedge]*dims+iDim];
			Su[endVertices[iedge]*dims+iDim]   += (fuci - fuce)*gamblend[0];
			Su[startVertices[iedge]*dims+iDim] -= (fuci - fuce)*gamblend[0];
		}
	}

}
