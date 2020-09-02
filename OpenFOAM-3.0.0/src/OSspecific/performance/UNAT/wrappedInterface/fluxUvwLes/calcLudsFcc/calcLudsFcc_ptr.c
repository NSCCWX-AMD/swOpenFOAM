#include "calcLudsFcc.h"
#include "swMacro.h"
#include <assert.h>

define_e2v_FunPtr(swCalcLudsFcc)
{
	swFloat *massFlux = accessArray(frontEdgeData,0);
	swFloat *cellx    = accessArray(data->vertexData,0);
	swFloat *facex    = accessArray(data->frontEdgeData,1);
	swFloat *fcc      = accessArray(data->frontEdgeData,0);
	swFloat *gamblend = accessArray(paraData, 0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim;
	swFloat facp, facn;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			facp = massFlux[iedge*dims+iDim] >= 0.0
				? massFlux[iedge*dims+iDim] : 0.0;
			facn = massFlux[iedge*dims+iDim] < 0.0
				? massFlux[iedge*dims+iDim] : 0.0;
			fcc[iedge*dims+iDim] = (facn*(facex[iedge*dims+iDim]
				- cellx[endVertices[iedge*dims+iDim]])
				+ facp*(facex[iedge*dims+iDim]
				- cellx[startVertices[iedge*dims+iDim]]))*gamblend[0];
		}
	}
}

define_e2v_FunPtr(swCalcLudsFcc_rface)
{
	swFloat *massFlux = accessArray(frontEdgeData,0);
	swFloat *rface0   = accessArray(data->frontEdgeData,0);
	swFloat *rface1   = accessArray(data->frontEdgeData,1);
	swFloat *paras    = accessArray(paraData, 0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim;
	swFloat facp, facn;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			facp = massFlux[iedge*dims+iDim] >= 0.0
				? massFlux[iedge*dims+iDim] : 0.0;
			facn = massFlux[iedge*dims+iDim] < 0.0
				? massFlux[iedge*dims+iDim] : 0.0;
			rface0[iedge*dims+iDim] -= facp;
			rface1[iedge*dims+iDim] += facn;
		}
	}
}

define_e2v_FunPtr(swCalcLudsFcc_rface_CDS)
{
	swFloat *facn     = accessArray(frontEdgeData,0);
	swFloat *rface0   = accessArray(data->frontEdgeData,0);
	swFloat *rface1   = accessArray(data->frontEdgeData,1);
	swFloat *paras    = accessArray(paraData, 0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			rface0[iedge*dims+iDim] -= 1.0-facn[iedge*dims+iDim];
			rface1[iedge*dims+iDim] += facn[iedge*dims+iDim];
		}
	}
}

