#include "integrate.h"
#include "swMacro.h"
#include "iterator.h"

//define 2 function pointers
define_e2v_FunPtr(swSurfaceIntegrate)
{
	swInt iDim,dims;
	swFloat* b = accessArray(vertexData,0);
	//frontEdge computation
	swFloat* upper	= accessArray(frontEdgeData, 0);
	swInt edgeNumber = getArraySize(frontEdgeData);
	dims = getArrayDims(frontEdgeData, 0);
	swInt iedge;
	for( iedge = 0; iedge < edgeNumber; iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			b[startVertices[iedge]*dims+iDim]
				+= upper[iedge*dims+iDim];
		}
	}

	//backEdge computation
	swFloat* lower	= accessArray(backEdgeData, 0);
	dims = getArrayDims(backEdgeData, 0);
	edgeNumber = getArraySize( backEdgeData );
	for( iedge = 0; iedge < edgeNumber; iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			b[endVertices[iedge]*dims+iDim]
				-= lower[iedge*dims+iDim];	
		}
	}
}


