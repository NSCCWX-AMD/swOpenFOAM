#include "interpolation.h"
#include "swMacro.h"
#include "iterator.h"

//define 2 function pointers
define_e2v_FunPtr(interpolation)
{
	swInt iDim,dims;
	swFloat* x = accessArray(vertexData,0);
	//frontEdge computation
	swFloat* ssf	= accessArray(data->frontEdgeData, 0);
//	swFloat* coef	= accessArray(frontEdgeData, 1);
	swInt edgeNumber = getArraySize(data->frontEdgeData);
	dims = getArrayDims(data->frontEdgeData, 0);
	swInt iedge;
	for( iedge = 0; iedge < edgeNumber; iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			ssf[iedge*dims+iDim]
//				= coef[iedge]
				= ssf[iedge*dims+iDim]
				* (x[endVertices[iedge]*dims+iDim]-x[startVertices[iedge]*dims+iDim]);
//if(iedge==7831) printf("%f,%f,%f\n",ssf[iedge],x[endVertices[iedge]],x[startVertices[iedge]]);
		}
	}
}


