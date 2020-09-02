#include "spMV.h"
#include "swMacro.h"
#include "iterator.h"

//define 2 function pointers
define_e2v_FunPtr(swSpMV)
{
	//selfConn computation
	swFloat* diag	= accessArray(data->selfConnData, 0);
	swFloat* x		= accessArray(data->vertexData, 0);
	swFloat* b		= accessArray(data->vertexData, 1);
	swInt iDim,dims;
	
	swInt vertexNum = getArraySize(data->selfConnData);
	dims = getArrayDims(data->selfConnData, 0);
	swInt ivertex;
	for( ivertex = 0; ivertex < vertexNum; ivertex++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
//if(ivertex==895) printf("diag: %f,%f,%f\n",b[ivertex],diag[ivertex],x[ivertex]);
			b[ivertex*dims+iDim]
				+= diag[ivertex*dims+iDim]*x[ivertex*dims+iDim];
//			b[ivertex]
//				+= diag[ivertex]*x[ivertex];
		}
	}
	
	//frontEdge computation
	swFloat* upper	= accessArray(data->frontEdgeData, 0);
	swInt edgeNumber = getArraySize(data->frontEdgeData);
	dims = getArrayDims(data->frontEdgeData, 0);
	swInt iedge;
	for( iedge = 0; iedge < edgeNumber; iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
//if(startVertices[iedge]==33971) printf("owner:%d,%d,%f,%f,%f\n",iedge,endVertices[iedge],b[startVertices[iedge]],upper[iedge],x[endVertices[iedge]]);
			b[startVertices[iedge]*dims+iDim]
				+= upper[iedge*dims+iDim]*x[endVertices[iedge]*dims+iDim];
//			b[startVertices[iedge]]
//				+= upper[iedge]*x[endVertices[iedge]];

		}
	}

	//backEdge computation
	edgeNumber = getArraySize( data->backEdgeData );
	if(edgeNumber<=0) return;
	swFloat* lower	= accessArray(data->backEdgeData, 0);
	dims = getArrayDims(data->backEdgeData, 0);
//	printf("edgeNumber: %d\n",edgeNumber);
	for( iedge = 0; iedge < edgeNumber; iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
//if(endVertices[iedge]==15) printf("neighbor:%d,%d,%f,%f,%f\n",iedge,startVertices[iedge],b[endVertices[iedge]],lower[iedge],x[startVertices[iedge]]);
			b[endVertices[iedge]*dims+iDim]
				+= lower[iedge*dims+iDim]*x[startVertices[iedge]*dims+iDim];	
//			b[endVertices[iedge]]
//				+= lower[iedge]*x[startVertices[iedge]];	

		}
	}
}


