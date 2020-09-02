#include "u_to_akh.h"
#include <math.h>
#include "swMacro.h"
#include "iterator.h"

//define 2 function pointers
define_e2v_FunPtr(u_to_akh)
{
	swInt iDim,dims;
	swFloat* b = accessArray(vertexData,1);
	swFloat* x = accessArray(vertexData,0);
	//frontEdge computation
	swFloat* upper1	= accessArray(frontEdgeData, 0);
	swFloat* upper2	= accessArray(frontEdgeData, 1);
	swInt edgeNumber = getArraySize(frontEdgeData);
	dims = getArrayDims(frontEdgeData, 0);
	swInt iedge;
	swFloat *u, *n, *v, *uo;
	swFloat un_a,vn,un;
	for( iedge = 0; iedge < edgeNumber; iedge++)
	{
		u = &b[startVertices[iedge]*dims];
		uo = &x[endVertices[iedge]*dims];
//		u[0] = 0.5*(u[0]+uo[0]);
//		u[1] = 0.5*(u[1]+uo[1]);
//		u[2] = 0.5*(u[2]+uo[2]);
		n = &upper1[iedge*dims];
		v = &upper2[iedge*dims];
		swFloat tmp = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
		u[0] = u[0]*tmp;
		u[1] = u[1]*tmp;
		u[2] = u[2]*tmp;
		un_a = u[0]*n[0]+u[1]*n[1]+u[2]*n[2];
		vn   = v[0]*n[0]+v[1]*n[1]+v[2]*n[2];
		un   = un_a-vn;
		un   = u[0]*un;
	}
}
