#ifndef U_TO_AKH_H
#define U_TO_AKH_H
#include "swMacro.h"
#include "iterator.h"

//using namespace UNAT;
#ifdef __cplusplus
extern "C"
{
#endif
	//define 2 function pointers
define_e2v_FunPtr(u_to_akh);
define_e2v_slaveFunPtr(u_to_akh);
//void slave_integrate(Arrays* backEdgeData, Arrays* frontEdgeData,
//			Arrays* selfConnData, Arrays* vertexData, swInt* startVertices,
//			swInt* endVertices);
//void integrate(Arrays* backEdgeData, Arrays* frontEdgeData,
//			Arrays* selfConnData, Arrays* vertexData, swInt* startVertices,
//			swInt* endVertices);
#ifdef __cplusplus
}
#endif

#endif
