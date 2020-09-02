
#ifndef interpolateParameter_H
#define interpolateParameter_H

#ifdef __cplusplus
extern "C"
{
#endif

#include "sunwayMacros.h"

//#include "rowSubsection.h"

typedef struct Parameter
{
	//the pointer of faces value
	char* fvPtr;

	//the pointer of cells value
	char* cvPtr;

	//the pointer of weights
	swFloat* weiPtr;

	//the pointer of 1-weights
	swFloat* _1wPtr;

	//the pointer of owner
	swInt* ownPtr;

	//the pointer of neighbor
	swInt* neiPtr;

	//the pointer of rowSubsection
	struct rowSubsection** subSecPtr;//

	//the number of sections in every sagment
	swInt secNum;

	//the number of element in a type
	swInt eleNum;

	//the size of type
	swInt typeSize;

	//some test information
	swFloat* usage;
}Parameter,*ParameterPtr;

#ifdef __cplusplus
}
#endif

#endif
