#ifndef SW_STRUCT_H
#define SW_STRUCT_H
#include "BlockOrdering.h"
//#define SCALAR double
//#define LABEL long int

typedef struct LDUMatrix{
	SCALAR *lower;
	SCALAR *upper;
	SCALAR *diag;
	LABEL *rowAddr;
	LABEL *colAddr;
	LABEL numCell;
	LABEL numEdge;
}Matrix;

#endif
