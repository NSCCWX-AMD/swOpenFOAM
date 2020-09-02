#ifndef SWRESTINTERSTRUCT_H
#define SWRESTINTERSTRUCT_H

#include "swMacro.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    const swInt*    mapPtr;  // restriction map
    const swFloat*  fPtr;    // fine data
    swFloat*  cPtr;		     // coarse data
    swInt**   localStartEnd; // local range of fine and coarse data in each slave core
    swInt     slaveCycles;   // slave cycles (how many cycles of slave cores to store the entire array)
}restStruct;


typedef struct
{
    swInt*    mapPtr;         // interpolation map
    swInt*    offsetMapPtr;   // offset map
    swFloat*  fPtr;           // fine data
    const swFloat*  cPtr;	  // coarse data
    swInt**   localStartEnd;  // local range of fine and coarse data in each slave core
    swInt     slaveCycles;    // slave cycles (how many cycles of slave cores to store the entire array)
}interStruct;



typedef struct
{
    const swInt*    mapPtr;  // restriction map
    const swFloat*  fPtr;    // fine upper data
    swFloat*  cUPtr;         // coarse upper data
    swFloat*  cDPtr;         // coarse diagonal data
    swInt**   localStartEnd; // local range of fine and coarse data in each slave core
    swInt     slaveCycles;   // slave cycles (how many cycles of slave cores to store the entire array)
}aggMatrixUpperStruct;


void restrictData_host(restStruct*);

void SLAVE_FUNC(restrictData_slave)(restStruct*);

void interpolateData_host(interStruct*);

void SLAVE_FUNC(interpolateData_slave)(interStruct*);

#ifdef __cplusplus
}
#endif

#endif
