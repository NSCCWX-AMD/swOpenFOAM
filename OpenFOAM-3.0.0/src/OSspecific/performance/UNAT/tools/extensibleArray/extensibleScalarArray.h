/*****************************************************************\
 * QDA represents Quasi-Diagonal Arrays, and this file contains
 * the macros to define subrutines on Sunway plat-forms to accele-
 * rate Caelus QDA operation kernels.
\* ***************************************************************/

#ifndef extensibleArray_H
#define extensibleArray_H

#include "stdlib.h"
#include "string.h"
#include "swMacro.h"

#ifdef __cplusplus
extern "C"
{
#endif

#define SCALAR swFloat
#define LABEL swInt

#define INITSIZE 32


typedef struct struct_extensibleSCALARArray
{
    SCALAR* data;
    LABEL size;
    LABEL maxSize;
} ExtensibleSCALARArray;

void extensibleSCALARArrayInit
(
 ExtensibleSCALARArray *array
);

void extensibleSCALARArrayAdd
(
 ExtensibleSCALARArray *array,
 SCALAR value
);

void extensibleSCALARArrayDestroy
(
 ExtensibleSCALARArray *array
);


#ifdef __cplusplus
}
#endif

#endif
