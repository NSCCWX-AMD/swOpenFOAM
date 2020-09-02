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

#define LABEL swInt

#define INITSIZE 32


typedef struct struct_extensibleLABELArray
{
    LABEL* data;
    LABEL size;
    LABEL maxSize;
} ExtensibleLABELArray;

void extensibleLABELArrayInit
(
 ExtensibleLABELArray *array
);

void extensibleLABELArrayAdd
(
 ExtensibleLABELArray *array,
 LABEL value
);

void extensibleLABELArrayDestroy
(
 ExtensibleLABELArray *array
);


#ifdef __cplusplus
}
#endif

#endif
