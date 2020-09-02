/*****************************************************************\
 * QDA represents Quasi-Diagonal Arrays, and this file contains
 * the macros to define subrutines on Sunway plat-forms to accele-
 * rate Caelus QDA operation kernels.
\* ***************************************************************/

#include "extensibleScalarArray.h"

void extensibleSCALARArrayInit
(
 ExtensibleSCALARArray *array
)
{
    //if(array->data != NULL) free(array->data);
    array->data = (SCALAR*) malloc( sizeof(SCALAR)*INITSIZE );
    array->size = 0;
    array->maxSize = INITSIZE;
}

void extensibleSCALARArrayAdd
(
 ExtensibleSCALARArray *array,
 SCALAR value
)
{
    if(array->size < array->maxSize)
    {
        array->data[array->size] = value;
        (array->size)++;
    }
    else
    {
        SCALAR* tmpData =
			(SCALAR*)malloc( sizeof(SCALAR)*2*(array->maxSize));
        array->maxSize *= 2;
        memcpy(tmpData, array->data, sizeof(SCALAR)*(array->size) );
        free(array->data);
        array->data = tmpData;
        array->data[array->size] = value;
        (array->size)++;
    }
}

void extensibleSCALARArrayDestroy
(
 ExtensibleSCALARArray *array
)
{
    free(array->data);
}
