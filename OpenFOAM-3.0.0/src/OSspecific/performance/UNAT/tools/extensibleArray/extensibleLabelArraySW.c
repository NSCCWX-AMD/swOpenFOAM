/*****************************************************************\
 * QDA represents Quasi-Diagonal Arrays, and this file contains
 * the macros to define subrutines on Sunway plat-forms to accele-
 * rate Caelus QDA operation kernels.
\* ***************************************************************/

#include "extensibleLabelArray.h"

void extensibleLABELArrayInit
(
 ExtensibleLABELArray *array
)
{
    //if(array->data != NULL) free(array->data);
    array->data = (LABEL*) malloc( sizeof(LABEL)*INITSIZE );
    array->size = 0;
    array->maxSize = INITSIZE;
}

void extensibleLABELArrayAdd
(
 ExtensibleLABELArray *array,
 LABEL value
)
{
    if(array->size < array->maxSize)
    {
        array->data[(array->size)] = value;
        (array->size)++;
    }
    else
    {
        LABEL* tmpData =
			(LABEL*)malloc( sizeof(LABEL)*2*(array->maxSize));
        array->maxSize *= 2;
        memcpy(tmpData, array->data, sizeof(LABEL)*(array->size) );
        free(array->data);
        array->data = tmpData;
        array->data[array->size] = value;
        (array->size)++;
    }
}

void extensibleLABELArrayDestroy
(
 ExtensibleLABELArray *array
)
{
    free(array->data);
}
