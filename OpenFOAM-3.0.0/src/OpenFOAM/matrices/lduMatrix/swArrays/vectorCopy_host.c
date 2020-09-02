#include <math.h>
#include <athread.h>
#include "userFunc_SAXPY.h"

void vectorCopy_host(MVM_Arrays *MVM_ArraysPtr)
{
    MVM_ParametersPack pack;
    pack.A1Ptr       = MVM_ArraysPtr->A1Ptr ;
    pack.A2Ptr       = MVM_ArraysPtr->A2Ptr ;
    pack.n1          = MVM_ArraysPtr->n1    ;

    athread_spawn(vectorCopy_slave, &pack);
    athread_join();
}
