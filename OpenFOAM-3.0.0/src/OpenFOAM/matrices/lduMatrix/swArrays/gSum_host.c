#include <math.h>
#include <athread.h>
#include "userFunc_SAXPY.h"

void gSum_host(MVM_Arrays *MVM_ArraysPtr, void (*userFuncPtr)(MVM_Arrays*))
{
    MVM_ParametersPack pack;
    pack.A1Ptr       = MVM_ArraysPtr->A1Ptr;
    pack.A2Ptr       = MVM_ArraysPtr->A2Ptr;
    pack.A3Ptr       = MVM_ArraysPtr->A3Ptr;
    pack.A4Ptr       = MVM_ArraysPtr->A4Ptr;
    pack.k1Ptr       = MVM_ArraysPtr->k1Ptr;
    pack.k2Ptr       = MVM_ArraysPtr->k2Ptr;
    pack.k1          = MVM_ArraysPtr->k1;
    pack.k2          = MVM_ArraysPtr->k2;
    pack.n1          = MVM_ArraysPtr->n1;
    pack.returnA2    = MVM_ArraysPtr->returnA2;
    pack.userFuncPtr = userFuncPtr;

    athread_spawn(gSum_slave, &pack);
    athread_join();
}
