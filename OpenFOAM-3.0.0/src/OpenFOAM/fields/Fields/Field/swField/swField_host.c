#include "athread.h"
#include "swField_host.h"
extern void SLAVE_FUN(swVectorDotTensor_slave)(swVectorDotTensor*);
void call_swVectorDotTensor_slave(swVectorDotTensor *data) {
    athread_spawn(swVectorDotTensor_slave, data);
    athread_join();
}

extern void SLAVE_FUN(swVectorDotVector_slave)(swVectorDotTensor*);
void call_swVectorDotVector_slave(swVectorDotTensor *data) {
    athread_spawn(swVectorDotVector_slave, data);
    athread_join();
}

extern void SLAVE_FUN(swVectorMulScalar_slave)(swVectorDotTensor*);
void call_swVectorMulScalar_slave(swVectorDotTensor *data)
{
    athread_spawn(swVectorMulScalar_slave, data);
    athread_join();
}
