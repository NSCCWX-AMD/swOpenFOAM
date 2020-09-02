#ifndef CALL_SWFIELD_H
#define CALL_SWFIELD_H
#include "sunwayMacros.h"
typedef swFloat scalar;
typedef swInt label;
#ifdef __cplusplus
extern "C" {
#endif

    typedef struct struct_swVectorDotTensor {
        scalar* igGrad;
        scalar* f1P; //vector results
        scalar* f2P; //tensor
        scalar* f3P; //vector
        scalar s1;
        label size;
        label vector_dim;
        label tensor_dim;
    } swVectorDotTensor;

    void call_swVectorDotTensor_slave(swVectorDotTensor *data);
    void call_swVectorDotVector_slave(swVectorDotTensor *data);

    void call_swVectorMulScalar_slave(swVectorDotTensor *data);
#ifdef __cplusplus
}
#endif

#endif /* CALL_SWFIELD_H */

