#ifndef CALL_SWFIELD_H
#define CALL_SWFIELD_H
#include "sunwayMacros.h"
#include "basicTypes.h"
typedef swFloat SCALAR;
typedef swInt LABEL;
#ifdef __cplusplus
extern "C" {
#endif

    typedef struct struct_swVectorDotTensor {
        SCALAR* igGrad;
        SCALAR* f1P; //vector results
        SCALAR* f2P; //tensor
        SCALAR* f3P; //vector
        SCALAR s1;
        LABEL size;
        LABEL vector_dim;
        LABEL tensor_dim;
    } swVectorDotTensor;

    void call_swVectorDotTensor_slave(swVectorDotTensor *data);
    void call_swVectorDotVector_slave(swVectorDotTensor *data);

    void call_swVectorMulScalar_slave(swVectorDotTensor *data);
#ifdef __cplusplus
}
#endif

#endif /* CALL_SWFIELD_H */

