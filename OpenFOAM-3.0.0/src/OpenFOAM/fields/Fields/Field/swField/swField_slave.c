#include "stdlib.h"
#include <slave.h>
#include <dma.h>
#include "swField_host.h"
#define ALIGNED(addr) ((((unsigned long)(addr-1)>>5)+1)<<5)

#define A_DMA_GET_SET(da,mode,len,re_addr)\
({  \
dma_set_op(&da,DMA_GET);    \
dma_set_mode(&da,mode); \
dma_set_size(&da,len);  \
dma_set_reply(&da,re_addr); \
})

#define A_DMA_GET_RUN(da,src,dest)  \
({  \
    dma(da,src,dest);   \
})

#define A_DMA_PUT_SET(da,mode,len,re_addr)\
({      \
dma_set_op(&da,DMA_PUT);        \
dma_set_mode(&da,mode); \
dma_set_size(&da,len);  \
dma_set_reply(&da,re_addr);     \
})

#define A_DMA_PUT_RUN(da,src,dest)      \
({      \
        dma(da,dest,src);       \
})

#define LDM_MALLOC(ptr,size) \
{ \
    if( get_allocatable_size() > size) \
        ptr = ldm_malloc(size); \
    else \
    { \
        printf("\n***Error: no enough LDM space to malloc\n"); \
        exit(-1); \
    } \
}

// can not be used in inbalance conditional branches
#define SERIAL_BEGIN \
{ \
    int serial_cpeid = athread_get_id(-1); \
    int serial_cpeIndex; \
    for(serial_cpeIndex = 0; \
        serial_cpeIndex < 64; \
        serial_cpeIndex++){\
        if(serial_cpeIndex == \
           serial_cpeid){

#define SERIAL_END \
        } \
        athread_syn(ARRAY_SCOPE, 0xFFFF); \
    } \
}

// can not be used in inbalance conditional branches
#define SPRINT(...) \
{ \
    int serial_cpeid = athread_get_id(-1); \
    int serial_cpeIndex; \
    for(serial_cpeIndex = 0; \
        serial_cpeIndex < 64; \
        serial_cpeIndex++){\
        if(serial_cpeIndex == \
           serial_cpeid){ \
            printf(__VA_ARGS__); \
        } \
        athread_syn(ARRAY_SCOPE, 0xFFFF); \
    } \
}

// can not be used in inbalance conditional branches
#define CPRINT(cpeID,...) \
{\
    if(athread_get_id(-1) == cpeID)\
        printf(__VA_ARGS__); \
}
//    athread_syn(ARRAY_SCOPE, 0xFFFF);\

// declare the global exception signal array
extern int CG_signals[];

#include "swField_slave.h"

void swVectorDotTensor_slave(swVectorDotTensor * data) {

    SCALAR *f1P = data->f1P; //res;
    SCALAR *f2P = data->f2P; //vector
    SCALAR *f3P = data->f3P; //tensor;
    const LABEL size = data->size; // how many vectors over all
    const LABEL vector_dim = data->vector_dim;
    const LABEL tensor_dim = data->tensor_dim;

    volatile LABEL nVector = (size + 63) / 64;
    volatile LABEL segStart = nVector * _MYID;

    // take care of the array tail
    volatile LABEL segEnd;
    if ((nVector * _MYID + nVector) < size) {
        segEnd = nVector * _MYID + nVector;
    } else {
        segEnd = size;
        if (segStart > segEnd) {
            nVector = 0;
        } else {
            nVector = segEnd - segStart;
        }
    }

    const LABEL nPerCycle0 = 512;

    SCALAR f1P_ldm[nPerCycle0 * 3];
    SCALAR f2P_ldm[nPerCycle0 * 3];
    SCALAR f3P_ldm[nPerCycle0 * 9];

    int nCycle0 = nVector / nPerCycle0;
    if (nCycle0 <= 0) {
        nCycle0 = 0;
    }
    int i, j;
    if (nVector > nPerCycle0)
        for (i = 0; i < nCycle0; i++) {
            const int vecStartPos = segStart * vector_dim + i * nPerCycle0*vector_dim;
            const int tensorStartPos = segStart * tensor_dim + i * nPerCycle0*tensor_dim;

            volatile unsigned long get_reply = 0;
            athread_get(PE_MODE, data->f2P + vecStartPos,
                    f2P_ldm, nPerCycle0 * vector_dim * sizeof (double), &get_reply, 0, 0, 0);

            athread_get(PE_MODE, data->f3P + tensorStartPos,
                    f3P_ldm, nPerCycle0 * tensor_dim * sizeof (double), &get_reply, 0, 0, 0);
            dma_wait(&get_reply, 2);

            for (j = 0; j < nPerCycle0; j++) {

                f1P_ldm[j * 3] = f2P_ldm[j * 3] * f3P_ldm[j * 9] +
                        f2P_ldm[j * 3 + 1] * f3P_ldm[j * 9 + 3] +
                        f2P_ldm[j * 3 + 2] * f3P_ldm[j * 9 + 6];

                f1P_ldm[j * 3 + 1] = f2P_ldm[j * 3] * f3P_ldm[j * 9 + 1] +
                        f2P_ldm[j * 3 + 1] * f3P_ldm[j * 9 + 4] +
                        f2P_ldm[j * 3 + 2] * f3P_ldm[j * 9 + 7];

                f1P_ldm[j * 3 + 2] = f2P_ldm[j * 3] * f3P_ldm[j * 9 + 2] +
                        f2P_ldm[j * 3 + 1] * f3P_ldm[j * 9 + 5] +
                        f2P_ldm[j * 3 + 2] * f3P_ldm[j * 9 + 8];
            }

            volatile unsigned long put_reply = 0;
            athread_put(PE_MODE, f1P_ldm, data->f1P + vecStartPos, nPerCycle0 * vector_dim * sizeof (double), &put_reply, 0, 0);
            dma_wait(&put_reply, 1);
        }

    const int n1Left = nVector - nCycle0*nPerCycle0;
    const int vecStartPos = segStart * vector_dim + nCycle0 * nPerCycle0 * vector_dim;
    const int tensorStartPos = segStart * tensor_dim + nCycle0 * nPerCycle0 * tensor_dim;


    if (n1Left > 0) {
        volatile unsigned long get_reply = 0;
        athread_get(PE_MODE, data->f2P + vecStartPos,
                f2P_ldm, n1Left * vector_dim * sizeof (double), &get_reply, 0, 0, 0);

        athread_get(PE_MODE, data->f3P + tensorStartPos,
                f3P_ldm, n1Left * tensor_dim * sizeof (double), &get_reply, 0, 0, 0);
        dma_wait(&get_reply, 2);
    }

    if (n1Left > 0) {
        for (j = 0; j < n1Left; j++) {

            f1P_ldm[j * 3] = f2P_ldm[j * 3] * f3P_ldm[j * 9] +
                    f2P_ldm[j * 3 + 1] * f3P_ldm[j * 9 + 3] +
                    f2P_ldm[j * 3 + 2] * f3P_ldm[j * 9 + 6];

            f1P_ldm[j * 3 + 1] = f2P_ldm[j * 3] * f3P_ldm[j * 9 + 1] +
                    f2P_ldm[j * 3 + 1] * f3P_ldm[j * 9 + 4] +
                    f2P_ldm[j * 3 + 2] * f3P_ldm[j * 9 + 7];

            f1P_ldm[j * 3 + 2] = f2P_ldm[j * 3] * f3P_ldm[j * 9 + 2] +
                    f2P_ldm[j * 3 + 1] * f3P_ldm[j * 9 + 5] +
                    f2P_ldm[j * 3 + 2] * f3P_ldm[j * 9 + 8];
        }

        if (n1Left > 0) {
            volatile unsigned long put_reply = 0;
            athread_put(PE_MODE, f1P_ldm, data->f1P + vecStartPos, n1Left * vector_dim * sizeof (double), &put_reply, 0, 0);
            dma_wait(&put_reply, 1);
        }
    }
}

void swVectorDotVector_slave(swVectorDotTensor * data) {

    SCALAR *f1P = data->f1P; //res;
    SCALAR *f2P = data->f2P; //vector
    SCALAR *f3P = data->f3P; //tensor;
    const LABEL size = data->size; // how many vectors over all
    const LABEL vector_dim = data->vector_dim;

    volatile LABEL nVector = (size + 63) / 64;
    volatile LABEL segStart = nVector * _MYID;

    // take care of the array tail
    volatile LABEL segEnd;
    if ((nVector * _MYID + nVector) < size) {
        segEnd = nVector * _MYID + nVector;
    } else {
        segEnd = size;
        if (segStart > segEnd) {
            nVector = 0;
        } else {
            nVector = segEnd - segStart;
        }
    }

    const LABEL nPerCycle0 = 512;

    SCALAR f1P_ldm[nPerCycle0 * 1];
    SCALAR f2P_ldm[nPerCycle0 * 3];
    SCALAR f3P_ldm[nPerCycle0 * 3];

    int nCycle0 = nVector / nPerCycle0;
    if (nCycle0 <= 0) {
        nCycle0 = 0;
    }
    int i, j;
    if (nVector > nPerCycle0)
        for (i = 0; i < nCycle0; i++) {
            const int vecStartPos = segStart * vector_dim + i * nPerCycle0*vector_dim;
            const int scaStartPos = segStart + i * nPerCycle0;

            volatile unsigned long get_reply = 0;
            athread_get(PE_MODE, data->f2P + vecStartPos,
                    f2P_ldm, nPerCycle0 * vector_dim * sizeof (double), &get_reply, 0, 0, 0);

            athread_get(PE_MODE, data->f3P + vecStartPos,
                    f3P_ldm, nPerCycle0 * vector_dim * sizeof (double), &get_reply, 0, 0, 0);
            dma_wait(&get_reply, 2);

            for (j = 0; j < nPerCycle0; j++) {

                f1P_ldm[j] = f2P_ldm[j * 3] * f3P_ldm[j * 3] +
                        f2P_ldm[j * 3 + 1] * f3P_ldm[j * 3 + 1] +
                        f2P_ldm[j * 3 + 2] * f3P_ldm[j * 3 + 2];
            }

            volatile unsigned long put_reply = 0;
            athread_put(PE_MODE, f1P_ldm, data->f1P + scaStartPos, nPerCycle0 * sizeof (double), &put_reply, 0, 0);
            dma_wait(&put_reply, 1);
        }

    const int n1Left = nVector - nCycle0*nPerCycle0;
    const int vecStartPos = segStart * vector_dim + nCycle0 * nPerCycle0 * vector_dim;
    const int scaStartPos = segStart + nCycle0 * nPerCycle0;


    if (n1Left > 0) {
        volatile unsigned long get_reply = 0;
        athread_get(PE_MODE, data->f2P + vecStartPos,
                f2P_ldm, n1Left * vector_dim * sizeof (double), &get_reply, 0, 0, 0);

        athread_get(PE_MODE, data->f3P + vecStartPos,
                f3P_ldm, n1Left * vector_dim * sizeof (double), &get_reply, 0, 0, 0);
        dma_wait(&get_reply, 2);
    }

    if (n1Left > 0) {
        for (j = 0; j < n1Left; j++) {

            f1P_ldm[j] = f2P_ldm[j * 3] * f3P_ldm[j * 3] +
                    f2P_ldm[j * 3 + 1] * f3P_ldm[j * 3 + 1] +
                    f2P_ldm[j * 3 + 2] * f3P_ldm[j * 3 + 2];
        }

        if (n1Left > 0) {
            volatile unsigned long put_reply = 0;
            athread_put(PE_MODE, f1P_ldm, data->f1P + scaStartPos, n1Left * sizeof (double), &put_reply, 0, 0);
            dma_wait(&put_reply, 1);
        }
    }
}

void swVectorMulScalar_slave(swVectorDotTensor * data) {

    SCALAR *f1P = data->f1P; //res;
    SCALAR *f2P = data->f2P; //vector
    const LABEL size = data->size; // how many vectors over all
    const SCALAR s1 = data->s1;

    volatile LABEL nVector = (size + 63) / 64;
    volatile LABEL segStart = nVector * _MYID;

    // take care of the array tail
    volatile LABEL segEnd;
    if ((nVector * _MYID + nVector) < size) {
        segEnd = nVector * _MYID + nVector;
    } else {
        segEnd = size;
        if (segStart > segEnd) {
            nVector = 0;
        } else {
            nVector = segEnd - segStart;
        }
    }

    const LABEL nPerCycle0 = 7168;


    SCALAR f2P_ldm[nPerCycle0];


    int nCycle0 = nVector / nPerCycle0;
    if (nCycle0 <= 0) {
        nCycle0 = 0;
    }
    int i, j;
    if (nVector > nPerCycle0)
        for (i = 0; i < nCycle0; i++) {

            const int scaStartPos = segStart + i * nPerCycle0;

            volatile unsigned long get_reply = 0;
            athread_get(PE_MODE, data->f2P + scaStartPos, f2P_ldm, nPerCycle0 * sizeof (double), &get_reply, 0, 0, 0);

            dma_wait(&get_reply, 1);

            for (j = 0; j < nPerCycle0; j++) {

                //                f1P_ldm[j] = f2P_ldm[j] * s1;
                f2P_ldm[j] *= s1;
            }

            volatile unsigned long put_reply = 0;
            athread_put(PE_MODE, f2P_ldm, data->f1P + scaStartPos, nPerCycle0 * sizeof (double), &put_reply, 0, 0);
            dma_wait(&put_reply, 1);
        }

    const int n1Left = nVector - nCycle0*nPerCycle0;
    const int scaStartPos = segStart + nCycle0 * nPerCycle0;

    if (n1Left > 0) {
        volatile unsigned long get_reply = 0;
        athread_get(PE_MODE, data->f2P + scaStartPos, f2P_ldm, n1Left * sizeof (double), &get_reply, 0, 0, 0);
        dma_wait(&get_reply, 1);
    }

    if (n1Left > 0) {
        for (j = 0; j < n1Left; j++) {
            //            f1P_ldm[j] = f2P_ldm[j] * s1;
            f2P_ldm[j] *= s1;
        }

        if (n1Left > 0) {
            volatile unsigned long put_reply = 0;
            athread_put(PE_MODE, f2P_ldm, data->f1P + scaStartPos, n1Left * sizeof (double), &put_reply, 0, 0);
            dma_wait(&put_reply, 1);
        }
    }
}
