#include <math.h>
#include "slave.h"
#include "vectorOps_struct.h"
#include "register_communication.h"

void residualNormFactor_slave(MVM_ParametersPack* MVM_ParametersPack_hostPtr)
{
	MVM_ParametersPack MVM_ParametersPack_slave;
    volatile LABEL get_reply, put_reply;
    volatile LABEL my_id;
    const    LABEL nThreads = 64;
    char     Array_slave[ArraySize];

    get_reply = 0;
    athread_get(PE_MODE,
                MVM_ParametersPack_hostPtr,
                &MVM_ParametersPack_slave,
                sizeof(MVM_ParametersPack),
                (LABEL*) (&get_reply),
                0,0,0);
    while (get_reply != 1);

    SCALAR* A1_slavePtr = NULL;
    SCALAR* A2_slavePtr = NULL;
    SCALAR* A3_slavePtr = NULL;
    SCALAR* A4_slavePtr = NULL;

    // user function pointer in slave core
    void (*userFunc_slavePtr)(MVM_Arrays*) = NULL;

    LABEL sizePerCycle;
    my_id = athread_get_id(-1);

    const LABEL n1      = MVM_ParametersPack_slave.n1;
    SCALAR *A1_hostPtr  = MVM_ParametersPack_slave.A1Ptr;
    SCALAR *A2_hostPtr  = MVM_ParametersPack_slave.A2Ptr;
    SCALAR *A3_hostPtr  = MVM_ParametersPack_slave.A3Ptr;
    SCALAR *A4_hostPtr  = MVM_ParametersPack_slave.A4Ptr;

    sizePerCycle = slaveDoubleNumber;

    // size: number of variables in each cycle
    LABEL remainder = n1 % nThreads;
    // two kinds of length
    LABEL lenShort  = n1 / nThreads;
    LABEL lenLong   = lenShort + 1;

    // lenLocal: local length of stride
    // sPosition: start position in each slave core
    LABEL lenLocal, sPosition;
    if (my_id < remainder)
    {
        lenLocal  = lenLong;
        sPosition = my_id * lenLocal;
    }
    else
    {
        lenLocal  = lenShort;
        sPosition = my_id * lenLocal + remainder;
    }

    LABEL i, j, remaining, offset_vector = sPosition;
    sizePerCycle = (lenLocal < sizePerCycle)? lenLocal : sizePerCycle;

    SCALAR k1Num_slave = 0.0, k2Denom_slave = 0.0;

    for(i=sizePerCycle; i<=lenLocal; i+=sizePerCycle)
    {
        A1_slavePtr = (SCALAR*) ALIGNED(Array_slave);
        A2_slavePtr = (SCALAR*) ALIGNED(A1_slavePtr + sizePerCycle);
        A3_slavePtr = (SCALAR*) ALIGNED(A2_slavePtr + sizePerCycle);

        get_reply = 0;
        athread_get(PE_MODE,
                    A1_hostPtr + offset_vector,
                    A1_slavePtr,
                    sizePerCycle * sizeof(SCALAR),
                    (LABEL*)&get_reply,
                    0,0,0);
        athread_get(PE_MODE,
                    A2_hostPtr + offset_vector,
                    A2_slavePtr,
                    sizePerCycle * sizeof(SCALAR),
                    (LABEL*)&get_reply,
                    0,0,0);
        athread_get(PE_MODE,
                    A3_hostPtr + offset_vector,
                    A3_slavePtr,
                    sizePerCycle * sizeof(SCALAR),
                    (LABEL*)&get_reply,
                    0,0,0);
        while(get_reply!=3);

        for(j=0; j<sizePerCycle; j++)
        {
            A1_slavePtr[j] = A2_slavePtr[j] - A3_slavePtr[j];
            k1Num_slave   += fabs(A2_slavePtr[j] - A3_slavePtr[j]);
            k2Denom_slave += fabs(A2_slavePtr[j]);
        }

        // return results
        put_reply = 0 ;
        athread_put(PE_MODE,
                    A1_slavePtr,
                    A1_hostPtr + offset_vector,
                    sizePerCycle * sizeof(SCALAR),
                    (LABEL*)&put_reply,
                    0,0);
        while(put_reply!=1);

        // update remainings and offset_vector
        offset_vector += sizePerCycle;
        remaining      = lenLocal - i;
        sizePerCycle   = ((remaining < sizePerCycle) && (remaining > 0))? remaining : sizePerCycle;
    }

    ALLSYN();
    doublev4 send;
    doublev4 recv;
    SCALAR *sendinfo = (SCALAR*)&send;
    SCALAR *recvinfo = (SCALAR*)&recv;
    const LABEL column_id = COL(my_id);
    const LABEL row_id = ROW(my_id);

    sendinfo[0] = k1Num_slave;
    sendinfo[1] = k2Denom_slave;

    recvinfo[0] = k1Num_slave;
    recvinfo[1] = k2Denom_slave;

    LABEL base = 0 ;
    for(i=1; i<8;  base|=i, i<<=1)
    {
        if((column_id & base) == base)
        {
            if(column_id & (i))
            {
                REG_GETR(recv);
                sendinfo[1] += recvinfo[1];
                sendinfo[0] += recvinfo[0];
            }
            else
            {
                REG_PUTR(send, column_id + i);
            }
        }
    }
    base = 0 ;

    for(i=1; i<8; base|=i, i<<=1)
    {
        if((row_id & base) == base)
        {
            if(row_id & i)
            {
                REG_GETC(recv);
                sendinfo[1] += recvinfo[1];
                sendinfo[0] += recvinfo[0];
            }
            else
            {
                REG_PUTC(send,row_id + i);
            }
        }
    }
    if(my_id == 63)
    {
        put_reply = 0 ;
        athread_put(PE_MODE,
                    &(sendinfo[0]),
                    MVM_ParametersPack_slave.k1Ptr,
                    2*sizeof(SCALAR),
                    (LABEL*)&put_reply,
                    0,0);
        while(put_reply!=1);
    }
}
