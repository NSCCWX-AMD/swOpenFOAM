#include <math.h>
#include "slave.h"
#include "vectorOps_struct.h"

void vectorCopy_slave(MVM_ParametersPack* MVM_ParametersPack_hostPtr)
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

    LABEL sizePerCycle;
    my_id = athread_get_id(-1);

    const LABEL n1      = MVM_ParametersPack_slave.n1;
    SCALAR *A1_hostPtr  = MVM_ParametersPack_slave.A1Ptr;
    SCALAR *A2_hostPtr  = MVM_ParametersPack_slave.A2Ptr;

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

    LABEL i, remaining, offset_vector = sPosition;
    sizePerCycle = (lenLocal < sizePerCycle)? lenLocal : sizePerCycle;

    for(i=sizePerCycle; i<=lenLocal; i+=sizePerCycle)
    {
        A1_slavePtr = (SCALAR*) ALIGNED(Array_slave);
        get_reply = 0;
        athread_get(PE_MODE,
                    A2_hostPtr + offset_vector,
                    A1_slavePtr,
                    sizePerCycle * sizeof(SCALAR),
                    (LABEL*)&get_reply,
                    0,0,0);
        while(get_reply!=1);

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
}
