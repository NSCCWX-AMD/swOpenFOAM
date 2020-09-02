#include <math.h>
#include "slave.h"
#include "vectorOps_struct.h"

void vectorOps_slave(MVM_ParametersPack* MVM_ParametersPack_hostPtr)
{
    MVM_ParametersPack MVM_ParametersPack_slave;
    MVM_Arrays MVM_Arrays_slave;
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
    const LABEL returnA2= MVM_ParametersPack_slave.returnA2;
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

    LABEL i, remaining, offset_vector = sPosition;
    sizePerCycle = (lenLocal < sizePerCycle)? lenLocal : sizePerCycle;

    MVM_Arrays_slave.k1 = MVM_ParametersPack_slave.k1;
    MVM_Arrays_slave.k2 = MVM_ParametersPack_slave.k2;
    userFunc_slavePtr   = MVM_ParametersPack_slave.userFuncPtr;

    for(i=sizePerCycle; i<=lenLocal; i+=sizePerCycle)
    {
        get_reply = 0;
        A1_slavePtr = (SCALAR*) ALIGNED(Array_slave);
        MVM_Arrays_slave.A1Ptr = A1_slavePtr;
        athread_get(PE_MODE,
                    A1_hostPtr + offset_vector,
                    A1_slavePtr,
                    sizePerCycle * sizeof(SCALAR),
                    (LABEL*)&get_reply,
                    0,0,0);
        while(get_reply!=1);

        if(A2_hostPtr)
        {
            get_reply = 0;
            A2_slavePtr = (SCALAR*) ALIGNED(A1_slavePtr + sizePerCycle);
            MVM_Arrays_slave.A2Ptr = A2_slavePtr;
            athread_get(PE_MODE,
                        A2_hostPtr + offset_vector,
                        A2_slavePtr,
                        sizePerCycle * sizeof(SCALAR),
                        (LABEL*)&get_reply,
                        0,0,0);
            while(get_reply!=1);
        }

        if(A3_hostPtr)
        {
            get_reply = 0;
            A3_slavePtr = (SCALAR*) ALIGNED(A2_slavePtr + sizePerCycle);
            MVM_Arrays_slave.A3Ptr = A3_slavePtr;
            athread_get(PE_MODE,
                        A3_hostPtr + offset_vector,
                        A3_slavePtr,
                        sizePerCycle * sizeof(SCALAR),
                        (LABEL*)&get_reply,
                        0,0,0);
            while(get_reply!=1);

            if(A4_hostPtr)
            {
                get_reply = 0;
                A4_slavePtr = (SCALAR*) ALIGNED(A3_slavePtr + sizePerCycle);
                MVM_Arrays_slave.A4Ptr = A4_slavePtr;
                athread_get(PE_MODE,
                            A4_hostPtr + offset_vector,
                            A4_slavePtr,
                            sizePerCycle * sizeof(SCALAR),
                            (LABEL*)&get_reply,
                            0,0,0);
                while(get_reply!=1);
            }
        }

        MVM_Arrays_slave.nCells = sizePerCycle;
        // calculate
        (*userFunc_slavePtr)(&MVM_Arrays_slave);

        // return results
        put_reply = 0 ;
        athread_put(PE_MODE,
                    A1_slavePtr,
                    A1_hostPtr + offset_vector,
                    sizePerCycle * sizeof(SCALAR),
                    (LABEL*)&put_reply,
                    0,0);
        while(put_reply!=1);

        if(returnA2)
        {
            put_reply = 0 ;
            athread_put(PE_MODE,
                        A2_slavePtr,
                        A2_hostPtr + offset_vector,
                        sizePerCycle * sizeof(SCALAR),
                        (LABEL*)&put_reply,
                        0,0);
            while(put_reply!=1);
        }

        // update remainings and offset_vector
        offset_vector += sizePerCycle;
        remaining      = lenLocal - i;
        sizePerCycle   = ((remaining < sizePerCycle) && (remaining > 0))? remaining : sizePerCycle;
    }
}
