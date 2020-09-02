#include "slave.h"
#include "swRestInterStruct.h"

void interpolateData_slave(interStruct* is)
{
	interStruct is_slave;
	volatile swInt get_reply, put_reply;
    volatile swInt my_id;
    const    swInt nThreads = 64;
    char     Array_slave[ArraySize];
    volatile swInt range_local[4];

    get_reply = 0;
    athread_get(PE_MODE,
                is,
                &is_slave,
                sizeof(interStruct),
                (swInt*) (&get_reply),
                0,0,0);
    while (get_reply != 1);

    swFloat* c_slavePtr   = NULL;
    swFloat* f_slavePtr   = NULL;
    swInt*   map_slavePtr = NULL;
    volatile swInt* offsetMap_slavePtr = NULL;

    swInt sizePerCycle;
    my_id = athread_get_id(-1);

    swInt*   map_hostPtr       = is_slave.mapPtr;
    swInt*   offsetMap_hostPtr = is_slave.offsetMapPtr;
    swFloat* f_hostPtr         = is_slave.fPtr;
    swFloat* c_hostPtr         = is_slave.cPtr;
    swInt**  range_hostPtr     = is_slave.localStartEnd;
    swInt    slaveCycles       = is_slave.slaveCycles;

    swInt cycleI;
    for(cycleI=0; cycleI<slaveCycles; ++cycleI)
    {
        get_reply = 0;
        athread_get(PE_MODE,
                    &range_hostPtr[_MYID + cycleI*64][0],
                    &range_local,
                    4*sizeof(swInt),
                    (swInt*) (&get_reply),
                    0,0,0);
        while (get_reply != 1);

        swInt fLenLocal, cLenLocal;
        fLenLocal = range_local[1] - range_local[0] + 1;
        cLenLocal = range_local[3] - range_local[2] + 1;

        if(fLenLocal > 5000)
        {
            sizePerCycle = 256;
        }
        else
        {
            sizePerCycle = 512;
        }

        f_slavePtr = (swFloat*) ALIGNED(Array_slave);

        swInt i, j, remaining, offset_vector = range_local[2];
        sizePerCycle = (cLenLocal < sizePerCycle)? cLenLocal : sizePerCycle;

        for(i=sizePerCycle; i<=cLenLocal; i+=sizePerCycle)
        {
            c_slavePtr         = (swFloat*) ALIGNED(f_slavePtr + fLenLocal);
            offsetMap_slavePtr = (swInt*)   ALIGNED(c_slavePtr + sizePerCycle);
            map_slavePtr       = (swInt*)   ALIGNED(offsetMap_slavePtr + sizePerCycle + 1);

            get_reply = 0;
            athread_get(PE_MODE,
                        c_hostPtr + offset_vector,
                        c_slavePtr,
                        sizePerCycle * sizeof(swFloat),
                        (swInt*)&get_reply,
                        0,0,0);
            athread_get(PE_MODE,
                        offsetMap_hostPtr + offset_vector,
                        offsetMap_slavePtr,
                        (sizePerCycle + 1) * sizeof(swInt),
                        (swInt*)&get_reply,
                        0,0,0);
            while(get_reply!=2);

            volatile swInt interMapSize = offsetMap_slavePtr[sizePerCycle] - offsetMap_slavePtr[0];
            get_reply = 0;
            athread_get(PE_MODE,
                        map_hostPtr + offsetMap_slavePtr[0],
                        map_slavePtr,
                        interMapSize * sizeof(swInt),
                        (swInt*)&get_reply,
                        0,0,0);
            while(get_reply!=1);

            //- compute
            for(j=0; j<sizePerCycle; ++j)
            {
                volatile swInt locSize = offsetMap_slavePtr[j+1] - offsetMap_slavePtr[j];
                swInt k = 0;

                for(k=0; k<locSize; ++k)
                {
                    volatile swInt fPos = map_slavePtr[offsetMap_slavePtr[j] + k - offsetMap_slavePtr[0]];
                    if(fPos >= range_local[0] && fPos <= range_local[1])
                    {
                        f_slavePtr[fPos - range_local[0]] = c_slavePtr[j];
                    }
                }
            }

            // update remainings and offset_vector for next loop
            offset_vector += sizePerCycle;
            remaining      = cLenLocal - i;
            sizePerCycle   = ((remaining < sizePerCycle) && (remaining > 0))? remaining : sizePerCycle;
        }

        // return results
        put_reply = 0 ;
        athread_put(PE_MODE,
                    f_slavePtr,
                    f_hostPtr+range_local[0],
                    fLenLocal*sizeof(swFloat),
                    (swInt*)&put_reply,
                    0,0);
        while(put_reply!=1);
    }
}
