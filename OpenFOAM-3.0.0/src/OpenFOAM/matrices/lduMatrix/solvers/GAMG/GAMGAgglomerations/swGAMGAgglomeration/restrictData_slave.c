#include "slave.h"
#include "swRestInterStruct.h"
/*#define athread_get(MODE, SRC, DST, SIZE, REPLY, x, y, z) {\
  memcpy(DST, SRC, SIZE);\
  (*(REPLY))++;\
}
#define athread_put(MODE, SRC, DST, SIZE, REPLY, x, y) {\
  memcpy(DST, SRC, SIZE);\
  (*(REPLY))++;\
}
__thread_local int fake_id;
#define athread_get_id(x) fake_id;*/
void restrictData_slave(restStruct* rs) //_original
{
	restStruct rs_slave;
	volatile swInt get_reply, put_reply;
    volatile swInt my_id;
    const    swInt nThreads = 64;
    char     Array_slave[ArraySize];
    volatile swInt range_local[4];

    get_reply = 0;
    athread_get(PE_MODE,
                rs,
                &rs_slave,
                sizeof(restStruct),
                (swInt*) (&get_reply),
                0,0,0);
    while (get_reply != 1);

    swFloat* c_slavePtr   = NULL;
    swFloat* f_slavePtr   = NULL;
    swInt*   map_slavePtr = NULL;

    swInt sizePerCycle;
    my_id = athread_get_id(-1);
    //if(my_id==0)printf("myid=%d\n",my_id);

    swFloat* f_hostPtr     = rs_slave.fPtr;
    swFloat* c_hostPtr     = rs_slave.cPtr;
    swInt*   map_hostPtr   = rs_slave.mapPtr;
    swInt**  range_hostPtr = rs_slave.localStartEnd;
    swInt    slaveCycles   = rs_slave.slaveCycles;

    swInt cycleI;
    for(cycleI=0; cycleI<slaveCycles; ++cycleI)
    {
        get_reply = 0;
        athread_get(PE_MODE,
                    &range_hostPtr[my_id + cycleI*64][0],
                    &range_local,
                    4*sizeof(swInt),
                    (volatile swInt*) (&get_reply),
                    0,0,0);
        while (get_reply != 1);

        //- range_local[0]: coarse left
        //- range_local[1]: coarse rightvolatile
        //- range_local[2]: fine left
        //- range_local[3]: fine right

        swInt fLenLocal, cLenLocal;
        fLenLocal = range_local[3] - range_local[2] + 1;
        cLenLocal = range_local[1] - range_local[0] + 1;

        if(cLenLocal > 5000)
        {
            sizePerCycle = 512;
        }
        else
        {
            sizePerCycle = 1024;
        }

        c_slavePtr = (swFloat*) ALIGNED(Array_slave);

        swInt i, j, remaining, offset_vector = range_local[2];
        sizePerCycle = (fLenLocal < sizePerCycle)? fLenLocal : sizePerCycle;

        for(i=0; i<cLenLocal; ++i)
        {
            c_slavePtr[i] = 0.0;
        }

        for(i=sizePerCycle; i<=fLenLocal; i+=sizePerCycle)
        {
            f_slavePtr = (swFloat*) ALIGNED(c_slavePtr + cLenLocal);
            map_slavePtr = (swInt*) ALIGNED(f_slavePtr + sizePerCycle);

            get_reply = 0;
            athread_get(PE_MODE,
                        f_hostPtr + offset_vector,
                        f_slavePtr,
                        sizePerCycle * sizeof(swFloat),
                        (swInt*)&get_reply,
                        0,0,0);
            athread_get(PE_MODE,
                        map_hostPtr + offset_vector,
                        map_slavePtr,
                        sizePerCycle * sizeof(swInt),
                        (swInt*)&get_reply,
                        0,0,0);
            while(get_reply!=2);

            //- compute
            for(j=0; j<sizePerCycle; ++j)
            {
                volatile swInt cPos = map_slavePtr[j];
                if(cPos >= range_local[0] && cPos <= range_local[1])
                {
                    c_slavePtr[cPos - range_local[0]] += f_slavePtr[j];
                    // c_hostPtr[cPos] += f_slavePtr[j];
                }
            }
            // update remainings and offset_vector for next loop
            offset_vector += sizePerCycle;
            remaining      = fLenLocal - i;
            sizePerCycle   = ((remaining < sizePerCycle) && (remaining > 0))? remaining : sizePerCycle;
        }

        // return results
        put_reply = 0 ;
        athread_put(PE_MODE,
                    c_slavePtr,
                    c_hostPtr+range_local[0],
                    cLenLocal*sizeof(swFloat),
                    (swInt*)&put_reply,
                    0,0);
        while(put_reply!=1);
    }
}
/*void restrictData_slave(restStruct* rs) {
    if (_MYID > 0) return;
    for (fake_id = 0; fake_id < 64; fake_id ++) {
        restrictData_slave_original(rs);
    }
}*/