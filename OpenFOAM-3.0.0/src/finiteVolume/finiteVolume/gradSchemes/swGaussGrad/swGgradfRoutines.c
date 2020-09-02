#include "athread.h"
#include "swGgradfRoutines.h"
#include "swGgradfRoutines_slave.h"

extern int SLAVE_FUN(swGgradfInner_scalar_slave)(struct swGgradf_scalar_paras *);
extern int SLAVE_FUN(swGgradfInner_vector_slave)(struct swGgradf_vector_paras *);
extern int SLAVE_FUN(swGgradfDivide_slave)(swGgradfDivide_paras *);

// global exception array
extern int  CG_signals[];
extern const int CG_signals_size;

// post exception check macro that
// takes a exception processing
// function
#define CG_CHECK(excptFun) \
{ \
    int i; \
    for(i=0; i<64; i++) \
    { \
        if(CG_signals[i]!=0) \
            excptFun(i, CG_signals[i]); \
    } \
}

// exception processing function that
// print massages corresponding to the
// slave function sig number
void swGgradfInner_excpt(label cpeID, label sigNum)
{
    printf("cpe %d exceptions :\n", cpeID);
    switch(sigNum)
    {
        case 1:
            printf("***SIGNAL_1: data size in ldm exceeds the bound\n");
            break;
        case 2:
            printf("***SIGNAL_2: data size in ldm exceeds the bound\n");
            break;
        case 3:
            printf("***SIGNAL_3: column is not local\n");
            break;
        case 4:
            printf("***SIGNAL_4: data size in ldm exceeds the bound\n");
            break;
        case 5:
            printf("***SIGNAL_5: data size in ldm exceeds the bound\n");
            break;
        case 6:
            printf("***SIGNAL_6: column is not local\n");
            break;
    }
    exit(-1);
}

void swGgradfInnerRoutine_scalar(
        const char* Sf_vptr,
        const scalar* issf_ptr,
        char* igGrad_vptr,
        const label* lPtr,
        const label* uPtr,
        const label fSize,
        const label cSize,
        const label vector_size,
        const label vector_offset,
        const struct rowSubsection** secs,
        const label secNumInseg,
        const label colRoundNum
        )
{
    struct swGgradf_scalar_paras paras =
    {
        Sf_vptr,
        issf_ptr,
        igGrad_vptr,
        lPtr,
        uPtr,
        fSize,
        cSize,
        vector_size,
        vector_offset,
        secs,
        secNumInseg,
        colRoundNum
    };

    memset(CG_signals, 0, CG_signals_size);
    athread_spawn(swGgradfInner_scalar_slave, &paras);
    athread_join();
    // post exception check
    CG_CHECK(swGgradfInner_excpt);
}


void swGgradfInnerRoutine_vector(
        const char* Sf_vptr,
        const char* issf_vptr,
        char* igGrad_tptr,
        const label* lPtr,
        const label* uPtr,
        const label fSize,
        const label cSize,
        const label vector_size,
        const label vector_offset,
        const label tensor_size,
        const label tensor_offset,
        const struct rowSubsection** secs,
        const label secNumInseg,
        const label colRoundNum
        )
{

    struct swGgradf_vector_paras paras =
    {
        Sf_vptr,
        issf_vptr,
        igGrad_tptr,
        lPtr,
        uPtr,
        fSize,
        cSize,
        vector_size,
        vector_offset,
        tensor_size,
        tensor_offset,
        secs,
        secNumInseg,
        colRoundNum
    };

    memset(CG_signals, 0, CG_signals_size);
    athread_spawn(swGgradfInner_vector_slave, &paras);
    athread_join();
    // post exception check
    CG_CHECK(swGgradfInner_excpt);
}

void swGgradfDivide_host
(
 scalar* igGrad,
 const scalar* volumn,
 const label cSize,
 const label size
)
{
    swGgradfDivide_paras para =
    {
        igGrad,
        volumn,
        cSize,
        size
    };

    athread_spawn(swGgradfDivide_slave, &para);
    athread_join();
}







