#include "athread.h"


void CG_init()
{
    athread_init();
}

void CG_halt()
{
    athread_halt();
}

// the array to hold exception signals
unsigned int CG_signals[64];
const int CG_signals_size = 64*sizeof(int);

