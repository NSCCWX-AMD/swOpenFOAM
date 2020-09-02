#include "swAmul_host.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "amulMacros.h"

#include "swAmul_slave.h"
#include <athread.h>

inline void amul_host(amul_para_ptr in,amul_translate_array_ptr in2)
{
    athread_spawn(ATmul_Diag, in);
    athread_join();

    athread_spawn(ATmul_Upper_Lowwer, in2);
    athread_join();
}

