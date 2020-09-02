#include <athread.h>
void SLAVE_FUN(interpolate1_)(struct Parameter* paramt);

void interpolate1_host(struct Parameter* paramt)
{
	athread_spawn(interpolate1_, paramt);
	athread_join();
}

