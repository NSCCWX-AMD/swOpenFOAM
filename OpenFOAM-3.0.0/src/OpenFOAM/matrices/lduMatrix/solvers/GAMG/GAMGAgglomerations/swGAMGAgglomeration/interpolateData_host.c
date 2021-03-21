#include <athread.h>
#include "swRestInterStruct.h"

void interpolateData_host(interStruct* is)
{
	athread_spawn(interpolateData_slave, is);
	athread_join();

	//printf("interpolateData_host is called\n");
}
