#include <athread.h>
#include "swRestInterStruct.h"

void restrictData_host(restStruct* rs)
{
	athread_spawn(restrictData_slave, rs);
	athread_join();

	//printf("restrictData_host is called\n");
	//exit(0);
}
