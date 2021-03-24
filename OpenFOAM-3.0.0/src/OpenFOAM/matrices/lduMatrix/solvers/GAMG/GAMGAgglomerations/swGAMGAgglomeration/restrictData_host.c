#include <athread.h>
#include "swRestInterStruct.h"
#include "mpi.h"
#include <time.h>
#include <stdio.h>

void restrictData_host(restStruct* rs)
{

#if 0  // check input
{
	// swInt*    mapPtr         
	// swFloat*  fPtr
	// swInt**   localStartEnd
	int fSize[2][4] = {
		{25000, 24800, 25240, 24800},
		{12500, 12440, 12620, 12440}};

	int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    char filename[256];
    sprintf(filename,"debug_intput/processor%d.dat",rank);
    FILE* debug_fp = fopen(filename,"a");
    static int ncall = 0;
    if (ncall==0)
    {
    	/*string dateString = clock::date();
        string timeString = clock::clockTime();
        fprintf(debug_fp,"Date: %s\n", dateString.c_str());
        fprintf(debug_fp,"Time: %s\n\n", timeString.c_str());*/
        time_t tt;
		time( &tt );
		tt = tt + 8*3600;
		struct tm* t= gmtime( &tt );
		fprintf(debug_fp,"time: %d-%02d-%02d %02d:%02d:%02d\n",
        	t->tm_year + 1900, t->tm_mon + 1, t->tm_mday,
        	t->tm_hour, t->tm_min, t->tm_sec);
    }
    swInt*    mapPtr         = rs->mapPtr;        
	swFloat*  fPtr           = rs->fPtr;
	swInt**   localStartEnd  = rs->localStartEnd;
    fprintf(debug_fp,"ncall=%d: \n", ncall+1);
    fprintf(debug_fp,"localStartEnd\n");
    for (int i = 0; i < 64; ++i)
    {
    	swInt* lsei = localStartEnd[i];
    	fprintf(debug_fp,"i=%d    ",i);
    	for (int j = 0; j < 4; ++j)
    	{
    		fprintf(debug_fp,"%d ", lsei[j]);
    	}
    	fprintf(debug_fp,"\n");
    }
    fprintf(debug_fp,"fSize=%d, nc  i  f2c   f, mapPtr=%p\n", fSize[ncall%2][rank], mapPtr);
    for(int i = 0 ; i < fSize[ncall%2][rank] ; i++) 
    {
    	if (mapPtr[i]>100000 || mapPtr[i]<0)
    	{
    		fprintf(debug_fp,"error:");
    	}
        fprintf(debug_fp,"nc=%d, i=%d:  %d, %15.14e\n", ncall+1, i, mapPtr[i], fPtr[i]);
    }
    ncall++;

    fclose(debug_fp);
}	
#endif

	athread_spawn(restrictData_slave, rs);
	athread_join();

	//printf("restrictData_host is called\n");
	//exit(0);
}
