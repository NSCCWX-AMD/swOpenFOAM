#ifndef SWMACRO_H
#define SWMACRO_H

#include "stdio.h"
#include "stdlib.h"

#define BLOCKNUM64K 64
#define SAFELDM 24000
#define EPS 1e-5
#define MAXSTEPS 6

// define the pattern of iterator
#define E2V 0
#define V2E 1
#define ARRAY 2
//#define DEBUG

#ifdef DEBUG
#define LOG(format,...) printf("File: "__FILE__",Line: %05d: "format"\n", __LINE__, ##__VA_ARGS__)
#else
#define LOG(format,...)
#endif


typedef int swInt;
typedef int swInt32;
typedef long swInt64;
typedef double swFloat;
typedef float swFloat32;
typedef double swFloat64;

// standard error dump to file
#define dumpErrorToFile(file, ...) \
{ \
	fprintf(file, "\n***Error in function \"%s\":\n",__FUNCTION__);\
	fprintf(file, __VA_ARGS__); \
	fprintf(file, "\n***at %s +%05d\n\n", __FILE__, __LINE__); \
	fflush(file); \
}

// standard error dump to stderr
#define dumpError(...) \
{ \
	dumpErrorToFile(stderr, __VA_ARGS__) \
}

// print array
#define printArray( unitFmt, array, length) \
{ \
	printf("%s(%d):{", #array, length); \
	size_t i; \
	for(i=0; i<(length); i++) \
	{ \
		printf(" ");\
		printf(unitFmt, (array)[i]); \
	} \
	printf("}\n"); \
}

#define MAX(x, y) ( ((x) > (y)) ? (x) : (y) )
#define MIN(x, y) ( ((x) < (y)) ? (x) : (y) )

// bisection table search
#define biSearch(posi, table, value, length) \
{ \
	size_t ubound, lbound, mbound;\
	lbound = 0; \
	ubound = (length)-1; \
	mbound = lbound; \
	while(lbound < ubound) \
	{ \
		if((table)[mbound] < (value)) \
		  lbound = mbound+1; \
		else if((table)[mbound] > (value)) \
		  ubound = mbound-1; \
		else \
			break; \
		mbound = (lbound+ubound)>>1; \
	} \
	if((value) == (table)[mbound]) \
		(posi) = mbound; \
	else \
	{ \
		(posi) = -1; \
	} \
}

// get precise time on unix
#define getTime(time) \
{ \
    struct timeval timer; \
    gettimeofday(&timer, 0); \
    time = ((double)(timer.tv_sec) + (double)(timer.tv_usec)*1.0e-6); \
}

// new, delete and resize an arrays
#define NEW( type, num ) \
( (type*) malloc((num)*sizeof(type)) )

#define DELETE( array ) \
{ free(array); }

#define RESIZE( type , array, numOld, numNew ) \
{ \
	if( numOld == 0 ) \
	{ \
		if(array != NULL) \
		{ \
			printf("***Warning: declare an un-null array with zero size!"); \
			DELETE( array );\
		} \
		(array) = NEW(type, numNew); \
	} \
	else if( numOld != 0 ) \
	{ \
		type* arrayTmp = NEW(type, numOld); \
		size_t i = numOld; \
		while(i--){ arrayTmp[i] = array[i]; } \
		DELETE( array ); \
		(array) = NEW(type, numNew); \
		i = numOld; \
		while(i--){ array[i] = arrayTmp[i]; } \
		DELETE( arrayTmp ); \
	} \
}

#endif

#define checkResult(array1, array2, count) \
{ \
	std::cout<<"check result..."<<std::endl; \
	for(int i=0;i<(count);i++) \
	{ \
		if(std::fabs(array1[i]-array2[i])>EPS) \
		{ \
			if(array1[i]==0) \
			{ \
				if(std::fabs(array2[i])>EPS) \
				{ \
					printf("Error on index[%d], %.8f, %.8f\n", \
								i, array1[i], array2[i]); \
					std::exit(-1); \
				} \
			} \
			else if(std::fabs((array1[i]-array2[i])/array1[i])>EPS) \
			{ \
				printf("Error on index[%d], %.8f, %.8f\n", \
							i, array1[i], array2[i]); \
				std::exit(-1); \
			} \
		} \
	} \
	std::cout<<"The result is correct!"<<std::endl; \
}


