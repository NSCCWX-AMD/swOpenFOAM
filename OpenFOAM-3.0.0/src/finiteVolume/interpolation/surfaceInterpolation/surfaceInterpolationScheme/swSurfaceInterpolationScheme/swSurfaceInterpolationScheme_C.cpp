#include <stdio.h>
#include "swSurfaceInterpolationScheme_C.H"
#include "interpolateParameter.h"
//#include "rowSubsection.h"
#include "rowSubsections.hpp"
#include "interpolate1_host.h"


//extern "C" {
//	#include <athread.h>
//	void SLAVE_FUN(interpolate1_)(struct Parameter* paramt);
//}

//static bool ifAthreadInit=false;

void interpolate1(
	swInt faceNum,
	swInt cellNum,
	const void* fvPtr,
	const void* cvPtr,
	const swFloat* weiPtr,
	const swFloat* _1wPtr,
	const swInt* ownPtr,
	const swInt* neiPtr,
	swInt typeSize)
{

    //the use rate of LDM
    swFloat usage[64][20];

    //subsection information
    static SWFoam::RowSubsections _1double( faceNum, cellNum, 64,4, ownPtr, neiPtr, 700, sizeof(swFloat));
    static SWFoam::RowSubsections _3double( faceNum, cellNum, 64,4, ownPtr, neiPtr, 300, 3*sizeof(swFloat));
    static SWFoam::RowSubsections _6double( faceNum, cellNum, 64,4, ownPtr, neiPtr, 200, 6*sizeof(swFloat));
    static SWFoam::RowSubsections _9double( faceNum, cellNum, 64,4, ownPtr, neiPtr, 90,  9*sizeof(swFloat));

    SWFoam::RowSubsections* subsecPtr;

    int eleNum;

    if(      typeSize ==  8 ){
    	subsecPtr = & _1double;
    	eleNum = 1;
    }
    else if( typeSize == 24 ){
    	subsecPtr = & _3double;
    	eleNum = 3;
    }
    else if( typeSize == 48 ){
    	subsecPtr = & _6double;
    	eleNum = 6;
    }
    else if( typeSize == 72 ){
    	subsecPtr = & _9double;
    	eleNum = 9;
    }
    else {
        printf("\n***Error: fail to select type in inpolate() surfaceInterpolationScheme!\n");
        //exit(1);
    }

	struct Parameter paramt;
	paramt.fvPtr = (char*)fvPtr;
	paramt.cvPtr = (char*)cvPtr;
	paramt.weiPtr = (swFloat*)weiPtr;
	paramt._1wPtr = (swFloat*)_1wPtr;
	paramt.ownPtr = (swInt*)ownPtr;
	paramt.neiPtr = (swInt*)neiPtr;
	paramt.subSecPtr = (rowSubsection**)(subsecPtr->getSubsections());
	paramt.secNum = subsecPtr->getSecNumInSeg();
	paramt.eleNum = eleNum;
	paramt.typeSize = typeSize;
	paramt.usage = &usage[0][0];

    interpolate1_host(&paramt);
	//void* tmp = &paramt;
	//__real_athread_spawn((void*)slave_interpolate1_, tmp);
	//athread_join();

}
