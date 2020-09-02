/***********slave***********/
#include <stdio.h>
#include <stdlib.h>
#include "slave.h"
#include "rowSubsection.h"
#include "interpolateParameter.h"
#include "slaveUtils.h"

void interpolate1_(struct Parameter* paramt)
{
	//接收传进来的参数
	char* fvPtr = paramt -> fvPtr;
	char* cvPtr = paramt -> cvPtr;
	swFloat* weiPtr = paramt -> weiPtr;
	swFloat* _1wPtr = paramt -> _1wPtr;
	swInt* ownPtr = paramt -> ownPtr;
	swInt* neiPtr = paramt -> neiPtr;
	const struct rowSubsection** subSecPtr = paramt -> subSecPtr;
	swInt secNum = paramt -> secNum;
	swInt typeSize = paramt -> typeSize;
	swInt eleNum = paramt -> eleNum;
	swFloat* usage = paramt -> usage;

	swInt myId, myNF, myStart, lDblVal, lIntVal;
	volatile int reply;
	struct rowSubsection subSec;
	swInt i,j,k,l,m,n,start,count;
	swInt NF, NC, NFN, EM, totalSize;

	myId = athread_get_id(-1);

	//按段数循环，每段内把该拷贝的数据拷贝进来后，按face循环计算
	for( n=0; n<secNum; n++)
	{
		//先把rowSubsection考进来，也可以不用athread_get
		reply = 0;
		athread_get(PE_MODE, &subSecPtr[myId][n], &subSec, sizeof(struct rowSubsection), &reply, 0, 0, 0);
		while (reply != 1);

		const label colSAndC[subSec.nSecs*2];
		dma_desc gv =0;
		volatile unsigned long gReply=0;
		A_DMA_GET_SET(gv,PE_MODE,subSec.nSecs*sizeof(label)*2,&gReply);
		A_DMA_GET_RUN(gv,subSec.colStartsAndCounts,&colSAndC[0]);
		DMA_Wait(&gReply,1);

		//把owner和neighbor数组拷贝进来
		myNF = subSec.nFaces;
		myStart = subSec.faceStart;
		lDblVal = myNF * sizeof(swFloat);
		lIntVal = myNF * sizeof(swInt);
		NF = myNF;
		swInt owner[NF], neighbor[NF];
		reply = 0;
		athread_get(PE_MODE,  ownPtr + myStart,    owner, lIntVal, &reply, 0, 0, 0);
		athread_get(PE_MODE,  neiPtr + myStart, neighbor, lIntVal, &reply, 0, 0, 0);
		while (reply != 2);

		//这时就能算出要拷贝的double行数据的数量，算一下是否超过了LDM的大小，超过则报错，调整maxFaces和subsectionSize，重新产生RowSubsections
		//计算通过owner索引的cell value的数量
		NC = owner[myNF-1] -owner[0] +1;
		//计算通过neighbor索引的cell value 的数量
		NFN = 0;
		for( m=0; m < subSec.nSecs; m++) NFN +=colSAndC[2*m+1];
		//模板形参中元素的个数
		EM = eleNum;
		totalSize = (NF + NC + NFN) * EM * sizeof(swFloat) //face value, cell value by owner and by neighbor
				  + 2 * NF * sizeof(swFloat)               //owner and neighbor
				  + 2 * NF * sizeof(swInt);                 //weight and 1-weght
	    //*(usage + myId*20 + n)=totalSize/(64*1024*1.0);
		//计算，把owner和neighbor站的空间也加起来
		if(totalSize > 62*1024)
		{
			printf("\n***Error: the size of arrays exceed that of LDM! Please decrease the \"max\" in rowSubSection.\n\n");
			printf("myId=%d, totalSize=%d\n",myId,totalSize);
		}
		else{
			//debug
			//usage[myId][n]=totalSize/(64*1024*1.0);
			/*if(myId==30)
			{
				char str[256];
				sprintf(str,"\nmyId=%d, secNum=%d, the usage of LDM is %.2lf.\n",myId,secNum,totalSize/(64*1024*1.0));//*100
				printf(str);
			}*/
			//debug
		}
		//声明数组
		//face value 的值，通过owner索引的cell value的值，通过neighbor索引的cell value 的值
		swFloat faceVal[NF*EM], cellValOw[NC*EM], cellValNe[NFN*EM];
		swFloat weight[NF], _1weight[NF];

		//开始拷贝weights， 和 1-weights
		reply = 0;
		athread_get(PE_MODE,  weiPtr + myStart,   weight, lDblVal, &reply, 0, 0, 0);
		athread_get(PE_MODE,  _1wPtr + myStart, _1weight, lDblVal, &reply, 0, 0, 0);
		while (reply != 2);

		//拷贝通过owner索引的cell value的值
		reply = 0;
		athread_get(PE_MODE, cvPtr + owner[0]*typeSize, cellValOw, (owner[myNF-1]-owner[0]+1)* typeSize, &reply, 0, 0, 0);//&cvPtr[owner[0]]
		while (reply != 1);

		//拷贝通过neighbor索引的cell value 的值
		char* cellValNePtr = (char*)(&cellValNe[0]);
        for( m=0; m < subSec.nSecs; m++)
		{
			start = colSAndC[2*m];
			count = colSAndC[2*m+1];

			reply = 0;
			athread_get(PE_MODE, cvPtr + start*typeSize, cellValNePtr, count* typeSize, &reply, 0, 0, 0);
			while (reply != 1);
			cellValNePtr += count* typeSize;
		}

		//计算
		swFloat  faceVal_, cellValOw_, cellValNe_, weight_, _1weight_;
		swFloat *cvOwPtr, *cvNePtr, *faVaPtr;
		swInt colNum, localN, _2j;
		for (i = 0; i < myNF; i++)
		{
			//找通过neighbor索引的cell value 的值
			colNum = neighbor[i];//int
			localN = 0;//int
			for(j = 0; j < subSec.nSecs; j++)
			{
				_2j = 2 * j;//int
				start = colSAndC[_2j];
				count = colSAndC[_2j+1];
				if( colNum>start+count-1 || colNum<start )localN += count;
				else break;
			}
			localN += colNum-start;

			//按元素个数循环计算
			weight_   = weight[i];
			_1weight_ = 1 - weight_;
			cvOwPtr = &cellValOw[0] +(owner[i]-owner[0])*eleNum;
			cvNePtr = &cellValNe[0] +localN*eleNum ;
			faVaPtr = &faceVal[0]   +i*eleNum ;
			for(l = 0; l < eleNum; l++)
			{
				cellValOw_ = *cvOwPtr;
				cellValNe_ = *cvNePtr;
				faceVal_   = weight_ * cellValOw_ + _1weight_ * cellValNe_;
				*faVaPtr   = faceVal_;

				cvOwPtr++;
				cvNePtr++;
				faVaPtr++;
			}
		}

		//发送结果
		reply = 0;
		athread_put(PE_MODE, &faceVal[0], fvPtr + myStart*typeSize, myNF*typeSize, &reply, 0, 0);//
		while (reply != 1);
	}
}
