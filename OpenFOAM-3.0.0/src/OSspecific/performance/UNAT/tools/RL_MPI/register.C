#include "register.H"
#include <iostream>
#include <math.h>
#include <assert.h>
#include "swMacro.h"
using namespace std;
extern "C"{
#include <athread.h>
// #include <algorithm>
void slave_initTable(RlmpiInfo *data);
}
//void initSendList(int *dataSendList, int mshBlockNum)
//{
//	printf("start Init SendList...\n");
//	schedule_data = new RlmpiInfo[mshBlockNum];;
//	for(int blockIdx=0;blockIdx<mshBlockNum;blockIdx++){
//		RlmpiInitializer reg;
//		vector<vector<int> > sendList(BLOCKNUM64K);
//		for(int i=0;i<BLOCKNUM64K;i++){
//			int row = blockIdx*BLOCKNUM64K+i;
//			for(int j=0;j<BLOCKNUM64K;j++){
//				int index = row*BLOCKNUM64K+j;
//				for(int k=0;k<dataSendList[index];k++){
//					sendList[i].push_back(j);
//				}
//			}
//		}
//		reg.init(sendList);
//		reg.copyRlmpiInfo(&schedule_data[blockIdx]);
//	}
//
//	printf("finish Init SendList...\n");
//}
//
//void initSendList(int *dataSendList, vector<vector<int> > dataList,
//			int mshBlockNum)
//{
//	printf("start Init SendList...\n");
//	schedule_data = new RlmpiInfo[mshBlockNum];;
//	for(int blockIdx=0;blockIdx<mshBlockNum;blockIdx++){
//		RlmpiInitializer reg;
//		vector<vector<int> > sendList(BLOCKNUM64K);
//		vector<vector<int> > nDataList(BLOCKNUM64K);
//		for(int i=0;i<BLOCKNUM64K;i++){
//			int row = blockIdx*BLOCKNUM64K+i;
//			int packIdx = 0;
//			for(int j=0;j<BLOCKNUM64K;j++){
//				int index = row*BLOCKNUM64K+j;
//				for(int k=0;k<dataSendList[index];k++){
////					printf("%d,%d,%d,%d,%d\n",row,j,packIdx,dataList[row][packIdx],dataSendList[index]);
//					sendList[i].push_back(j);
//					nDataList[i].push_back(dataList[row][packIdx]);
//					packIdx++;
//				}
//			}
//		}
//		reg.init(sendList, nDataList);
//		reg.copyRlmpiInfo(&schedule_data[blockIdx]);
//	}
//
//	printf("finish Init SendList...\n");
//}

void initTable(int index){
	__real_athread_spawn((void*)slave_initRlmpiInfo,&schedule_data[index]);
	athread_join();
}

void destroyTable(int index){
	__real_athread_spawn((void*)slave_destroyRlmpiInfo,&schedule_data[index]);
	athread_join();
}
