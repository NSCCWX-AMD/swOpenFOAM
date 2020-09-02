#include "swNegSumDiag_reg.h"
//#include <stdlib.h>
//#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <assert.h>
//#include "./RegMsgPas/RegisterLevelMsgPass.hxx"
using namespace std;
extern "C"{
#include "athread.h"
#include <algorithm>
void slave_initTable(Schedule *data);
}
void SMVM_reg(int *dataSendList, int mshBlockNum)
{
	printf("start Init SendList...\n");
	schedule_data = new Schedule[mshBlockNum];;
	for(int blockIdx=0;blockIdx<mshBlockNum;blockIdx++){
		RlmpiInitializer reg;
		vector<vector<ThreadID> > sendList(SPNUMs);
		for(int i=0;i<SPNUMs;i++){
			int row = blockIdx*SPNUMs+i;
			for(int j=i+1;j<SPNUMs;j++){
				int index = row*SPNUMs+j;
				for(int k=0;k<dataSendList[index];k++){
					sendList[i].push_back(j);
				}
			}
		}
		reg.init(sendList);
		reg.copyinfo(&schedule_data[blockIdx]);
		printf("%d\n",schedule_data[blockIdx].nCycle);
		for(int i=0;i<SPNUMs;i++){
//			printf("%d,%d,%d\n",blockIdx,i,schedule_data[blockIdx].putr_schedules[i][0]);
		}
	}

	printf("finish Init SendList...\n");
}

void initTable(int index){
	schedule_data[index].destroy = 0;
	__real_athread_spawn((void*)slave_initTable,&schedule_data[index]);
	athread_join();
}

void destroyTable(int index){
	schedule_data[index].destroy = 1;
	__real_athread_spawn((void*)slave_initTable,&schedule_data[index]);
	athread_join();
}
