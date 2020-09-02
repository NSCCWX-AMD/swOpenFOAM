#ifndef SMVM_reg_H
#define SMVM_reg_H

#include "RlmpiInitializer.hxx"
#define SPNUMs 64
#ifdef __cplusplus
extern "C" {
#endif
	void SMVM_reg(int *dataSendList, int mshBlockNum);
	void initTable(int index);
	void destroyTable(int index);
	Schedule *schedule_data;
#ifdef __cplusplus
}
#endif
//extern "C" void SMVM_reg();
#endif
