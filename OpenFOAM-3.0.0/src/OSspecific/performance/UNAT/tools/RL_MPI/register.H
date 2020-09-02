#ifndef REGISTER_H
#define REGISTER_H

#include "RlmpiInitializer.hxx"
// #include <vector>
// #include <iostream>
#include <map>
#include <cstring>
//#ifdef __cplusplus
//extern "C"{
//#endif
// void initSendList(int *dataSendList, int mshBlockNum);
// void initSendList(int *dataSendList, vector<vector<int> > dataList,
			// int mshBlockNum);
void initTable(int index);
void destroyTable(int index);
RlmpiInfo *schedule_data;
//#ifdef __cplusplus
//}
//#endif
#endif
