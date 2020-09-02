#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include "./negSumDiag_host.h"
//extern "C"{
//#include "RegisterLevelMsgPass.hxx"
//}

int* readFile(char *name);
int main(){
	int i,j,k,row,col;

	//构造存储形式为LDU格式的矩阵Matrix
	struct LDUMatrix Matrix;
	char owner[] = "owner_512000";
	char neighbour[] = "neighbor_512000";
	Matrix.rowAddr = readFile(owner);
	Matrix.colAddr = readFile(neighbour);
	Matrix.lower = (SCALAR*)malloc(sizeof(SCALAR)*NONZERONUM);
	Matrix.upper = (SCALAR*)malloc(sizeof(SCALAR)*NONZERONUM);
	j = 0 ;
	for(i=0;i<NONZERONUM;i++){
		Matrix.upper[i] = (double)(Matrix.rowAddr[i]+1)/(Matrix.colAddr[i]+1);
		Matrix.lower[i] = (double)(Matrix.colAddr[i]+1)/(Matrix.rowAddr[i]+1);
	}
	Matrix.diag = (SCALAR*)malloc(sizeof(SCALAR)*CELLNUM);
	for(i=0;i<CELLNUM;i++){
		Matrix.diag[i] = i;
	}
	Matrix.numCell = CELLNUM;
	Matrix.numEdge = NONZERONUM;
	
	//构造自变量x
	SCALAR *x = (SCALAR*)malloc(sizeof(SCALAR)*CELLNUM);
	for(i=0;i<CELLNUM;i++){
		x[i] = (double)(i+1)/(i+2);
	}
	SCALAR *b_MLBSlaveReordered = (SCALAR*)malloc(sizeof(SCALAR)*CELLNUM);
	SCALAR *b_MLBReordered = (SCALAR*)malloc(sizeof(SCALAR)*CELLNUM);
	SCALAR* b = (SCALAR*)malloc(sizeof(SCALAR)*CELLNUM);
	for(i=0;i<CELLNUM;i++){b_MLBSlaveReordered[i]=Matrix.diag[i];}
	for(i=0;i<CELLNUM;i++){b[i]=Matrix.diag[i];}

	BlockOrderingMLB(Matrix);
	struct timeval start,end;

	if(TEST==0){	
		//输出结果（调用主核计算BlockOrdering）
//		gettimeofday(&start,NULL);
//		SMVM_MP_BlockOrder(Matrix,x,b_MLBReordered);
//		gettimeofday(&end,NULL);
//		int timeuse = 1000000*(end.tv_sec-start.tv_sec)+end.tv_usec-start.tv_usec;
//		printf("Master Processor Time(MLB): %f us\n",(double)timeuse);
	}else{
		//输出结果（调用从核计算）
		//gettimeofday(&start,NULL);
		SMVM_SP_BlockOrder(Matrix);
		//gettimeofday(&end,NULL);
		//int timeuse = 1000000*(end.tv_sec-start.tv_sec)+end.tv_usec-start.tv_usec;
		//printf("Slave Time: %f us\n",(double)timeuse);
	}
	//计算结果b（单核版本）
	//memset大坑，千万别用
	//memset(b,0,CELLNUM);
	gettimeofday(&start,NULL);
	for(i=0;i<NONZERONUM;i++){
		b[Matrix.rowAddr[i]] -= Matrix.upper[i];
//		if(Matrix.rowAddr[i]==6400){
//			printf("row: %d,%f,%f\n",i,Matrix.upper[i],Matrix.diag[Matrix.rowAddr[i]]);
//		}
		b[Matrix.colAddr[i]] -= Matrix.lower[i];
//		if(Matrix.colAddr[i]==6400){
//			printf("col: %d,%f,%f\n",i,Matrix.lower[i],Matrix.diag[Matrix.colAddr[i]]);
//		}
//		b[Matrix.rowAddr[i]] += Matrix.upper[i]*x[Matrix.colAddr[i]];
//		b[Matrix.colAddr[i]] += Matrix.lower[i]*x[Matrix.rowAddr[i]];
	}
//	for(i=0;i<CELLNUM;i++){
//		b[i] += Matrix.diag[i]*x[i];
//	}
	gettimeofday(&end,NULL);
	int timeuse = 1000000*(end.tv_sec-start.tv_sec)+end.tv_usec-start.tv_usec;
	printf("Master Processor Time: %f us\n",(double)timeuse);

	//检查结果正确性
	int flag =0;
//	for(i=0;i<CELLNUM;i++){
//		if(fabs((b[i]-b_MLBReordered[i])/b[i])>EPSILON){
//			flag = 1;
//			printf("The MLB_MP result is wrong: i=%d,b[i]=%.9f,b_MLB[i]=%.9f\n",i,b[i],b_MLBReordered[i]);
//			printf("The index of reordered cell[%d] is %d\n",i,postCellOrder[i]);
//			break;
//		}
//	}
	for(i=0;i<CELLNUM;i++){
		if(fabs((Matrix.diag[i]-b[i])/Matrix.diag[i])>EPSILON){
			flag = 1;
			printf("The MLB_SP result is wrong: i=%d,b[i]=%.9f,b_MLB[i]=%.9f\n",i,Matrix.diag[i],b[i]);
			printf("The index of reordered cell[%d] is %d\n",i,postCellOrder[i]);
			break;
		}
	}

#if DEBUG
	//打印输出矩阵
	for(i=0;i<20;i++){
		printf("%5.2f",b_MLBReordered[i]);
	}
	printf("\n");
	SCALAR test[20][20] = {0};
	for(i=0;i<20;i++){
		if(Matrix.rowAddr[i]<20 && Matrix.colAddr[i]<20){
			test[Matrix.rowAddr[i]][Matrix.colAddr[i]] = Matrix.upper[i];
			test[Matrix.colAddr[i]][Matrix.rowAddr[i]] = Matrix.lower[i];
		}
	}
	for(i=0;i<20;i++){
		test[i][i] = Matrix.diag[i];
	}
	printf("A[0~20][0~20]\n");
	for(i=0;i<20;i++){
		for(j=0;j<20;j++){
			printf("%5.2f",test[i][j]);
		}
		printf("\n");
	}
#endif

	//释放内存空间
	free(Matrix.lower);
	free(Matrix.upper);
	free(Matrix.rowAddr);
	free(Matrix.colAddr);
	free(Matrix.diag);
	free(x);
	free(b);
	free(b_MLB);
	free(postCellOrder);
	free(postEdgeOrder);
	free(blockStarts);
	free(cellStarts);
	free(ownReordered);
	free(neiReordered);
	free(upperReordered);
	free(lowerReordered);
	free(diagReordered);
	free(xReordered);
	free(b_MLBReordered);
	return 0;
}

int* readFile(char* name){
	FILE *fp = fopen(name,"r");
	if(fp==NULL){
		return 0;
	}
	char c1[100],c2[100];
	int cellNum;
	while(1)
	{
		fgets(c1,100,fp);
		if(c1[0]=='('){
			cellNum = atoi(c2);
			break;
		}
		fgets(c2,100,fp);
		if(c2[0]=='('){
			cellNum = atoi(c1);
			break;
		}else{
		}
	}
	int *buff = (int*)malloc(sizeof(int)*NONZERONUM);
	int i=0;
	while(i<NONZERONUM)
	{
		fgets(c1,100,fp);
		buff[i] = atoi(c1);
		i++;
	}
	printf("InternalFaceNum = %d\n",i);
	return buff;
}


