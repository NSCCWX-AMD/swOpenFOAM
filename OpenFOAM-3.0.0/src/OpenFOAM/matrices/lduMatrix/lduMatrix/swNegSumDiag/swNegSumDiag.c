#include "swNegSumDiag.h"
//#include "negSumDiag_slave.h"
#include "swNegSumDiag_host.h"
#include "athread.h"
#include <assert.h>

void negSumDiag_host(const SCALAR *lower,const SCALAR *upper,SCALAR *Diag,const LABEL *l,const LABEL *u,LABEL nCells,LABEL nFaces){
	LABEL face;
//	printf("nCells:%d,nFaces:%d\n",nCells,nFaces);
//	for(face=0;face<nFaces;face++){
//		//printf("%d,%d,%d\n",face,l[face],u[face]);
//		//assert(u[face]<=face);
//		//assert(l[face]<=face);
//		Diag[l[face]] -= lower[face];
//		Diag[u[face]] -= upper[face];
//	}
	struct LDUMatrix Matrix;
	Matrix.lower = lower;
	Matrix.upper = upper;
	Matrix.diag  = Diag;
	Matrix.rowAddr     = l;
	Matrix.colAddr     = u;
	Matrix.numCell= nCells;
	Matrix.numEdge= nFaces;
	//athread_init();
//	athread_spawn(negSumDiag_slave,&Para);
//	athread_join();
	static int flag=0;
	if(flag==0){
		printf("start BlockOrdering...\n");
		BlockOrderingMLB(Matrix);
		initOwnNeiSendList();
		SMVM_reg(ownNeiSendList,mshBlockNum);
		flag=1;
	}
	printf("start slave compute...\n");
	SMVM_SP_BlockOrder(Matrix);
	printf("finish slave compute...\n");
//	athread_halt();
}
