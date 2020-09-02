#include<stdio.h>
#include "surfaceIntegrate_link.h"
//#include <string.h>
//#include<time.h>
#include"surfaceIntegrate_struct.h"
#include"rowSubsection.h"
#include"rowSubsections.hpp"

namespace Foam{
    namespace fvc{
 void surfaceIntegrate_host
 (const int  *owner,
  const int  *neighbour,
  const double *issf,
  double  *ivf,
  int issfSize,
  int ivfSize,
  int vectorSize,
  int owner_size,
  int neighbour_size
  )

 {
     //--------the data of computing owner and neighbour-------

    struct surfaceIntegrate_para para;
 	para.issf_Ptr=issf;
 	para.ivf_Ptr=ivf;
 	para.owner_Ptr=owner;
    para.neighbour_Ptr=neighbour;

    const int nFaces=issfSize;
    const int nCells=ivfSize;
	para.nFaces=nFaces;
	para.nCells=nCells;


    const int CEPs=64;
	int secNumInSeg;
	const rowSubsection **seg;
	int colRoundNum;


	if(vectorSize==1)
	{
	static const int faceSize1 =(int)(ArraySize/(sizeof(double)*2*3+sizeof(double)+2*sizeof(int)));
    static SWFoam::RowSubsections subsections_double(nFaces, nCells, CEPs, 4,owner, neighbour, faceSize1, 1);
	secNumInSeg=subsections_double.getSecNumInSeg();
	seg=subsections_double.getSubsections();
    colRoundNum=subsections_double.getColRoundNum();
	}
	else if(vectorSize==3)
	{
	static const int faceSize3 =(int)(ArraySize/(sizeof(double)*3*3*2+sizeof(double)+ 2*sizeof(int)));
	static SWFoam::RowSubsections subsections_Vector(nFaces, nCells, CEPs, 4,owner, neighbour, faceSize3, 3);
	secNumInSeg=subsections_Vector.getSecNumInSeg();
	seg=subsections_Vector.getSubsections();
    colRoundNum=subsections_Vector.getColRoundNum();
	}
	else if(vectorSize==6)
	{
	static const int faceSize6 =(int)(ArraySize/(sizeof(double)*3*6*2 + sizeof(double)+2*sizeof(int)));
	static SWFoam::RowSubsections subsections_Tensor(nFaces, nCells, CEPs, 4,owner, neighbour, faceSize6, 6);
	secNumInSeg=subsections_Tensor.getSecNumInSeg();
	seg=subsections_Tensor.getSubsections();
    colRoundNum=subsections_Tensor.getColRoundNum();
	}
	else if(vectorSize==9)
	{
	static const int faceSize9 =(int)(ArraySize/(sizeof(double)*3*9*2 +sizeof(double)+2*sizeof(int)));
	static SWFoam::RowSubsections subsections_SymmTensor(nFaces, nCells, CEPs, 4,owner, neighbour, faceSize9, 9);
	secNumInSeg=subsections_SymmTensor.getSecNumInSeg();
	seg=subsections_SymmTensor.getSubsections();
    colRoundNum=subsections_SymmTensor.getColRoundNum();
	}

	para.secs=seg;
	para.secNumInSeg=secNumInSeg;
    para.colRoundNum=colRoundNum;
    para.CEPs=CEPs;
    para.vector_size=vectorSize;

	 //--------call spe algorithm-------
	SWlink_slave(&para);
    }
}
}
