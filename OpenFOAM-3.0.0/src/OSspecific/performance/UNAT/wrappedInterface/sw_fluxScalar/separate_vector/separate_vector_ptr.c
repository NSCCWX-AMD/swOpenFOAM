#include "separate_vector.h"
#include "swMacro.h"
#include <assert.h>

define_array_FunPtr(swSeparate_vector)
{
	swFloat *fu        = accessArray(data->frontEdgeData,0);
	swFloat *face_x0   = accessArray(data->frontEdgeData,1);
	swFloat *XU0       = accessArray(data->frontEdgeData,2);
	swFloat *dPhidXU0  = accessArray(data->frontEdgeData,3);
	swFloat *face_x1   = accessArray(data->frontEdgeData,4);
	swFloat *XU1       = accessArray(data->frontEdgeData,5);
	swFloat *dPhidXU1  = accessArray(data->frontEdgeData,6);
	swFloat *face_x2   = accessArray(data->frontEdgeData,7);
	swFloat *XU2       = accessArray(data->frontEdgeData,8);
	swFloat *dPhidXU2  = accessArray(data->frontEdgeData,9);
	swFloat *mf        = accessArray(data->frontEdgeData,10);
	swFloat *phiUDS    = accessArray(data->frontEdgeData,11);
	swFloat *visac     = accessArray(data->frontEdgeData,12);
	swFloat *face_D    = accessArray(data->frontEdgeData,13);
	swFloat *Xpn0      = accessArray(data->frontEdgeData,14);
	swFloat *dPhidXac0 = accessArray(data->frontEdgeData,15);
	swFloat *Xpn1      = accessArray(data->frontEdgeData,16);
	swFloat *dPhidXac1 = accessArray(data->frontEdgeData,17);
	swFloat *Xpn2      = accessArray(data->frontEdgeData,18);
	swFloat *dPhidXac2 = accessArray(data->frontEdgeData,19);
	swFloat *face_n0   = accessArray(data->frontEdgeData,20);
	swFloat *face_n1   = accessArray(data->frontEdgeData,21);
	swFloat *face_n2   = accessArray(data->frontEdgeData,22);
	swFloat *rface_1   = accessArray(data->frontEdgeData,23);
	swFloat *rface_2   = accessArray(data->frontEdgeData,24);
	swFloat *rcpac     = accessArray(data->frontEdgeData,25);
	swFloat *denac     = accessArray(data->frontEdgeData,26);
	swFloat *Xac0      = accessArray(data->frontEdgeData,27);
	swFloat *Xac1      = accessArray(data->frontEdgeData,28);
	swFloat *Xac2      = accessArray(data->frontEdgeData,29);
	swFloat *phiCDS    = accessArray(data->frontEdgeData,30);
	swFloat *phiCDS2   = accessArray(data->frontEdgeData,31);
	swFloat *paras     = accessArray(paraData, 0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat phiDelta, fce, fci, blend, fdi, fde;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			if(paras[7]==1)
			{
				fce = mf[idx] *phiUDS[idx];
			} else if(paras[7]==2)
			{
				phiDelta = (face_x0[idx]-Xac0[idx])*dPhidXac0[idx]
					+ (face_x1[idx]-Xac1[idx])*dPhidXac1[idx]
					+ (face_x2[idx]-Xac2[idx])*dPhidXac2[idx];
				fce = mf[idx]*(phiCDS[idx]+phiDelta*paras[8]);
			} else if(paras[7]==3)
			{
				fce = mf[idx]*phiCDS2[idx];
			} else if(paras[7]==4)
			{
				phiDelta = (face_x0[idx]-XU0[idx])*dPhidXU0[idx]
					+ (face_x1[idx]-XU1[idx])*dPhidXU1[idx]
					+ (face_x2[idx]-XU2[idx])*dPhidXU2[idx];
				fce = mf[idx]*(phiUDS[idx]+phiDelta);
			}
			fci = mf[idx]*phiUDS[idx];
			blend = paras[0]*(fce-fci);
			if(paras[6]==1)
			{
				visac[idx] = paras[4]*rcpac[idx]+(visac[idx]-paras[1])/paras[2];
			} else if(paras[6]==2)
			{
				visac[idx] = paras[5]*denac[idx]+(visac[idx]-paras[1]/paras[2]);
			} else
			{
				visac[idx] = paras[1]/paras[3]+(visac[idx]-paras[1])/paras[2];
			}
			fdi = visac[idx]*face_D[idx]*(dPhidXac0[idx]*Xpn0[idx]
						+ dPhidXac1[idx]*Xpn1[idx]
						+ dPhidXac2[idx]*Xpn2[idx]);
			fde = visac[idx]*(dPhidXac0[idx]*face_n0[idx]
						+ dPhidXac1[idx]*face_n1[idx]
						+ dPhidXac2[idx]*face_n2[idx]);
			fu[idx] = -blend+fde-fdi;
			rface_1[idx] = -(visac[idx]*face_D[idx])
				- MAX(mf[idx], 0.0);
			rface_2[idx] = -(visac[idx]*face_D[idx])
				+ MIN(mf[idx], 0.0);
		}
	}
}


