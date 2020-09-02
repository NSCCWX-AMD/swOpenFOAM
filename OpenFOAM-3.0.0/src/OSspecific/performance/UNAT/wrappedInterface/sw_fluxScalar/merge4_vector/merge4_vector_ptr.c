#include "merge4_vector.h"
#include "swMacro.h"
#include <assert.h>

define_array_FunPtr(swMerge4_vector_1)
{
	swFloat *phiDelta  = accessArray(data->frontEdgeData,0);
	swFloat *face_x0   = accessArray(data->frontEdgeData,1);
	swFloat *XU0       = accessArray(data->frontEdgeData,2);
	swFloat *dPhidXU0  = accessArray(data->frontEdgeData,3);

	swFloat *paras  = accessArray(paraData, 0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			phiDelta[idx] = (face_x0[idx]-XU0[idx])*dPhidXU0[idx];
		}
	}
}

define_array_FunPtr(swMerge4_vector_2)
{
	swFloat *fde       = accessArray(data->frontEdgeData,0);
	swFloat *visac     = accessArray(data->frontEdgeData,1);
	swFloat *face_n0   = accessArray(data->frontEdgeData,2);
	swFloat *dPhidXac0 = accessArray(data->frontEdgeData,3);

	swFloat *paras  = accessArray(paraData, 0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			visac[idx] = paras[3]/paras[6]+(visac[idx]-paras[3])/paras[4];
			fde[idx] = visac[idx]*face_n0[idx]*dPhidXac0[idx];
		}
	}
}

define_array_FunPtr(swMerge4_vector_3)
{
	swFloat *fde       = accessArray(data->frontEdgeData,0);
	swFloat *visac     = accessArray(data->frontEdgeData,1);
	swFloat *face_n0   = accessArray(data->frontEdgeData,2);
	swFloat *dPhidXac0 = accessArray(data->frontEdgeData,3);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			fde[idx] = visac[idx]*face_n0[idx]*dPhidXac0[idx];
		}
	}
}

define_array_FunPtr(swMerge4_vector_4)
{
	swFloat *fce       = accessArray(data->frontEdgeData,0);
	swFloat *mf        = accessArray(data->frontEdgeData,1);
	swFloat *phiUDS    = accessArray(data->frontEdgeData,2);
	swFloat *phiDelta  = accessArray(data->frontEdgeData,3);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			fce[idx] = mf[idx]*(phiUDS[idx]+phiDelta[idx]);
		}
	}
}

define_array_FunPtr(swMerge4_vector_5)
{
	swFloat *fce       = accessArray(data->frontEdgeData,0);
	swFloat *mf        = accessArray(data->frontEdgeData,1);
	swFloat *phiUDS    = accessArray(data->frontEdgeData,2);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			fce[idx] = mf[idx]*phiUDS[idx];
		}
	}
}
