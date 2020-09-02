#include "merge5.h"
#include "swMacro.h"
#include <assert.h>

define_e2v_FunPtr(swMerge5_mf_1)
{
	swFloat *mf = accessArray(frontEdgeData,0);
	swFloat *cell_x = accessArray(data->vertexData,0);
	swFloat *XacU   = accessArray(data->frontEdgeData,0);
	swFloat *face_x = accessArray(data->frontEdgeData,1);

	swInt edgeNumber = getArraySize(frontEdgeData);
	swInt dims = getArrayDims(frontEdgeData,0);
	int iedge,iDim;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			faceVal = mf[iedge*dims+iDim] >= 0.0 
				? cell_x[endVertices[iedge]*dims+iDim]
				: cell_x[startVertices[iedge]*dims+iDim];
			XacU[iedge*dims+iDim] = face_x[iedge*dims+iDim] - faceVal;
		}
	}
}

define_e2v_FunPtr(swMerge5_mf_2)
{
	swFloat *mf = accessArray(frontEdgeData,0);
	swFloat *gradPhi  = accessArray(data->vertexData,0);
	swFloat *XacU     = accessArray(data->frontEdgeData,0);
	swFloat *phiDelta = accessArray(frontEdgeData,1);

	swInt edgeNumber = getArraySize(frontEdgeData);
	swInt dims = getArrayDims(frontEdgeData,0);
	int iedge,iDim;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			phiDelta[iedge*dims+iDim] += mf[iedge*dims+iDim] >= 0.0 
				? gradPhi[endVertices[iedge]*dims+iDim]
				: XacU[startVertices[iedge]*dims+iDim];
		}
	}
}

define_e2v_FunPtr(swMerge5_mf_3)
{
	swFloat *mf = accessArray(frontEdgeData,0);
	swFloat *phi  = accessArray(data->vertexData,0);
	swFloat *su   = accessArray(data->vertexData,1);
	swFloat *phiDelta = accessArray(frontEdgeData,1);
	swFloat *paras  = accessArray(paraData,0);

	swInt edgeNumber = getArraySize(frontEdgeData);
	swInt dims = getArrayDims(frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat fci,fce,blend,phiUDS;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			phiUDS = mf[idx] >= 0.0
				? phi[endVertices[iedge]*dims+iDim]
				: phi[startVertices[iedge]*dims+iDim];
			fci = MIN(mf[idx], 0.0)*phi[startVertices[iedge]*dims+iDim]
				+ MAX(mf[idx], 0.0)*phi[endVertices[iedge]*dims+iDim];
			fce = mf[idx] * (phiUDS+phiDelta[idx]);
			blend = paras[3]*(fce - fci);
//if(startVertices[iedge]==405||endVertices[iedge]==405) printf("%d,%f,%f,%f\n",iedge, paras[3],fce,fci);
			su[endVertices[iedge]*dims+iDim] -= blend;
			su[startVertices[iedge]*dims+iDim] += blend;
		}
	}
}

define_e2v_FunPtr(swMerge5_lam_1)
{
	swFloat *face_lam = accessArray(frontEdgeData,0);
	swFloat *vis      = accessArray(vertexData,0);
	swFloat *face_D   = accessArray(data->frontEdgeData,0);
	swFloat *visFace  = accessArray(data->frontEdgeData,1);
	swFloat *paras    = accessArray(paraData, 0);

	swInt edgeNumber = getArraySize(frontEdgeData);
	swInt dims = getArrayDims(frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			faceVal = face_lam[idx]*vis[startVertices[iedge]*dims+iDim]
				+ (1-face_lam[idx])*vis[endVertices[iedge]*dims+iDim];
			faceVal = paras[0]/paras[1]+(faceVal-paras[0])/paras[2];
			visFace[idx] = faceVal * face_D[idx];
		}
	}
}

define_e2v_FunPtr(swMerge5_lam_2)
{
	swFloat *face_lam = accessArray(frontEdgeData,0);
	swFloat *vis      = accessArray(vertexData,0);
	swFloat *visFace  = accessArray(data->frontEdgeData,0);
	swFloat *paras    = accessArray(paraData, 0);

	swInt edgeNumber = getArraySize(frontEdgeData);
	swInt dims = getArrayDims(frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			faceVal = face_lam[idx]*vis[startVertices[iedge]*dims+iDim]
				+ (1-face_lam[idx])*vis[endVertices[iedge]*dims+iDim];
			visFace[idx] = paras[0]/paras[1]+(faceVal-paras[0])/paras[2];
		}
	}
}

define_e2v_FunPtr(swMerge5_lam_3)
{
	swFloat *face_lam = accessArray(frontEdgeData,0);
	swFloat *cell_x   = accessArray(data->vertexData,0);
	swFloat *gradPhi  = accessArray(data->vertexData,1);
	swFloat *dPhi     = accessArray(frontEdgeData,1);

	swInt edgeNumber = getArraySize(frontEdgeData);
	swInt dims = getArrayDims(frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat dPhidXac,xpn;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			xpn = cell_x[startVertices[iedge]*dims+iDim]
				- cell_x[endVertices[iedge]*dims+iDim];
			dPhidXac = face_lam[idx]*gradPhi[startVertices[iedge]*dims+iDim]
				+ (1-face_lam[idx])*gradPhi[endVertices[iedge]*dims+iDim];
			dPhi[idx] = xpn*dPhidXac;
		}
	}
}

define_e2v_FunPtr(swMerge5_lam_4)
{
	swFloat *face_lam = accessArray(frontEdgeData,0);
	swFloat *gradPhi  = accessArray(data->vertexData,0);
	swFloat *dPhidFn  = accessArray(frontEdgeData,1);
	swFloat *face_n   = accessArray(data->frontEdgeData,0);

	swInt edgeNumber = getArraySize(frontEdgeData);
	swInt dims = getArrayDims(frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat dPhidXac;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			dPhidXac = face_lam[idx]*gradPhi[startVertices[iedge]*dims+iDim]
				+ (1-face_lam[idx])*gradPhi[endVertices[iedge]*dims+iDim];
			dPhidFn[idx] = dPhidXac*face_n[idx];
		}
	}
}

define_e2v_FunPtr(swMerge5_su_1)
{
	swFloat *vis  = accessArray(frontEdgeData,0);
	swFloat *dPhi = accessArray(data->frontEdgeData,0);
	swFloat *su   = accessArray(vertexData,0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			faceVal = vis[idx] * dPhi[idx];
			su[endVertices[iedge]*dims+iDim] -= faceVal;
			su[startVertices[iedge]*dims+iDim] += faceVal;
		}
	}
}

define_e2v_FunPtr(swMerge5_su_2)
{
	swFloat *vis  = accessArray(data->frontEdgeData,0);
	swFloat *dPhi = accessArray(data->frontEdgeData,1);
	swFloat *su   = accessArray(vertexData,0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			faceVal = vis[idx] *dPhi[idx];
			su[endVertices[iedge]*dims+iDim] += faceVal;
			su[startVertices[iedge]*dims+iDim] -= faceVal;
		}
	}
}

define_e2v_FunPtr(swMerge5_su_3)
{
	swFloat *visFace = accessArray(frontEdgeData,0);
	swFloat *mf      = accessArray(frontEdgeData,1);
	swFloat *rface   = accessArray(data->frontEdgeData,0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			rface[idx] = -visFace[idx]-MAX(mf[idx], 0.0);
		}
	}
}

define_e2v_FunPtr(swMerge5_su_4)
{
	swFloat *visFace = accessArray(data->frontEdgeData,0);
	swFloat *dPhi    = accessArray(data->frontEdgeData,1);
	swFloat *fd      = accessArray(data->frontEdgeData,2);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims = getArrayDims(data->frontEdgeData,0);
	int iedge,iDim,idx;
	swFloat faceVal;
	for(iedge=0;iedge<edgeNumber;iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			idx = iedge*dims+iDim;
			fd[idx] = visFace[idx]*dPhi[idx];
		}
	}
}
