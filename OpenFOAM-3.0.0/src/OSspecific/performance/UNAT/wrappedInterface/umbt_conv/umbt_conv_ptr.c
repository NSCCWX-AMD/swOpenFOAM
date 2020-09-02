#include "umbt_conv.h"
#include "swMacro.h"
#include "iterator.h"
#include "math.h"

#define rhs_conv_roe() \
{ \
	/* 273 flops */ \
	swFloat uL[5],uR[5],duc[5],u[5],eig[3],e[3],ad[5],tv[5]; \
	swFloat aL,aR,hL,hR,unL,unR,d,h,a,ke,un,vn,area,uan,eiga,eigb,c1,c2; \
	swFloat tmp1,tmp2; \
	/* 31 flops */ \
	tmp1 = 0.5*(uLf[1]*uLf[1]+uLf[2]*uLf[2]+uLf[3]*uLf[3]); \
	tmp2 = 0.5*(uRf[1]*uRf[1]+uRf[2]*uRf[2]+uRf[3]*uRf[3]); \
	duc[0] = uRf[0]-uLf[0]; \
	duc[1] = uRf[0]*uRf[1]-uLf[0]*uLf[1]; \
	duc[2] = uRf[0]*uRf[2]-uLf[0]*uLf[2]; \
	duc[3] = uRf[0]*uRf[3]-uLf[0]*uLf[3]; \
	duc[4] = (uRf[4]-uLf[4])/0.4+uRf[0]*tmp2-uLf[0]*tmp1; \
	/* 34 flops*/ \
	aL = sqrt(1.4*uLf[4]/uLf[0]); \
	hL = aL*aL/0.4+tmp1; \
	aR = sqrt(1.4*uRf[4]/uRf[0]); \
	hR = aR*aR/0.4+tmp2; \
	/* 58 flops*/ \
	u[0] = sqrt(uLf[0]*uRf[0]); \
	d    = sqrt(uRf[0]/uLf[0]); \
	tmp1 = 1.0/(1.0+d); \
	u[1] = (uLf[1]+d*uRf[1])*tmp1; \
	u[2] = (uLf[2]+d*uRf[2])*tmp1; \
	u[3] = (uLf[3]+d*uRf[3])*tmp1; \
	h    = (hL+d*hR)*tmp1; \
	ke   = 0.5*(u[1]*u[1]+u[2]*u[2]+u[3]*u[3]); \
	u[4] = (h-ke)*0.4*u[0]/1.4; \
	a    = sqrt(abs(1.4*u[4]/u[0])); \
	/* 17 flops*/\
	e[0] = fac_1d[0]; \
	e[1] = fac_1d[1]; \
	e[2] = fac_1d[2]; \
	area = fac_1d[3]; \
	vn   = vgf[0]*e[0]+vgf[1]*e[1]+vgf[2]*e[2]; \
	unL = e[0]*uLf[1]+e[1]*uLf[2]+e[2]*uLf[3]-vn; \
	unR = e[0]*uRf[1]+e[1]*uRf[2]+e[2]*uRf[3]-vn; \
	/* 27 flops */\
	if(1) \
	{ \
		tmp1  = uLf[0]*unL; \
		tmp2  = uRf[0]*unR; \
		tv[0] = tmp1+tmp2; \
		tv[1] = tmp1*uLf[1]+tmp2*uRf[2]+e[0]*(uLf[4]+uRf[4]); \
		tv[2] = tmp1*uLf[2]+tmp2*uRf[3]+e[1]*(uLf[4]+uRf[4]); \
		tv[3] = tmp1*uLf[3]+tmp2*uRf[4]+e[2]*(uLf[4]+uRf[4]); \
		tv[4] = tmp1*hL    +tmp2*hR    +vn  *(uLf[4]+uRf[4]); \
	} else \
	{ \
	} \
	/* 35 flops */\
	uan    = u[1]*e[0]+u[2]*e[1]+u[3]*e[2]; \
	un     = uan-vn; \
	eig[0] = abs(un); \
	eig[1] = abs(un+a); \
	eig[2] = abs(un-a); \
	tmp1   = 2.0*abs(unR+aR-unL-aL); \
	if(eig[1]<tmp1) eig[1] = 0.5*(eig[1]*eig[1]/tmp1+tmp1); \
	tmp1   = 2.0*abs(unR-aR-unL+aL); \
	if(eig[2]<tmp1) eig[2] = 0.5*(eig[2]*eig[2]/tmp1+tmp1); \
	/* 16 flops */\
	eiga = 0.5*(eig[1]+eig[2])-eig[0]; \
	eigb = 0.5*(eig[1]-eig[2]); \
	tmp1 = -uan*duc[0]+e[0]*duc[1]+e[1]*duc[2]+e[2]*duc[3]; \
	tmp2 = ke*duc[0]-u[1]*duc[1]-u[2]*duc[2]-u[3]*duc[3]+duc[4]; \
	tmp2 = tmp2*0.4; \
	/* 55 flops */\
	a  = 1.0/a; \
	c1 = (eiga*tmp2*a+eigb*tmp1)*a; \
	c2 = eigb*tmp2*a+eiga*tmp1; \
	ad[0] = eig[0]*duc[0]+c1; \
	ad[1] = eig[0]*duc[1]+c1*u[1]+c2*e[0]; \
	ad[2] = eig[0]*duc[2]+c1*u[2]+c2*e[1]; \
	ad[3] = eig[0]*duc[3]+c1*u[3]+c2*e[2]; \
	ad[4] = eig[0]*duc[4]+c1*h   +c2*uan; \
	rhsl[0] = MAX(0.1, 0.5*area*(tv[0]-ad[0])); \
	rhsl[1] = MAX(0.1, 0.5*area*(tv[1]-ad[1])); \
	rhsl[2] = MAX(0.1, 0.5*area*(tv[2]-ad[2])); \
	rhsl[3] = MAX(0.1, 0.5*area*(tv[3]-ad[3])); \
	rhsl[4] = MAX(0.1, 0.5*area*(tv[4]-ad[4])); \
}

//define 2 function pointers
define_e2v_FunPtr(umbt_conv)
{
	//selfConn computation
	int i,j;
	swFloat* u	 = accessArray(data->vertexData, 1);
	swFloat* rhs = accessArray(data->vertexData, 0);
	swFloat* gra = accessArray(data->vertexData, 2);
	swInt iDim,dims;
	swInt dims_rhs = getArrayDims(data->vertexData, 0);
	
	
	//frontEdge computation
	swFloat* n_vg	 = accessArray(data->frontEdgeData, 0);
	swFloat* cen	 = accessArray(data->frontEdgeData, 1);
	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt dims_vg    = getArrayDims(data->frontEdgeData, 0);
	swInt dims_cen   = getArrayDims(data->frontEdgeData, 1);
	swInt dims_u     = getArrayDims(data->vertexData, 1);
	swInt dims_gra   = getArrayDims(data->vertexData, 2);
	swInt iedge;
	swFloat fac_1d[4], vgf[3], uLf[5],uRf[5], dL[3], dR[3], duL[5], duR[5];
	swFloat limiter_L=1.0, limiter_R=1.0;
	swFloat rhsl[5],rhsD[5];
	for( iedge = 0; iedge < edgeNumber; iedge++)
	{	
		// 36 flops
		fac_1d[0] = n_vg[iedge*dims_vg+0];
		fac_1d[1] = n_vg[iedge*dims_vg+1];
		fac_1d[2] = n_vg[iedge*dims_vg+2];
		fac_1d[3] = n_vg[iedge*dims_vg+3];
		vgf[0]    = n_vg[iedge*dims_vg+5];
		vgf[1]    = n_vg[iedge*dims_vg+6];
		vgf[2]    = n_vg[iedge*dims_vg+7];

		uLf[0] = u[startVertices[iedge]*dims_u+0];
		uLf[1] = u[startVertices[iedge]*dims_u+1];
		uLf[2] = u[startVertices[iedge]*dims_u+2];
		uLf[3] = u[startVertices[iedge]*dims_u+3];
		uLf[4] = u[startVertices[iedge]*dims_u+4];

		uRf[0] = u[endVertices[iedge]*dims_u+0];
		uRf[1] = u[endVertices[iedge]*dims_u+1];
		uRf[2] = u[endVertices[iedge]*dims_u+2];
		uRf[3] = u[endVertices[iedge]*dims_u+3];
		uRf[4] = u[endVertices[iedge]*dims_u+4];

		dL[0] = cen[iedge*dims_cen+0];
		dL[1] = cen[iedge*dims_cen+1];
		dL[2] = cen[iedge*dims_cen+2];

		dR[0] = cen[iedge*dims_cen+0];
		dR[1] = cen[iedge*dims_cen+1];
		dR[2] = cen[iedge*dims_cen+2];

		// 190 flops
		for(j=0;j<5;j++)
		{
			duL[j] = gra[startVertices[iedge]*dims_gra+3*j+0]*dL[0]
				+ gra[startVertices[iedge]*dims_gra+3*j+1]*dL[1]
				+ gra[startVertices[iedge]*dims_gra+3*j+2]*dL[2];
			duR[j] = gra[endVertices[iedge]*dims_gra+3*j+0]*dL[0]
				+ gra[endVertices[iedge]*dims_gra+3*j+1]*dL[1]
				+ gra[endVertices[iedge]*dims_gra+3*j+2]*dL[2];
			if(0)
			{
				limiter_L = gra[startVertices[iedge]*dims_gra+18];
				limiter_R = gra[endVertices[iedge]*dims_gra+18];
			}
			if(1)
			{
				uLf[j] += limiter_L*duL[j];
				uRf[j] += limiter_R*duR[j];
			}
			if(0)
			{
				uLf[j] = u[startVertices[iedge]*dims_u+j];
				uRf[j] = u[endVertices[iedge]*dims_u+j];
			}
		}
		// 273 flops
		if(1)
		{
			rhs_conv_roe();
		}
		// 15 flops
		for(j=0;j<5;j++)
		{
			//rhs[startVertices[iedge]*dims_rhs+j] += duR[j];
			rhs[startVertices[iedge]*dims_rhs+j] += rhsl[j];
			//rhs[endVertices[iedge]*dims_rhs+j]   -= 1.4;
			rhs[endVertices[iedge]*dims_rhs+j]   -= rhsl[j];
		}

	}
}


