#ifndef USERFUNC_SAXPY_H
#define USERFUNC_SAXPY_H
#include "vectorOps_struct.h"

#define SLAVE_FUNC(funcname) slave_##funcname

extern void SLAVE_FUNC(vectorOps_slave)(MVM_ParametersPack*);
extern void SLAVE_FUNC(vectorCopy_slave)(MVM_ParametersPack*);
extern void SLAVE_FUNC(scalingFactor_slave)(MVM_ParametersPack*);
extern void SLAVE_FUNC(gSum_slave)(MVM_ParametersPack*);
extern void SLAVE_FUNC(residualNormFactor_slave)(MVM_ParametersPack*);

void SLAVE_FUNC(userFunc_aEbPk1Mua)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_aEaMik1Mub)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_aEaPk1Mub)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_aEbMic)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_aEk1MuaPk2MubMuScMidS)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_aE1Db)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_aEbMuc)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_aEbMuScMiaSMuk1)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_aEaPb)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_jacobi)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_aEk1Mua)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_aEbPk1MuSaMik2MucS)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_aEbMik1Muc)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_aEaPk1MubPk2Muc)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_residualSum)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_residualSumK)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_digPrecondSum)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_sumProd)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_sumSqr)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_aEcPk1Mua_bEdPk1Mub)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_aEaPk1Muc_bEbMik1Mud)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_aEcMuSdMiaSMuk1_bEbPa)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_sumFabs)(MVM_Arrays*);
void SLAVE_FUNC(userFunc_aEbDc)(MVM_Arrays*);

#endif
