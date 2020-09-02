#include "sunwayMacros.h"
#ifdef __cplusplus
extern "C" {
#endif

void interpolate1(
	swInt faceNum,
	swInt cellNum,
	const void* fvPtr,
	const void* cvPtr,
	const swFloat* weiPtr,
	const swFloat* _1wPtr,
	const swInt* ownPtr,
	const swInt* neiPtr,
	swInt typeSize);

#ifdef __cplusplus
}
#endif
