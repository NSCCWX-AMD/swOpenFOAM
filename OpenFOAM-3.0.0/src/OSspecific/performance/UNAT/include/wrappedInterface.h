#ifndef WRAPPEDINTERFACE_H
#define WRAPPEDINTERFACE_H

//#ifdef __cplusplus
//extern "C"
//{
//#endif

#include "topology.hpp"
#include <sys/time.h>

#include "integrate.h"
#include "interpolation.h"
#include "spMV.h"
#include "spMV_host.hpp"
#include "u_to_akh.h"
#include "negSumDiag.h"

#include "calculateVisFlux.h"
#include "calculateVisFlux_host.hpp"
#include "interpolateViscosity.h"
#include "interpolateViscosity_host.hpp"
#include "calculateUvwFlux.h"
#include "calculateUvwFlux_host.hpp"
#include "calcLudsFcc.h"
#include "calcLudsFcc_host.hpp"
#include "calculateFccFlux.h"
#include "calculateFccFlux_host.hpp"

#include "merge5.h"

#include "merge4_interpolation.h"
#include "merge4_vector.h"

#include "separate_vector.h"
#include "separate_vector_host.hpp"
#include "separate_interpolation.h"
#include "separate_interpolation_host.hpp"
#include "separate_integrate.h"
#include "separate_integrate_host.hpp"

#include "umbt_conv.h"
#include "umbt_conv_host.hpp"

//#ifdef __cplusplus
//}
//#endif

#endif
