/*=======================================================================*/
/*
 * The row oriented subsection decomposer
 * Hu Ren
 * renhu@nsccwx.mail.cn
 * latest modification on 2017-08-02
 */
/*========================================================================*/

#ifndef rowSubsection_H
#define rowSubsection_H

#include "swMacro.h"

#ifdef FOAM_LABEL64
    typedef swInt64 label;
#else
    typedef swInt32 label;
#endif

#ifdef WM_DP
    typedef swFloat64 scalar;
#else
    typedef swFloat32 scalar;
#endif

#ifdef __cplusplus
extern "C"
{
#endif

#define MAX_COL_SUBSECTIONS 64
struct rowSubsection
{
    // The lower bound for face field DMA
    label faceStart;
    // The number of data for face field DMA
    label nFaces;
	// The lower bound for cell field DMA
	label cellStart;
	// The number of data for cell field DMA
	label nCells;
    // the round at which to write columes
    label colRound;
    // The number of subsections of columns(neighbors)
    label nSecs;
    // The array holds the subsection "lower bound and data count" pairs of columns(neighbors)
    //label colStartsAndCounts[2*MAX_COL_SUBSECTIONS];
    label *colStartsAndCounts;
};
//#undef MAX_COL_SUBSECTIONS

#ifdef __cplusplus
}
#endif



#endif
