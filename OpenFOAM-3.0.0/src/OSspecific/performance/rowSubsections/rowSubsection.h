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

#include "sunwayMacros.h"
#include "basicTypes.h"i

typedef swInt LABEL;
typedef swFloat SCALAR;

#ifdef __cplusplus
extern "C"
{
#endif

#define MAX_COL_SUBSECTIONS 64
struct rowSubsection
{
    // The lower bound for face field DMA
    LABEL faceStart;
    // The number of data for face field DMA
    LABEL nFaces;
    // the round at which to write columes
    LABEL colRound;
    // The number of subsections of columns(neighbors)
    LABEL nSecs;
    // The array holds the subsection "lower bound and data count" pairs of columns(neighbors)
    //LABEL colStartsAndCounts[2*MAX_COL_SUBSECTIONS];
    LABEL *colStartsAndCounts;
};
//#undef MAX_COL_SUBSECTIONS

#ifdef __cplusplus
}
#endif



#endif
