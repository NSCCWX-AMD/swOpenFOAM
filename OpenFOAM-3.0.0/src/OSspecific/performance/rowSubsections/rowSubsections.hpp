/*=======================================================================*/
/*
 * The row oriented subsection decomposer
 * Hu Ren
 * renhu@nsccwx.mail.cn
 * latest modification on 2017-08-02
 */
/*========================================================================*/

// Algorithm:
//  Often, the unstructured mesh faces or sparse matrix coefficients
//  are ROW oriented COO or CSR format. This class is aimed at
//  decomposing the faces to segment acording the core numbers in
//  many-core processors and then subsections acording the max faces
//  data you expect to copy into device memory at each computing round.
//  This class uses the face row projections to locate the column
//  sections and use the grouping strategy to avoid column conflicts
//  while writing back arrays.
//
// Usage:
//  parameters for constructor
//      const label nFace:              The face field length
//      const label nCell:              The cell field length
//      const label nCPE:               The computing core number
//      const label nMC:                The memory channel number
//      const label * lPtr:             The row(owner) index pointer
//      const label * uPtr:             The column(neighbor) index pointer
//      const label maxFaces:           The max face counts for LDM to include all data array needed
//      const label subsectionSize:     The subsection length to do agglomeration


#ifndef rowSubsections_HPP
#define rowSubsections_HPP

#include "stdlib.h"
#include "rowSubsection.h"


namespace SWFoam
{

class RowSubsections
{
    // the rowSubsection pointer
    rowSubsection ** _subsections;
    // the memory array holds subsections
    rowSubsection* _subsectionsHolder;

    label _segNum; // total segment number
    label _secNumInSeg; // subsection number in each segment

    // the matrix band width
    label _bandwidth;
    label _colRoundNum;

    // disable the copy constructors
    RowSubsections(const RowSubsections&);

    void operator=(const RowSubsections&);

 public:

    // open data from constructors
    const label _nFace; // The face field length
    const label _nCell; // The cell field length
    const label _nCPE; // The computing core number
    const label _nMC; // The memory channel number
    const label * _lPtr; // The row(owner) index pointer
    const label * _uPtr; // The column(neighbor) index pointer
    const label _maxFaces; // The max face counts for LDM to include all data array needed
    const label _dataSize; // The subsection length to do agglomeration

 public:

    // constructors
    RowSubsections(
                      const label nFace, // The face field length
                      const label nCell, // The cell field length
                      const label nCPE, // The computing core number
                      const label nMC, // The memory channel number
                      const label * lPtr, // The row(owner) index pointer
                      const label * uPtr, // The column(neighbor) index pointer
                      const label maxFaces, // The max face counts for LDM to include all data array needed
                      const label dataSize// The subsection length to do agglomeration
                  );
    // deconstructors
    virtual ~RowSubsections()
    {
        free(_subsectionsHolder);
        free(_subsections);
    }

    // print subsections on standard i/o terminal
    virtual void stdPrint() const;

    const rowSubsection** getSubsections() const
    { return (const rowSubsection**) _subsections; }

    const rowSubsection* getSubsections( const label iSeg, const label iSec) const
    {
        return &_subsections[iSeg][iSec];
    }

    // return the subsection number in each CPE segment
    const label& getSecNumInSeg() const { return _secNumInSeg; }
    // return the matrix bandwidth
    const label& getBandwidth() const { return _bandwidth; }
    // return the column write round number
    const label& getColRoundNum() const { return _colRoundNum; }

    // verify this rowSubsections
    void verify();

 private:
    // decompose the face set into CPE odd/even segments, and then subsections
    virtual void decomposeFaces();
    // project the columns(neighbors) for each subsection
    virtual void projectCols();
    virtual void projectCols(rowSubsection **);
    // check the column subsection overlap of odd/even segments
    virtual void checkColOverlap();
};


} // end namespace SWFoam


#endif
