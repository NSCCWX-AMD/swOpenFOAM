/*=======================================================================*/
/*
 * The row oriented subsection decomposer
 * Hu Ren
 * renhu@nsccwx.mail.cn
 * latest modification on 2017-08-02
 */
/*========================================================================*/





#include "rowSubsections.hpp"
// #include "messageStream.H"
#include "stdio.h"
#include "stdlib.h"
#include "mpi.h"
#include <iostream>
#include <vector>
#include <set>

#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define SUBSECTIONSSIZE 512
#define MIN_EFFECT_RATIO 0.3
#define MIN_COLSEC_SIZE 32

//using namespace std;

SWFoam::RowSubsections::RowSubsections(
        const label nFace, // The face field length
        const label nCell, // The cell field length
        const label nCPE, // The computing core number
        const label nMC, // The memory channel number
        const label * lPtr, // The row(owner) index pointer
        const label * uPtr, // The column(neighbor) index pointer
        const label maxFaces, // The max face counts for LDM to include all data array needed
        const label dataSize // The subsection length to do agglomeration
        )
:
    _nFace(nFace),
    _nCell(nCell),
    _nCPE(nCPE),
    _nMC(nMC),
    _lPtr(lPtr),
    _uPtr(uPtr),
    _maxFaces(maxFaces),
    _dataSize(dataSize)
{
    if(_nCPE % _nMC != 0)
    {
        printf("\n***Error: CPE number should be integer times of MC number\n\n");
        exit(1);
    }
    decomposeFaces();
    projectCols();
    checkColOverlap();
#if FULL_DEBUG
    //* self-check
    verify();
#endif
}

void SWFoam::RowSubsections::verify()
{
	//* arrays to record the lower and upper
	//* bound of the column indices in each
	//* subsection
	label colStart[_segNum][_secNumInSeg];
	label colEnd[_segNum][_secNumInSeg];

//* std error output macro
#define CERR std::cout<<"***Error(seg " \
			<<segIter<<", sec " \
			<<secIter<<"): "
#define ENDL "***"<<std::endl

	//* check the inner subsections
	for(label segIter=0; segIter<_segNum; segIter++)
		for(label secIter=0; secIter<_secNumInSeg; secIter++){
			//** take care of the first subsection
			label preFaceEnd;
			if((secIter==0 && segIter==0))
				preFaceEnd = 0;
			else if(secIter == 0)
				preFaceEnd =
					_subsections[segIter-1][_secNumInSeg-1].faceStart
					+_subsections[segIter-1][_secNumInSeg-1].nFaces;
			else
				preFaceEnd =
					_subsections[segIter][secIter-1].faceStart
					+_subsections[segIter][secIter-1].nFaces;

			rowSubsection* sec = &_subsections[segIter][secIter];
			//** check the face start, subsection length
			//** and row decomposition
			if(sec->faceStart != preFaceEnd)
				CERR<<"subsection face is not continuous"
					<<ENDL;
			if( sec->nFaces <= 0 )
				CERR<<"subsection length is smaller than 1"
					<<ENDL;
			if( sec->faceStart != 0 &&
                    secIter == 0 &&
					_lPtr[sec->faceStart] == _lPtr[sec->faceStart-1])
				CERR<<"segment is not decomposed by rows"
					<<ENDL;
			//** check the column subsections
			label* colSandC = sec->colStartsAndCounts;
			label minColStart=_nFace, maxColEnd=0;
			for(label subSecIter=0;
				subSecIter< sec->nSecs;
				subSecIter++)
			{
				//* check the overlap
				label i = subSecIter;
				while(i--)
					if(!(colSandC[2*subSecIter] >= colSandC[2*i]+colSandC[2*i+1] ||
							colSandC[2*subSecIter]+colSandC[2*subSecIter+1] <= colSandC[2*i]))
						CERR<<"column subsection "<<subSecIter
							<<" overlaped with "<<i
							<<ENDL;
				//* check the length
				if(colSandC[2*subSecIter+1] <= 0)
					CERR<<"column subsection "<<subSecIter
						<<" has a length no larger than 0"
						<<ENDL;
				//* update the column subsection span
				minColStart = MIN(minColStart, colSandC[2*subSecIter]);
				maxColEnd   = MAX(maxColEnd, colSandC[2*subSecIter]+colSandC[2*subSecIter+1]);
			}
			//** record the column subsection span
			colStart[segIter][secIter] = minColStart;
			colEnd[segIter][secIter]   = maxColEnd;
			//** check the column index for completeness
			for(label faceIter=sec->faceStart;
				faceIter<sec->faceStart+sec->nFaces;
				faceIter++)
			{
				label colIndex = _uPtr[faceIter];
				label subSecIter=0;
				while(subSecIter<sec->nSecs)
					if(!(colIndex>=colSandC[2*subSecIter] &&
							colIndex<colSandC[2*subSecIter]+colSandC[2*subSecIter+1]))
						subSecIter++;
                    else
                        break;

				if(subSecIter == sec->nSecs)
					CERR<<"column index was not in column subsections"<<ENDL;
			}
		}

	//* check the last subsection
	//** check the face end of the last subsection
	rowSubsection* lastSec = &_subsections[_segNum-1][_secNumInSeg-1];
	if(lastSec->faceStart+lastSec->nFaces != _nFace)
		std::cout<<"subsection is not ended at last face"
			<<std::endl;

	//* check the overlap between rounds
	for(label secIter=0; secIter<_secNumInSeg; secIter++)
		for(label segIter=0; segIter<_segNum; segIter++){
			label i = segIter;
			while(i-_colRoundNum>=0){
				i -= _colRoundNum;
				if(!(colStart[segIter][secIter]>=colEnd[i][secIter] ||
						colEnd[segIter][secIter]<=colStart[i][secIter]))
					CERR<<"segment "<<segIter
                        <<"("<<colStart[segIter][secIter]<<", "
                        <<colEnd[segIter][secIter]<<")"
						<<" overlapped with segment "<<i
                        <<"("<<colStart[i][secIter]<<", "
                        <<colEnd[i][secIter]<<")"<<ENDL;
			}
		}
}


void SWFoam::RowSubsections::decomposeFaces()
{
    _segNum  = 1 * _nCPE;
    if( _segNum * 10 >  _nCell )
    {
        printf("\n***Error: too few cells(rows) to creat sufficient segments !\n\n");
        exit(1);
    }

    // average length of  segment
    const label avgSegLength = label( (scalar)_nFace/(scalar)_segNum );
    label maxSegLength = 0; // max length of segment
#if FULL_DEBUG
    printf("avgSegLength is %d\n", avgSegLength);
#endif
    label segStart[_segNum+1]; // start point to decompose face set

    /* decompose the face set into segment by rows */
    segStart[0] = 0;
    for(int i = 1; i < _segNum; i++)
    {
        label tmpStart = segStart[i-1] + avgSegLength;
        if( tmpStart >= _nFace - 1 )
        {
            printf("\n***Error: segment start exceed the face number!\n\n");
            exit(1);
        }
        if( _lPtr[tmpStart] == _nCell-1 && _lPtr[tmpStart] == _lPtr[tmpStart-1] )
        {
            printf("\n***Error: segment start exceed the face number!\n\n");
            exit(1);
        }
        while( _lPtr[tmpStart] == _lPtr[tmpStart-1] ) tmpStart++;
        segStart[i] = tmpStart;
        // update the max length of segment
        maxSegLength = (segStart[i] - segStart[i-1] > maxSegLength) ?
                        segStart[i] - segStart[i-1] : maxSegLength;
    }
    segStart[_segNum] = _nFace;
    maxSegLength = (segStart[_segNum] - segStart[_segNum-1] > maxSegLength) ?
                    segStart[_segNum] - segStart[_segNum-1] : maxSegLength;

#if FULL_DEBUG
    printf("the final segStarts are:\n");
    for(int i = 0; i <= _segNum; i++)
    {
        printf("%-8d", segStart[i]);
    }
    printf("\n");
#endif

    /* initialize subsection arrays */
    // up round the number of subsections in each segment with
    // the max segment length
    _secNumInSeg = (maxSegLength+_maxFaces-1)/_maxFaces;

    const label memory4Seg = _secNumInSeg * sizeof(rowSubsection);
    const label memorySize = _segNum * memory4Seg;

    _subsections  = (rowSubsection **) malloc(_segNum * sizeof(void*));
    _subsectionsHolder  = (rowSubsection *)  malloc(memorySize);

    for(int i = 0; i < _segNum; i++)
    {
        _subsections[i]  = &_subsectionsHolder[i * _secNumInSeg];
    }

    /* calculate face start points and face counts for each subsection */
    for(int i = 0; i < _segNum; i++ )
    {
        _subsections[i][0].faceStart  = segStart[i];

        for(int j = 1; j < _secNumInSeg; j++)
        {
            // up round the average section length
            label avgSecSpan = ( segStart[i+1] - _subsections[i][j-1].faceStart + _secNumInSeg -j )
                                / (_secNumInSeg - j + 1);

            label tmpStart   = _subsections[i][j-1].faceStart  + avgSecSpan;
            if( tmpStart >= segStart[i+1] )
            {
                printf("\n***Error: subsection face starts exceed length of segment %d !\n\n", i);
                exit(1);
            }
            if( _lPtr[tmpStart] == _lPtr[segStart[i+1]-1] && _lPtr[tmpStart] == _lPtr[tmpStart-1] )
            {
                printf("\n***Error: subsection face starts exceed rows in segment %d !\n\n", i);
                exit(1);
            }
            //while( _lPtr[tmpStart]  == _lPtr[tmpStart-1] )  tmpStart++;
            _subsections[i][j].faceStart  = tmpStart;
            _subsections[i][j-1].nFaces   = tmpStart - _subsections[i][j-1].faceStart;
        }
        _subsections[i][_secNumInSeg-1].nFaces = segStart[i+1] - _subsections[i][_secNumInSeg - 1].faceStart;
    }

#if FULL_DEBUG
    for(int i = 0; i < _segNum; i++ )
        for(int j = 0; j < _secNumInSeg; j++)
            std::cout<<"subsection "<<j<<" in segment "<<i<<", faceStart is "
                <<_subsections[i][j].faceStart<<", nFaces is "<<_subsections[i][j].nFaces<<std::endl;
#endif
}


void SWFoam::RowSubsections::projectCols()
{
    projectCols(_subsections);
}

void SWFoam::RowSubsections::projectCols( rowSubsection ** subsection )
{
    // initialize the bandWidth
    _bandwidth = 1;

    /* project the column(neighbor) subsection */
    for(int i = 0; i < _segNum; i++)
	{
		int subSectionSizeTmp = SUBSECTIONSSIZE;
        for(int j = 0; j < _secNumInSeg; j++)
        {
			std::set<label>* effectiveColumns = new std::set<label>;
//            label lAndUBounds[sizeof(subsection[0][0].colStartsAndCounts)/sizeof(label)/2][2];
			std::vector<std::vector<int> > lAndUBounds(2);
//			std::vector<std::vector<int> > colStartsAndCounts(2);

            label faceStart = subsection[i][j].faceStart;
            label nFaces    = subsection[i][j].nFaces;
            label maxSecNum;


//            lAndUBounds[0][0] = _uPtr[faceStart];
//            lAndUBounds[0][1] = lAndUBounds[0][0];
			lAndUBounds[0].push_back(_uPtr[faceStart]);
			lAndUBounds[1].push_back(lAndUBounds[0][0]);
            maxSecNum = 1;
#define ABS(expr) \
    (((expr) > 0) ? (expr) : (-(expr)))
            // update the bandWidth
            label tmpBandwidth = ABS(lAndUBounds[0][0] - _lPtr[faceStart]);
            _bandwidth = ( _bandwidth < tmpBandwidth ) ? tmpBandwidth : _bandwidth;


            for(int faceDelta = 1; faceDelta < nFaces; faceDelta++)
            {
                label colIndex = _uPtr[faceStart+faceDelta];

				effectiveColumns->insert(colIndex);

                // update the bandWidth
                tmpBandwidth = ABS(colIndex - _lPtr[faceStart+faceDelta]);
                _bandwidth = ( _bandwidth < tmpBandwidth ) ? tmpBandwidth : _bandwidth;

                label secIndex = 0;
                label distance = ABS( colIndex - (lAndUBounds[0][0]+lAndUBounds[1][0]) / 2 );
                for(int k = 0; k < maxSecNum; k++)
                {
                    if( colIndex <= lAndUBounds[1][k] && colIndex >= lAndUBounds[0][k])
                    {
                        secIndex = k;
                        distance = 0;
                        break;
                    }
                    label tmpDist = MIN( ABS( colIndex - (lAndUBounds[0][k]+lAndUBounds[1][k]) / 2 ), distance );
                    if (  tmpDist < distance )
                    {
                        secIndex = k;
                        distance = tmpDist;
                    }
                }
#undef ABS
                if( distance < subSectionSizeTmp/_dataSize )
                {
                    if( colIndex > lAndUBounds[1][secIndex] )
                        lAndUBounds[1][secIndex] = colIndex;
                    else if(colIndex < lAndUBounds[0][secIndex])
                        lAndUBounds[0][secIndex] = colIndex;
                }
                else
                {
//                    if ( (unsigned)maxSecNum + 1  > sizeof(subsection[0][0].colStartsAndCounts)/sizeof(label)/2 )
//                    {
//                        printf("\n***Error: column subsection number exceeds the max value in segment %d subsection %d \n\n", i, j);
//                        exit(1);
//                    }
                    lAndUBounds[0].push_back(colIndex);
                    lAndUBounds[1].push_back(colIndex);
//                    lAndUBounds[maxSecNum][1] = colIndex;
                    maxSecNum++;
                }
            }

            subsection[i][j].nSecs = maxSecNum;
			subsection[i][j].colStartsAndCounts = (label*)malloc(maxSecNum*sizeof(label)*2);
			label totalColNumber = 0;
            for(int k = 0; k < maxSecNum; k++)
            {
                subsection[i][j].colStartsAndCounts[2*k]   = lAndUBounds[0][k];
                //subsection[i][j].colSAndC[2*k]   = lAndUBounds[0][k];
                subsection[i][j].colStartsAndCounts[2*k+1] =
                //subsection[i][j].colSAndC[2*k+1] =
                    lAndUBounds[1][k] - lAndUBounds[0][k] + 1;
				totalColNumber += subsection[i][j].colStartsAndCounts[2*k+1];
            }
			label colCounts = effectiveColumns->size();
			effectiveColumns->clear();
			delete effectiveColumns;	
			// if the total columns in column section is bigger than face number
			// reproject the present subsection with the half subSectionSize
			if( MIN_EFFECT_RATIO*totalColNumber > colCounts && totalColNumber > nFaces) 
			{
//printf("too large subSectionSizeTmp %d, in seg %d in sec %d, half!\n", totalColNumber, i, j);	
//printf("SUBSECTIONSIZE %d, nFaces %d, _dataSize %d, maxSecNum %d\n", SUBSECTIONSSIZE, nFaces, _dataSize, maxSecNum);	
				j--;
				subSectionSizeTmp /= 2;
				if( subSectionSizeTmp < MIN_COLSEC_SIZE )
				{
					printf("***Error: too sparse column subsection distribution!\n");
					exit(-1);
				}
				//break;
			}
			else
				subSectionSizeTmp = SUBSECTIONSSIZE;
        }
	}
#if FULL_DEBUG
        for(int j = 0; j < _secNumInSeg; j++)
            for(int i = 0; i < _segNum; i++)
            {
                printf("In seg %d: ", i);
                for(int k = 0; k < _subsections[i][j].nSecs; k++)
                    printf("sec%d: %d, %dn; ", k, _subsections[i][j].colStartsAndCounts[2*k], _subsections[i][j].colStartsAndCounts[2*k+1]);
                printf("\n");
            }
        printf("Max column subsection size is %d\n\n", SUBSECTIONSSIZE);
#endif
}

void SWFoam::RowSubsections::checkColOverlap()
{
    // calculate the min segment span and max round number
    label segSpan  = label((scalar)_bandwidth / ((scalar)_nCell/(scalar)_segNum) + 0.5);
    label groupNum = (_segNum%segSpan == 0) ? (_segNum/segSpan) : (_segNum/segSpan + 1);
#if FULL_DEBUG
    printf("The matrix bandwidth is %d, and the average rows for each segment is %d\n",
            _bandwidth, label((scalar)_nCell/(scalar)_segNum));
#endif

    // round groupNum to be integer times of memory channel number, if possible
    if( groupNum > _nMC )
    {
        groupNum = (groupNum/_nMC) * _nMC;
        segSpan  = (_segNum%groupNum == 0) ? (_segNum/groupNum) : (_segNum/groupNum + 1);
    }

    // check overlaps
    bool overlapped;
    do
    {
        overlapped = false;
        // check each subsections in each segment to make sure no overlaps
        for(int j = 0; j < _secNumInSeg; j++)
        {
            for(int i = 0; i < _segNum; i++)
            {
                // get the subjective subsection infomation
                const label & secNum = _subsections[i][j].nSecs;
                const label * secs = _subsections[i][j].colStartsAndCounts;
                // check with the subsections after the subjective one
                for(int si = i+segSpan; si < _segNum; si+=segSpan)
                {
                    // get the objective subsection information
                    const label & psecNum = _subsections[si][j].nSecs;
                    const label * psecs = _subsections[si][j].colStartsAndCounts;
                    // check each column sections in subjective subsection
                    // with the ones in objective subsections
                    for(int k = 0; k < secNum; k++)
                    {
                        for(int sk = 0; sk < psecNum; sk++)
                        {
                            if( !(secs[2*k] >= psecs[2*sk]+psecs[2*sk+1] ||
                                  secs[2*k]+secs[2*k+1] <= psecs[2*sk])
                              )
                            {
                                overlapped = true;
                                break;
                            }
                        }
                        if(overlapped)
                        {
                            break;
                        }
                    }
                    if(overlapped) break;
                }
                if(overlapped) break;
            }
            if(overlapped) break;
        }
        if(overlapped)
        {
            if(groupNum <= _nMC)
            {
                segSpan++;
                groupNum = (_segNum%segSpan == 0) ? (_segNum/segSpan) : (_segNum/segSpan + 1);
            }
            else
            {
                groupNum = groupNum - _nMC;
                segSpan  = (_segNum%groupNum == 0) ? (_segNum/groupNum) : (_segNum/groupNum + 1);
            }
        }
    } while(overlapped);
#if FULL_DEBUG
    printf("The final sagment span is %d, final group number is %d\n", segSpan, groupNum);
#endif
    if(groupNum < _nMC) printf(
            "\n***Warning: too big matrix bandwidth to utilize all %d momery channel\n\n",
            _nMC);

    _colRoundNum = segSpan;
    for(int i = 0; i < _segNum; i++)
    {
        label round = i%segSpan;
        for(int j = 0; j < _secNumInSeg; j++)
            _subsections[i][j].colRound = round;
    }
}

void SWFoam::RowSubsections::stdPrint() const
{
    if( _nCell > 200)
    {
        printf("\n***Error: stdPrint not suitable for more than 200 cells!\n\n");
    }
    else{
        printf("==================================================================\n");
        printf("                    The subsetions decoposition                   \n");
        printf("==================================================================\n");
        printf("FORMAT:\n");
        printf("subsection ID\n    non-zero positions\n    column subsections\n");
        printf("intersected segments in each round\n\n");

        for(int i = 0; i < _segNum; i++)
            for(int j = 0; j < _secNumInSeg; j++)
            {
                printf("subsection %d in segment %d:\n", j, i);
                label nFaces    = _subsections[i][j].nFaces;
                label faceStart = _subsections[i][j].faceStart;

                // Print face(coefficients) positions
                printf("\r%3d|\033[%dCX", _lPtr[faceStart], _uPtr[faceStart]);
                for(int k = 1; k < nFaces; k++)
                {
                    label rowIndex = _lPtr[ faceStart+k ];
                    label colIndex = _uPtr[ faceStart+k ];
                    if( rowIndex != _lPtr[faceStart+k-1] )
                        printf("\n%3d|", rowIndex);
                    printf("\r\033[%dCX", colIndex + 4);
                }
                printf("\n");


                // Print the column subsections
                const label & secNum = _subsections[i][j].nSecs;
                const label * secs   = _subsections[i][j].colStartsAndCounts;
                for(int k = 0; k < secNum; k++)
                {
                    printf("\r\033[%dC", secs[2*k]+4);
                    printf("\b\b\b%3d", secs[2*k]);
                    for( int sk = 0; sk < secs[2*k+1] ; sk++)
                        printf("-");
                    printf("%-3d", secs[2*k]+secs[2*k+1]-1);
                    printf("\n");
                }
            }
    }
    // Print the colume write priority
    printf("\n------------------------------------------------------------------\n");
    printf("column write priority for each segment: \n");
    for(int i = 0; i < _segNum; i++)
    {
        printf("(%d,%d) ", i,_subsections[i][0].colRound);
    }
    printf("\n");
}

#undef MIN
#undef MAX












