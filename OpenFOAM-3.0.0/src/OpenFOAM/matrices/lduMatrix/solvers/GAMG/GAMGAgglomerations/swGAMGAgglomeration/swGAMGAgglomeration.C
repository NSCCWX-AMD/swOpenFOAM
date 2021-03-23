#include "GAMGAgglomeration.H"
#include "swGAMGAgglomeration.H"
#include "mapDistribute.H"
#include "globalIndex.H"
#include "label.H"
#include "mpi.h"


#include <fstream>

#define swForAll(i, size) for(int i=0; i<size; ++i)
namespace Foam
{
	//- restrict
	bool* swRestInterMap::restFirstUse_ = NULL;
	restStruct* swRestInterMap::restStructLevels_ = NULL;

	//- face restrict
	bool* swRestInterMap::faceRestFirstUse_ = NULL;
	restStruct* swRestInterMap::faceRestStructLevels_ = NULL;

	//- interpolate
	bool* swRestInterMap::interFirstUse_ = NULL;
	interStruct* swRestInterMap::interStructLevels_ = NULL;

	//- agglomerate matrix upper
	bool* swRestInterMap::aggMatrixUpperFirstUse_ = NULL;
	aggMatrixUpperStruct* swRestInterMap::aggMatrixUpperStructLevels_ = NULL;


	label swRestInterMap::bandSize_ = 6250;
	label swRestInterMap::minCellsUsingSW_ = 10000;
}

Foam::swRestInterMap::swRestInterMap(const GAMGAgglomeration& aggl)
	:
		aggl_(aggl)
{
	// creat std::vector<std::vector<int>> restrictAddressing_int32 
	// from PtrList<labelField> restrictAddressing_
	restrictAddressing_int32.resize(100);
	
}

void Foam::swRestInterMap::initRestInterMap()
{
	const label levelSize = aggl_.size();

	if(!restFirstUse_)
	{
		restFirstUse_ = new bool [levelSize];
		restStructLevels_ = new restStruct [levelSize];

		swForAll(i, levelSize)
		{
			restFirstUse_[i] = true;
		}
	}

	if(!faceRestFirstUse_)
	{
		faceRestFirstUse_ = new bool [levelSize];
		faceRestStructLevels_ = new restStruct [levelSize];

		swForAll(i, levelSize)
		{
			faceRestFirstUse_[i] = true;
		}
	}


	if(!interFirstUse_)
	{
		interFirstUse_ = new bool [levelSize];
		interStructLevels_ = new interStruct [levelSize];

		swForAll(i, levelSize)
		{
			interFirstUse_[i]  = true;
		}
	}

	if(!aggMatrixUpperFirstUse_)
	{
		aggMatrixUpperFirstUse_ = new bool [levelSize];
		aggMatrixUpperStructLevels_ = new aggMatrixUpperStruct [levelSize];

		swForAll(i, levelSize)
		{
			aggMatrixUpperFirstUse_[i]  = true;
		}
	}
}


void Foam::swRestInterMap::initRestStruct
(
	Field<scalar>& cf,
	const Field<scalar>& ff,
	const label fineLevelIndex
)
{
	if(restFirstUse_[fineLevelIndex])
	{
		//const swInt* restrictMap  = aggl_.restrictAddressing(fineLevelIndex).begin();
		const labelField& raField = aggl_.restrictAddressing(fineLevelIndex);
		restrictAddressing_int32[fineLevelIndex].clear();
		restrictAddressing_int32[fineLevelIndex].resize( raField.size() );

		for (int ilf = 0; ilf < raField.size(); ++ilf)
		{
			restrictAddressing_int32[fineLevelIndex][ilf] = raField[ilf];
		}
		const swInt* restrictMap  = &restrictAddressing_int32[fineLevelIndex][0];
#if 0
{
		static int icheck = 0;
		icheck++;
		int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		//printf("check restrictAddressing_int32\n");
        char filename[256];
        sprintf(filename,"debug/processor%d.dat",rank);
        //FILE* debug_fp=fopen(filename,"w");
        //fprintf(debug_fp,"date: %s\n",__DATE__);
        //fprintf(debug_fp,"time: %s\n",__TIME__);
        int errorNum = 0;
        for(int i = 0 ; i < restrictAddressing_int32[fineLevelIndex].size() ; i++) 
        {
        	//if(restrictAddressing_int32[fineLevelIndex][i]!=raField[i])
        	if(*(restrictMap+i)!=raField[i])
        	{
        		errorNum++;
        		//fprintf(debug_fp, "error\n");
        	}
            /*fprintf(debug_fp,"i=%d:  %d, %ld\n",
                i,restrictAddressing_int32[fineLevelIndex][i], raField[i]);*/
        }
        //fclose(debug_fp);
        printf("rank%d: checkinitRestStruct, icheck=%d, errorNum=%d\n",rank, icheck, errorNum);
        //printf("end check restrictAddressing_int32\n");
}
#endif

		const swInt fSize = ff.size();
		const swInt cSize = cf.size();
		swInt slaveCycles = cSize / (bandSize_*64);
	    if(cSize % (bandSize_*64)) ++slaveCycles;
	    swInt slaveCores  = 64 * slaveCycles;

	    swInt remainder = cSize % slaveCores;
	    swInt lenShort  = cSize / slaveCores;
	    swInt lenLong   = lenShort + 1;

	    swInt fRemainder = fSize % slaveCores;
	    swInt fLenShort  = fSize / slaveCores;
	    swInt fLenLong   = fLenShort + 1;

	    //- allocate range
	    swInt** range = new swInt*[slaveCores];
	    swForAll(i, slaveCores)
	    {
	        range[i] = new swInt[4];
	    }

	    swForAll(i, slaveCores)
	    {
	        if (i < remainder)
	        {
	            range[i][0] = i * lenLong;
	            range[i][1] = range[i][0] + lenLong - 1;
	        }
	        else
	        {
	            range[i][0] = i * lenShort + remainder;
	            range[i][1] = range[i][0] + lenShort - 1;
	        }

	        //- range of fine data estimated
	        if (i < fRemainder)
	        {
	            range[i][2] = i * fLenLong;
	            range[i][3] = range[i][2] + fLenLong - 1;
	        }
	        else
	        {
	            range[i][2] = i * fLenShort + fRemainder;
	            range[i][3] = range[i][2] + fLenShort - 1;
	        }
	    }

	    //- check how many points will not be computed in slave cores
	    swForAll(i, fSize)
	    {
	    	swInt cPos = restrictMap[i];
	        swInt cSlavePos = -1;
	        if(cPos < remainder*lenLong)
	        {
	            cSlavePos = cPos / lenLong;
	        }
	        else
	        {
	            cSlavePos = (cPos - remainder*lenLong) / lenShort + remainder;
	        }

	        if(i < range[cSlavePos][2])
	        {
	            range[cSlavePos][2] = range[cSlavePos][2] < i? range[cSlavePos][2] : i;
	        }
	        else if(i > range[cSlavePos][3])
	        {
	            range[cSlavePos][3] = range[cSlavePos][3] > i? range[cSlavePos][3] : i;
	        }
	    }

	    restStructLevels_[fineLevelIndex].mapPtr = restrictMap;
	    restStructLevels_[fineLevelIndex].localStartEnd = range;
	    restStructLevels_[fineLevelIndex].slaveCycles   = slaveCycles;

		restFirstUse_[fineLevelIndex] = false;

#if 0
{
	int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	printf("rank%d: mapPtrsize %ld, rangsize %dx4, slaveCycle %d, ff %ld, cf %ld\n", rank,
		restrictAddressing_int32[fineLevelIndex].size(), 
		slaveCores, slaveCycles, 
		ff.size(),cf.size());
}
#endif

	}

	restStructLevels_[fineLevelIndex].fPtr = ff.begin();
	restStructLevels_[fineLevelIndex].cPtr = cf.begin();

}


void Foam::swRestInterMap::initInterStruct
(
	Field<scalar>& ff,
    const Field<scalar>& cf,
    const label levelIndex
)
{
	if(interFirstUse_[levelIndex])
	{
		//const swInt* restrictMap  = aggl_.restrictAddressing(levelIndex).begin();
		const labelField& raField = aggl_.restrictAddressing(levelIndex);
		restrictAddressing_int32[levelIndex].clear();
		restrictAddressing_int32[levelIndex].resize( raField.size() );

		for (int ilf = 0; ilf < raField.size(); ++ilf)
		{
			restrictAddressing_int32[levelIndex][ilf] = raField[ilf];
		}
#if 0
{
		static int icheck = 0;
		icheck++;
		int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		//printf("check restrictAddressing_int32\n");
        char filename[256];
        sprintf(filename,"debug/processor%d.dat",rank);
        //FILE* debug_fp=fopen(filename,"w");
        //fprintf(debug_fp,"date: %s\n",__DATE__);
        //fprintf(debug_fp,"time: %s\n",__TIME__);
        int errorNum = 0;
        for(int i = 0 ; i < restrictAddressing_int32[levelIndex].size() ; i++) 
        {
        	if(restrictAddressing_int32[levelIndex][i]!=raField[i])
        	{
        		errorNum++;
        		//fprintf(debug_fp, "error\n");
        	}
            /*fprintf(debug_fp,"i=%d:  %d, %ld\n",
                i,restrictAddressing_int32[levelIndex][i], raField[i]);*/
        }
        //fclose(debug_fp);
        printf("rank%d: check initInterStruct, icheck=%d, errorNum=%d\n",rank, icheck, errorNum);
        //printf("end check restrictAddressing_int32\n");
}
#endif
		const swInt* restrictMap  = &restrictAddressing_int32[levelIndex][0];
		const swInt fSize = ff.size();
		const swInt cSize = cf.size();

		//- translate restrictMap
		swInt* interpolateMap = new swInt[fSize];
		swInt* interMapOffset = new swInt[cSize+1];

		swInt* offsetTemp = new swInt[cSize];

		swForAll(i, cSize)
		{
			offsetTemp[i] = 0;
		}

		swForAll(i, fSize)
		{
			swInt cPos = restrictMap[i];
			offsetTemp[cPos]++;
		}

		interMapOffset[0] = 0;
		swForAll(i, cSize)
		{
			interMapOffset[i+1] = interMapOffset[i] + offsetTemp[i];
			offsetTemp[i] = 0;
		}

		swForAll(i, fSize)
		{
			swInt cPos = restrictMap[i];
			interpolateMap[interMapOffset[cPos]+offsetTemp[cPos]] = i;
			offsetTemp[cPos]++;
		}

		delete []offsetTemp;

		swInt slaveCycles = fSize / (bandSize_*64);
		if(fSize % (bandSize_*64)) ++slaveCycles;
		swInt slaveCores  = 64 * slaveCycles;

	    swInt remainder = fSize % slaveCores;
	    swInt lenShort  = fSize / slaveCores;
	    swInt lenLong   = lenShort + 1;

	    swInt cRemainder = cSize % slaveCores;
	    swInt cLenShort  = cSize / slaveCores;
	    swInt cLenLong   = cLenShort + 1;

	    //- allocate range
	    swInt** range = new swInt*[slaveCores];
	    swForAll(i, slaveCores)
	    {
	        range[i] = new swInt[4];
	    }

	    swForAll(i, slaveCores)
	    {
	        //- range of fine data
	        if (i < remainder)
	        {
	            range[i][0] = i * lenLong;
	            range[i][1] = range[i][0] + lenLong - 1;
	        }
	        else
	        {
	            range[i][0] = i * lenShort + remainder;
	            range[i][1] = range[i][0] + lenShort - 1;
	        }

	        //- range of coarse data estimated
	        if(i < cRemainder)
	        {
	        	range[i][2] = i * cLenLong;
	        	range[i][3] = range[i][2] + cLenLong - 1;
	        }
	        else
	        {
	        	range[i][2] = i * cLenShort + cRemainder;
	            range[i][3] = range[i][2] + cLenShort - 1;
	        }
	    }



	    //- check how many points will not be computed in slave cores
	    swForAll(i, cSize)
	    {
	    	swInt fCells = interMapOffset[i+1] - interMapOffset[i];
	    	swInt fMax = -1;
	    	swInt fMin = 1e+8;

	    	swForAll(j, fCells)
	    	{
	    		swInt fPos = interpolateMap[interMapOffset[i]+j];
	    		fMax = fMax < fPos? fPos : fMax;
	    		fMin = fMin > fPos? fPos : fMin;
	    	}

	    	swInt fMaxSlavePos = -1, fMinSlavePos = -1;
	    	if(fMax < remainder*lenLong)
	    	{
	    		fMaxSlavePos = fMax / lenLong;
	    	}
	    	else
	    	{
	    		fMaxSlavePos = (fMax - remainder*lenLong) / lenShort + remainder;
	    	}

	    	if(fMin < remainder*lenLong)
	    	{
	    		fMinSlavePos = fMin / lenLong;
	    	}
	    	else
	    	{
	    		fMinSlavePos = (fMin - remainder*lenLong) / lenShort + remainder;
	    	}


	    	if(i < range[fMaxSlavePos][2])
	    	{
	    		range[fMaxSlavePos][2] = range[fMaxSlavePos][2] < i? range[fMaxSlavePos][2] : i;
	    	}
	    	else if(i > range[fMaxSlavePos][3])
	    	{
	    		range[fMaxSlavePos][3] = range[fMaxSlavePos][3] > i? range[fMaxSlavePos][3] : i;
	    	}


	    	if(i < range[fMinSlavePos][2])
	    	{
	    		range[fMinSlavePos][2] = range[fMinSlavePos][2] < i? range[fMinSlavePos][2] : i;
	    	}
	    	else if(i > range[fMinSlavePos][3])
	    	{
	    		range[fMinSlavePos][3] = range[fMinSlavePos][3] > i? range[fMinSlavePos][3] : i;
	    	}
	    }

	    interStructLevels_[levelIndex].mapPtr = interpolateMap;
    	interStructLevels_[levelIndex].offsetMapPtr = interMapOffset;
    	interStructLevels_[levelIndex].localStartEnd = range;
    	interStructLevels_[levelIndex].slaveCycles   = slaveCycles;

    	interFirstUse_[levelIndex] = false;
	}

	interStructLevels_[levelIndex].fPtr = ff.begin();
	interStructLevels_[levelIndex].cPtr = cf.begin();
}


template<>
void Foam::GAMGAgglomeration::restrictField
(
    Field<scalar>& cf,
    const Field<scalar>& ff,
    const label fineLevelIndex
) const
{
    cf = pTraits<scalar>::zero;


    if(ff.size() > swRestInterMap::minCellsUsingSW_)
    {
    	swRestInterMap restMap(*this);
    	restMap.initRestStruct(cf, ff, fineLevelIndex);
   
    	restrictData_host(&swRestInterMap::restStructLevels_[fineLevelIndex]);

#if 0
    	// 输出从核计算结果
    	int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        char filename[256];
        sprintf(filename,"debug/processor%d.dat",rank);
        FILE* debug_fp=fopen(filename,"w");
        fprintf(debug_fp,"date: %s\n",__DATE__);
        fprintf(debug_fp,"time: %s\n",__TIME__);
        fprintf(debug_fp,"cf.size()=%ld, i, cf[i]\n", cf.size());
        for(int i = 0 ; i < cf.size() ; i++) 
        {
            fprintf(debug_fp,"%d, %g\n", i,cf[i]);
        }
        fclose(debug_fp);
        std::exit(0);
#endif 
    }
    else
    {
    	const labelList& fineToCoarse = restrictAddressing_[fineLevelIndex];
    	forAll(ff, i)
	    {
	        cf[fineToCoarse[i]] += ff[i];
	    }
/*#if 1
    	// 输出从核计算结果
if(ff.size() > swRestInterMap::minCellsUsingSW_)
{
    	int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        char filename[256];
        sprintf(filename,"debug_host/processor%d.dat",rank);
        FILE* debug_fp=fopen(filename,"w");
        fprintf(debug_fp,"date: %s\n",__DATE__);
        fprintf(debug_fp,"time: %s\n",__TIME__);
        fprintf(debug_fp,"cf.size()=%ld, i, cf[i]\n", cf.size());
        for(int i = 0 ; i < cf.size() ; i++) 
        {
            fprintf(debug_fp,"%d, %g\n", i,cf[i]);
        }
        fclose(debug_fp);
        std::exit(0);
}
#endif*/

    }
}


template<>
void Foam::GAMGAgglomeration::restrictField
(
    Field<scalar>& cf,
    const Field<scalar>& ff,
    const labelList& fineToCoarse
) const
{
    cf = pTraits<scalar>::zero;

    forAll(ff, i)
    {
        cf[fineToCoarse[i]] += ff[i];
    }
}


template<>
void Foam::GAMGAgglomeration::restrictField
(
    Field<scalar>& cf,
    const Field<scalar>& ff,
    const label fineLevelIndex,
    const bool procAgglom
) const
{
    const labelList& fineToCoarse = restrictAddressing_[fineLevelIndex];

    if (!procAgglom && ff.size() != fineToCoarse.size())
    {
        FatalErrorIn
        (
            "void GAMGAgglomeration::restrictField"
            "(Field<scalar>& cf, const Field<scalar>& ff, "
            "const label fineLevelIndex) const"
        )   << "field does not correspond to level " << fineLevelIndex
            << " sizes: field = " << ff.size()
            << " level = " << fineToCoarse.size()
            << abort(FatalError);
    }

    restrictField(cf, ff, fineLevelIndex);

    label coarseLevelIndex = fineLevelIndex+1;

    if (procAgglom && hasProcMesh(coarseLevelIndex))
    {
        label fineComm = UPstream::parent(procCommunicator_[coarseLevelIndex]);

        const List<label>& procIDs = agglomProcIDs(coarseLevelIndex);
        const labelList& offsets = cellOffsets(coarseLevelIndex);

        globalIndex::gather
        (
            offsets,
            fineComm,
            procIDs,
            cf,
            UPstream::msgType(),
            Pstream::nonBlocking    //Pstream::scheduled
        );
    }
}


template<>
void Foam::GAMGAgglomeration::restrictFaceField
(
    Field<scalar>& cf,
    const Field<scalar>& ff,
    const label fineLevelIndex
) const
{
    const labelList& fineToCoarse = faceRestrictAddressing_[fineLevelIndex];

    if (ff.size() != fineToCoarse.size())
    {
        FatalErrorIn
        (
            "void GAMGAgglomeration::restrictFaceField"
            "(Field<scalar>& cf, const Field<scalar>& ff, "
            "const label fineLevelIndex) const"
        )   << "field does not correspond to level " << fineLevelIndex
            << " sizes: field = " << ff.size()
            << " level = " << fineToCoarse.size()
            << abort(FatalError);
    }

    cf = pTraits<scalar>::zero;

    forAll(fineToCoarse, ffacei)
    {
        label cFace = fineToCoarse[ffacei];

        if (cFace >= 0)
        {
            cf[cFace] += ff[ffacei];
        }
    }
}


template<>
void Foam::GAMGAgglomeration::prolongField
(
    Field<scalar>& ff,
    const Field<scalar>& cf,
    const label levelIndex,
    const bool procAgglom
) const
{
    const labelList& fineToCoarse = restrictAddressing_[levelIndex];

    label coarseLevelIndex = levelIndex+1;

    if (procAgglom && hasProcMesh(coarseLevelIndex))
    {
        label coarseComm = UPstream::parent
        (
            procCommunicator_[coarseLevelIndex]
        );

        const List<label>& procIDs = agglomProcIDs(coarseLevelIndex);
        const labelList& offsets = cellOffsets(coarseLevelIndex);

        label localSize = nCells_[levelIndex];

        Field<scalar> allCf(localSize);
        globalIndex::scatter
        (
            offsets,
            coarseComm,
            procIDs,
            cf,
            allCf,
            UPstream::msgType(),
            Pstream::nonBlocking    //Pstream::scheduled
        );

        forAll(fineToCoarse, i)
        {
            ff[i] = allCf[fineToCoarse[i]];
        }
    }
    else
    {
        /*if(ff.size() > swRestInterMap::minCellsUsingSW_)
        {
        	swRestInterMap interMap(*this);
        	interMap.initInterStruct(ff, cf, levelIndex);
   
        	interpolateData_host(&swRestInterMap::interStructLevels_[levelIndex]);
        }
        else
        {*/
        	forAll(fineToCoarse, i)
	        {
	            ff[i] = cf[fineToCoarse[i]];
	        }
        //}
    }
}


// void Foam::swRestInterMap::initAggMatrixUpperStruct
// (
// 	scalarField& coarseUpper,
// 	scalarField& coarseDiag,
// 	const scalarField& fineUpper,
// 	const label& fineLevelIndex
// )
// {
// 	if(aggMatrixUpperFirstUse_[fineLevelIndex])
// 	{
// 		const swInt* restrictMap  = aggl_.faceRestrictAddressing(fineLevelIndex).begin();
// 		const swInt fSize  = fineUpper.size();
// 		const swInt cUSize = coarseUpper.size();
// 		const swInt cDSize = coarseDiag.size();
//
// 		swInt slaveCycles = cSize / (bandSize_*64);
// 	    if(cSize % (bandSize_*64)) ++slaveCycles;
// 	    swInt slaveCores  = 64 * slaveCycles;
//
// 	    swInt remainder = cSize % slaveCores;
// 	    swInt lenShort  = cSize / slaveCores;
// 	    swInt lenLong   = lenShort + 1;
//
// 	    //- allocate range
// 	    swInt** range = new swInt*[slaveCores];
// 	    swForAll(i, slaveCores)
// 	    {
// 	        range[i] = new swInt[4];
// 	    }
//
// 	    swForAll(i, slaveCores)
// 	    {
// 	        if (i < remainder)
// 	        {
// 	            range[i][0] = i * lenLong;
// 	            range[i][1] = range[i][0] + lenLong - 1;
// 	        }
// 	        else
// 	        {
// 	            range[i][0] = i * lenShort + remainder;
// 	            range[i][1] = range[i][0] + lenShort - 1;
// 	        }
//
// 	        //- range of fine data estimated
// 	        range[i][2] = 2 * range[i][0];
// 	        range[i][3] = 2 * range[i][1] + 1;
// 	    }
//
// 	    //- check how many points will not be computed in slave cores
// 	    swForAll(i, fSize)
// 	    {
// 	    	swInt cPos = restrictMap[i];
// 	        swInt cSlavePos = -1;
// 	        if(cPos < remainder*lenLong)
// 	        {
// 	            cSlavePos = cPos / lenLong;
// 	        }
// 	        else
// 	        {
// 	            cSlavePos = (cPos - remainder*lenLong) / lenShort + remainder;
// 	        }
//
// 	        if(i < range[cSlavePos][2])
// 	        {
// 	            range[cSlavePos][2] = range[cSlavePos][2] < i? range[cSlavePos][2] : i;
// 	        }
// 	        else if(i > range[cSlavePos][3])
// 	        {
// 	            range[cSlavePos][3] = range[cSlavePos][3] > i? range[cSlavePos][3] : i;
// 	        }
// 	    }
//
// 	    aggMatrixUpperStructLevels_[fineLevelIndex].mapPtr = restrictMap;
// 	    aggMatrixUpperStructLevels_[fineLevelIndex].localStartEnd = range;
// 	    aggMatrixUpperStructLevels_[fineLevelIndex].slaveCycles   = slaveCycles;
//
// 		aggMatrixUpperFirstUse_[fineLevelIndex] = false;
// 	}
//
// 	aggMatrixUpperStructLevels_[fineLevelIndex].fPtr = ff.begin();
// 	aggMatrixUpperStructLevels_[fineLevelIndex].cPtr = cf.begin();
// }


void  Foam::swRestInterMap::agglomerateMatrixUpper
(
	scalarField& coarseUpper,
	scalarField& coarseDiag,
	const scalarField& fineUpper,
	const label fineLevelIndex
)
{
    // if(fineUpper.size() > minCellsUsingSW_)
    // {
    // 	initAggMatrixUpperStruct(coarseUpper, coarseDiag, fineUpper, fineLevelIndex);

    // 	agglomerateMatrixUpper_host(&aggMatrixUpperStructLevels_[fineLevelIndex]);
    // }
    // else
    {
    	// Get face restriction map for current level
	    const labelList& faceRestrictAddr =
	        aggl_.faceRestrictAddressing(fineLevelIndex);

	    /*if(faceRestrictAddr.size() < 1000000)
        {
            std::ofstream fout;
            label size = faceRestrictAddr.size();
            label cDSize = coarseDiag.size();
            label cUSize = coarseUpper.size();
            fout.open("faceRestrictAddressing.dat");
            fout << "fine size: " << size << std::endl;
            fout << "coarse diag size: " << cDSize << std::endl;
            fout << "coarse upper size: " << cUSize << std::endl;
            for(int i=0; i<size; i++)
            {
                fout << i << " " << faceRestrictAddr[i] << " " << std::endl;
            }
            fout.close();
            Info << "finish writing" << endl;
            int myi = 0;
            while (0 == myi)
            {
                Foam::sleep(5);
            }
        }*/

	    forAll(faceRestrictAddr, fineFacei)
	    {
	        label cFace = faceRestrictAddr[fineFacei];

	        if (cFace >= 0)
	        {
	            coarseUpper[cFace] += fineUpper[fineFacei];
	        }
	        else
	        {
	            // Add the fine face coefficient into the diagonal.
	            coarseDiag[-1 - cFace] += 2*fineUpper[fineFacei];
	        }
	    }
    }
}


