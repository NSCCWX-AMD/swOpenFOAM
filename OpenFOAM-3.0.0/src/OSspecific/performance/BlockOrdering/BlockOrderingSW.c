/*#########################
	将数据结构进行到底

	MLB means Multilevel_Blocking
#########################*/
#include "BlockOrdering.h"
#define LOG (printf("%s (%d) - <%s>\n",__FILE__,__LINE__,__FUNCTION__),printf)

void MLB_Multilevel_ordering
(
 MLB_graph graph,
 LABEL levels,		// number of levels
 LABEL* blockNums,	// number of blocks for each level
 // blockStarts: Four-LABEL element array.
 // The first is block-row index, and second block-column.
 // The third is the start(include) edge index for each block.
 // The fourth is the end(exclude) edge index for each block.
 // Only upper triagle stored.
 LABEL* blockStarts, 
 LABEL* cellStarts, // start cell index for each block
 // postCellOrder & postEdgeOrder:
 // Index mirror from the old order to the new one.
 // The position of the element corresponds to the old order,
 // and the element value the new one.
 // The negative value in postEdgeOrder means the face is 
 // swapped, since the new owner is bigger tha new neighbor.
 LABEL* postCellOrder,
 LABEL* postEdgeOrder
)
{
	LABEL* owner		= graph.owner;
 	LABEL* neighbor		= graph.neighbor;
 	LABEL* cellWeights	= graph.cellWeights;
 	LABEL* edgeWeights	= graph.edgeWeights;
 	LABEL cellNum		= graph.cellNum;
 	LABEL edgeNum		= graph.edgeNum;
	// construct the cell index array
	LABEL* cellID =
		(LABEL*)malloc(sizeof(LABEL)*cellNum);
	MLB_generateCellID
	( owner,
	  neighbor,
	  cellID,
	  cellNum,
	  edgeNum
	);
	// initiallize multilevel Block infomation
	// initiallize the original cell order
	LABEL icell;
	for(icell=0; icell<cellNum; icell++)
		postCellOrder[icell] = icell;
	// get the total blocks number	
	LABEL blockSpan = 1;
	LABEL leveli;
	for(leveli=0; leveli<levels; leveli++)
		blockSpan *= blockNums[leveli];
	// construct a graph array that holds all the subgraphs
	MLB_graph block_graph[blockSpan/blockNums[levels-1]];

	// initialize the first graph
	block_graph[0].owner = 
		(LABEL*)malloc(sizeof(LABEL)*edgeNum);
	block_graph[0].neighbor = 
		(LABEL*)malloc(sizeof(LABEL)*edgeNum);
	block_graph[0].cellWeights =
		(LABEL*)malloc(sizeof(LABEL)*cellNum);
	block_graph[0].edgeWeights = 
		(LABEL*)malloc(sizeof(LABEL)*edgeNum);
	// offset edge to avoid discontinuous cell index
	MLB_offsetEdges
	(
	 owner,
	 neighbor,
	 cellID,
	 block_graph[0].owner,
	 block_graph[0].neighbor,
	 cellNum,
	 edgeNum
	);
	memcpy
	(
	 block_graph[0].cellWeights, 
	 cellWeights, 
	 sizeof(LABEL)*cellNum
	);
	memcpy
	(
	 block_graph[0].edgeWeights, 
	 edgeWeights, 
	 sizeof(LABEL)*edgeNum
	);
	block_graph[0].cellNum = cellNum;
	block_graph[0].edgeNum = edgeNum;
	// set the initial number of blocks at
	// set the first cell start
	cellStarts[0] = 0;
	// present level, "1" for first level
	LABEL blockNum  = 1;

	// multilevel decomposing
	for(leveli=0; leveli<levels; leveli++)
	{
		LABEL blocki;
		for(blocki=0; blocki<blockNum; blocki++)
		{
			// get initial block index
			LABEL blockIndex = blocki*blockSpan;
			// get the father graph
			MLB_graph father_graph = block_graph[blockIndex/blockNums[levels-1]];
			// initial ordering result arrays
			LABEL* block_postCellOrder =
				(LABEL*)malloc(sizeof(LABEL)*father_graph.cellNum);
			LABEL icell;
			for(icell=0; icell<father_graph.cellNum; icell++)
				block_postCellOrder[icell] = icell;
			LABEL* block_cellStarts = 
				(LABEL*)malloc(sizeof(LABEL)*(blockNums[leveli]+1));
			// decompose present block
			MLB_ordering
			(
			 father_graph,
			 blockNums[leveli],
			 block_postCellOrder,
			 block_cellStarts
			);
			// update global arrays
			// update the global cell order
			LABEL* postCellOrder_ptr = &postCellOrder[cellStarts[blockIndex]];
			MLB_postLABEL
			(
			 block_postCellOrder,
			 postCellOrder_ptr,
			 father_graph.cellNum,
			 1
			);
			// update the global cell starts
			LABEL jblock;
			for(jblock=0; jblock<blockNums[leveli]+1; jblock++)
				cellStarts[ blockIndex+jblock*blockSpan/blockNums[leveli] ]
					= cellStarts[blockIndex]+block_cellStarts[jblock];
			// reorder the father graph cell weight
			MLB_postLABEL
			(
			 block_postCellOrder,
			 father_graph.cellWeights,
			 father_graph.cellNum,
			 1
			);
			// construct the next level graphs
			if(leveli != levels-1)
			{
				for(jblock=blockNums[leveli]-1; jblock>=0; jblock--)
				{
					MLB_graph* presentGraph = 
						&block_graph
						[ 
						 (blockIndex+jblock*blockSpan/blockNums[leveli])
						 /blockNums[levels-1] 
						];

					ExtensibleLABELArray presentOwner;
					ExtensibleLABELArray presentNeighbor;
					ExtensibleLABELArray presentEdgeWeight;
					extensibleLABELArrayInit(&presentOwner);
					extensibleLABELArrayInit(&presentNeighbor);
					extensibleLABELArrayInit(&presentEdgeWeight);
					LABEL iface;
					for(iface=0; iface<father_graph.edgeNum; iface++ )
					{
						LABEL new_ownerOrder = 
							block_postCellOrder[father_graph.owner[iface]];
						LABEL new_neighOrder = 
							block_postCellOrder[father_graph.neighbor[iface]];
						if(
							new_ownerOrder < block_cellStarts[jblock+1]
							&& new_neighOrder < block_cellStarts[jblock+1]
							&& new_ownerOrder >= block_cellStarts[jblock]
							&& new_neighOrder >= block_cellStarts[jblock]
						  )
						{
							extensibleLABELArrayAdd(
								&presentOwner, father_graph.owner[iface] );
							extensibleLABELArrayAdd(
								&presentNeighbor, father_graph.neighbor[iface] );
							extensibleLABELArrayAdd
							( 
							 &presentEdgeWeight, 
							 father_graph.edgeWeights[iface] 
							);
						}	
					}

					LABEL presentCellNum = block_cellStarts[jblock+1] 
										   - block_cellStarts[jblock];
					LABEL* presentCellWeight =
						(LABEL*)malloc(sizeof(LABEL)*presentCellNum);

					memcpy( presentCellWeight, 
							father_graph.cellWeights + block_cellStarts[jblock], 
							sizeof(LABEL)*presentCellNum );
					
					if(jblock==0)
					{
						free(father_graph.owner);
						free(father_graph.neighbor);
						free(father_graph.cellWeights);
						free(father_graph.edgeWeights);
					}

					presentGraph->cellNum = presentCellNum;
					presentGraph->cellWeights =	
						(LABEL*)malloc(sizeof(LABEL)*presentCellNum);
					memcpy( presentGraph->cellWeights, 
							presentCellWeight, 
							sizeof(LABEL)*presentCellNum );
					free(presentCellWeight);

					presentGraph->edgeNum = presentOwner.size;
					presentGraph->owner =
						(LABEL*)malloc(sizeof(LABEL)*presentOwner.size);
					presentGraph->neighbor =
						(LABEL*)malloc(sizeof(LABEL)*presentOwner.size);
					presentGraph->edgeWeights =
						(LABEL*)malloc(sizeof(LABEL)*presentOwner.size);
					LABEL* present_cellID =
						(LABEL*)malloc(sizeof(LABEL)*presentCellNum);
					MLB_generateCellID
					(
					 presentOwner.data,
					 presentNeighbor.data,
					 present_cellID,
					 presentCellNum,
					 presentOwner.size
					);
					MLB_offsetEdges
					(
					 presentOwner.data,
					 presentNeighbor.data,
					 present_cellID,
					 presentGraph->owner,
					 presentGraph->neighbor,
					 presentCellNum,
					 presentOwner.size
					);
					free(present_cellID);
					//memcpy
					//( 
					// presentGraph->owner, 
					// presentOwner.data, 
					// sizeof(LABEL)*presentOwner.size
					//);
					//memcpy
					//( 
					// presentGraph->neighbor, 
					// presentNeighbor.data, 
					// sizeof(LABEL)*presentOwner.size
					//);
					memcpy
					( 
					 presentGraph->edgeWeights, 
					 presentEdgeWeight.data, 
					 sizeof(LABEL)*presentOwner.size
					);
					extensibleLABELArrayDestroy(&presentOwner);
					extensibleLABELArrayDestroy(&presentNeighbor);
					extensibleLABELArrayDestroy(&presentEdgeWeight);
				}
			}
			else
			{
				free(father_graph.owner);
				free(father_graph.neighbor);
				free(father_graph.cellWeights);
				free(father_graph.edgeWeights);
			}
			free(block_postCellOrder);
			free(block_cellStarts);
  	}
  	// update the block number and span for the next level 
  	blockNum*=blockNums[leveli];
  	blockSpan/=blockNums[leveli];
  }

  // generate the global edge order array and generate blockStarts
  // swap the global cell order array
	LABEL* postCellOrder_tmp = (LABEL*)malloc(sizeof(LABEL)*cellNum);
	memcpy(postCellOrder_tmp, postCellOrder, sizeof(LABEL)*cellNum);
	for(icell=0; icell<cellNum; icell++)
		postCellOrder[postCellOrder_tmp[icell]] = icell;
	free(postCellOrder_tmp);
	// recount the block num
	blockSpan = 1;
	for(leveli=0; leveli<levels; leveli++)                                
		blockSpan *= blockNums[leveli];
	// construct new local edge
	LABEL* newLocalOwner = 
		(LABEL*)malloc(sizeof(LABEL)*edgeNum);	
	LABEL* newLocalNeigh = 
		(LABEL*)malloc(sizeof(LABEL)*edgeNum);
	LABEL iedge;
	for(iedge=0; iedge<edgeNum; iedge++)
	{
		newLocalOwner[iedge] = postCellOrder
							   [
					    	    MLB_find( owner[iedge], 
							  		      cellID, 
							  		      cellNum )
							   ];
		newLocalNeigh[iedge] = postCellOrder
							   [
					    	    MLB_find( neighbor[iedge], 
							  	 	      cellID, 
							  		      cellNum )
							   ];
	}
	free(cellID);
	// pick up the edges
	LABEL edgeNumPicked = 0;
	LABEL brow, bcol, blockCount;
	blockCount = 0;
    //==================flat version=================
    LABEL matrixBlockNum = (blockSpan+1)*blockSpan/2;
    ExtensibleLABELArray* edgeInBlock =
        (ExtensibleLABELArray*)malloc(
                sizeof(ExtensibleLABELArray)*matrixBlockNum);
    for(blockCount=0; blockCount<matrixBlockNum; blockCount++)
        extensibleLABELArrayInit(&edgeInBlock[blockCount]);
    for(iedge=0; iedge<edgeNum; iedge++)
    {
        LABEL NLOwner = newLocalOwner[iedge];
        LABEL NLNeigh = newLocalNeigh[iedge];

        brow = 0, bcol = 0;
        LABEL uBound = blockSpan;
        LABEL lBound = 0;
        brow = (uBound+lBound)/2;
		while(1)
        {
            if(NLOwner < cellStarts[brow])
                uBound=brow-1;
            else if(NLOwner >= cellStarts[brow+1])
                lBound=brow+1;
            else
                break;
            brow = (uBound+lBound)/2;
            if( (NLOwner < cellStarts[brow]
                || NLOwner >= cellStarts[brow+1])
                && uBound == lBound )
            {
                CERR("edge not linked any local Cell!\n");
                exit(-1);
            }
        }
        uBound = blockSpan;
        lBound = 0;
        bcol = (uBound+lBound)/2;
		while(1)
        {
            if(NLNeigh < cellStarts[bcol])
                uBound=bcol-1;
            else if(NLNeigh >= cellStarts[bcol+1])
                lBound=bcol+1;
            else
                break;
            bcol = (uBound+lBound)/2;
            if( (NLNeigh < cellStarts[bcol]
                || NLNeigh >= cellStarts[bcol+1])
                && uBound == lBound )
            {
                CERR("edge not linked any local Cell!\n");
                exit(-1);
            }
        }
        if(brow<=bcol)
        {
            blockCount = (2*blockSpan-brow+1)*brow/2+(bcol-brow);
            if(blockCount<0||blockCount>matrixBlockNum)
                CERR("block location out of range\n");
            extensibleLABELArrayAdd(
                    &edgeInBlock[blockCount],
                    iedge);
        }
        if(brow>bcol)
        {

            blockCount = (2*blockSpan-bcol+1)*bcol/2+(brow-bcol);
            if(blockCount<0||blockCount>matrixBlockNum)
                CERR("block location out of range\n");
            extensibleLABELArrayAdd(
                    &edgeInBlock[blockCount],
                    -1*iedge);
        }
    }

	for(brow=0; brow<blockSpan; brow++)
	for(bcol=brow; bcol<blockSpan; bcol++)
	{
        blockCount = (2*blockSpan-brow+1)*brow/2+(bcol-brow);
        if(blockCount<0||blockCount>matrixBlockNum)
            CERR("block location out of range\n");
		
        blockStarts[blockCount*4] = brow;
		blockStarts[blockCount*4+1] = bcol;
		blockStarts[blockCount*4+2] = edgeNumPicked;

        memcpy(
            &postEdgeOrder[edgeNumPicked],
            edgeInBlock[blockCount].data,
            edgeInBlock[blockCount].size*sizeof(LABEL)
            );
        edgeNumPicked += edgeInBlock[blockCount].size;
		
        blockStarts[blockCount*4+3] = edgeNumPicked;

        extensibleLABELArrayDestroy(
                &edgeInBlock[blockCount]);
    }
    free(edgeInBlock);

    //============multilevel loop version============
	//for(brow=0; brow<blockSpan; brow++)
	//{
	//	for(bcol=brow; bcol<blockSpan; bcol++)
	//	{
	//		blockStarts[blockCount*4] = brow;
	//		blockStarts[blockCount*4+1] = bcol;
	//		blockStarts[blockCount*4+2] = edgeNumPicked;
	//		for(iedge=0; iedge<edgeNum; iedge++)
	//		{
	//			LABEL NLOwner = newLocalOwner[iedge];
	//			LABEL NLNeigh = newLocalNeigh[iedge];

	//			if( NLOwner >= cellStarts[brow] 
	//				&& NLOwner < cellStarts[brow+1] 
	//				&& NLNeigh >= cellStarts[bcol] 
	//				&& NLNeigh < cellStarts[bcol+1] )
	//			{
	//				postEdgeOrder[edgeNumPicked] = iedge;
	//				edgeNumPicked++;
	//			}
	//			else 
	//			if(
	//				NLNeigh >= cellStarts[brow] 
	//				&& NLNeigh < cellStarts[brow+1] 
	//				&& NLOwner >= cellStarts[bcol] 
	//				&& NLOwner < cellStarts[bcol+1] )
	//			{
	//				postEdgeOrder[edgeNumPicked] = -1*iedge;
	//				edgeNumPicked++;
	//			}
	//		}
    //        blockStarts[blockCount*4+3] = edgeNumPicked;
    //        blockCount++;
	//	}
	//}

	if(edgeNumPicked != edgeNum )
	{
		CERR("more or less edge !!!\n");
		exit(-1);
	}
	free(newLocalOwner);
	free(newLocalNeigh);
	// swap the global edge order array
	LABEL* postEdgeOrder_tmp = (LABEL*)malloc(sizeof(LABEL)*edgeNum);
	memcpy(postEdgeOrder_tmp, postEdgeOrder, sizeof(LABEL)*edgeNum);
	for(iedge=0; iedge<edgeNum; iedge++)
	{
		if(postEdgeOrder_tmp[iedge]>=0)
			postEdgeOrder[postEdgeOrder_tmp[iedge]] = iedge;
		else
			postEdgeOrder[-1*postEdgeOrder_tmp[iedge]] = -1*iedge;
	}
	free(postEdgeOrder_tmp);
}

void MLB_ordering
(
 MLB_graph graph,
 LABEL blockNum,	// number of blocks for each level
 // postCellOrder:
 // Index mirror from the old order to the new one.
 // The position of the element corresponds to the old order,
 // and the element value the new one.
 LABEL* postCellOrder,	// starts from zero
 LABEL* cellStarts	// cells in each block
)
{
	LABEL* owner		= graph.owner;
 	LABEL* neighbor		= graph.neighbor;
 	LABEL* cellWeights	= graph.cellWeights;
 	LABEL* edgeWeights	= graph.edgeWeights;
 	LABEL cellNum		= graph.cellNum;
 	LABEL edgeNum		= graph.edgeNum;

	idx_t* xadj = 
		(idx_t*)malloc(sizeof(idx_t)*(cellNum+1));
	idx_t* adjncy = 
		(idx_t*)malloc(sizeof(idx_t)*2*edgeNum);
	idx_t* edge_wgt = 
		(idx_t*)malloc(sizeof(idx_t)*2*edgeNum);
	// construct metis arrays
	MLB_constructMetisCSR
	(
	 owner,
	 neighbor,
	 edgeWeights,
	 xadj,
	 adjncy,
	 edge_wgt,
	 cellNum,
	 edgeNum
	);
    //=============metis options setup================
	//set the options
	idx_t options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	// explicitly minimize subdomain degree
	options[METIS_OPTION_MINCONN] = 1;
	// set the method: "k-way" prior to k_way
	char method[] = "k_way"; 
    //=================================================
	idx_t* blockCells = 
		(idx_t*)malloc(sizeof(idx_t)*cellNum);
	LABEL edgeCut = 0;	
	// decompose graph with metis
	MLB_metis_decompose
	(
	 cellNum,
	 edgeNum*2,
	 xadj,
	 adjncy,
	 cellWeights,
	 edge_wgt,
	 blockNum,
	 options,
	 &edgeCut,
	 blockCells,
	 method
	);
	// update the cell order
	LABEL blockj, cellj;
	LABEL newCellPosi = 0;
    //==============flat version=====================
	//ExtensibleLABELArray* cellsInBlock =
	//	(ExtensibleLABELArray*)malloc( 
    //            sizeof(ExtensibleLABELArray)*blockNum
    //            );
	//for(blockj=0; blockj<blockNum; blockj++)
	//  extensibleLABELArrayInit(&cellsInBlock[blockj]);

	//for(cellj=0; cellj<cellNum; cellj++)
	//{
	//	extensibleLABELArrayAdd(
	//		&cellsInBlock[blockCells[cellj]],
	//		cellj 
    //        );
	//}

	//for(blockj=0; blockj<blockNum; blockj++)
	//{
	//	cellStarts[blockj] = newCellPosi;
	//	memcpy( 
	//		&postCellOrder[newCellPosi],
	//		cellsInBlock[blockj].data,
	//		sizeof(LABEL)*cellsInBlock[blockj].size
	//		);
	//	newCellPosi+=cellsInBlock[blockj].size;
	//	extensibleLABELArrayDestroy(
	//		&cellsInBlock[blockj] );
	//}
	//cellStarts[blockNum] = newCellPosi;
	//if(newCellPosi != cellNum)
	//	CERR("Unexpected number of decomposed cells!!\n");
	//free(cellsInBlock);
    //// swap postCellOrder
	//LABEL* postCellOrder_tmp =
	//	(LABEL*)malloc(sizeof(LABEL)*cellNum);
	//memcpy(postCellOrder_tmp, 
	//			postCellOrder, sizeof(LABEL)*cellNum);
	//for(cellj=0; cellj<cellNum; cellj++)
	//  postCellOrder[postCellOrder_tmp[cellj]]
	//	  = cellj;
    //free(postCellOrder_tmp);

	//============two level loop version=================
	for(blockj=0; blockj<blockNum; blockj++)
	{
		cellStarts[blockj] = newCellPosi;
		for(cellj=0; cellj<cellNum; cellj++)
		{
			if(blockCells[cellj]==blockj)
			{
			  postCellOrder[cellj] = newCellPosi;
			  newCellPosi++;
			}
		}
	}
	cellStarts[blockNum] = newCellPosi;
	
	if(newCellPosi != cellNum)
		CERR("Unexpected number of decomposed cells!!\n");

	// release memory
	free(xadj);
	free(adjncy);
	free(edge_wgt);
	free(blockCells);
}


