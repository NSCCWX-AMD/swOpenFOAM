#ifndef UNAT_MULTILEVELBLOCKITERATOR_HPP
#define UNAT_MULTILEVELBLOCKITERATOR_HPP

// #include <iostream>
#include "swMacro.h"
#include "iterator.hpp"
#include "iterator.h"
#include "RlmpiInitializer.hxx"
//#include "topology.H"
//using namespace std;
namespace UNAT
{

// ------------------------------------------------------------------------
//                     class MultiLevelBlockIterator
// ------------------------------------------------------------------------

class MultiLevelBlockIterator : public Iterator
{
	private:
		swInt  _cpeBlockNum;
		swInt  _mshBlockNum;
		swInt  _mtxBlockNum;
		swInt* _blockStarts;
		swInt* _blockStartsUnsymm;
		swInt* _vertexStarts;
		swInt  _maxXNum;
		swInt  _maxCells;
		swInt  _maxEdges;
		swInt  _maxEdgesUnsymm;
		swInt* _owner;
		swInt* _neighbor;
		swInt* _postEdgeOrder;
		swInt* _postVertexOrder;
		// 此处保存一份未重排的COO形式拓扑拷贝
		// 重排过的的拓扑保存在Topology结构体中
		swInt* _firstEdgeVertices;
		swInt* _vertexNeighbours;

		swInt* _recvStarts;

		// Register communication topology
		RlmpiInfo* _schedule_data;

		void reorderVertexArray(swFloat* array);
		void reorderEdgeArrayUnsymm(swFloat* array);
		void MLBReorder(Topology &topo, swInt ref, swInt lastBlockNum, swInt iter);

		void initOwnNeiSendList();
		void writeTopology();
	public:
		// Constructors
		MultiLevelBlockIterator(Topology &topo, swInt *vertexWeights,
					swInt* edgeWeights, bool duplicate = true);
		// Deconstructors
		~MultiLevelBlockIterator(){};

		virtual void reorderEdges(swInt* startVertices, swInt* endVertices,
					swInt edgeNumber, swInt vertexNumber);
		virtual void reorderNeighbor(swInt* firstEdgeVertices,
					swInt* vertexNeighbours,
					swInt edgeNumber, swInt vertexNumber);
		virtual void edge2VertexIteration(Arrays* backEdgeData,
					Arrays* frontEdgeData, Arrays* selfConnData,
					Arrays* vertexData,
					e2v_hostFunPtr fun_host, e2v_slaveFunPtr fun_slave){};
		void        edge2VertexIteration(Arrays* paraData,
					coupledOperator* cOpt, int optNum);
		virtual void arrayIteration(Arrays* paraData,
                    coupledArrayOperator* cOpt, int optNum);
		virtual void vertex2EdgeIteration(Arrays* paraData,
					coupledOperator *cOpt, swInt optNum);
		virtual void reorderEdgeData(Arrays* backEdgeData,
					Arrays* frontEdgeData);
		void reorderEdgeDataUnsymm(Arrays* edgeData);
		virtual void reorderVertexData(Arrays* vertexData);
		void restoreVertexData(Arrays* vertexData);
		void restoreEdgeData(Arrays* backEdgeData, Arrays* frontEdgeData);

		swInt getCpeBlockNum(){return this->_cpeBlockNum;};
		swInt getMshBlockNum(){return this->_mshBlockNum;};
		swInt getMtxBlockNum(){return this->_mtxBlockNum;};
		swInt getMaxXNum(){return this->_maxXNum;};
		swInt getMaxCells(){return this->_maxCells;};
		swInt getMaxEdges(){return this->_maxEdges;};
		swInt getMaxEdgesUnsymm(){return this->_maxEdgesUnsymm;};
		swInt* getBlockStarts(){return this->_blockStarts;};
		swInt* getBlockStartsUnsymm(){return this->_blockStartsUnsymm;};
		swInt* getVertexStarts(){return this->_vertexStarts;};
		swInt* getPostEdgeOrder(){return this->_postEdgeOrder;};
		swInt* getPostVertexOrder(){return this->_postVertexOrder;};
		swInt* getFirstEdgeVertices() {return this->_firstEdgeVertices;};
		swInt* getVertexNeighbours(){return this->_vertexNeighbours;};
};

}

#endif
