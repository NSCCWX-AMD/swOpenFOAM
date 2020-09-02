#ifndef ITERATOR_HPP
#define ITERATOR_HPP

// #include <iostream>
#include <map>
#include <cstring>
#include "swMacro.h"
#include "topology.hpp"
#include "iterator.h"

// ------------------------------------------------------------------------
//                             class iterator
// ------------------------------------------------------------------------
namespace UNAT
{

class Iterator
{
	private:
		Topology* _topo;
		map<swInt, swInt> _edgeMap;
		map<swInt, swInt> _vertexMap;
		swInt* _vertexWeights;
		swInt* _edgeWeights;

		// duplicate indicator
		bool duplicate_;

	public:

		// Constructors
		Iterator(){}
		Iterator(Topology &topo, swInt* vertexWeights, swInt* edgeWeights, bool duplicate = false)
		{
			if(duplicate)
			{
				this->_topo = topo.clone();
				swInt vertexNum = _topo->getVertexNumber();
				swInt edgeNum = _topo->getEdgeNumber();
				_vertexWeights = new swInt[vertexNum];
				memcpy(_vertexWeights, vertexWeights, sizeof(swInt)*vertexNum);
				_edgeWeights = new swInt[edgeNum];
				memcpy(_edgeWeights, edgeWeights, sizeof(swInt)*edgeNum);
				duplicate_ = true;
			}
			else
			{
				_topo = &topo;
				_vertexWeights = vertexWeights;
				_edgeWeights = edgeWeights;
				duplicate_ = false;
			}
		};

		Iterator(Topology &topo, bool duplicate = false)
		{
			if(duplicate)
			{
				this->_topo = topo.clone();
				duplicate_ = true;
			} else 
			{
				_topo = &topo;
				duplicate_ = false;
			}
			_vertexWeights = NULL;
			_edgeWeights = NULL;
		}

		// Deconstructors
		virtual ~Iterator()
		{
			if(duplicate_)
			{
				delete _topo;
				delete _vertexWeights;
				delete _edgeWeights;
			}
		}
		

		void reformInnerTopology()
		{
			reorderEdges(this->_topo->getStartVertices(),
						this->_topo->getEndVertices(),
						this->_topo->getEdgeNumber(),
						this->_topo->getVertexNumber());
			reorderNeighbor(this->_topo->getFirstEdgeVertices(),
						this->_topo->getVertexNeighbours(),
						this->_topo->getEdgeNumber(),
						this->_topo->getVertexNumber());
		}

		map<swInt, swInt>& getEdgeMap(){return this->_edgeMap;}
		map<swInt, swInt>& getVertexMap(){return this->_vertexMap;}
 
		virtual void reorderEdges(swInt* startVertices,
					swInt* endVertices,
					swInt edgeNumber, swInt vertexNumber){
			// To do
		}
		virtual void reorderNeighbor(swInt* firstEdgeVertices,
					swInt* vertexNeighbours,
					swInt edgeNumber, swInt vertexNumber) = 0;
		virtual void reorderEdgeData(Arrays* backEdgeData,
					Arrays* frontEdgeData){
			// To do
		}
		virtual void reorderNeighborData(Arrays* edgeData){
			// To do	
		}
		virtual void reorderVertexData(Arrays* edgeData){
			// To do
		}

		virtual void edge2VertexIteration(Arrays* paraData,
					coupledOperator *cOpt, swInt optNum) = 0;
		virtual void arrayIteration(Arrays* paraData, 
					coupledOperator* cOpt, int optNum) = 0;
		virtual void vertex2EdgeIteration(Arrays* paraData, 
					coupledOperator *cOpt, swInt optNum) = 0;
		
		Topology* getTopology(){return this->_topo;}
		swInt* getVertexWeights(){return this->_vertexWeights;}
		swInt* getEdgeWeights(){return this->_edgeWeights;}
};

} // namespace UNAT

#endif
