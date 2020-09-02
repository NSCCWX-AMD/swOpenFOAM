#ifndef UNAT_ROWSUBSECTIONITERATOR_HPP
#define UNAT_ROWSUBSECTIONITERATOR_HPP

// #include <iostream>
#include "swMacro.h"
#include "iterator.hpp"
#include "iterator.h"
#include "rowSubsections.hpp"
//#include "topology.H"
//using namespace std;
namespace UNAT
{

// ------------------------------------------------------------------------
//                     class RowSubsectionIterator
// ------------------------------------------------------------------------

class RowSubsectionIterator : public Iterator
{
	private:
		SWFoam::RowSubsections *secs_;
		swInt maxFaces_;
		swInt dataSize_;
		swInt *maxFacesInSeg_;
		swInt *maxCellsInSeg_;
		swInt *maxColsInSeg_;
		swInt maxColNum_;

	public:
		// Constructors
		RowSubsectionIterator(coupledOperator* opt, Topology &topo, swInt* cellWeights, 
					swInt* edgeWeights, int pattern, bool duplicate = true);
		// Deconstructors
		~RowSubsectionIterator(){};

		virtual void reorderEdges(swInt* startVertices, swInt* endVertices,
					swInt edgeNumber, swInt vertexNumber){};
		virtual void reorderNeighbor(swInt* firstEdgeVertices,
					swInt* vertexNeighbours,
					swInt edgeNumber, swInt vertexNumber){};
		virtual void edge2VertexIteration(Arrays* backEdgeData,
					Arrays* frontEdgeData, Arrays* selfConnData,
					Arrays* vertexData,
					e2v_hostFunPtr fun_host, e2v_slaveFunPtr fun_slave){};
		virtual void edge2VertexIteration(Arrays* paraData,
					coupledOperator* cOpt, int optNum);
		virtual void arrayIteration(Arrays* paraData,
                    coupledOperator* cOpt, int optNum);
		virtual void vertex2EdgeIteration(Arrays* paraData,
					coupledOperator *cOpt, swInt optNum);

		swInt getMaxFaces(){return this->maxFaces_;};
		swInt getDataSize(){return this->dataSize_;};
};

}

#endif
