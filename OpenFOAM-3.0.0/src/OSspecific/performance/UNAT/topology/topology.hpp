#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

using namespace std;
#include "swMacro.h"

// ------------------------------------------------------------------------
//                         class Topology
// ------------------------------------------------------------------------
namespace UNAT
{

class Topology{
	private:
		swInt  _vertexNumber;
		swInt  _edgeNumber;
		swInt* _startVertices;
		swInt* _endVertices;
		swInt* _startVertexNumbers;
		swInt* _accuStartVertexNumbers;
		swInt* _firstEdgeVertices;
		swInt* _vertexNeighbours;
		swInt* _vertexEdgeNumbers;
		swInt* _accuVertexEdgeNumbers;
		void EdgeBasedInit();
		void VertexBasedInit();
		void copy(const Topology &topo);

		// duplicate indicator
		// 0 for constructFromEdge don't duplicate
		// 1 for constructFromEdge duplicate
		// 2 for constructFromEdge don't duplicate
		// 3 for constructFromEdge duplicate
		swInt duplicate_;

		//friend Topology* constructFromEdge(swInt *startVertices,
		//			swInt *endVertices, swInt edgeNumber, bool duplicate);
		//friend Topology* constructFromVertex(swInt *accuVertexEdgeNumbers,
		//			swInt *vertexNeighbours,swInt vertexNumber, bool duplicate);
	public:
		// Constructors
		Topology();
		// Factory constructors
		static Topology* constructFromEdge(swInt *startVertices, 
					swInt *endVertices, swInt edgeNumber, bool duplicate = true);
		static Topology* constructFromVertex(swInt *accuVertexEdgeNumbers,
					swInt *vertexNeighbours,swInt vertexNumber, bool duplicate = true);
		// Copy constructors
		Topology(const Topology &topo);
		Topology& operator=(const Topology& topo);
		Topology* clone() const {return new Topology(*this);};

		// Deconstructors
		~Topology();

		// Data accessers
		swInt  getVertexNumber();
		swInt  getEdgeNumber();
		swInt* getStartVertices();
		swInt* getEndVertices();
		swInt* getStartVertexNumbers();
		swInt* getAccuStartVertexNumbers();
		swInt* getFirstEdgeVertices();
		swInt* getVertexNeighbours();
		swInt* getVertexEdgeNumbers();
		swInt* getAccuVertexEdgeNumbers();

		void setStartVertices(swInt* startVertices)
		{ _startVertices = startVertices; }
		void setEndVertices(swInt* endVertices)
		{ _endVertices = endVertices; }
		void setEdgeNumber(swInt edgeNumber)
		{ _edgeNumber = edgeNumber; }
		void setAccuVertexEdgeNumbers(swInt* accuVertexEdgeNumbers)
		{ _accuVertexEdgeNumbers = accuVertexEdgeNumbers; }
		void setVertexNeighbours(swInt* vertexNeighbours)
		{ _vertexNeighbours = vertexNeighbours; }
		void setVertexNumber(swInt vertexNumber)
		{ _vertexNumber = vertexNumber; }
		void setDuplicate( swInt duplicate)
		{ duplicate_ = duplicate; }

		// TODO: 对于未按顺序排列的LDU方式存储拓扑，提供排序接口，
		// 返回类型为LDUMatrix，CSR同理
		// 提供LDU拓扑与CSR拓扑之间的相互转换接口
		void addEdge();
		void addVertex();
		void removeEdge();
		void removeVertex();
		void transpose();

		void sortAndCompress();
		void edgeBasedToVertexBased();
		void vertexBasedToEdgeBased();
};


} // namespace UNAT

//#include "topology.C"

#endif
