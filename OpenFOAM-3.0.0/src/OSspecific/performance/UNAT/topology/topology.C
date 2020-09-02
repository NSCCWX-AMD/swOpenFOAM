#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <assert.h>
#include "topology.hpp"

using namespace std;

namespace UNAT
{

// Dummy constructor
Topology::Topology()
{
	// cout<<"default constructor"<<endl;
	this->_vertexNumber=0;
	this->_edgeNumber=0;
	this->_startVertices=NULL;
	this->_endVertices=NULL;
	this->_startVertexNumbers=NULL;
	this->_accuStartVertexNumbers=NULL;
	this->_firstEdgeVertices=NULL;
	this->_vertexNeighbours=NULL;
	this->_vertexEdgeNumbers=NULL;
	this->_accuVertexEdgeNumbers=NULL;
}

// Factory Constructor
Topology* Topology::constructFromEdge(swInt *startVertices,
			swInt *endVertices, swInt edgeNumber, bool duplicate)
{
	Topology* topo = new Topology();
	topo->setStartVertices(startVertices);
	topo->setEndVertices(endVertices);
	topo->setEdgeNumber(edgeNumber);
//printf("edgeNumber is %d\n", topo->_edgeNumber);
//printArray("%d", topo->_startVertices, edgeNumber);
//printArray("%d", topo->_endVertices, edgeNumber);
	// cout<<"Constructor based on LDU"<<endl;
	topo->setDuplicate(0);
	topo->EdgeBasedInit();
	if( ! duplicate )
	{
		return topo;
	}
	else
	{
		Topology *topoCopy = new Topology();
		topoCopy->copy(*topo);
		delete topo;
		topoCopy->setDuplicate(1);
//printf("edgeNumber is %d\n", topoCopy->_edgeNumber);
//printArray("%d", topoCopy->_startVertices, edgeNumber);
//printArray("%d", topoCopy->_endVertices, edgeNumber);
		return topoCopy;
	}
}

Topology* Topology::constructFromVertex(swInt *accuVertexEdgeNumbers,
			swInt *vertexNeighbours, swInt vertexNumber, bool duplicate)
{
	Topology* topo = new Topology();
	topo->setAccuVertexEdgeNumbers(accuVertexEdgeNumbers);
	topo->setVertexNeighbours(vertexNeighbours);
	topo->setVertexNumber(vertexNumber);
	cout<<"Constructor based on CSR"<<endl;
	topo->VertexBasedInit();
	topo->setDuplicate(2);
	if( ! duplicate )
	{
		return topo;
	}
	else
	{
		Topology *topoCopy = new Topology();
		topoCopy->copy(*topo);
		delete topo;
		topoCopy->setDuplicate(3);
		return topoCopy;
	}
}

// Copy constructors
Topology::Topology(const Topology &topo)
{
	// cout<<"Copy constructor"<<endl;
	this->copy(topo);
}

Topology& Topology::operator=(const Topology& topo)
{
	cout<<"operator="<<endl;
	if(this == &topo) return *this;
	this->~Topology();
	this->copy(topo);
	return *this;
}

void Topology::copy(const Topology &topo)
{
	this->_edgeNumber = topo._edgeNumber;
	this->_vertexNumber = topo._vertexNumber;
	this->_startVertices=
		(swInt*)malloc(sizeof(swInt)*this->_edgeNumber);
	this->_endVertices=
		(swInt*)malloc(sizeof(swInt)*this->_edgeNumber);
	this->_firstEdgeVertices=
		(swInt*)malloc(sizeof(swInt)*this->_edgeNumber*2);
	this->_vertexNeighbours=
		(swInt*)malloc(sizeof(swInt)*this->_edgeNumber*2);
	for(int i=0;i<this->_edgeNumber;i++)
	{
		this->_startVertices[i] = topo._startVertices[i];
		this->_endVertices[i]   = topo._endVertices[i];
	}
	for(int i=0;i<this->_edgeNumber*2;i++)
	{
		this->_firstEdgeVertices[i] = topo._firstEdgeVertices[i];
		this->_vertexNeighbours[i]  = topo._vertexNeighbours[i];
	}

	this->_startVertexNumbers=
		(swInt*)malloc(sizeof(swInt)*this->_vertexNumber);
	this->_vertexEdgeNumbers=
		(swInt*)malloc(sizeof(swInt)*this->_vertexNumber);
	for(int i=0;i<this->_vertexNumber;i++)
	{
		this->_startVertexNumbers[i] = topo._startVertexNumbers[i];
		this->_vertexEdgeNumbers[i] = topo._vertexEdgeNumbers[i];
	}

	this->_accuStartVertexNumbers=
		(swInt*)malloc(sizeof(swInt)*this->_vertexNumber+1);
	this->_accuVertexEdgeNumbers=
		(swInt*)malloc(sizeof(swInt)*this->_vertexNumber+1);
	for(int i=0;i<this->_vertexNumber+1;i++)
	{
		this->_accuStartVertexNumbers[i] = topo._accuStartVertexNumbers[i];
		this->_accuVertexEdgeNumbers[i] = topo._accuVertexEdgeNumbers[i];
	}
}

// Deconstructors
Topology::~Topology()
{
	if(duplicate_ == 0)
	{
		free(_startVertexNumbers);
		free(_accuStartVertexNumbers);
		free(_firstEdgeVertices);
		free(_vertexNeighbours);
		free(_vertexEdgeNumbers);
		free(_accuVertexEdgeNumbers);
	}
	else if(duplicate_ == 2)
	{
		free(_startVertices);
		free(_endVertices);
		free(_startVertexNumbers);
		free(_accuStartVertexNumbers);
		free(_firstEdgeVertices);
		free(_vertexEdgeNumbers);
	}
	else
	{
		free(_startVertices);
		free(_endVertices);
		free(_startVertexNumbers);
		free(_accuStartVertexNumbers);
		free(_firstEdgeVertices);
		free(_vertexNeighbours);
		free(_vertexEdgeNumbers);
		free(_accuVertexEdgeNumbers);
	}
}

void Topology::sortAndCompress()
{
}

void Topology::EdgeBasedInit()
{
	for(int i=0;i<this->_edgeNumber;i++)
	{
//		cout<<i<<","<<this->_endVertices[i]<<endl;
		this->_vertexNumber=max(this->_vertexNumber,this->_endVertices[i]);
		this->_vertexNumber=max(this->_vertexNumber,this->_startVertices[i]);
	}
	this->_vertexNumber++;

	// record the face numbers of each row
	this->_startVertexNumbers=
		(swInt*)malloc(sizeof(swInt)*this->_vertexNumber);
	this->_vertexEdgeNumbers=
		(swInt*)malloc(sizeof(swInt)*this->_vertexNumber);
	for(int i=0;i<this->_vertexNumber;i++)
	{
		this->_startVertexNumbers[i]=0;
		this->_vertexEdgeNumbers[i]=0;
	}
	for(int i=0;i<this->_edgeNumber;i++)
	{
		this->_startVertexNumbers[this->_startVertices[i]]++;
		this->_vertexEdgeNumbers[this->_startVertices[i]]++;
		this->_vertexEdgeNumbers[this->_endVertices[i]]++;
	}

	// record the row offsets
	this->_accuStartVertexNumbers=
		(swInt*)malloc(sizeof(swInt)*(this->_vertexNumber+1));
	this->_accuVertexEdgeNumbers=
		(swInt*)malloc(sizeof(swInt)*(this->_vertexNumber+1));
	this->_accuStartVertexNumbers[0] = 0;
	this->_accuVertexEdgeNumbers[0] = 0;
	this->_accuVertexEdgeNumbers[1] = 0;
	for(int i=1;i<this->_vertexNumber;i++)
	{
		this->_accuVertexEdgeNumbers[i+1]
			= this->_accuVertexEdgeNumbers[i]
			+ this->_vertexEdgeNumbers[i-1];
	}
	for(int i=0;i<this->_vertexNumber;i++)
	{
		this->_accuStartVertexNumbers[i+1]
			= this->_accuStartVertexNumbers[i]
			+ this->_startVertexNumbers[i];
	}


	// edgeBase to vertexBase: lower triangle
	this->_vertexNeighbours =
		(swInt*)malloc(sizeof(swInt)*(this->_edgeNumber*2));
	this->_firstEdgeVertices =
		(swInt*)malloc(sizeof(swInt)*(this->_edgeNumber*2));
	for(int i=0;i<this->_edgeNumber;i++)
	{
		swInt row = this->_endVertices[i]+1;
		swInt col = this->_startVertices[i];
		this->_vertexNeighbours[this->_accuVertexEdgeNumbers[row]] = col;
		this->_firstEdgeVertices[this->_accuVertexEdgeNumbers[row]] = row-1;
		this->_accuVertexEdgeNumbers[row]++;
	}
	// vertexBase to edgeBase: upper triangle
	for(int i=0;i<this->_edgeNumber;i++)
	{
		swInt row = this->_startVertices[i]+1;
		swInt col = this->_endVertices[i];
		this->_vertexNeighbours[this->_accuVertexEdgeNumbers[row]] = col;
		this->_firstEdgeVertices[this->_accuVertexEdgeNumbers[row]] = row-1;
		this->_accuVertexEdgeNumbers[row]++;
	}
}

void Topology::VertexBasedInit()
{
	this->_edgeNumber = this->_accuVertexEdgeNumbers[this->_vertexNumber]/2;
	int edgeIndex = 0;

	// vertexBase to edgeBase
	this->_startVertices =
		(swInt*)malloc(sizeof(swInt)*(this->_edgeNumber));
	this->_endVertices =
		(swInt*)malloc(sizeof(swInt)*(this->_edgeNumber));
	this->_firstEdgeVertices =
		(swInt*)malloc(sizeof(swInt)*(this->_edgeNumber*2));
	for(int i=0;i<this->_vertexNumber;i++)
	{
		for(int j=this->_accuVertexEdgeNumbers[i];
					j<this->_accuVertexEdgeNumbers[i+1];j++)
		{
			this->_firstEdgeVertices[j] = i;
			if(this->_vertexNeighbours[j] > i)
			{
				this->_startVertices[edgeIndex] = i;
				this->_endVertices[edgeIndex] = this->_vertexNeighbours[j];
				edgeIndex++;
			}
		}
	}
	// record the face numbers of each row
	this->_startVertexNumbers=
		(swInt*)malloc(sizeof(swInt)*this->_vertexNumber);
	this->_vertexEdgeNumbers=
		(swInt*)malloc(sizeof(swInt)*this->_vertexNumber);
	for(int i=0;i<this->_vertexNumber;i++)
	{
		this->_startVertexNumbers[i]=0;
	}
	for(int i=0;i<this->_edgeNumber;i++)
	{
		this->_startVertexNumbers[this->_startVertices[i]]++;
		this->_vertexEdgeNumbers[this->_startVertices[i]]++;
		this->_vertexEdgeNumbers[this->_endVertices[i]]++;
	}

	this->_accuStartVertexNumbers=
		(swInt*)malloc(sizeof(swInt)*(this->_vertexNumber+1));
	this->_accuStartVertexNumbers[0] = 0;
	for(int i=0;i<this->_vertexNumber;i++)
	{
		this->_accuStartVertexNumbers[i+1]
			= this->_accuStartVertexNumbers[i]
			+ this->_startVertexNumbers[i];
	}
}

swInt Topology::getVertexNumber()
{
	return this->_vertexNumber;
}

swInt Topology::getEdgeNumber()
{
	return this->_edgeNumber;
}

swInt* Topology::getStartVertices()
{
	return this->_startVertices;
}

swInt* Topology::getEndVertices()
{
	return this->_endVertices;
}

swInt* Topology::getStartVertexNumbers()
{
	return this->_startVertexNumbers;
}

swInt* Topology::getAccuStartVertexNumbers()
{
	return this->_accuStartVertexNumbers;
}

swInt* Topology::getFirstEdgeVertices()
{
	return this->_firstEdgeVertices;
}

swInt* Topology::getVertexNeighbours()
{
	return this->_vertexNeighbours;
}

swInt* Topology::getVertexEdgeNumbers()
{
	return this->_vertexEdgeNumbers;
}

swInt* Topology::getAccuVertexEdgeNumbers()
{
	return this->_accuVertexEdgeNumbers;
}


} // namespace UNAT


