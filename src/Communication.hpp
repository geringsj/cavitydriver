
#ifndef Communication_hpp
#define Communication_hpp

#include "GridFunction.hpp"
#include "Domain.hpp"

/** 
 * Communication via MPI.
 * 
 * @author becherml, friesfn, geringsj
 * @date 11/2014
 */
class Communication
{
private:
	int m_numProcs;

	Dimension m_procsGrid_dim;
	Index m_procsGrid_myPosition;

	/* rank 0 process is manager (and also worker), others are workers */
	int m_myRank;
	int m_downRank;
	int m_upRank;
	int m_leftRank;
	int m_rightRank;

	Dimension m_globalDomain_dim;
	Dimension m_myDomain_dim;
	Color m_myDomainFirstCellColor;
	Index m_myOffsetToGlobalDomain;

	uint m_bufferSize;
	Real* m_sendBuffer;
	Real* m_recvBuffer;

	/* use these to send/receive the filled send/receive buffer */
	void m_sendToOne(int count, int dest, int tag);
	int m_recvFromOne(int source, int tag);

	/* maybe we won't need those later, because we can tell MPI directly to sum over all received residdums and stuff
	 * so maybe use an MPI call on check*SORCycle / getGlobalTimeStep */
	void m_sendToAll(void *buffer, int count);
	void m_recvFromAll(/*void* sendbuf, int sendcount,*/void* recvbuf, int recvcount);

	void ExchangeTwoGridFunctions(GridFunction& one, GridFunction& two, Domain domain);

	void ExchangeOneGridFunction(GridFunction& one, Domain domain);

public:

	bool getValid() const { return m_myRank>=0; }
	Dimension getLocalDimensions() const { return m_myDomain_dim; }

	enum class Handle {
		Pressure,
		Velocities,
		PreliminaryVelocities
	};

	Communication(Dimension globalDomainDim);
	~Communication();

	void exchangeGridBoundaryValues(Domain domain, Handle grid, Color handleColorCells=Color::All);

	void exchangeGridInnerValues(Domain domain, Handle grid);

	bool checkForAnotherSORCycle(Real mySubResiduum);

	Real getGlobalTimeStep(Delta myMaxValues);
};


#endif
