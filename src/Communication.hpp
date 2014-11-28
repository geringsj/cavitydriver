
#ifndef Communication_hpp
#define Communication_hpp

#include "GridFunction.hpp"
#include "Domain.hpp"


/** 
 * Communication via MPI 
 * 
 * @author becherml, friesfn, geringsj
 * @date 11/2014
 */
class Communication
{
	/* MPI will be asked to construct a cartesian grid of processes */
private:
	Dimension m_procsGrid_dim;
	Dimension m_procsGrid_myPosition;
	Dimension m_globalDomain_dim;
	Dimension m_myDomain_dim;
	Color m_myDomainFirstCellColor;
	Dimension m_myOffsetToGlobalDomain;

	int m_numProcs;

	/* rank 0 process is manager (and also worker), others are workers */
	int m_myRank;
	int m_downRank;
	int m_upRank;
	int m_leftRank;
	int m_rightRank;

	void m_sendToOne(/* ... */);
	void m_recvFromOne(/* ... */);

	void m_sendToAll(/* ... */);
	void m_recvFromAll(/* ... */);

	bool SqrtIsEven(int number);

public:

	enum class Handle {
		Pressure,
		Velocities,
		PreliminaryVelocities
	};

	Communication(Dimension globalDomainDim, /* MPI_Init needs argc and argv */int argc, char** argv);
	~Communication();

	void exchangeGridBoundaryValues(Domain domain, Handle grid, Color handleColorCells=Color::All);

	void exchangeGridInnerValues(Domain domain, Handle grid);

	bool checkForAnotherSORCycle(Real mySubResiduum);

	Real getGlobalTimeStep(Delta myMaxValues);
};


#endif
