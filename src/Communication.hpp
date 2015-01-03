
#ifndef Communication_hpp
#define Communication_hpp

#include "GridFunction.hpp"
#include "Domain.hpp"
#include "Boundary.hpp"

/** 
 * Communication via MPI.
 * 
 * @author becherml, friesfn, geringsj
 * @date 11/2014
 */
class Communication
{
private:
	void* mycom;

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

	Range m_globalInnerRange;
	Range m_localInnerRange;

	uint m_bufferSize;
	Real* m_sendBuffer;
	Real* m_recvBuffer;

	/* use these to send/receive the filled send/receive buffer */
	void sendBufferTo(int count, int dest, int tag=0);
	void recvBufferFrom(int source, int tag=0);

	void exchangeGridBoundaryValues(
			GridFunction& gf,
			Range inner);

public:

	enum class Handle {
		Pressure,
		Velocities,
		PreliminaryVelocities
	};

	Communication(Dimension globalDomainDim);
	~Communication();

	void exchangeGridBoundaryValues(Domain& domain, Handle grid);
	Real getGlobalResidual(Real mySubResidual);
	Delta getGlobalMaxVelocities(Delta myMaxValues);
	bool checkGlobalFinishSOR(bool myLoopCondition);
	Color getFirstCellColor() const { return m_myDomainFirstCellColor; }
	//void exchangeGridInnerValues(Domain domain, Handle grid);
	Boundary::Competence getBoundaryCompetence() const 
	{ return Boundary::Competence(m_upRank<0, m_downRank<0, m_leftRank<0, m_rightRank<0); }

	int getRank() const { return m_myRank; }
	int getProcsCount() const { return m_numProcs; }

	Dimension getProcsGridDim() const { return m_procsGrid_dim; }
	Index getProcsGridPosition() const { return m_procsGrid_myPosition; }

	Dimension getGlobalDimension() const { return m_globalDomain_dim; }
	Range getGlobalInnerRange() const { return m_globalInnerRange; }

	Dimension getLocalDimension() const { return m_myDomain_dim; }
	Range getLocalInnerRange() const { return m_localInnerRange; }
	Range getProcLocalInnerRange(int procrank) const;
};


#endif
