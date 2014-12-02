
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
	void sendBufferTo(int count, int dest, int tag=0);
	void recvBufferFrom(int source, int tag=0);

	void exchangeGridBoundaryValues(
			GridFunction& gf,
			Index ibegin, 
			Index iend);

	void ExchangeVelocityBoundaryValues(GridFunction& u, GridFunction& v,
		Index u_begin, Index u_end, Index v_begin, Index v_end);

	void ExchangePressureBoundaryValues(GridFunction& p,
		Index dim_one_begin, Index dim_one_end);

public:

	enum class Handle {
		Pressure,
		Velocities,
		PreliminaryVelocities
	};

	Communication(Dimension globalDomainDim);
	~Communication();

	void exchangeGridBoundaryValues(Domain domain, Handle grid);

	//void exchangeGridInnerValues(Domain domain, Handle grid);

	Real getGlobalResidual(Real mySubResidual);

	Delta getGlobalMaxVelocities(Delta myMaxValues);

	bool getProcessValid() const { return m_myRank>=0; }
	Dimension getLocalDimensions() const { return m_myDomain_dim; }
	Color getFirstCellColor() const { return m_myDomainFirstCellColor; }
};


#endif
