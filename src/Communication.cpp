#include "Communication.hpp"
#include "Debug.hpp"

//#define WITHMPI

#ifdef WITHMPI
	#include <mpi.h>
#endif
#include <cmath>

Communication::Communication(Dimension globalDomainDim)
{
#ifdef WITHMPI
	/* init MPI */
	int rc = MPI_Init(0, NULL);
	if (rc != MPI_SUCCESS) 
	{
		log_err("Error starting MPI program. Terminating.");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	MPI_Comm_size(MPI_COMM_WORLD, &m_numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &m_myRank);
#else
	m_numProcs = 1;
	m_myRank = 0;
#endif
	//if(m_myRank==0)
	//	log_info("number of tasks=%i - my rank=%i", m_numProcs, m_myRank);

	/* set mpi cartesian grid with dimension infos */
	m_globalDomain_dim = globalDomainDim;
	m_globalInnerRange = Range(Index(1,1,1), m_globalDomain_dim);

	int x_procs = (int)( floor(sqrt((Real)(m_numProcs))) );
	int y_procs = m_numProcs/x_procs;
	/* y_procs >= x_procs, so switch dimensions if there are more cells in x direction */
	if(globalDomainDim.i > globalDomainDim.j)
	{
		int tmp = x_procs; x_procs = y_procs; y_procs = tmp;
	}
	m_procsGrid_dim.i = x_procs;
	m_procsGrid_dim.j = y_procs;
	//if(m_myRank == 0)
	//	log_info("going for a %ix%i processors grid", x_procs, y_procs);
	/* column-major order of processes in grid 
	 * (only because we use this in GridFunction too)*/
	m_procsGrid_myPosition.i = m_myRank / y_procs;
	m_procsGrid_myPosition.j = m_myRank % y_procs;
	/* disable processes which don't fit into grid */
	if(m_myRank >= x_procs*y_procs)
	{
		log_info("[P%i] no place in processes grid. terminating. (at position [%i,%i])", 
				m_myRank, m_procsGrid_myPosition.i, m_procsGrid_myPosition.j);
		//m_myRank = -1;
		/* return; ? */
	}

	/* make a new MPI communication domain for processes we actually will use, 
	 * let the rest finalize out. */
#ifdef WITHMPI
	mycom = new MPI_Comm;
	*(MPI_Comm*)mycom = MPI_COMM_WORLD;
	if(m_numProcs > x_procs*y_procs)
	{
		/* the following is taken from: 
		 * http://stackoverflow.com/questions/13774968/mpi-kill-unwanted-processes */
		MPI_Group world_group;
		MPI_Comm_group(MPI_COMM_WORLD, &world_group);
		MPI_Group new_group;
		int ranges[3] = { x_procs*y_procs, m_numProcs-1, 1 };
		MPI_Group_range_excl(world_group, 1, &ranges, &new_group);
		/* create new communicator for non-dead processes */
		MPI_Comm_create(MPI_COMM_WORLD, new_group, (MPI_Comm*)mycom);
		if (*(MPI_Comm*)mycom == MPI_COMM_NULL)
		{
		   //log_info("Bye bye cruel world. (%i)", m_myRank);
		   MPI_Finalize();
		   exit(0);
		}
	}
#else
	mycom = NULL;
#endif

	/* compute direct neighbours */
	m_downRank = (m_procsGrid_myPosition.j-1 >= 0) 
		? (m_myRank-1) : (-1);
	m_upRank = (m_procsGrid_myPosition.j+1 < y_procs) 
		? (m_myRank+1) : (-1);
	m_leftRank = (m_procsGrid_myPosition.i-1 >= 0) 
		? (m_myRank-y_procs) : (-1);
	m_rightRank = (m_procsGrid_myPosition.i+1 < x_procs) 
		? (m_myRank+y_procs) : (-1);
	//log_info("rank=%i, right=%i, down=%i, left=%i, up=%i", 
	//		m_myRank, m_rightRank, m_downRank, m_leftRank, m_upRank);

	/* compute global/local dimensions, checkerboard colors and stuff */
	int x_pcells = m_globalDomain_dim.i / m_procsGrid_dim.i;
	int y_pcells = m_globalDomain_dim.j / m_procsGrid_dim.j;

	Index myOffsetToGlobalDomain = Dimension(
		m_procsGrid_myPosition.i * x_pcells ,
		m_procsGrid_myPosition.j * y_pcells);

	m_myDomain_dim = Dimension(
			x_pcells + ((m_procsGrid_myPosition.i < m_procsGrid_dim.i-1) ?
				(0) : (m_globalDomain_dim.i % m_procsGrid_dim.i)) ,
			y_pcells + ((m_procsGrid_myPosition.j < m_procsGrid_dim.j-1) ?
				(0) : (m_globalDomain_dim.j % m_procsGrid_dim.j)) );

	m_localInnerRange = Range(
		Index(
			myOffsetToGlobalDomain.i + 1,
			myOffsetToGlobalDomain.j + 1,
			myOffsetToGlobalDomain.k + 1
			), 
		Index(
			myOffsetToGlobalDomain.i + m_myDomain_dim.i,
			myOffsetToGlobalDomain.j + m_myDomain_dim.j,
			myOffsetToGlobalDomain.k + m_myDomain_dim.k
			) );

	/* global (1,1)  cell is Red */
	m_myDomainFirstCellColor = 
		(/* have fun reversing this */
		 !(((myOffsetToGlobalDomain.i % 2)+(myOffsetToGlobalDomain.j % 2))%2)
		 ) 
			? (Color::Red) : (Color::Black);

	/* set send/recv buffers. the longest side of domain will be buffers size. */
	m_bufferSize = (m_myDomain_dim.i > m_myDomain_dim.j) 
		? (m_myDomain_dim.i) : (m_myDomain_dim.j);
	m_bufferSize += 3;
#ifdef WITHMPI
	m_sendBuffer = new Real[m_bufferSize];
	m_recvBuffer = new Real[m_bufferSize];
#else
	m_sendBuffer = NULL;
	m_recvBuffer = NULL;
#endif

	//log_info("rank=%i has offset=(%i,%i), mydims=(%i,%i), firstColor=%s, bufferSize=%i",
	//		m_myRank, myOffsetToGlobalDomain.i, myOffsetToGlobalDomain.j,
	//		m_myDomain_dim.i, m_myDomain_dim.j, 
	//		(m_myDomainFirstCellColor==Color::Red)?("Red"):("Black"), m_bufferSize);
}

Communication::~Communication()
{
#ifdef WITHMPI
	delete[] m_recvBuffer;
	delete[] m_sendBuffer;
	delete (MPI_Comm*)mycom;
	MPI_Finalize();
#endif
}

void Communication::sendBufferTo(int count, int dest, int tag)
{
#ifdef WITHMPI
	MPI_Send(m_sendBuffer, count, 
			MPI_DOUBLE, dest, tag, *(MPI_Comm*)mycom);
#else
	count = count*dest*tag; // this is to get rid of warnings
#endif
}

void Communication::recvBufferFrom(int source, int tag)
{
#ifdef WITHMPI
	MPI_Recv(m_recvBuffer, m_bufferSize, 
			MPI_DOUBLE, source, tag, *(MPI_Comm*)mycom, MPI_STATUS_IGNORE);
#else
	source = source*tag; // this is to get rid of warnings
#endif
}

Real Communication::getGlobalResidual(Real mySubResidual)
{
#ifdef WITHMPI
	MPI_Allreduce(&mySubResidual, m_recvBuffer, 1, MPI_DOUBLE, 
			MPI_SUM, *(MPI_Comm*)mycom);
	return sqrt(m_recvBuffer[0]);
#else
	return sqrt(mySubResidual);
#endif
}

int Communication::getGlobalFluidCellsCount(int local_FluidCells)
{
#ifdef WITHMPI
	int recv;
	MPI_Allreduce(&local_FluidCells, &recv, 1, MPI_INT, 
			MPI_SUM, *(MPI_Comm*)mycom);
	return recv;
#else
	return local_FluidCells;
#endif
}

bool Communication::checkGlobalFinishSOR(bool myLoopCondition)
{
#ifdef WITHMPI
	MPI_Allreduce(MPI_IN_PLACE, &myLoopCondition, 1, MPI_C_BOOL, 
			MPI_LOR, *(MPI_Comm*)mycom);
	return myLoopCondition;
#else
	return myLoopCondition;
#endif
}

Delta Communication::getGlobalMaxVelocities(Delta myMaxValues)
{
#ifdef WITHMPI
	m_sendBuffer[0] = myMaxValues.x;
	m_sendBuffer[1] = myMaxValues.y;
	m_sendBuffer[2] = 0.0;//myMaxValues.z; TODO 3D

	MPI_Allreduce(m_sendBuffer, m_recvBuffer, 3, MPI_DOUBLE, 
			MPI_MAX, *(MPI_Comm*)mycom);

	return Delta(m_recvBuffer[0], m_recvBuffer[1], m_recvBuffer[2]);
#else
	return myMaxValues;
#endif
}

void Communication::exchangeGridBoundaryValues(
		Domain& domain, Handle grid) 
{
#ifdef WITHMPI
	switch (grid)
	{
	case Communication::Handle::Pressure:
		exchangeGridBoundaryValues(domain.p(), domain.getWholeInnerRange());
		break;
	case Communication::Handle::Velocities:
		exchangeGridBoundaryValues(domain.u(), domain.getWholeInnerRange());
		exchangeGridBoundaryValues(domain.v(), domain.getWholeInnerRange());
		break;
	case Communication::Handle::PreliminaryVelocities:
		exchangeGridBoundaryValues(domain.F(), domain.getWholeInnerRange());
		exchangeGridBoundaryValues(domain.G(), domain.getWholeInnerRange());
		break;
	default:
		break;
	}
#else
	grid = (Handle)(domain.getDimension().i * (int)grid * 0);
	/* just to get rid of warnings. hopefully nothing bad happpens. */
#endif
}

void Communication::exchangeGridBoundaryValues(
		GridFunction& gf,
		Range inner)
{
	/* 
	 *handle all neighbour boundaries 
	 */

	/* a) send left - receive right */
	if(m_leftRank >= 0)
	{
		// copy data to buffer
		for (int j = inner.begin[1], n=0; j <= inner.end[1]; j++, n++)
			m_sendBuffer[n] = (gf(inner.begin[0], j));

		// send buffer to m_leftRank
		sendBufferTo(inner.end[1]-inner.begin[1]+1, m_leftRank, m_myRank);
	}
	if(m_rightRank >= 0)
	{
		// receive buffer from m_rightRank
		recvBufferFrom(m_rightRank, m_rightRank);

		// copy data from buffer to grid
		for (int j = inner.begin[1], n=0; j <= inner.end[1]; j++, n++)
			 gf(inner.end[0]+1, j) = m_recvBuffer[n];
	}

	/* b) send right - receive left */
	if(m_rightRank >= 0)
	{
		// copy data to buffer
		for (int j = inner.begin[1], n=0; j <= inner.end[1]; j++, n++)
			m_sendBuffer[n] = (gf(inner.end[0], j));

		// send buffer to m_rightRank
		sendBufferTo(inner.end[1]-inner.begin[1]+1, m_rightRank, m_myRank);
	}
	if(m_leftRank >= 0)
	{
		// receive buffer from m_leftRank
		recvBufferFrom(m_leftRank, m_leftRank);

		// copy data from buffer to grid
		for (int j = inner.begin[1], n=0; j <= inner.end[1]; j++, n++)
			 gf(inner.begin[0]-1, j) = m_recvBuffer[n];
	}

	/* c) send up - receive down */
	if(m_upRank >= 0)
	{
		// copy data to buffer
		for (int i = inner.begin[0], n=0; i <= inner.end[0]; i++, n++)
			m_sendBuffer[n] = (gf(i, inner.end[1]));

		// send buffer to m_upRank
		sendBufferTo(inner.end[0]-inner.begin[0]+1, m_upRank, m_myRank);
	}
	if(m_downRank >= 0)
	{
		// receive buffer from m_downRank
		recvBufferFrom(m_downRank, m_downRank);

		// copy data from buffer to grid
		for (int i = inner.begin[0], n=0; i <= inner.end[0]; i++, n++)
			gf(i, inner.begin[1]-1) = m_recvBuffer[n];
	}

	/* d) send down - receive up */
	if(m_downRank >= 0)
	{
		// copy data to buffer
		for (int i = inner.begin[0], n=0; i <= inner.end[0]; i++, n++)
			m_sendBuffer[n] = (gf(i, inner.begin[1]));

		// send buffer to m_downRank
		sendBufferTo(inner.end[0]-inner.begin[0]+1, m_downRank, m_myRank);
	}
	if(m_upRank >= 0)
	{
		// receive buffer from m_upRank
		recvBufferFrom(m_upRank, m_upRank);

		// copy data from buffer to grid
		for (int i = inner.begin[0], n=0; i <= inner.end[0]; i++, n++)
			gf(i, inner.end[1]+1) = m_recvBuffer[n];
	}

	/* attention: domain handles real boundaries somewhere else, not in this function */
}

Range Communication::getProcLocalInnerRange(int procrank) const 
{ 
	int x_pcells = m_globalDomain_dim.i / m_procsGrid_dim.i;
	int y_pcells = m_globalDomain_dim.j / m_procsGrid_dim.j;

	Index localProcGridPosition
		(procrank / m_procsGrid_dim.j, 
		 procrank % m_procsGrid_dim.j);

	Index procOffsetToGlobalDomain = Dimension(
		localProcGridPosition.i * x_pcells ,
		localProcGridPosition.j * y_pcells);

	Dimension localDomain_dim = Dimension(
			x_pcells + ((localProcGridPosition.i < m_procsGrid_dim.i-1) ?
				(0) : (m_globalDomain_dim.i % m_procsGrid_dim.i)) ,
			y_pcells + ((localProcGridPosition.j < m_procsGrid_dim.j-1) ?
				(0) : (m_globalDomain_dim.j % m_procsGrid_dim.j)) );

	Range ret(
		Index(
			procOffsetToGlobalDomain.i + 1,
			procOffsetToGlobalDomain.j + 1,
			procOffsetToGlobalDomain.k + 1
			), 
		Index(
			procOffsetToGlobalDomain.i + localDomain_dim.i,
			procOffsetToGlobalDomain.j + localDomain_dim.j,
			procOffsetToGlobalDomain.k + localDomain_dim.k
			) );
	return ret; 
}

