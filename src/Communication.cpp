#include "Communication.hpp"
#include "Debug.hpp"

#include <mpi.h>
#include <cmath>


Communication::Communication(Dimension globalDomainDim)
{
	/* init MPI */
	int rc = MPI_Init(0, NULL);
	if (rc != MPI_SUCCESS) 
	{
		log_err("Error starting MPI program. Terminating.");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	MPI_Comm_size(MPI_COMM_WORLD, &m_numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &m_myRank);
	log_info("number of tasks=%i - my rank=%i", m_numProcs, m_myRank);

	/* set mpi cartesian grid with dimension infos */
	m_globalDomain_dim = globalDomainDim;

	int x_procs = (int)( floor(sqrt((Real)(m_numProcs))) );
	int y_procs = m_numProcs/x_procs;
	/* y_procs >= x_procs, so switch dimensions if there are more cells in x direction */
	if(globalDomainDim.i > globalDomainDim.j)
	{
		int tmp = x_procs; x_procs = y_procs; y_procs = tmp;
	}
	m_procsGrid_dim.i = x_procs;
	m_procsGrid_dim.j = y_procs;
	if(m_myRank == 0)
		debug("going for a %ix%i processors grid", x_procs, y_procs);
	/* column-major order of processes in grid 
	 * (only because we use this in GridFunction too)*/
	m_procsGrid_myPosition.i = m_myRank / y_procs;
	m_procsGrid_myPosition.j = m_myRank % y_procs;
	/* disable processes which don't fit into grid */
	if(m_myRank >= x_procs*y_procs)
	{
		log_info("process %i has no place in grid (at position [%i,%i])", 
				m_myRank, m_procsGrid_myPosition.i, m_procsGrid_myPosition.j);
		m_myRank = -1;
		/* return; ? */
	}

	/* compute direct neighbours */
	m_downRank = (m_procsGrid_myPosition.j-1 >= 0) 
		? (m_myRank-1) : (-1);
	m_upRank = (m_procsGrid_myPosition.j+1 < y_procs) 
		? (m_myRank+1) : (-1);
	m_leftRank = (m_procsGrid_myPosition.i-1 >= 0) 
		? (m_myRank-y_procs) : (-1);
	m_rightRank = (m_procsGrid_myPosition.i+1 < x_procs) 
		? (m_myRank+y_procs) : (-1);
	log_info("rank=%i, right=%i, down=%i, left=%i, up=%i", 
			m_myRank, m_rightRank, m_downRank, m_leftRank, m_upRank);

	/* compute global/local dimensions, checkerboard colors and stuff */
	int x_pcells = m_globalDomain_dim.i / x_procs;
	int y_pcells = m_globalDomain_dim.j / y_procs;

	m_myOffsetToGlobalDomain = Dimension(
			m_procsGrid_myPosition.i * x_pcells ,
			m_procsGrid_myPosition.j * y_pcells);

	m_myDomain_dim = Dimension(
			x_pcells + ((m_procsGrid_myPosition.i < x_procs-1) ?
				(0) : (m_globalDomain_dim.i % x_procs)) ,
			y_pcells + ((m_procsGrid_myPosition.j < y_procs-1) ?
				(0) : (m_globalDomain_dim.j % y_procs)) );

	/* global (1,1)  cell is Red */
	m_myDomainFirstCellColor = 
		(/* have fun reversing this */
		 !(((m_myOffsetToGlobalDomain.i % 2)+(m_myOffsetToGlobalDomain.j % 2))%2)
		 ) 
			? (Color::Red) : (Color::Black);

	/* set send/recv buffers. the longest side of domain will be buffers size. */
	m_bufferSize = (m_myDomain_dim.i > m_myDomain_dim.j) 
		? (m_myDomain_dim.i) : (m_myDomain_dim.j);
	m_bufferSize = (m_bufferSize > 3) ? (m_bufferSize) : (3); /* TODO also numProcs ?? */
	m_sendBuffer = new Real[m_bufferSize];
	m_recvBuffer = new Real[m_bufferSize];

	log_info("rank=%i has offset=(%i,%i), mydims=(%i,%i), firstColor=%s, bufferSize=%i",
			m_myRank, m_myOffsetToGlobalDomain.i, m_myOffsetToGlobalDomain.j,
			m_myDomain_dim.i, m_myDomain_dim.j, 
			(m_myDomainFirstCellColor==Color::Red)?("Red"):("Black"), m_bufferSize);

}

Communication::~Communication()
{
	delete[] m_recvBuffer;
	delete[] m_sendBuffer;
	MPI_Finalize();
}

void Communication::sendBufferTo(int count, int dest, int tag)
{
	MPI_Send(m_sendBuffer, count, 
			/*TODO switch?*/MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
}

void Communication::recvBufferFrom(int source, int tag)
{
	MPI_Recv(m_recvBuffer, m_bufferSize, 
			/*TODO switch?*/MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

Real Communication::getGlobalResidual(Real mySubResidual)
{
	MPI_Allreduce(&mySubResidual, m_recvBuffer, 1, /*TODO switch?*/MPI_DOUBLE, 
			MPI_SUM, MPI_COMM_WORLD);
	return sqrt(m_recvBuffer[0]);
}

Delta Communication::getGlobalMaxVelocities(Delta myMaxValues)
{
	m_sendBuffer[0] = myMaxValues.x;
	m_sendBuffer[1] = myMaxValues.y;
	m_sendBuffer[2] = myMaxValues.z;

	MPI_Allreduce(m_sendBuffer, m_recvBuffer, 3, /*TODO switch?*/MPI_DOUBLE, 
			MPI_MAX, MPI_COMM_WORLD);

	myMaxValues.x = m_recvBuffer[0];
	myMaxValues.y = m_recvBuffer[1];
	myMaxValues.z = m_recvBuffer[2];
	return myMaxValues;
}

void Communication::exchangeGridBoundaryValues(
		Domain domain, Handle grid) //, Color handleColorCells)
{
	switch (grid)
	{
	case Communication::Handle::Pressure:
		ExchangePressureBoundaryValues(domain.p(),
			domain.getBeginInnerDomains(), domain.getEndInnerDomainP());
		break;
	case Communication::Handle::Velocities:
		ExchangeVelocityBoundaryValues(domain.u(), domain.v(),
			domain.getBeginInnerDomains(), domain.getEndInnerDomainU(), 
			domain.getBeginInnerDomains(), domain.getEndInnerDomainV());
		break;
	case Communication::Handle::PreliminaryVelocities:
		ExchangeVelocityBoundaryValues(domain.F(), domain.G(),
			domain.getBeginInnerDomains(), domain.getEndInnerDomainU(),
			domain.getBeginInnerDomains(), domain.getEndInnerDomainV());
		break;
	default:
		break;
	}
}

/*
void Communication::exchangeGridInnerValues(Domain domain, Handle grid)
{
	switch (grid)
	{
	case Communication::Handle::Pressure:
		break;
	case Communication::Handle::Velocities:
		break;
	case Communication::Handle::PreliminaryVelocities:
		break;
	default:
		break;
	}
}
*/

void Communication::exchangeGridBoundaryValues(
		GridFunction& gf,
		Index ibegin, 
		Index iend)
{
	/* 
	 *handle all neighbour boundaries 
	 */

	/* a) send left - receive right */
	if(m_leftRank >= 0)
	{
		// copy data to buffer
		for (int j = ibegin[1], n=0; j <= iend[1]; j++, n++)
			m_sendBuffer[n] = (gf(ibegin[0], j));

		// send buffer to m_leftRank
		sendBufferTo(iend[1]-ibegin[1]+1, m_leftRank, m_myRank);
	}
	if(m_rightRank >= 0)
	{
		// receive buffer from m_rightRank
		recvBufferFrom(m_rightRank, m_myRank);

		// copy data from buffer to grid
		for (int j = ibegin[1], n=0; j <= iend[1]; j++, n++)
			 gf(iend[0]+1, j) = m_recvBuffer[n];
	}

	/* b) send right - receive left */
	if(m_rightRank >= 0)
	{
		// copy data to buffer
		for (int j = ibegin[1], n=0; j <= iend[1]; j++, n++)
			m_sendBuffer[n] = (gf(iend[0], j));

		// send buffer to m_rightRank
		sendBufferTo(iend[1]-ibegin[1]+1, m_rightRank, m_myRank);
	}
	if(m_leftRank >= 0)
	{
		// receive buffer from m_leftRank
		recvBufferFrom(m_leftRank, m_myRank);

		// copy data from buffer to grid
		for (int j = ibegin[1], n=0; j <= iend[1]; j++, n++)
			 gf(ibegin[0]-1, j) = m_recvBuffer[n];
	}

	/* c) send up - receive down */
	if(m_upRank >= 0)
	{
		// copy data to buffer
		for (int i = ibegin[0], n=0; i <= iend[0]; i++, n++)
			m_sendBuffer[n] = (gf(i, iend[1]));

		// send buffer to m_upRank
		sendBufferTo(iend[0]-ibegin[0]+1, m_upRank, m_myRank);
	}
	if(m_downRank >= 0)
	{
		// receive buffer from m_downRank
		recvBufferFrom(m_downRank, m_myRank);

		// copy data from buffer to grid
		for (int i = ibegin[0], n=0; i <= iend[0]; i++, n++)
			gf(i, ibegin[1]-1) = m_recvBuffer[n];
	}

	/* d) send down - receive up */
	if(m_downRank >= 0)
	{
		// copy data to buffer
		for (int i = ibegin[0], n=0; i <= iend[0]; i++, n++)
			m_sendBuffer[n] = (gf(i, ibegin[1]));

		// send buffer to m_downRank
		sendBufferTo(iend[0]-ibegin[0]+1, m_downRank, m_myRank);
	}
	if(m_upRank >= 0)
	{
		// receive buffer from m_upRank
		recvBufferFrom(m_upRank, m_myRank);

		// copy data from buffer to grid
		for (int i = ibegin[0], n=0; i <= iend[0]; i++, n++)
			gf(i, iend[1]+1) = m_recvBuffer[n];
	}

	/* attention: Domain should handle real boundaries somewhere else, 
	 * not in this function 
	 * maybe give back to Domain class ? */
}

void Communication::ExchangeVelocityBoundaryValues(GridFunction& u, GridFunction& v,
		Index u_begin, Index u_end, Index v_begin, Index v_end)
{
	// a)
	if (m_leftRank >= 0)
	{
		// SERGEJ CHECK/FIX THIS
		// fill buffer with u values from left border
		for(int i = u_begin[1]; i <= u_end[1]; i++) 
			//DOUBLE-TODO!!!! check if this size is correct
			m_sendBuffer[i /* <- this is dangerous TODO */ ] = (u(u_begin[0], i));

		// send buffer to m_leftRank
		sendBufferTo(u_end[1]-u_begin[1]+1,m_leftRank,m_myRank); 
		// TODO check size

		// TODO receive buffer from m_rightRank
		recvBufferFrom(m_rightRank, m_myRank);
		// TODO do something with the values in the buffer

		// fill buffer with v values from left border
		for (int i = v_begin[1]; i <= v_end[1]; i++)
			m_sendBuffer[i] = (v(v_begin[0], i));

		// send buffer to m_leftRank
		sendBufferTo(v_end[1]-v_begin[1]+1, m_leftRank, m_myRank); 
		// TODO check size

		// TODO receive buffer from m_rightRank
		recvBufferFrom(m_rightRank, m_myRank);
		// TODO do something with the values in the buffer
	}
	else
	{
		// TODO ask domain to set left boundary values
	}

	// b)
	if (m_rightRank != -1)
	{
		// fill buffer with u values from left border
		for (int i = u_begin[1]; i <= u_end[1]; i++) 
			//DOUBLE-TODO!!!! check if this size is correct
			m_sendBuffer[i /* <- this is STILL dangerous , TODO */] = (u(u_end[0], i));

		// send buffer to m_rightRank
		sendBufferTo(u_end[1]-u_begin[1]+1, m_rightRank, m_myRank); 
		// TODO check size

		// TODO receive buffer from m_leftRank
		recvBufferFrom(m_leftRank, m_myRank);
		// TODO do something with the values in the buffer

		// fill buffer with v values from right border
		for (int i = v_begin[1]; i <= v_end[1]; i++)
			m_sendBuffer[i] = (v(v_end[0], i));

		// send buffer to m_rightRank
		sendBufferTo(v_end[1]-v_begin[1]+1, m_rightRank, m_myRank); 
		// TODO check size

		// TODO receive buffer from m_leftRank
		recvBufferFrom(m_leftRank, m_myRank);
		// TODO do something with the values in the buffer
	}
	else
	{
		// TODO ask domain to set right boundary values
	}

	// c)
	if (m_upRank != -1)
	{
		// fill buffer with values from upper border
		for (int i = u_begin[0]; i <= u_end[0]; i++) //DOUBLE-TODO!!!! check if this size is correct
			m_sendBuffer[i] = (u(i, u_end[1]));

		// send buffer to m_upRank
		sendBufferTo(u_end[0]-u_begin[0]+1, m_upRank, m_myRank); // TODO check size

		// TODO receive buffer from m_downRank
		recvBufferFrom(m_downRank, m_myRank);
		// do something with the values in the buffer

		// fill buffer with values from upper border
		for (int i = v_begin[0]; i <= v_end[0]; i++) //DOUBLE-TODO!!!! check if this size is correct
			m_sendBuffer[i] = (v(i, v_end[1]));

		// send buffer to m_upRank
		sendBufferTo(v_end[0]-v_begin[0]+1, m_upRank, m_myRank); // TODO check size

		// TODO receive buffer from m_downRank
		recvBufferFrom(m_downRank, m_myRank);
		// do something with the values in the buffer
	}
	else
	{
		// TODO ask domain to set upper boundary values
	}

	// d)
	if (m_downRank != -1)
	{
		// fill buffer with values from lower border
		for (int i = u_begin[0]; i <= u_end[0]; i++) //DOUBLE-TODO!!!! check if this size is correct
			m_sendBuffer[i] = (u(i,u_begin[1]));

		// send buffer to m_downRank
		sendBufferTo(u_end[0]-u_begin[0]+1, m_downRank, m_myRank); // TODO check size

		// TODO receive buffer from m_upRank
		recvBufferFrom(m_upRank, m_myRank);
		// do something with the values in the buffer

		// fill buffer with values from upper border
		for (int i = v_begin[0]; i <= v_end[0]; i++) //DOUBLE-TODO!!!! check if this size is correct
			m_sendBuffer[i] = (v(i,v_begin[1]));

		// send buffer to m_downRank
		sendBufferTo(v_end[0]-v_begin[0]+1, m_downRank, m_myRank); // TODO check size

		// TODO receive buffer from m_upRank
		recvBufferFrom(m_upRank, m_myRank);
		// do something with the values in the buffer
	}
	else
	{
		// TODO ask domain to set lower boundary values
	}
}

void Communication::ExchangePressureBoundaryValues(GridFunction& p,
		Index p_begin, Index p_end)
{
	// a)
	if (m_leftRank != -1)
	{
		// SERGEJ CHECK/FIX THIS
		// fill buffer with p values from left border
		for (int i = p_begin[1]; i <= p_end[1]; i++) //DOUBLE-TODO!!!! check if this size is correct
			m_sendBuffer[i] = (p(p_begin[0], i));

		// send buffer to m_leftRank
		sendBufferTo(p_end[1]-p_begin[1]+1, m_leftRank, m_myRank); // TODO check size

		// TODO receive buffer from m_rightRank
		recvBufferFrom(m_rightRank, m_myRank);
		// do something with the values in the buffer
	}
	else
	{
		// TODO ask domain to set left boundary values
	}

	// b)
	if (m_rightRank != -1)
	{
		// fill buffer with u values from left border
		for (int i = p_begin[1]; i <= p_end[1]; i++) //DOUBLE-TODO!!!! check if this size is correct
			m_sendBuffer[i] = (p(p_end[0], i));

		// send buffer to m_rightRank
		sendBufferTo(p_end[1]-p_begin[1]+1, m_rightRank, m_myRank); // TODO check size

		// TODO receive buffer from m_leftRank
		recvBufferFrom(m_leftRank, m_myRank);
		// do something with the values in the buffer
	}
	else
	{
		// TODO ask domain to set right boundary values
	}

	// c)
	if (m_upRank != -1)
	{
		// fill buffer with values from upper border
		for (int i = p_begin[0]; i <= p_end[0]; i++) //DOUBLE-TODO!!!! check if this size is correct
			m_sendBuffer[i] = (p(i, p_end[1]));

		// send buffer to m_upRank
		sendBufferTo(p_end[0]-p_begin[0]+1, m_upRank, m_myRank); // TODO check size

		// TODO receive buffer from m_downRank
		recvBufferFrom(m_downRank, m_myRank);
		// do something with the values in the buffer
	}
	else
	{
		// TODO ask domain to set upper boundary values
	}

	// d)
	if (m_downRank != -1)
	{
		// fill buffer with values from lower border
		for (int i = p_begin[0]; i <= p_end[0]; i++) //DOUBLE-TODO!!!! check if this size is correct
			m_sendBuffer[i] = (p(i, p_begin[1]));

		// send buffer to m_downRank
		sendBufferTo(p_end[0]-p_begin[0]+1, m_downRank, m_myRank); // TODO check size

		// receive buffer from m_upRank
		recvBufferFrom(m_upRank, m_myRank);
		// do something with the values in the buffer
	}
	else
	{
		// TODO ask domain to set lower boundary values
	}
}

