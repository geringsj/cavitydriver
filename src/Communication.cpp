#include "Communication.hpp"
#include "Debug.hpp"


#include <mpi.h>
#include <vector>
#include <cmath>



/**
 * There was no other way to find out the thread causing the error
 * so I had to use this "madness" to solve the problem.
 * 
 * We can define our own error handling function for MPI errors.
 * The problem ist this is not allowed to be a member function of
 * Communication. This leads us to the fact that we cann only use
 * static variables to "communicate" between this function and
 * our Communication object.
 * 
 * There is an bool array of a certain length and if a thread
 * throws an error the respective value in the array is set to
 * true.
 *
 * One could think of making this a std::pair, and pair the
 * rank with the error but for now this works.
 */
#define MAX_NUMBER_OF_THREADS 500
static bool error[MAX_NUMBER_OF_THREADS];

// void my_errhandler(MPI_Comm* Comm, int* errcode)
// {
// 	int myid;
// 	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
// 	error[myid] = true;
// 
// 	/**
// 	 * In case we need a detailed error description just uncomment
// 	 * the next lines.
// 	 */
// 	/**
// 	 * char      hostname[80];
// 	 * char      errstring[MPI_MAX_ERROR_STRING];
// 	 * int	      resultlen, myid;
// 	 * 
// 	 * MPI_Comm_rank(MPI_COMM_WORLD, &myid);
// 	 * MPI_Error_string(*errcode, errstring, &resultlen);
// 	 * printf("Proc %d on %-15.12s: UserErrorhandler - ErrorCode %d Message %s\n",
// 	 * 	myid, hostname, *errcode, errstring);
// 	 */
// 
// }

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
	/* TODO */
	m_myDomain_dim = Dimension(0,0);
	m_myDomainFirstCellColor = Color::Red;
	m_myOffsetToGlobalDomain = Dimension(0,0);

	/* set send/recv buffers. the longest side of domain will be buffers size. */
	/* TODO */



//	/**
//	 * Create the Grid of threads that compute the solution.
//	 * If we have more threads than equal size patches we
//	 * "kill" the threads that are left over, by letting them
//	 * call the error handling function.
//	 */
//	int dim[2]; // dim[0] number of patches in x direction, dim[1] number of patches in y direction
//	int periods[2] = {false,false}; // our grid isn't periodic isn't it?	
//	if (SqrtIsEven(m_numProcs))
//	{
//		/**
//		 * We can simply divide the grid into equal size
//		 * patches by setting the number of patches along
//		 * the x and y axes to sqrt(m_numProcs).
//		 */
//		dim[0] = (int)sqrt(m_numProcs);
//		dim[1] = (int)sqrt(m_numProcs);
//	}
//	else
//	{ //if size = 2 then dims = 2, 1; size = 4 then 2,2; 8 = 4, 2...
//		dim[1] = (int)sqrt(m_numProcs + m_numProcs);
//		dim[0] = dim[1] / 2;
//	}
//	printf("%d | %d \n", dim[0], dim[1]);
//
//	/**
//	 * Error handling
//	 * We create a handle to the function my_errhandler and a handle
//	 * for the MPI_COMM_WORLD. If there is an error in the MPI_Cart_create
//	 * function our function is called and the thread is terminated.
//	 */
//	MPI_Errhandler errh, store_errh;
//	MPI_Comm_get_errhandler(MPI_COMM_WORLD, &store_errh);
//	MPI_Comm_create_errhandler((MPI_Handler_function *)(&my_errhandler), &errh);
//	MPI_Comm_set_errhandler(MPI_COMM_WORLD, errh);
//
//	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, periods, 1, &m_comm);
//
//	MPI_Errhandler_free(&errh);
//	MPI_Comm_set_errhandler(m_comm, store_errh);
//
//	/**
//	 * We have found an error so set the Communication to
//	 * not valid and return. The Thread is than terminated
//	 * in the main.cpp.
//	 */
//	if(error[m_myRank])
//	{
//		m_valid = false;
//		return;
//	}
//
//	/**
//	 * This is just here for testing!
//	 */
//	int coordinates[2];
//	MPI_Cart_coords(m_comm, m_myRank, 2, coordinates);
//	MPI_Cart_shift(m_comm, 1, -1, &m_rightRank, &m_leftRank);
//	MPI_Cart_shift(m_comm, 0, -1, &m_downRank, &m_upRank);
//
}

void Communication::exchangeGridBoundaryValues(Domain domain, Handle grid, Color handleColorCells)
{
	/* We need an empty GridFunction constructor :)*/
	GridFunction& one = domain.u();
	GridFunction& two = domain.u();
	switch (grid)
	{
	case Communication::Handle::Pressure:
		one = domain.p();
		ExchangeOneGridFunction(one, domain,
			domain.getBeginInnerDomains(), domain.getEndInnerDomainP());
		break;
	case Communication::Handle::Velocities:
		one = domain.u();
		two = domain.v();
		ExchangeTwoGridFunctions(one, two, domain,
			domain.getBeginInnerDomains(), domain.getEndInnerDomainU(), 
			domain.getBeginInnerDomains(), domain.getEndInnerDomainV());
		break;
	case Communication::Handle::PreliminaryVelocities:
		one = domain.F();
		two = domain.G();
		ExchangeTwoGridFunctions(one, two, domain,
			domain.getBeginInnerDomains(), domain.getEndInnerDomainU(),
			domain.getBeginInnerDomains(), domain.getEndInnerDomainV());
		break;
	default:
		break;
	}
}

void Communication::exchangeGridInnerValues(Domain domain, Handle grid)
{
	/* We need an empty GridFunction constructor :)*/
	GridFunction& one = domain.u();
	GridFunction& two = domain.u();
	switch (grid)
	{
	case Communication::Handle::Pressure:
		one = domain.p();
		ExchangeOneGridFunction(one, domain,
			domain.getBeginInnerDomains(), domain.getEndInnerDomainP());
		break;
	case Communication::Handle::Velocities:
		one = domain.u();
		two = domain.v();
		ExchangeTwoGridFunctions(one, two, domain,
			domain.getBeginInnerDomains(), domain.getEndInnerDomainU(),
			domain.getBeginInnerDomains(), domain.getEndInnerDomainV());
		break;
	case Communication::Handle::PreliminaryVelocities:
		one = domain.F();
		two = domain.G();
		ExchangeTwoGridFunctions(one, two, domain,
			domain.getBeginInnerDomains(), domain.getEndInnerDomainU(),
			domain.getBeginInnerDomains(), domain.getEndInnerDomainV());
		break;
	default:
		break;
	}
}

bool Communication::checkForAnotherSORCycle(Real mySubResiduum)
{
	if (m_myRank != 0)
	{
		Real buffer[1];
		buffer[0] = mySubResiduum;
		m_sendToOne(1, 0, 0);
	}
	else
	{
		std::vector<Real> recvbuf;
		m_recvFromAll(&recvbuf, 1);
		std::vector<Real> Residuuen;
		Real tmp_residuum;
		for (unsigned int i = 0; i < recvbuf.size(); i++)
		{
			tmp_residuum = recvbuf.at(i);
		}
		// do something with the data
	}
}

Real Communication::getGlobalTimeStep(Delta myMaxValues)
{
	if (m_myRank != 0)
	{
		Real buffer[3];
		buffer[0] = myMaxValues.x;
		buffer[1] = myMaxValues.y;
		buffer[2] = myMaxValues.z;
		m_sendToOne(3, 0, 0);
	}
	else
	{
		std::vector<Real> recvbuf;
		m_recvFromAll(&recvbuf, 3);
		std::vector<Delta> maxValues;
		Delta tmp_maxValue;
		for (unsigned int i = 0; i < recvbuf.size(); i = i + 3)
		{
			tmp_maxValue.x = recvbuf.at(i);
			tmp_maxValue.y = recvbuf.at(i + 1);
			tmp_maxValue.z = recvbuf.at(i + 2);
		}
		// todo: do something with the data.
	}
}

void Communication::m_sendToOne(int count, int dest, int tag)
{
	MPI_Send(m_sendBuffer, count, /*TODO switch*/MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
}

int Communication::m_recvFromOne(int source, int tag)
{
	MPI_Recv(m_recvBuffer, m_bufferSize, /*TODO switch*/MPI_DOUBLE, source, tag, MPI_COMM_WORLD, NULL);

	return 0; /*???*/
}

void Communication::m_sendToAll(void* buf, int count)
{
	MPI_Bcast(buf, count, MPI_DOUBLE, m_myRank, MPI_COMM_WORLD);
}

void Communication::m_recvFromAll(/*void* sendbuf, int sendcount,*/ void* recvbuf, int recvcount)
{
	std::vector<Real> data;
	for (int i = 1; i < m_numProcs; i++)
	{
		m_recvFromOne(i,m_myRank);
		for (int j = 0; j < recvcount; j++)
		{
			data.push_back(m_recvBuffer[j]);
		}
	}
	recvbuf = &data;
	
	///**
	// * Gather the data from all threads in m_comm and distribute the combined data
	// * to all threads.
	// * We could use a foor loop and loop over all other threads in m_comm and just
	// * recieve the data but i think this function could prove usefull.
	// */
	//MPI_Allgather(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcount, MPI_DOUBLE, m_comm);
}

void Communication::ExchangeTwoGridFunctions(GridFunction& one, GridFunction& two, Domain domain,
	Dimension dim_one_begin, Dimension dim_one_end, Dimension dim_two_begin, Dimension dim_two_end)
{
	// a)
	if (m_leftRank != -1)
	{
		// SERGEJ CHECK/FIX THIS
		std::vector<Real> buffer;
		buffer.reserve(dim_one_end[1] - dim_one_begin[1] + 1); //TODO!!!! check if this size is correct
		// TODO fill buffer with u values from left border
		for (int i = dim_one_begin[1]; i <= dim_one_end[1]; i++) //DOUBLE-TODO!!!! check if this size is correct
			buffer.push_back(one(dim_one_begin[1], i));

		// TODO send buffer to m_leftRank
		m_sendBuffer = buffer.data();
		m_sendToOne((int)buffer.size(), m_leftRank,m_myRank);

		// TODO receive buffer from m_rightRank
		m_recvFromOne(m_rightRank, m_myRank);
		// do something with the values in the buffer
		buffer.clear();
		buffer.reserve(dim_two_end[1] - dim_two_begin[1] + 1); //TODO!!!! check if this size is correct

		// TODO fill buffer with v values from left border
		for (int i = dim_two_begin[1]; i <= dim_two_end[1]; i++)
			buffer.push_back(two(dim_two_begin[1], i));

		// TODO send buffer to m_leftRank
		m_sendBuffer = buffer.data();
		m_sendToOne((int)buffer.size(), m_leftRank, m_myRank);

		// TODO receive buffer from m_rightRank
		m_recvFromOne(m_rightRank, m_myRank);
		// do something with the values in the buffer
	}
	else
	{
		// TODO ask domain to set left boundary values
	}

	// b)
	if (m_rightRank != -1)
	{
		std::vector<Real> buffer;
		buffer.reserve(dim_one_end[1] - dim_one_begin[1] + 1); //TODO!!!! check if this size is correct
		// TODO fill buffer with u values from left border
		for (int i = dim_one_begin[1]; i <= dim_one_end[1]; i++) //DOUBLE-TODO!!!! check if this size is correct
			buffer.push_back(one(dim_one_end[1], i));

		// TODO send buffer to m_rightRank
		m_sendBuffer = buffer.data();
		m_sendToOne((int)buffer.size(), m_rightRank, m_myRank);

		// TODO receive buffer from m_leftRank
		m_recvFromOne(m_leftRank, m_myRank);
		// do something with the values in the buffer
		buffer.clear();
		buffer.reserve(dim_two_end[1] - dim_two_begin[1] + 1); //TODO!!!! check if this size is correct

		// TODO fill buffer with v values from right border
		for (int i = dim_two_begin[1]; i <= dim_two_end[1]; i++)
			buffer.push_back(two(dim_two_end[1], i));

		// TODO send buffer to m_rightRank
		m_sendBuffer = buffer.data();
		m_sendToOne((int)buffer.size(), m_rightRank, m_myRank);

		// TODO receive buffer from m_leftRank
		m_recvFromOne(m_leftRank, m_myRank);
		// do something with the values in the buffer
	}
	else
	{
		// TODO ask domain to set right boundary values
	}

	// c)
	if (m_upRank != -1)
	{
		std::vector<Real> buffer;
		buffer.reserve(dim_one_end[0] - dim_one_begin[0] + 1); //TODO!!!! check if this size is correct
		// TODO fill buffer with values from upper border
		for (int i = dim_one_begin[0]; i <= dim_one_end[0]; i++) //DOUBLE-TODO!!!! check if this size is correct
			buffer.push_back(one(i, dim_one_end[0]));

		// TODO send buffer to m_upRank
		m_sendBuffer = buffer.data();
		m_sendToOne((int)buffer.size(), m_upRank, m_myRank);

		// TODO receive buffer from m_downRank
		m_recvFromOne(m_downRank, m_myRank);
		// do something with the values in the buffer
		buffer.clear();
		buffer.reserve(dim_two_end[0] - dim_two_begin[0] + 1); //TODO!!!! check if this size is correct

		// TODO fill buffer with values from upper border
		for (int i = dim_two_begin[0]; i <= dim_two_end[0]; i++) //DOUBLE-TODO!!!! check if this size is correct
			buffer.push_back(two(i, dim_two_end[0]));

		// TODO send buffer to m_upRank
		m_sendBuffer = buffer.data();
		m_sendToOne((int)buffer.size(), m_upRank, m_myRank);

		// TODO receive buffer from m_downRank
		m_recvFromOne(m_downRank, m_myRank);
		// do something with the values in the buffer
	}
	else
	{
		// TODO ask domain to set upper boundary values
	}

	// d)
	if (m_downRank != -1)
	{
		std::vector<Real> buffer;
		buffer.reserve(dim_one_end[0] - dim_one_begin[0] + 1); //TODO!!!! check if this size is correct
		// TODO fill buffer with values from upper border
		for (int i = dim_one_begin[0]; i <= dim_one_end[0]; i++) //DOUBLE-TODO!!!! check if this size is correct
			buffer.push_back(one(dim_one_end[0], i));

		// TODO send buffer to m_downRank
		m_sendBuffer = buffer.data();
		m_sendToOne((int)buffer.size(), m_downRank, m_myRank);

		// TODO receive buffer from m_upRank
		m_recvFromOne(m_upRank, m_myRank);
		// do something with the values in the buffer
		buffer.clear();
		buffer.reserve(dim_two_end[0] - dim_two_begin[0] + 1); //TODO!!!! check if this size is correct

		// TODO fill buffer with values from upper border
		for (int i = dim_two_begin[0]; i <= dim_two_end[0]; i++) //DOUBLE-TODO!!!! check if this size is correct
			buffer.push_back(two(dim_two_end[0], i));

		// TODO send buffer to m_downRank
		m_sendToOne((int)buffer.size(), m_downRank, m_myRank);

		// TODO receive buffer from m_upRank
		m_recvFromOne(m_upRank, m_myRank);
		// do something with the values in the buffer
	}
	else
	{
		// TODO ask domain to set lower boundary values
	}
}

void Communication::ExchangeOneGridFunction(GridFunction& one, Domain domain,
	Dimension dim_one_begin, Dimension dim_one_end)
{
	// a)
	if (m_leftRank != -1)
	{
		// SERGEJ CHECK/FIX THIS
		std::vector<Real> buffer;
		buffer.reserve(dim_one_end[1] - dim_one_begin[1] + 1); //TODO!!!! check if this size is correct
		// TODO fill buffer with u values from left border
		for (int i = dim_one_begin[1]; i <= dim_one_end[1]; i++) //DOUBLE-TODO!!!! check if this size is correct
			buffer.push_back(one(dim_one_begin[1], i));

		// TODO send buffer to m_leftRank
		m_sendBuffer = buffer.data();
		m_sendToOne((int)buffer.size(), m_leftRank, m_myRank);

		// TODO receive buffer from m_rightRank
		m_recvFromOne(m_rightRank, m_myRank);
		// do something with the values in the buffer
	}
	else
	{
		// TODO ask domain to set left boundary values
	}

	// b)
	if (m_rightRank != -1)
	{
		std::vector<Real> buffer;
		buffer.reserve(dim_one_end[1] - dim_one_begin[1] + 1); //TODO!!!! check if this size is correct
		// TODO fill buffer with u values from left border
		for (int i = dim_one_begin[1]; i <= dim_one_end[1]; i++) //DOUBLE-TODO!!!! check if this size is correct
			buffer.push_back(one(dim_one_end[1], i));

		// TODO send buffer to m_rightRank
		m_sendBuffer = buffer.data();
		m_sendToOne((int)buffer.size(), m_rightRank, m_myRank);

		// TODO receive buffer from m_leftRank
		m_recvFromOne(m_leftRank, m_myRank);
		// do something with the values in the buffer
	}
	else
	{
		// TODO ask domain to set right boundary values
	}

	// c)
	if (m_upRank != -1)
	{
		std::vector<Real> buffer;
		buffer.reserve(dim_one_end[0] - dim_one_begin[0] + 1); //TODO!!!! check if this size is correct
		// TODO fill buffer with values from upper border
		for (int i = dim_one_begin[0]; i <= dim_one_end[0]; i++) //DOUBLE-TODO!!!! check if this size is correct
			buffer.push_back(one(i, dim_one_end[0]));

		// TODO send buffer to m_upRank
		m_sendBuffer = buffer.data();
		m_sendToOne((int)buffer.size(), m_upRank, m_myRank);

		// TODO receive buffer from m_downRank
		m_recvFromOne(m_downRank, m_myRank);
		// do something with the values in the buffer
	}
	else
	{
		// TODO ask domain to set upper boundary values
	}

	// d)
	if (m_downRank != -1)
	{
		std::vector<Real> buffer;
		buffer.reserve(dim_one_end[0] - dim_one_begin[0] + 1); //TODO!!!! check if this size is correct
		// TODO fill buffer with values from upper border
		for (int i = dim_one_begin[0]; i <= dim_one_end[0]; i++) //DOUBLE-TODO!!!! check if this size is correct
			buffer.push_back(one(dim_one_end[0], i));

		// TODO send buffer to m_downRank
		m_sendBuffer = buffer.data();
		m_sendToOne((int)buffer.size(), m_downRank, m_myRank);

		// TODO receive buffer from m_upRank
		m_recvFromOne(m_upRank, m_myRank);
		// do something with the values in the buffer
	}
	else
	{
		// TODO ask domain to set lower boundary values
	}
}

Communication::~Communication()
{
	MPI_Finalize();
}







//if (m_numProcs == 1)
//{
//	/**
//	* We just have one patch, so we only have to create
//	* one Domain with the given dimensions.
//	*/
//	dim[0] = 1;
//	dim[1] = 1;
//}
//else if (m_numProcs == 2)
//{
//	/**
//	* We create 2 patchen with size x/2, y.
//	*/
//	dim[0] = 2;
//	dim[1] = 1;
//}
//else if (SqrtIsEven(m_numProcs))
//{
//	/**
//	* We can simply divide the grid into equal size
//	* patches by setting the number of patches along
//	* the x and y axes to sqrt(m_numProcs).
//	*/
//	dim[0] = (int)sqrt(m_numProcs);
//	dim[1] = (int)sqrt(m_numProcs);
//}
//else if ((m_numProcs % 2) == 1)
//{
//	/**
//	* We can only create one row of patches because
//	* the pathces have to be of the same size.
//	*/
//	dim[0] = m_numProcs;
//	dim[1] = 1;
//}
//else if ((m_numProcs % 2) == 0)
//{
//	/**
//	* We can at least create 2 rows of patches of
//	* the same size.
//	*/
//	if ((m_numProcs % 4) == 0)
//	{
//		/**
//		* We can create 4 rows.
//		*/
//		dim[0] = 4;
//		dim[1] = (int)(m_numProcs / 4);
//	}
//	else
//	{
//		/**
//		* We can create 2 rows.
//		*/
//		dim[0] = 2;
//		dim[1] = (int)(m_numProcs / 2);
//	}
//}











// /**
//  * Todo:
//  * I don't know what to to with the size, reorder, dims and periods
//  * variables (they where in the UML of the Communication class).
//  * If we don't need them Communication can be turned into a namespace.
//  * 
//  * I used the MPI_Send and MPI_Recv functions as described in the 
//  * "MPI lecture". But i'm not sure if the data is automatically copied
//  * into the given GridFunction (Send and Recv expect a pointer to the
//  * data). If yes this implementation should work. If not there is more
//  * work to be done.
//  *
//  * The fifth argument in MPI_Send and MPI_Recv could be set to the rank
//  * of the corresponding process. But so far i don't know how the ranks
//  * are assigned, therefore the value is 0.
//  */
// 
// Communication::Communication(int size, int reorder, MPI_Comm comm_cart)
// {
// 	this->m_size = size;
// 	this->m_reorder = reorder;
// 	this->m_comm_cart = comm_cart;
// 	this->m_dims = new int[this->m_size];
// 	this->m_periods = new int[this->m_size];
// }
// 
// void Communication::ExchangePValues(int rank_to, int rank_from, 
// 	int count_p, GridFunction& p, GridFunction& p_fromleft)
// {
// 	MPI_Send(&p, count_p, MPI_DOUBLE, rank_to, 0, m_comm_cart);
// 
// 	MPI_Recv(&p_fromleft, count_p, MPI_DOUBLE, rank_from, 0, m_comm_cart, &status);
// }
// 
// void Communication::ExchangeUVValues(int rank_to, int rank_from, 
// 	int count_u, GridFunction& u, int count_v, GridFunction& v,
// 	GridFunction& u_fromleft, GridFunction& v_fromleft)
// {
// 	MPI_Send(&u, count_u, MPI_DOUBLE, rank_to, 0, m_comm_cart);
// 	MPI_Send(&v, count_v, MPI_DOUBLE, rank_to, 0, m_comm_cart);
// 
// 	MPI_Recv(&u_fromleft, count_u, MPI_DOUBLE, rank_from, 0, m_comm_cart, &status);
// 	MPI_Recv(&u_fromleft, count_v, MPI_DOUBLE, rank_from, 0, m_comm_cart, &status);
// }
// 
// void Communication::ExchangeFGValues(int rank_to, int rank_from, 
// 	int count_F, GridFunction& F, int count_G, GridFunction& G,
// 	GridFunction& F_fromleft, GridFunction& G_fromleft)
// {
// 	MPI_Send(&F, count_F, MPI_DOUBLE, rank_to, 0, m_comm_cart);
// 	MPI_Send(&G, count_G, MPI_DOUBLE, rank_to, 0, m_comm_cart);
// 
// 	MPI_Recv(&F_fromleft, count_F, MPI_DOUBLE, rank_from, 0, m_comm_cart, &status);
// 	MPI_Recv(&G_fromleft, count_G, MPI_DOUBLE, rank_from, 0, m_comm_cart, &status);
// }
// 
// void Communication::ExchangeRHSValues(int rank_to, int rank_from, 
// 	int count_rhs, GridFunction& rhs, GridFunction& rhs_fromleft)
// {
// 	MPI_Send(&rhs, count_rhs, MPI_DOUBLE, rank_to, 0, m_comm_cart);
// 
// 	MPI_Recv(&rhs_fromleft, count_rhs, MPI_DOUBLE, rank_from, 0, m_comm_cart, &status);
// }
