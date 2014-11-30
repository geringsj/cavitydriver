#include "Communication.hpp"

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

void my_errhandler(MPI_Comm* Comm, int* errcode)
{
	int myid;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	error[myid] = true;

	/**
	 * In case we need a detailed error description just uncomment
	 * the next lines.
	 */
	/**
	 * char      hostname[80];
	 * char      errstring[MPI_MAX_ERROR_STRING];
	 * int	      resultlen, myid;
	 * 
	 * MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	 * MPI_Error_string(*errcode, errstring, &resultlen);
	 * printf("Proc %d on %-15.12s: UserErrorhandler - ErrorCode %d Message %s\n",
	 * 	myid, hostname, *errcode, errstring);
	 */

}

Communication::Communication(Dimension globalDomainDim, int argc, char** argv)
{
	m_valid = true;
	m_globalDomain_dim = globalDomainDim;
	/**
	 * MPI initialization
	 * We need argc and argv here so i had to pass them.
	 * rc stores the errors that can occur during int.
	 * m_numProcs stores the number of threads we create.
	 * m_myRank stores the rank of each thread so we 
	 *			know who the communication belongs to.
	 * rc stores the return value from MPI_Init.
	 *
	 * Don't forget to call the destructor from Communication
	 * because it calls MPI_Finalize to end the MPI programm.
	 */
	int rc = MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS) {
		printf("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	MPI_Comm_size(MPI_COMM_WORLD, &m_numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &m_myRank);
	printf("Number of tasks= %d My rank= %d \n", m_numProcs, m_myRank);

	/**
	 * Create the Grid of threads that compute the solution.
	 * If we have more threads than equal size patches we
	 * "kill" the threads that are left over, by letting them
	 * call the error handling function.
	 */
	int dim[2]; // dim[0] number of patches in x direction, dim[1] number of patches in y direction
	int periods[2] = {false,false}; // our grid isn't periodic isn't it?	
	if (SqrtIsEven(m_numProcs))
	{
		/**
		 * We can simply divide the grid into equal size
		 * patches by setting the number of patches along
		 * the x and y axes to sqrt(m_numProcs).
		 */
		dim[0] = (int)sqrt(m_numProcs);
		dim[1] = (int)sqrt(m_numProcs);
	}
	else
	{ //if size = 2 then dims = 2, 1; size = 4 then 2,2; 8 = 4, 2...
		dim[1] = (int)sqrt(m_numProcs + m_numProcs);
		dim[0] = dim[1] / 2;
	}
	printf("%d | %d \n", dim[0], dim[1]);

	/**
	 * Error handling
	 * We create a handle to the function my_errhandler and a handle
	 * for the MPI_COMM_WORLD. If there is an error in the MPI_Cart_create
	 * function our function is called and the thread is terminated.
	 */
	MPI_Errhandler errh, store_errh;
	MPI_Comm_get_errhandler(MPI_COMM_WORLD, &store_errh);
	MPI_Comm_create_errhandler((MPI_Handler_function *)(&my_errhandler), &errh);
	MPI_Comm_set_errhandler(MPI_COMM_WORLD, errh);

	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, periods, 1, &m_comm);

	MPI_Errhandler_free(&errh);
	MPI_Comm_set_errhandler(m_comm, store_errh);

	/**
	 * We have found an error so set the Communication to
	 * not valid and return. The Thread is than terminated
	 * in the main.cpp.
	 */
	if(error[m_myRank])
	{
		m_valid = false;
		return;
	}

	/**
	 * This is just here for testing!
	 */
	int coordinates[2];
	MPI_Cart_coords(m_comm, m_myRank, 2, coordinates);
	MPI_Cart_shift(m_comm, 1, -1, &m_rightRank, &m_leftRank);
	MPI_Cart_shift(m_comm, 0, -1, &m_downRank, &m_upRank);

	printf("rank = %d, rightrank = %d, downrank = %d, leftrank = %d, uprank = %d\n", m_myRank, m_rightRank, m_downRank, m_leftRank, m_upRank);
}

bool Communication::SqrtIsEven(int number){
	double root = sqrt(number);
	int root_int = (int)floor(root);
	root_int = (int)pow(root_int, 2);
	if (number == root_int)	return true;
	else return false;
}

void Communication::exchangeGridBoundaryValues(Domain domain, Handle grid, Color handleColorCells)
{
	int coordinates[2];
	MPI_Cart_coords(m_comm, m_myRank, 2, coordinates);
	MPI_Cart_shift(m_comm, 1, -1, &m_rightRank, &m_leftRank);
	MPI_Cart_shift(m_comm, 0, -1, &m_downRank, &m_upRank);
	
	printf("rank = %d, rightrank = %d, downrank = %d, leftrank = %d, uprank = %d\n", m_myRank, m_rightRank, m_downRank, m_leftRank, m_upRank);

	switch (grid)
	{
	case Communication::Handle::Pressure:
		break;
	case Communication::Handle::Velocities:
	{
		GridFunction& u = domain.u();
		GridFunction& v = domain.v();
		MPI_Status status;

		// a)
		if(m_leftRank != -1)
		{
			// SERGEJ CHECK/FIX THIS
			std::vector<double> buffer;
			buffer.reserve(domain.getBeginInnerDomains()[1]-domain.getEndInnerDomainU()[1]+1); //TODO!!!! check if this size is correct
			// TODO fill buffer with u values from left border
			for(int i=domain.getBeginInnerDomains()[1]; i <= domain.getEndInnerDomainU()[1]; i++) //DOUBLE-TODO!!!! check if this size is correct
				buffer.push_back(u(domain.getBeginInnerDomains()[0],i));

			// TODO send buffer to m_leftRank
			m_sendToOne(&buffer, (int)buffer.size(), m_leftRank, m_myRank);

			// TODO receive buffer from m_rightRank
			m_recvFromOne(&buffer, (int)buffer.size(), m_rightRank,m_rightRank, &status);
			// do something with the values in the buffer
			buffer.clear();
			buffer.reserve(domain.getBeginInnerDomains()[1] - domain.getEndInnerDomainV()[1] + 1); //TODO!!!! check if this size is correct

			// TODO fill buffer with v values from left border
			for(int i=domain.getBeginInnerDomains()[1]; i <= domain.getEndInnerDomainV()[1]; i++)
				buffer.push_back(v(domain.getBeginInnerDomains()[1],i));

			// TODO send buffer to m_leftRank
			m_sendToOne(&buffer,(int)buffer.size(),m_leftRank,m_myRank);

			// TODO receive buffer from m_rightRank
			m_recvFromOne(&buffer, (int)buffer.size(), m_rightRank, m_rightRank, &status);
			// do something with the values in the buffer
		}
		else
		{
			// TODO ask domain to set left boundary values
		}

		// b)
		if(m_rightRank != -1)
		{
			std::vector<double> buffer;
			buffer.reserve(domain.getBeginInnerDomains()[1]-domain.getEndInnerDomainU()[1]+1); //TODO!!!! check if this size is correct
			// TODO fill buffer with u values from left border
			for(int i=domain.getBeginInnerDomains()[1]; i <= domain.getEndInnerDomainU()[1]; i++) //DOUBLE-TODO!!!! check if this size is correct
				buffer.push_back(u(domain.getEndInnerDomainU()[0],i));

			// TODO send buffer to m_rightRank
			m_sendToOne(&buffer, (int)buffer.size(), m_rightRank, m_myRank);

			// TODO receive buffer from m_leftRank
			m_recvFromOne(&buffer, (int)buffer.size(), m_leftRank, m_leftRank, &status);
			// do something with the values in the buffer
			buffer.clear();
			buffer.reserve(domain.getBeginInnerDomains()[1] - domain.getEndInnerDomainV()[1] + 1); //TODO!!!! check if this size is correct

			// TODO fill buffer with v values from left border
			for(int i=domain.getBeginInnerDomains()[1]; i <= domain.getEndInnerDomainV()[1]; i++)
				buffer.push_back(v(domain.getEndInnerDomainV()[0],i));

			// TODO send buffer to m_rightRank
			m_sendToOne(&buffer,(int)buffer.size(),m_rightRank,m_myRank);

			// TODO receive buffer from m_leftRank
			m_recvFromOne(&buffer, (int)buffer.size(), m_leftRank, m_leftRank, &status);
			// do something with the values in the buffer
		}
		else
		{
			// TODO ask domain to set right boundary values
		}

		// c)
		if(m_upRank != -1)
		{
			std::vector<double> buffer;
			buffer.reserve(domain.getBeginInnerDomains()[0] - domain.getEndInnerDomainU()[0] + 1); //TODO!!!! check if this size is correct
			// TODO fill buffer with values from upper border
			for (int i = domain.getBeginInnerDomains()[0]; i <= domain.getEndInnerDomainU()[0]; i++) //DOUBLE-TODO!!!! check if this size is correct
				buffer.push_back(u(i,domain.getEndInnerDomainU()[0]));

			// TODO send buffer to m_upRank
			m_sendToOne(&buffer, (int)buffer.size(), m_upRank, m_myRank);

			// TODO receive buffer from m_downRank
			m_recvFromOne(&buffer, (int)buffer.size(), m_downRank, m_downRank, &status);
			// do something with the values in the buffer
			buffer.clear();
			buffer.reserve(domain.getBeginInnerDomains()[0] - domain.getEndInnerDomainV()[0] + 1); //TODO!!!! check if this size is correct

			// TODO fill buffer with values from upper border
			for (int i = domain.getBeginInnerDomains()[0]; i <= domain.getEndInnerDomainV()[0]; i++) //DOUBLE-TODO!!!! check if this size is correct
				buffer.push_back(v(i, domain.getEndInnerDomainV()[0]));

			// TODO send buffer to m_upRank
			m_sendToOne(&buffer, (int)buffer.size(), m_upRank, m_myRank);

			// TODO receive buffer from m_downRank
			m_recvFromOne(&buffer, (int)buffer.size(), m_downRank, m_downRank, &status);
			// do something with the values in the buffer
		}
		else
		{
			// TODO ask domain to set upper boundary values
		}

		// d)
		if(m_downRank != -1)
		{
			std::vector<double> buffer;
			buffer.reserve(domain.getBeginInnerDomains()[0] - domain.getEndInnerDomainU()[0] + 1); //TODO!!!! check if this size is correct
			// TODO fill buffer with values from upper border
			for (int i = domain.getBeginInnerDomains()[0]; i <= domain.getEndInnerDomainU()[0]; i++) //DOUBLE-TODO!!!! check if this size is correct
				buffer.push_back(u(domain.getEndInnerDomainU()[0],i));

			// TODO send buffer to m_downRank
			m_sendToOne(&buffer, (int)buffer.size(), m_downRank, m_myRank);

			// TODO receive buffer from m_upRank
			m_recvFromOne(&buffer, (int)buffer.size(), m_upRank, m_upRank, &status);
			// do something with the values in the buffer
			buffer.clear();
			buffer.reserve(domain.getBeginInnerDomains()[0] - domain.getEndInnerDomainV()[0] + 1); //TODO!!!! check if this size is correct

			// TODO fill buffer with values from upper border
			for (int i = domain.getBeginInnerDomains()[0]; i <= domain.getEndInnerDomainU()[0]; i++) //DOUBLE-TODO!!!! check if this size is correct
				buffer.push_back(v(domain.getEndInnerDomainV()[0], i));

			// TODO send buffer to m_downRank
			m_sendToOne(&buffer, (int)buffer.size(), m_downRank, m_myRank);

			// TODO receive buffer from m_upRank
			m_recvFromOne(&buffer, (int)buffer.size(), m_upRank, m_upRank, &status);
			// do something with the values in the buffer
		}
		else
		{
			// TODO ask domain to set lower boundary values
		}


		break;
	}
	case Communication::Handle::PreliminaryVelocities:
		break;
	default:
		break;
	}
}

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

bool Communication::checkForAnotherSORCycle(Real mySubResiduum)
{

}

Real Communication::getGlobalTimeStep(Delta myMaxValues)
{

}

void Communication::m_sendToOne(void* buf, int count, int dest, int tag)
{
	MPI_Send(&buf, count, MPI_DOUBLE, dest, tag, m_comm);
}

void Communication::m_recvFromOne(void* buf, int count, int source, int tag, MPI_Status* status)
{
	MPI_Recv(&buf, count, MPI_DOUBLE, source, tag, m_comm, status);
}

void Communication::m_sendToAll(void* buf, int count)
{
	MPI_Bcast(buf, count, MPI_DOUBLE, m_myRank, m_comm);
}

void Communication::m_recvFromAll(void* sendbuf, int sendcount, void* recvbuf, int recvcount)
{
	/**
	 * Gather the data from all threads in m_comm and distribute the combined data
	 * to all threads.
	 * We could use a foor loop and loop over all other threads in m_comm and just
	 * recieve the data but i think this function could prove usefull.
	 */
	MPI_Allgather(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcount, MPI_DOUBLE, m_comm);
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
