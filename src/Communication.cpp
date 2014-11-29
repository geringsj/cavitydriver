#include "Communication.hpp"


Communication::Communication(Dimension globalDomainDim, int argc, char** argv)
{
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
		std::cout << "good\n";
	}
	else
	{ //if size = 2 then dims = 2, 1; size = 4 then 2,2; 8 = 4, 2...
		dim[1] = (int)sqrt(m_numProcs + m_numProcs);
		dim[0] = dim[1] / 2;
	}

	printf("%d | %d \n", dim[0], dim[1]);

	MPI_Dims_create(m_numProcs, 2, dim);
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, periods, 1, &comm);
	MPI_Barrier(MPI_COMM_WORLD);

	int coordinates[2];
	MPI_Cart_coords(comm, m_myRank, 2, coordinates);
	
	MPI_Cart_shift(comm, 1, -1, &m_rightRank, &m_leftRank);
	MPI_Cart_shift(comm, 0, -1, &m_downRank, &m_upRank);
	
	
	
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
	//int coordinates[2];
	//int shiftsource, shiftdest;
	//int ix, swap;
	//MPI_Status status;

	//MPI_Cart_coords(comm, m_myRank, 2, coordinates);
	//
	//MPI_Cart_shift(comm, 1, -1, &m_rightRank, &m_leftRank);
	//MPI_Cart_shift(comm, 0, -1, &m_downRank, &m_upRank);
	//
	//
	//
	//printf("rank = %d, rightrank = %d, downrank = %d, leftrank = %d, uprank = %d\n", m_myRank, m_rightRank, m_downRank, m_leftRank, m_upRank);



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
