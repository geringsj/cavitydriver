
#include "Communication.hpp"

#include <mpi.h>


Communication::Communication(Dimension globalDomainDim)
{
}

~Communication::Communication()
{
}

















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
