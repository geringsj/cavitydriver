#include "GridFunction.hpp"

#include <mpi.h>

class Communication
{
public:
	Communication(int size, int reorder, MPI_Comm comm_cart);
	void ExchangePValues(int rank_to, int rank_from, int count_p, GridFunction& p,
		GridFunction& p_fromleft);
	void ExchangeUVValues(int rank_to, int rank_from, int count_u, GridFunction& u, 
		int count_v, GridFunction& v, GridFunction& u_fromleft, GridFunction& v_fromleft);
	void ExchangeFGValues(int rank_to, int rank_from, int count_F, GridFunction& F, 
		int count_G, GridFunction& G, GridFunction& F_fromleft, GridFunction& G_fromleft);
	void ExchangeRHSValues(int rank_to, int rank_from, int count_rhs, GridFunction& rhs, 
		GridFunction& rhs_fromleft);
private:
	MPI_Comm m_comm_cart;
	MPI_Status status;
	int m_size;
	int* m_dims;
	int* m_periods;
	int m_reorder;
};