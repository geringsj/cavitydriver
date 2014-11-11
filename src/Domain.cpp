#include "Domain.hpp"


Domain::Domain(Dimension dimension, Point delta,

	std::function<Real(Index,GridFunction&,Dimension)> in_u,
	std::function<Real(Index,GridFunction&,Dimension)> in_v,
	std::function<Real(Index, GridFunction&, Dimension)> in_w,
	std::function<Real(Index, GridFunction&, Dimension)> in_p,

	std::function<Real(Point)> in_gx,
	std::function<Real(Point)> in_gy,
	std::function<Real(Point)> in_gz
	) : 
		m_dimension(dimension),
		m_inner_begin(1,1,1), 
		m_delta(delta),
		m_p(dimension[0]+2,dimension[1]+2), 
		m_p_rhs(dimension[0]+2,dimension[1]+2),
		m_velocities(dimension),
		m_preliminary_velocities_FGH(dimension),
		m_infunc_gx(in_gx), m_infunc_gy(in_gy), m_infunc_gz(in_gz)
{
	for(uint d=0; d<DIMENSIONS; d++)
		for(uint d2=0; d2<DIMENSIONS; d2++)
			m_inner_end[d][d2] = dimension[d2] + ((d2==d)?(-1):(0)); 
	m_inner_end[3] = dimension; /* pressure */

	for(uint d=0; d<DIMENSIONS; d++)
	{
		Dimension max(m_inner_end[d][0]+1, m_inner_end[d][1]+1, m_inner_end[d][2]+1);

		if(d==0)
			m_infunc_u = std::bind(in_u,std::placeholders::_1, std::ref(u()), max);
		if(d==1)
			m_infunc_v = std::bind(in_v,std::placeholders::_1, std::ref(v()), max);
		if(d==2)
			m_infunc_w = std::bind(in_w,std::placeholders::_1, std::ref(w()), max);
	}

	/* init all to start */
	for(uint d=0; d<DIMENSIONS; d++)
	{
		for(int i=m_inner_begin[0]; i<=m_inner_end[d][0]; i++)
			for(int j=m_inner_begin[1]; j<=m_inner_end[d][1]; j++)
				//for(int k=m_inner_begin[2]; k<=m_inner_end[d][2]; k++)
		{
			Dimension current(i, j);//,k);
			if(d==0)
				u()(i,j)/*,k)*/ = m_infunc_u(current);
			if(d==1)
				v()(i, j)/*,k)*/ = m_infunc_v(current);
			if(d==3)
				w()(i, j)/*,k)*/ = m_infunc_w(current);
		}
	}
	for (int i = m_inner_begin[0]; i <= m_inner_end[3][0]; i++)
		for (int j = m_inner_begin[1]; j <= m_inner_end[3][1]; j++)
			//for (int k = m_inner_begin[2]; k <= m_inner_end[3][2]; k++)
				p()(i,j) = in_p(Index(i, j), p(), Dimension(dimension.i + 2, dimension.j + 2, dimension.k + 2));
	/* and boundaries */
	this->setVelocitiesBoundaries();
	this->setPreliminaryVelocitiesBoundaries();
	this->setPressureBoundaries();
}

Domain::~Domain()
{
}
