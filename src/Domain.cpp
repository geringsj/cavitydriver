#include "Domain.hpp"


Domain::Domain(Dimension dimension, Delta delta,

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
		m_p(Dimension(dimension[0]+2,dimension[1]+2)), 
		m_p_rhs(Dimension(dimension[0]+2,dimension[1]+2)),
		m_velocities(dimension),
		m_preliminary_velocities_FGH(dimension),
		m_forcefunc_gx(in_gx), m_forcefunc_gy(in_gy), m_forcefunc_gz(in_gz)
{
	/* set indices for end of inner of grids u, v, w and p individually */
	/* for now, everything < inner_begin and everything > inner_end 
	 * is border of the grid */
	for(uint d=0; d<DIMENSIONS; d++)
		for(uint e=0; e<DIMENSIONS; e++)
			m_inner_end[d][e] = dimension[e] + ((e==d)?(-1):(0)); 

	/* the pressure is fourth entry of the array and has symmetric dimensions */
	m_inner_end[3] = dimension; 

	/* bind given border functions to operate on gridfunctions u(), v(), w()
	 * and to automatically use the right (upper bounds) dimensions of the grid */
	for(uint d=0; d<DIMENSIONS; d++)
	{
		Dimension max(m_inner_end[d][0]+1, m_inner_end[d][1]+1, m_inner_end[d][2]+1);
		if(d==0)
			m_borderfunc_u = std::bind(in_u,std::placeholders::_1, std::ref(u()), max);
		if(d==1)
			m_borderfunc_v = std::bind(in_v,std::placeholders::_1, std::ref(v()), max);
		if(d==2)
			m_borderfunc_w = std::bind(in_w,std::placeholders::_1, std::ref(w()), max);
	}

	/* init all grids to start values given from outside
	 * yes, the borderfunctions shall return us initial values if not 
	 * evaluated on the border. somebody should document this demand. */
	for(uint d=0; d<DIMENSIONS; d++)
	{
		forall(i,j,m_inner_begin,m_inner_end[d])
		{
			Dimension current(i, j);//,k);
			if(d==0)
				u()(i,j)/*,k)*/ = m_borderfunc_u(current);
			if(d==1)
				v()(i, j)/*,k)*/ = m_borderfunc_v(current);
			if(d==3)
				w()(i, j)/*,k)*/ = m_borderfunc_w(current);
		}
	}

	/* init pressure and rhs 
	 * we don't use the pressure init function in_p after this,
	 * maybe kick it and just init with zero? would the script be against us? */
	forall(i,j,m_inner_begin,m_inner_end[3])
	{
		p()(i,j) = in_p(Index(i,j), p(), 
				Dimension(dimension.i+2, dimension.j+2, dimension.k+2));
		rhs()(i,j) = 0.0;
	}

	/* and don't forget the boundaries */
	this->setVelocitiesBoundaries();
	this->setPreliminaryVelocitiesBoundaries();
	this->setPressureBoundaries();
}

Domain::~Domain()
{
}

#define LEFT(start,end)	for(int i = start[0]-1, j = start[1]; j <= end[1]; j++)
#define RIGHT(start,end) for(int i = end[0] + 1, j = start[1]; j <= end[1]; j++)
#define TOP(start,end) for(int j = end[1] + 1, i = start[0]; i <= end[0]; i++)
#define BOTTOM(start,end) for(int j = start[1]-1, i = start[0]; i <= end[0]; i++)

#define VELOCITIESBOUNDARIES(i,j,d) do{\
	;\
	if(d == 0)\
		u()(i, j)/*,k)*/ = m_borderfunc_u(Dimension(i,j));\
	if(d == 1)\
		v()(i, j)/*,k)*/ = m_borderfunc_v(Dimension(i,j));\
	if(d == 2)\
		w()(i, j)/*,k)*/ = m_borderfunc_w(Dimension(i,j));\
}while(0)

#define PRELIMINARYVELOCITIESBOUNDARIES(i,j,d) do{\
	if (d == 0)\
		F()(i, j) = u()(i, j);\
	if (d == 1)\
		G()(i, j) = v()(i, j);\
	if (d == 2)\
		H()(i, j) = w()(i, j);\
}while(0)

#define PRESSUREBOUNDARIES(i,j) do{\
	if (i == 0) p()(i, j) = p()(i + 1, j);\
	if (i == m_dimension[0] + 1) p()(i, j) = p()(i - 1, j);\
	if (j == m_dimension[1] + 1) p()(i, j) = p()(i, j - 1);\
	if (j == 0) p()(i, j) = p()(i, j + 1);\
}while(0)


void Domain::setVelocitiesBoundaries()
{
	Dimension current;
	for(uint d=0; d<DIMENSIONS; d++)
	{
		LEFT(m_inner_begin, m_inner_end[d])
		{
			VELOCITIESBOUNDARIES(i, j, d);
		}
		RIGHT(m_inner_begin, m_inner_end[d]) 
		{
			VELOCITIESBOUNDARIES(i, j, d);
		}
		TOP(m_inner_begin, m_inner_end[d])
		{
			VELOCITIESBOUNDARIES(i, j, d);
		}
		BOTTOM(m_inner_begin, m_inner_end[d])
		{
			VELOCITIESBOUNDARIES(i, j, d);
		}
	}
}

void Domain::setPreliminaryVelocitiesBoundaries()
{
	for (uint d = 0; d<DIMENSIONS; d++)
	{
		LEFT(m_inner_begin, m_inner_end[d])
		{
			PRELIMINARYVELOCITIESBOUNDARIES(i, j, d);
		}
		RIGHT(m_inner_begin, m_inner_end[d])
		{
			PRELIMINARYVELOCITIESBOUNDARIES(i, j, d);
		}
		TOP(m_inner_begin, m_inner_end[d])
		{
			PRELIMINARYVELOCITIESBOUNDARIES(i, j, d);
		}
		BOTTOM(m_inner_begin, m_inner_end[d])
		{
			PRELIMINARYVELOCITIESBOUNDARIES(i, j, d);
		}
	}
}

void Domain::setPressureBoundaries()
{
	LEFT(m_inner_begin, m_dimension)
	{
		PRESSUREBOUNDARIES(i, j);
	}
	RIGHT(m_inner_begin, m_dimension)
	{
		PRESSUREBOUNDARIES(i, j);
	}
	TOP(m_inner_begin, m_dimension)
	{
		PRESSUREBOUNDARIES(i, j);
	}
	BOTTOM(m_inner_begin, m_dimension)
	{
		PRESSUREBOUNDARIES(i, j);
	}
}
