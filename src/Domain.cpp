#include "Domain.hpp"
#include "Debug.hpp"


Domain::Domain(
		Dimension dimension, Delta delta,

		Boundary boundary,

		/* field forces */
		Real in_gx,
		Real in_gy,
		//Real in_gz,

		/* initial grid values */
		Real in_uinit,
		Real in_vinit,
		//Real in_winit,
		Real in_pinit,

		/* color for SOR Red/Black pattern */
		Color firstCellColor
	) : 
		m_dimension(dimension),
		m_boundary(boundary),
		m_FirstCellColor(firstCellColor),
		m_delta(delta),

		m_p(Dimension(dimension[0]+2,dimension[1]+2)), 
		m_p_rhs(Dimension(dimension[0]+2,dimension[1]+2)),

		m_velocities(dimension,m_boundary.getCompetence()),
		m_preliminary_velocities_FGH(dimension,m_boundary.getCompetence()),

		m_force_gx(in_gx), m_force_gy(in_gy) // , m_force_gz(in_gz)
{
	m_inner_ranges[0] = m_boundary.getInnerRanges(Boundary::Grid::U);
	m_inner_ranges[1] = m_boundary.getInnerRanges(Boundary::Grid::V);
	//m_inner_ranges[2] = m_boundary.getInnerRanges(Boundary::Grid::W);
	m_inner_ranges[3] = m_boundary.getInnerRanges(Boundary::Grid::P);
	m_whole_inner_range = m_boundary.getWholeInnerRange();

	/* init all grids to start values given from outside */
	for(uint d=0; d<DIMENSIONS; d++)
	{
		for_vecrange(i,j,m_inner_ranges[d])
		{
			if(d==0)
				u()(i,j) = in_uinit;
			if(d==1)
				v()(i, j) = in_vinit;
			//if(d==3)
			//	w()(i, j)/*,k)*/ = ;
		}
	}

	/* init pressure and rhs */
	for_vecrange(i,j,m_inner_ranges[3])
	{
		p()(i,j) = in_pinit;
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

void Domain::setVelocitiesBoundaries()
{
	this->m_boundary.setBoundary(Boundary::Grid::U, this->u());
	this->m_boundary.setBoundary(Boundary::Grid::V, this->v());
}

void Domain::setPreliminaryVelocitiesBoundaries()
{
	this->m_boundary.copyGridBoundary(Boundary::Grid::F, this->u(), this->F());
	this->m_boundary.copyGridBoundary(Boundary::Grid::G, this->v(), this->G());
}

void Domain::setPressureBoundaries()
{
	this->m_boundary.setBoundary(Boundary::Grid::P, this->p());
}

