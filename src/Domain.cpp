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
		Real in_tinit,

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

		m_force_gx(in_gx), m_force_gy(in_gy), // m_force_gz(in_gz),

		m_temperature(Dimension(dimension[0]+2,dimension[0]+2))
{
	m_inner_ranges[0] = m_boundary.getInnerRanges(Boundary::Grid::U);
	m_inner_ranges[1] = m_boundary.getInnerRanges(Boundary::Grid::V);
	//m_inner_ranges[2] = m_boundary.getInnerRanges(Boundary::Grid::W);
	m_inner_ranges[3] = m_boundary.getInnerRanges(Boundary::Grid::P);
	m_inner_ranges[4] = m_boundary.getInnerRanges(Boundary::Grid::T);
	m_whole_inner_range = m_boundary.getWholeInnerRange();

	/* init all grids to start values given from outside */
	for(uint d=0; d<DIMENSIONS; d++)
	{
		for_vecrange(i,j,m_inner_ranges[d])
		{
			if(d==0)
				u()(i,j) = in_uinit;
			if(d==1)
				v()(i,j) = in_vinit;
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

	/* init temperature */
	for_vecrange(i,j,m_inner_ranges[4])
		t()(i,j) = in_tinit;

	/* and don't forget the boundaries */
	this->setVelocitiesBoundaries();
	this->setPreliminaryVelocitiesBoundaries();
	this->setPressureBoundaries();
	this->setTemperatureBoundaries();
}

Domain::~Domain()
{
}

void Domain::setVelocitiesBoundaries()
{
	this->m_boundary.setBoundary(Boundary::Grid::U, this->u());
	this->m_boundary.setBoundary(Boundary::Grid::V, this->v());

//	/* U */
//	Range RangeU = m_inner_ranges[0];//[0];
//	Range upRangeU(Index(RangeU.begin.i, RangeU.end.j), RangeU.end);
//	Range downRangeU(RangeU.begin, Index(RangeU.end.i,RangeU.begin.j));
//	Range leftRangeU(RangeU.begin, Index(RangeU.begin.i,RangeU.end.j));
//	Range rightRangeU(Index(RangeU.end.i,RangeU.begin.j), RangeU.end);
//	/* Up */
//	for_range(i,j,upRangeU)
//	{
//		Index bindex(i,j+1);
//		Index iindex(i,j);
//		/* NOSLIP */
//		u()( bindex ) = - u()( iindex );
//	}
//	/* Down */
//	for_range(i,j,downRangeU)
//	{
//		Index bindex(i,j-1);
//		Index iindex(i,j);
//		/* NOSLIP */
//		u()( bindex ) = - u()( iindex );
//	}
//	/* Left */
//	for_range(i,j,leftRangeU)
//	{
//		Index bindex(i-1,j);
//		Index iindex(i,j);
//		/* OUTFLOW */
//		u()( bindex ) = u()( iindex );
//	}
//	/* Right */
//	for_range(i,j,rightRangeU)
//	{
//		Index bindex(i+1,j);
//		Index iindex(i,j);
//		/* OUTFLOW */
//		u()( bindex ) = u()( iindex );
//	}
//
//	/* V */
//	Range RangeV = m_inner_ranges[1];//[0];
//	Range upRangeV(Index(RangeV.begin.i, RangeV.end.j), RangeV.end);
//	Range downRangeV(RangeV.begin, Index(RangeV.end.i,RangeV.begin.j));
//	Range leftRangeV(RangeV.begin, Index(RangeV.begin.i,RangeV.end.j));
//	Range rightRangeV(Index(RangeV.end.i,RangeV.begin.j), RangeV.end);
//	/* Up */
//	for_range(i,j,upRangeV)
//	{
//		Index bindex(i,j+1);
//		Index iindex(i,j);
//		/* NOSLIP */
//		v()( bindex ) = 0.0;//- v()( iindex );
//	}
//	/* Down */
//	for_range(i,j,downRangeV)
//	{
//		Index bindex(i,j-1);
//		Index iindex(i,j);
//		/* NOSLIP */
//		v()( bindex ) = 0.0;//- v()( iindex );
//	}
//	/* Left */
//	for_range(i,j,leftRangeV)
//	{
//		Index bindex(i-1,j);
//		Index iindex(i,j);
//		/* OUTFLOW */
//		v()( bindex ) = v()( iindex );
//	}
//	/* Right */
//	for_range(i,j,rightRangeV)
//	{
//		Index bindex(i+1,j);
//		Index iindex(i,j);
//		/* OUTFLOW */
//		v()( bindex ) = v()( iindex );
//	}
}

void Domain::setPreliminaryVelocitiesBoundaries()
{
	this->m_boundary.copyGridBoundary(Boundary::Grid::F, this->u(), this->F());
	this->m_boundary.copyGridBoundary(Boundary::Grid::G, this->v(), this->G());

//	/* U / F */
//	Range RangeU = m_inner_ranges[0];//[0];
//	Range upRangeU(Index(RangeU.begin.i, RangeU.end.j), RangeU.end);
//	Range downRangeU(RangeU.begin, Index(RangeU.end.i,RangeU.begin.j));
//	Range leftRangeU(RangeU.begin, Index(RangeU.begin.i,RangeU.end.j));
//	Range rightRangeU(Index(RangeU.end.i,RangeU.begin.j), RangeU.end);
//	/* Up */
//	for_range(i,j,upRangeU)
//	{
//		Index bindex(i,j+1);
//		F()( bindex ) = u()( bindex );
//	}
//	/* Down */
//	for_range(i,j,downRangeU)
//	{
//		Index bindex(i,j-1);
//		F()( bindex ) = u()( bindex );
//	}
//	/* Left */
//	for_range(i,j,leftRangeU)
//	{
//		Index bindex(i-1,j);
//		F()( bindex ) = u()( bindex );
//	}
//	/* Right */
//	for_range(i,j,rightRangeU)
//	{
//		Index bindex(i+1,j);
//		F()( bindex ) = u()( bindex );
//	}
//
//	/* V */
//	Range RangeV = m_inner_ranges[1];//[0];
//	Range upRangeV(Index(RangeV.begin.i, RangeV.end.j), RangeV.end);
//	Range downRangeV(RangeV.begin, Index(RangeV.end.i,RangeV.begin.j));
//	Range leftRangeV(RangeV.begin, Index(RangeV.begin.i,RangeV.end.j));
//	Range rightRangeV(Index(RangeV.end.i,RangeV.begin.j), RangeV.end);
//	/* Up */
//	for_range(i,j,upRangeV)
//	{
//		Index bindex(i,j+1);
//		G()( bindex ) = v()( bindex );
//	}
//	/* Down */
//	for_range(i,j,downRangeV)
//	{
//		Index bindex(i,j-1);
//		G()( bindex ) = v()( bindex );
//	}
//	/* Left */
//	for_range(i,j,leftRangeV)
//	{
//		Index bindex(i-1,j);
//		G()( bindex ) = v()( bindex );
//	}
//	/* Right */
//	for_range(i,j,rightRangeV)
//	{
//		Index bindex(i+1,j);
//		G()( bindex ) = v()( bindex );
//	}
}

void Domain::setPressureBoundaries()
{
	this->m_boundary.setBoundary(Boundary::Grid::P, this->p());

//	Real inflowval = 0.1;
//
//	Range RangeP = m_inner_ranges[3];//[0];
//	Range upRangeP(Index(RangeP.begin.i, RangeP.end.j), RangeP.end);
//	Range downRangeP(RangeP.begin, Index(RangeP.end.i,RangeP.begin.j));
//	Range leftRangeP(RangeP.begin, Index(RangeP.begin.i,RangeP.end.j));
//	Range rightRangeP(Index(RangeP.end.i,RangeP.begin.j), RangeP.end);
//	/* Up */
//	for_range(i,j,upRangeP)
//	{
//		Index bindex(i,j+1);
//		Index iindex(i,j);
//		/* OUTFLOW */
//		p()( bindex ) = p()( iindex );
//	}
//	/* Down */
//	for_range(i,j,downRangeP)
//	{
//		Index bindex(i,j-1);
//		Index iindex(i,j);
//		/* OUTFLOW */
//		p()( bindex ) = p()( iindex );
//	}
//	/* Left */
//	for_range(i,j,leftRangeP)
//	{
//		Index bindex(i-1,j);
//		Index iindex(i,j);
//		/* INFLOW */
//		p()( bindex ) = 2.0*inflowval - p()( iindex );
//	}
//	/* Right */
//	for_range(i,j,rightRangeP)
//	{
//		Index bindex(i+1,j);
//		Index iindex(i,j);
//		/* NOSLIP */
//		p()( bindex ) = - p()( iindex );
//	}
}

void Domain::setTemperatureBoundaries()
{
	this->m_boundary.setBoundary(Boundary::Grid::T, this->t());
}
