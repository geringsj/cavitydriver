
#include "Boundary.hpp"
#include "Debug.hpp"

#include <algorithm>
#include <functional>

void Boundary::applyBoundaries(
		Grid grid, 
		GridFunction& gf,
		std::vector<Entry> boundary) const
{
	for(auto& b : boundary)
	{
		switch(b.condition)
		{
			case Condition::NOSLIP:
				computeNOSLIP(gf, b.bposition, b.iposition, b.direction, grid);
				break;
			case Condition::INFLOW:
				computeINFLOW(gf, b.bposition, b.iposition, b.direction, grid, 
					b.condition_value);
				break;
			case Condition::OUTFLOW:
				computeOUTFLOW(gf, b.bposition, b.iposition);
				break;
			case Condition::SLIP:
				computeSLIP(gf, b.bposition, b.iposition, b.direction, grid);
				break;
			default:
				break;
		}
	}
}

void Boundary::copyBoundaries(
		std::vector<Entry> boundary, 
		const GridFunction& source, 
		GridFunction& target) const
{
	for(auto& b : boundary)
	{
		target(b.bposition) = source(b.bposition);
	}
}

void Boundary::computeNOSLIP(
		GridFunction& g, 
		const Index bindex, 
		const Index iindex, 
		const Direction dir,
		const Grid grid) const
{
	if(
			(grid==Grid::U && (dir==Direction::Right || dir==Direction::Left) ) ||
			(grid==Grid::V && (dir==Direction::Up || dir==Direction::Down) ) 
	)
	{
		g(bindex) = 0.0;
	}
	else
	{
		g(bindex) = -g(iindex);
	}
}

/* TODO: this is not correct yet! */
void Boundary::computeSLIP(
		GridFunction& g, 
		const Index bindex, 
		const Index iindex, 
		const Direction dir,
		Grid grid) const
{
	if(grid==Grid::V && 1==(int)dir /* ...? */) /* flow in U direction and we want to set V */
	{
		g(bindex) = 0;
	}
	if(grid==Grid::U && 1 /* ...? */) /* flow in V direction and we want to set U */
	{
		g(iindex) = 0;
	}
	if(grid==Grid::V)
	{
	}
	if(grid==Grid::U)
	{
	}
}

void Boundary::computeINFLOW(
		GridFunction& g, 
		const Index bindex, 
		const Index iindex, 
		const Direction dir,
		const Grid grid,
		const Real value) const
{
	if(
			(grid==Grid::U && (dir==Direction::Right || dir==Direction::Left) ) ||
			(grid==Grid::V && (dir==Direction::Up || dir==Direction::Down) ) )
	{
		g(bindex) = value;
	}
	else
	{
		g(bindex) = 2.0*value - g(iindex);
	}
}

void Boundary::computeOUTFLOW(
		GridFunction& g, 
		const Index bindex, 
		const Index iindex) const
{
	/* boundary cell must be set to value of shifted cell, 
	 * so that as a result we get:
	 *    (f_b - f_i)/delta = 0
	 *    (OUTFLOW Condition)
	 */
	g(bindex) = g(iindex);
}

bool IndexIsInRange(Index index, Range range)
{
	return (
		index.i >= range.begin.i &&
		index.i <= range.end.i &&
		index.j >= range.begin.j &&
		index.j <= range.end.j
		);
}

std::vector<Range> Boundary::getInnerRanges(
		const Grid grid) const
{
	switch(grid)
	{
		case Grid::P:
			return getInnerRanges(m_inner_extent, m_boundaries_P[0]);
			break;
		case Grid::U:
		case Grid::F:
			return this->getInnerRanges(m_inner_extent, m_boundaries_U[0]);
			break;
		case Grid::V:
		case Grid::G:
			return getInnerRanges(m_inner_extent, m_boundaries_V[0]);
			break;
		default:
			return std::vector<Range>(); /* empty... */
			break;
	}
}

std::vector<Range> Boundary::getInnerRanges(
		const Range inner_extent,
		const std::vector<Entry> boundary) const
{
	std::vector<Range> ranges;
	GridFunction pattern(Dimension(inner_extent.end.i+2, inner_extent.end.j+2));

	/* mark inner cells as not reachable */
	for_range(i,j,inner_extent)
		pattern(i,j) = -1.0;
	/* mark all boundary cells as such */
	for(auto cell : boundary)
		pattern(cell.bposition.i, cell.bposition.j) = 0.0;
	/* mark inner cells next to boundary as reachable */
	for(auto cell : boundary)
	{
		int i = cell.iposition.i; int j = cell.iposition.j;
		pattern(i,j) = 1.0;
	}
	/* propagate 1s */
	bool haveReachableBegin = false;
	for(int j=inner_extent.begin.j; j<=inner_extent.end.j; j++)
	{
		haveReachableBegin = false;
		for(int i=inner_extent.begin.i; i<=inner_extent.end.i; i++)
		{
			if(pattern(i,j) > 0.0) /* at marked inner */
			{
				if(!haveReachableBegin) /* beginning after boundary */
				{
					haveReachableBegin = true;
				}
				else /* arriving at boundary */
				{
					haveReachableBegin = false;
				}
			}
			else
			if(pattern(i,j) == 0.0) /* at boundary */
			{
				haveReachableBegin = false;
			}
			else /* at blacklisted inner */
			{
				if(haveReachableBegin) /* have way to this blacklisted cell */
				{
					pattern(i,j) = 1.0;
				}
			}
		}
	}

	/* collect ranges of reachable inner cells */
	haveReachableBegin = false;
	Range lineRange;
	for(int j=inner_extent.begin.j; j<=inner_extent.end.j; j++)
	{
		for(int i=inner_extent.begin.i; i<=inner_extent.end.i; i++)
		{
			if(pattern(i,j) == 1.0) /* on reachable cell */
			{
				if(! haveReachableBegin) /* starting new range */
				{
					haveReachableBegin = true;
					lineRange.begin = Index(i,j);
					lineRange.end = Index(i,j);
				}
				else /* extend current range by one */
				{
					lineRange.end = Index(i,j);
				}
			}
			else /* non-reachable cell */
			{
				if(haveReachableBegin) /* write out and reset */
				{
					haveReachableBegin = false;
					ranges.push_back(lineRange);
				}
				else /* just go on */
				{
				}
			}
		} 
		/* line done */
		if(haveReachableBegin) /* write out and reset */
		{
			haveReachableBegin = false;
			ranges.push_back(lineRange);
		}
	}

	/* TODO: merge good Ranges to bigger cells ? */

	return ranges;
}

void Boundary::initBoundaries(
		Range localSubInnerPRange,
		std::vector<BoundaryPiece>& boundary_conditions)
{
	/* transfoprm boundaries to local Range [(1,1),(iMax,jMax)] */
	Index negOffsetTo11
		(1- localSubInnerPRange.begin.i, 1- localSubInnerPRange.begin.j);
	Index nego = negOffsetTo11;
	localSubInnerPRange.begin.i += nego.i;
	localSubInnerPRange.begin.j += nego.j;
	localSubInnerPRange.end.i += nego.i;
	localSubInnerPRange.end.j += nego.j;
	m_inner_extent = localSubInnerPRange;
	for(auto& bpiece : boundary_conditions)
	{
		/* transform boundary to coordinates of local domain */
		bpiece.range.begin.i += nego.i;
		bpiece.range.begin.j += nego.j;
		bpiece.range.end.i += nego.i;
		bpiece.range.end.j += nego.j;
	}

	/* split boundary pieces to cell-wise portions 
	 * and keep those for local domain */
	for(auto bpiece : boundary_conditions)
	for_range(i,j,bpiece.range)
	if( IndexIsInRange(Index(i,j), localSubInnerPRange) )
	{
		Boundary::Entry te(
			bpiece.direction, bpiece.condition,
			Index(i,j), Index(i,j),
			bpiece.condition_value );
		switch(bpiece.gridtype)
		{
			case Grid::P:
				this->m_boundaries_P[0].push_back(te); break;
			case Grid::U:
				this->m_boundaries_U[0].push_back(te); break;
			case Grid::V:
				this->m_boundaries_V[0].push_back(te); break;
			default: break;
		}
	}

	/* correct boundaries for U/V by shiftig right/upper boundaries one left/down.
	 * we do this because the boundaries input is given with respect to P and we
	 * need to account for the shift of U and V values on the grid */
	for(auto& bound : this->m_boundaries_U[0])
		if(bound.direction == Direction::Right)
		{
			bound.bposition.i--;
			bound.iposition.i--;
		}
	for(auto& bound : this->m_boundaries_V[0])
		if(bound.direction == Direction::Up)
		{
			bound.bposition.j--;
			bound.iposition.j--;
		}

	/* move boundary indices to be at boundary. 
	 * up to this point boundary and inner indices are at the boundary cell */
	for(auto& blist : 
		{&this->m_boundaries_U[0], &this->m_boundaries_V[0], &this->m_boundaries_P[0]} )
		for(auto& entry : *blist)
		{
			switch(entry.direction)
			{
				case Direction::Up:
					entry.bposition.j ++; break;
				case Direction::Down:
					entry.bposition.j --; break;
				case Direction::Left:
					entry.bposition.i --; break;
				case Direction::Right:
					entry.bposition.i ++; break;
				default: break;
			}
		}

		/* sort boundary cells in column major order */
	for(auto& blist : 
		{&(this->m_boundaries_U[0]),&(this->m_boundaries_V[0]),&(this->m_boundaries_P[0])})
	if(! blist->empty())
	std::sort(blist->begin(), blist->end(),
			[](const Entry& a, const Entry& b)
			{ 
				return
				(a.bposition.j < b.bposition.j) ? (1)/* a lower than b */
				: ((a.bposition.j > b.bposition.j) ? (0) /* a above b */
					: (a.bposition.i <= b.bposition.i) );/* at same line, a comes before b? */
			} );
}

Boundary::Boundary(
		Range localSubInnerPRange, 
		Boundary::Competence competence)
{
	/* sadly, here we have to compute the boundary conditions ourself. 
	 * for the driven cavity, as a fallback */
	std::vector<BoundaryPiece> boundary_conditions;

	/* add pressure boundary conditions: "copy inner to boundary" */
	if(competence.Up)
	boundary_conditions.push_back(
		BoundaryPiece(Direction::Up, Condition::OUTFLOW, Grid::P, 0.0, Range(
			Index(localSubInnerPRange.begin.i, localSubInnerPRange.end.j),
			Index(localSubInnerPRange.end.i, localSubInnerPRange.end.j) ) ));
	if(competence.Down)
	boundary_conditions.push_back(
		BoundaryPiece(Direction::Down, Condition::OUTFLOW, Grid::P, 0.0, Range(
			Index(localSubInnerPRange.begin.i, localSubInnerPRange.begin.j),
			Index(localSubInnerPRange.end.i, localSubInnerPRange.begin.j) ) ));
	if(competence.Left)
	boundary_conditions.push_back(
		BoundaryPiece(Direction::Left, Condition::OUTFLOW, Grid::P, 0.0, Range(
			Index(localSubInnerPRange.begin.i, localSubInnerPRange.begin.j),
			Index(localSubInnerPRange.begin.i, localSubInnerPRange.end.j) ) ));
	if(competence.Right)
	boundary_conditions.push_back(
		BoundaryPiece(Direction::Right, Condition::OUTFLOW, Grid::P, 0.0, Range(
			Index(localSubInnerPRange.end.i, localSubInnerPRange.begin.j),
			Index(localSubInnerPRange.end.i, localSubInnerPRange.end.j) ) ));

	/* add v boundaries: 0 everywhere */
	if(competence.Up)
	boundary_conditions.push_back(
		BoundaryPiece(Direction::Up, Condition::INFLOW, Grid::V, 0.0, Range(
			Index(localSubInnerPRange.begin.i, localSubInnerPRange.end.j),
			Index(localSubInnerPRange.end.i, localSubInnerPRange.end.j) ) ));
	if(competence.Down)
	boundary_conditions.push_back(
		BoundaryPiece(Direction::Down, Condition::NOSLIP, Grid::V, 0.0, Range(
			Index(localSubInnerPRange.begin.i, localSubInnerPRange.begin.j),
			Index(localSubInnerPRange.end.i, localSubInnerPRange.begin.j) ) ));
	if(competence.Left)
	boundary_conditions.push_back(
		BoundaryPiece(Direction::Left, Condition::NOSLIP, Grid::V, 0.0, Range(
			Index(localSubInnerPRange.begin.i, localSubInnerPRange.begin.j),
			Index(
				localSubInnerPRange.begin.i, 
				localSubInnerPRange.end.j -(competence.Up)) ) ));
	if(competence.Right)
	boundary_conditions.push_back(
		BoundaryPiece(Direction::Right, Condition::NOSLIP, Grid::V, 0.0, Range(
			Index(localSubInnerPRange.end.i, localSubInnerPRange.begin.j),
			Index(
				localSubInnerPRange.end.i, 
				localSubInnerPRange.end.j -(competence.Up)) ) ));

	/* add u boundaries: 1 at top, 0 everywhere else */
	if(competence.Up)
	boundary_conditions.push_back(
		BoundaryPiece(Direction::Up, Condition::INFLOW, Grid::U, 1.0, Range(
			Index(localSubInnerPRange.begin.i, localSubInnerPRange.end.j),
			Index(
				localSubInnerPRange.end.i -(competence.Right), 
				localSubInnerPRange.end.j) ) ));
	if(competence.Down)
	boundary_conditions.push_back(
		BoundaryPiece(Direction::Down, Condition::NOSLIP, Grid::U, 0.0, Range(
			Index(localSubInnerPRange.begin.i, localSubInnerPRange.begin.j),
			Index(
				localSubInnerPRange.end.i -(competence.Right), 
				localSubInnerPRange.begin.j) ) ));
	if(competence.Left)
	boundary_conditions.push_back(
		BoundaryPiece(Direction::Left, Condition::NOSLIP, Grid::U, 0.0, Range(
			Index(localSubInnerPRange.begin.i, localSubInnerPRange.begin.j),
			Index(localSubInnerPRange.begin.i, localSubInnerPRange.end.j) ) ));
	if(competence.Right)
	boundary_conditions.push_back(
		BoundaryPiece(Direction::Right, Condition::NOSLIP, Grid::U, 0.0, Range(
			Index(localSubInnerPRange.end.i, localSubInnerPRange.begin.j),
			Index(localSubInnerPRange.end.i, localSubInnerPRange.end.j) ) ));
	this->m_competence = competence;

	/* *** */
	initBoundaries(localSubInnerPRange, boundary_conditions);

	/* TODO: sort by boundary type ? */
}

Boundary::Boundary(
		Range localSubInnerPRange, 
		std::vector<BoundaryPiece> boundary_conditions)
{
	this->m_competence = Boundary::Competence(0,0,0,0);

	initBoundaries(localSubInnerPRange, boundary_conditions);

	/* TODO: sort by boundary type ? */
}

Boundary::~Boundary()
{
}

void Boundary::setBoundary(
		Grid grid, 
		GridFunction& gf) const
{
	switch(grid)
	{
		case Grid::U:
			applyBoundaries(grid, gf, m_boundaries_U[0]);
			break;
		case Grid::V:
			applyBoundaries(grid, gf, m_boundaries_V[0]);
			break;
		case Grid::P:
			applyBoundaries(grid, gf, m_boundaries_P[0]);
			break;
		default:
			break;
	}
}

/* for F and G */
void Boundary::copyGridBoundary(
		Grid grid, 
		const GridFunction& source, 
		GridFunction& target) const
{
	switch(grid)
	{
		case Grid::F:
		case Grid::U:
			copyBoundaries(m_boundaries_U[0], source, target);
			break;
		case Grid::G:
		case Grid::V:
			copyBoundaries(m_boundaries_V[0], source, target);
			break;
		default:
			break;
	}
}

