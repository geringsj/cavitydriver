
#include "Boundary.hpp"
#include "Debug.hpp"

#include <algorithm>
#include <functional>
#include <list>


void Boundary::copyBoundaries(
		const std::vector<Entry> boundary, 
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

void Boundary::computeSLIP(
		GridFunction& g,
		const Index bindex,
		const Index iindex, 
		const Direction dir,
		const Grid grid) const
{
	/* the dot-product condition of SLIP: (U,V) dot normalDirection = 0 */
	if(grid==Grid::V && (dir==Direction::Up || dir==Direction::Down)) 
	{
		g(bindex) = 0;
		return;
	}
	if(grid==Grid::U && (dir==Direction::Right || dir==Direction::Left))
	{
		g(bindex) = 0;
		return;
	}
	/* derivative in normal direction must be zero => boundary == inner cell */
	g(bindex) = g(iindex);
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

static bool IndexIsInRange(Index index, Range range)
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
			return getInnerRanges(m_inner_extent, m_boundaries_P);
			break;
		case Grid::U:
		case Grid::F:
			return getInnerRanges(m_inner_extent, m_boundaries_U);
			break;
		case Grid::V:
		case Grid::G:
			return getInnerRanges(m_inner_extent, m_boundaries_V);
			break;
		default:
			return std::vector<Range>(); /* empty... */
			break;
	}
}

std::vector<Range> Boundary::getInnerRanges(
		const Range inner_extent,
		const std::vector<Entry> boundary[]) const
{
	Dimension pattern_dim(inner_extent.end.i+2, inner_extent.end.j+2);
	GridFunction pattern(pattern_dim);
	std::list<Index> cells;

	/* mark all cells as not reachable */
	for_range(i,j,Range(Index(0,0), Index(pattern_dim.i-1,pattern_dim.j-1)))
		pattern(i,j) = -1.0;
	/* mark all boundary and neighbour inner cells as such */
	for(int i : {0,1,2,3} )
	for(auto cell : boundary[i])
	{
		pattern(cell.iposition) = 1.0;
		cells.push_back(cell.iposition);
	}
	for(int i : {0,1,2,3} )
	for(auto cell : boundary[i])
		pattern(cell.bposition) = 0.0;

	/* propagate 1s */
	while(! cells.empty())
	{
		Index ci = cells.front(); cells.pop_front();
		for(int i : {1, -1})
		for(int j : {1, -1})
		{
			Index ni(ci.i+i, ci.j+j);
			if( pattern(ni) == -1.0 && IndexIsInRange(ni,inner_extent)) 
			/* cell reachable from current inner cell, 
			 * not outside domain and not marked yet: add to list */
			{
				pattern(ni) = 1.0;
				cells.push_back(ni);
			}
		}
	}

	/* collect ranges of reachable inner cells */
	std::list<Range> ranges;
	bool haveReachableBegin = false;
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

	/* merge fitting neighbour ranges to bigger ones */
	std::vector<Range> vranges;
	if(! ranges.empty())
	{
		vranges.push_back( ranges.front() );
		ranges.pop_front();
	}
	while(! ranges.empty())
	{
		if(
			vranges.back().begin.i == ranges.front().begin.i &&
			vranges.back().end.i == ranges.front().end.i &&
			vranges.back().begin.j < ranges.front().begin.j &&
			vranges.back().end.j+1 == ranges.front().end.j
		  )
			/* merge ranges */
		{
			vranges.back().end.j++;
			ranges.pop_front();
		}
		else /* start new range to maybe append lines on */
		{
			vranges.push_back( ranges.front() );
			ranges.pop_front();
		}
	}

	return vranges;
}

void Boundary::initBoundaries(
		Range localSubInnerPRange,
		std::vector<BoundaryPiece>& boundary_conditions)
{
	/* transfoprm boundaries to local Range [(1,1),(iMax,jMax)] */
	Index negOffsetTo11(1-localSubInnerPRange.begin.i,1-localSubInnerPRange.begin.j);
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
			{&this->m_boundaries_U[0], &this->m_boundaries_V[0], &this->m_boundaries_P[0]} 
		)
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
		Boundary::Competence competence,
		std::vector<BoundaryPiece>& boundary_conditions)
{
	/* sadly, here we have to compute the boundary conditions ourself. 
	 * for the driven cavity, as a fallback */
	// std::vector<BoundaryPiece> boundary_conditions;

	this->m_competence = competence;

	if(boundary_conditions.size() == 0)
	{
		/* get missing boundary conditions for driven cavity */
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
					localSubInnerPRange.end.j /*-(competence.Up)*/ ) ) ));
		if(competence.Right)
		boundary_conditions.push_back(
			BoundaryPiece(Direction::Right, Condition::NOSLIP, Grid::V, 0.0, Range(
				Index(localSubInnerPRange.end.i, localSubInnerPRange.begin.j),
				Index(
					localSubInnerPRange.end.i, 
					localSubInnerPRange.end.j /*-(competence.Up)*/ ) ) ));

		/* add u boundaries: 1 at top, 0 everywhere else */
		if(competence.Up)
		boundary_conditions.push_back(
			BoundaryPiece(Direction::Up, Condition::INFLOW, Grid::U, 1.0, Range(
				Index(localSubInnerPRange.begin.i, localSubInnerPRange.end.j),
				Index(
					localSubInnerPRange.end.i /*-(competence.Right)*/,
					localSubInnerPRange.end.j) ) ));
		if(competence.Down)
		boundary_conditions.push_back(
			BoundaryPiece(Direction::Down, Condition::NOSLIP, Grid::U, 0.0, Range(
				Index(localSubInnerPRange.begin.i, localSubInnerPRange.begin.j),
				Index(
					localSubInnerPRange.end.i /*-(competence.Right)*/,
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
	}

	initBoundaries(localSubInnerPRange, boundary_conditions);
	splitBoundaryTypes();
}

void Boundary::splitBoundaryTypes()
{
	for(auto& bnds : {m_boundaries_U, m_boundaries_V, m_boundaries_P} )
	{
	for(auto e : bnds[0])
		for(auto cond : {Condition::INFLOW, Condition::OUTFLOW, Condition::SLIP} )
				if(e.condition == cond)
				bnds[static_cast<int>(cond)].push_back(e);

		bnds[0].erase( 
				std::remove_if(bnds[0].begin(), bnds[0].end(), 
					[](Boundary::Entry e)
					{
						return (e.condition != Condition::NOSLIP);
					} ),
				bnds[0].end() );
	}
}

Boundary::~Boundary()
{
}

void Boundary::setBoundary(
		const Grid grid, 
		GridFunction& gf) const
{
	const std::vector<Entry>* boundaries;
	switch(grid)
	{
		case Grid::U:
				boundaries = &m_boundaries_U[0];
				break;
		case Grid::V:
				boundaries = &m_boundaries_V[0];
				break;
			break;
		case Grid::P:
				boundaries = &m_boundaries_P[0];
			break;
		default:
				boundaries = NULL;
			break;
	}

	for(auto& b : boundaries[static_cast<int>(Condition::NOSLIP)])
				computeNOSLIP(gf, b.bposition, b.iposition, b.direction, grid);

	for(auto& b : boundaries[static_cast<int>(Condition::INFLOW)])
				computeINFLOW(gf, b.bposition, b.iposition, b.direction, grid, 
					b.condition_value);

	for(auto& b : boundaries[static_cast<int>(Condition::OUTFLOW)])
				computeOUTFLOW(gf, b.bposition, b.iposition);

	for(auto& b : boundaries[static_cast<int>(Condition::SLIP)])
				computeSLIP(gf, b.bposition, b.iposition, b.direction, grid);
}

/* for F and G */
void Boundary::copyGridBoundary(
		const Grid grid, 
		const GridFunction& source, 
		GridFunction& target) const
{
	switch(grid)
	{
		case Grid::F:
		case Grid::U:
			copyBoundaries(m_boundaries_U[0], source, target);
			copyBoundaries(m_boundaries_U[1], source, target);
			copyBoundaries(m_boundaries_U[2], source, target);
			copyBoundaries(m_boundaries_U[3], source, target);
			break;
		case Grid::G:
		case Grid::V:
			copyBoundaries(m_boundaries_V[0], source, target);
			copyBoundaries(m_boundaries_V[1], source, target);
			copyBoundaries(m_boundaries_V[2], source, target);
			copyBoundaries(m_boundaries_V[3], source, target);
			break;
		default:
			break;
	}
}

