
#include "Boundary.hpp"


void Boundary::computeNOSLIP(
		GridFunction& g, 
		const Range range, 
		const Direction dir) const
{
}
void Boundary::computeSLIP(
		GridFunction& g, 
		const Range range, 
		const Direction dir) const
{
}
void Boundary::computeINFLOW(
		GridFunction& g, 
		const Range range, 
		const Real value, 
		const Direction dir) const
{
}
void Boundary::computeOUTFLOW(
		GridFunction& g, 
		const Range  range, 
		const Direction dir) const
{
}

Boundary::Boundary(
		Range localSubInnerPRange, 
		Boundary::Competence competence)
{
	/* sadly, here we have to compute the boundary conditions ourself. 
	 * for the driven cavity, as a fallback */
	std::vector<BoundaryPiece> boundary_conditions;

	/* add pressure boundary conditions: "copy inner to boundary! */
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
					Index(localSubInnerPRange.begin.i, localSubInnerPRange.end.j) ) ));
	if(competence.Right)
	boundary_conditions.push_back(
			BoundaryPiece(Direction::Right, Condition::NOSLIP, Grid::V, 0.0, Range(
					Index(localSubInnerPRange.end.i, localSubInnerPRange.begin.j),
					Index(localSubInnerPRange.end.i, localSubInnerPRange.end.j) ) ));

	/* add u boundaries: 1 at top, 0 everywhere else */
	if(competence.Up)
	boundary_conditions.push_back(
			BoundaryPiece(Direction::Up, Condition::INFLOW, Grid::U, 1.0, Range(
					Index(localSubInnerPRange.begin.i, localSubInnerPRange.end.j),
					Index(localSubInnerPRange.end.i, localSubInnerPRange.end.j) ) ));
	if(competence.Down)
	boundary_conditions.push_back(
			BoundaryPiece(Direction::Down, Condition::NOSLIP, Grid::U, 0.0, Range(
					Index(localSubInnerPRange.begin.i, localSubInnerPRange.begin.j),
					Index(localSubInnerPRange.end.i, localSubInnerPRange.begin.j) ) ));
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

	/* *** */

	/* split boundary pieces to cell-wise portions */
	
	/* get local inner extents for p, u and v (starting with (1,1) of local Domain) */
}
Boundary::Boundary(
		Range localSubInnerPRange, 
		std::vector<BoundaryPiece> boundary_conditions)

	/* add pressure boudnaries */
{
	this->competence = Boundary::Competence(0,0,0,0);
	/* compare given boundary conditions to given inner range
	 * and keep the boudnaries for my inner of domain. */

	/* translate global ranges of boundary pieces to local ranges */

	/* *** */
	/* do from above constructor to handle boundary pieces */
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
			break;
		case Grid::V:
			break;
		case Grid::P:
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
			break;
		case Grid::G:
		case Grid::V:
			break;
		default:
			break;
	}
}

