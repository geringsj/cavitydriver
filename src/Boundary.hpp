
#ifndef Boundary_hpp
#define Boundary_hpp

#include "Structs.hpp"
#include "GridFunction.hpp"

#include <vector>

class Boundary {
public:
	struct Competence {
		/* to be used as booleans ! */
		int Up;
		int Down;
		int Left;
		int Right;

		Competence() : Up(1), Down(1), Left(1), Right(1) {}
		Competence(int u, int d, int l, int r) : Up(u), Down(d), Left(l), Right(r) {}

		bool operator[](int i)
		{
			switch(i)
			{
				case 0: /* U */
					return Right; break;
				case 1: /* V */
					return Up; break;
				case 2: /* W */
					return 0; break; /* TODO Front/Back */
				default: 
					return 0;  break;
			}
		}
	};

	enum class Direction {
		Up, Down, Left, Right
	};
	enum class Condition {
		NOSLIP=0, INFLOW=1, SLIP=2, OUTFLOW=3
	};
	enum class Grid {
		U, V, P, F, G /* , W, H */
	};

	struct BoundaryPiece {
		Direction direction;
		Condition condition;
		Grid gridtype;
		Real condition_value;
		Range range;

		BoundaryPiece(
			Direction direction, Condition condition,
			Grid gridtype, Real condition_value, Range range) :
			direction(direction), condition(condition),
			gridtype(gridtype), condition_value(condition_value), range(range) {}
	};

private:
	Competence competence;

	struct Entry {
		Direction direction;
		Condition condition;
		Index position;
		Real condition_value;
	};
	std::vector<Entry> m_boundaries_U;//[4];
	std::vector<Entry> m_boundaries_V;//[4];
	std::vector<Entry> m_boundaries_P;//[4];

	void computeNOSLIP(GridFunction& g, const Range range, const Direction dir) const;
	void computeSLIP(GridFunction& g, const Range range, const Direction dir) const;
	void computeINFLOW(
			GridFunction& g, const Range range, 
			const Real value, const Direction dir) const;
	void computeOUTFLOW(GridFunction& g, const Range  range, const Direction dir) const;

public:
	Boundary(Range localSubInnerPRange, Boundary::Competence competence=Competence());
	Boundary(Range localSubInnerPRange, std::vector<BoundaryPiece> boundary_conditions);
	~Boundary();

	Boundary::Competence getCompetence(){ return competence; }

	void setBoundary(Grid grid, GridFunction& gf) const;

	/* for F and G */
	void copyGridBoundary(Grid grid, const GridFunction& source, GridFunction& target) const;
};


#endif

