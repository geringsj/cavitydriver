
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
		NOSLIP=0, INFLOW=1, OUTFLOW=2,SLIP=3
	};
	enum class Grid {
		U, V, P, F, G, T /* , W, H */
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
	struct Entry {
		Direction direction;
		Condition condition;
		Index bposition;
		Index iposition;
		Real condition_value;
		Entry(
			Direction direction, Condition condition,
			Index bposition, Index iposition, Real condition_value) :
			direction(direction), condition(condition),
			bposition(bposition), iposition(iposition), 
			condition_value(condition_value) {}
	};

	Competence m_competence;
	Range m_inner_extent;

	std::vector<Entry> m_boundaries_U[4];
	std::vector<Entry> m_boundaries_V[4];
	std::vector<Entry> m_boundaries_P[4];

	void computeNOSLIP(GridFunction& g, const Index bindex, const Index iindex, const Direction dir, const Grid grid) const;
	void computeSLIP(GridFunction& g, const Index bindex, const Index iindex, const Direction dir, const Grid grid) const;
	void computeINFLOW(GridFunction& g, const Index bindex, const Index iindex, const Direction dir, const Grid grid, const Real value) const;
	void computeOUTFLOW(GridFunction& g, const Index bindex, const Index iindex) const;

	std::vector<Range> getInnerRanges(
		const Range inner_extent,
		const std::vector<Entry>* boundary) const;

	void initBoundaries(
		Range localSubInnerPRange, 
		std::vector<BoundaryPiece>& boundary_conditions);
	void copyBoundaries(
		const std::vector<Entry>& boundary,
		const GridFunction& source, GridFunction& target) const;
	void workBoundaries(const std::vector<Entry>* boundaries, const Grid grid, GridFunction& gf) const;
	void splitBoundaryTypes();

public:
	Boundary(
		Range localSubInnerPRange, 
		Boundary::Competence competence,
		std::vector<BoundaryPiece>& boundary_conditions);

	~Boundary();

	Boundary::Competence getCompetence() const { return m_competence; }

	void setBoundary(const Grid grid, GridFunction& gf) const;

	/* for F and G */
	void copyGridBoundary(
			const Grid grid, const GridFunction& source, GridFunction& target) const;

	std::vector<Range> getInnerRanges(const Grid grid) const;
	Range getWholeInnerRange() const { return m_inner_extent; }
};


#endif

