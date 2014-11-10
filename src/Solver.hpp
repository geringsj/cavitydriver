#ifndef Solver_hpp
#define Solver_hpp

#include "GridFunction.hpp"
#include "Structs.hpp"

namespace Solver
{

	Real computeResidual(
			GridFunction& p,
			GridFunction& rhs,
			const Point& delta,
			Dimension inner_begin, Dimension inner_end);
	
	void SORCycle(
			GridFunction& p,
			GridFunction& rhs,
			const Point& delta,
			Dimension inner_begin, Dimension inner_end,
			Real omega);
}

#endif
