#ifndef Solver_hpp
#define Solver_hpp

#include "GridFunction.hpp"
#include "Structs.hpp"

#include <math.h>

namespace Solver
{

	Real computeResidual(
			GridFunction& sourcegridfunction,
			GridFunction& rhs,
			const Point& h,
			int iMax, int jMax);
	
	void SORCycle(
			GridFunction& gridfunction,
			GridFunction& rhs,
			const Point& delta,
			Real omega);
	
}

#endif
