#ifndef Solver_hpp
#define Solver_hpp

#include "GridFunction.hpp"
#include "typedef.hpp"

#include <math.h>

namespace Solver
{

	RealType computeResidual(GridFunctionType& sourcegridfunction,
											GridFunctionType& rhs,
											const PointType& h,
											int iMax, int jMax);
	
	void SORCycle(GridFunction* gridfunction,
							GridFunctionType& rhs,
							const PointType& delta,
							RealType omega);
	
}

#endif
