#ifndef Solver_hpp
#define Solver_hpp

#include "typedef.hpp"

namespace Solver
{

	RealType computeResidual( GridFunctionType& sourcegridfunction, GridFunctionType& rhs, const PointType& h );
	
	void SORCycle( GridFunction* gridfunction, GridFunctionType& rhs, const PointType& delta, RealType omega );
	
}

#endif