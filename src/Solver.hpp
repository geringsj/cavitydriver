
/** 
 * We compute the residual and SOR cycle here.
 * 
 * @author becherml, friesfn, geringsj
 * @date 11/2014
 *
 */


#ifndef Solver_hpp
#define Solver_hpp

#include "GridFunction.hpp"

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
