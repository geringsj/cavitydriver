//! We compute the residual and SOR Cycle here
/*!
 * @author becherml, friesfn, geringsj
 * @date 2014
 */

#ifndef Solver_hpp
#define Solver_hpp

#include "GridFunction.hpp"

/* TODO: document !! */

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
