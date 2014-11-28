
/** 
 * We compute the residual and SOR cycle here.
 * 
 * @author becherml, friesfn, geringsj
 * @date 11/2014
 *
 */


#ifndef Solver_hpp
#define Solver_hpp

#include "Domain.hpp"
#include "GridFunction.hpp"

namespace Solver
{


	/**
	* Compute the residual res that is to be used as an exit condition.
	*
	* @return Computed residual res.
	*
	* Implements chapter "3.2.5 Time Step - Stability Conditions" formulas from
	* the script (page 25).
	*/
	Real computeResidual(
			GridFunction& p,
			GridFunction& rhs,
			const Point& delta,
			Dimension inner_begin, Dimension inner_end);
	
	/**
	* Compute one SOR iteration on the pressure p.
	*
	* Implements formula (4.1) from chapter "4.2 Improving Gauss-Seidel - SOR" from the script (page 28).
	*/
	void SORCycle(
			GridFunction& p,
			GridFunction& rhs,
			const Point& delta,
			Dimension inner_begin, Dimension inner_end,
			Real omega);

	/**
	* Compute one SOR iterator on every 'red' cell of the pressure field
	*/
	void SORCycleRedBlack(
			Domain& domain,
			const Point& delta,
			Real omega,
			Color color);

}

#endif
