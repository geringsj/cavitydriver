
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
			const GridFunction& p,
			const GridFunction& rhs,
			const Point delta,
			const Ranges& inner_range,
			const Dimension globalDims);

	/**
	* Compute the squared residual res that is to be used as an exit condition.
	*
	* @return Computed squared residual res.
	*
	* Implements chapter "3.2.5 Time Step - Stability Conditions" formulas from
	* the script (page 25).
	*/
	Real computeSquaredResidual(
			const GridFunction& p,
			const GridFunction& rhs,
			const Point& delta,
			const Ranges& inner_range,
			const Dimension globalDims);
	
	/**
	* Compute one SOR iteration on the pressure p.
	*
	* Implements formula (4.1) from chapter "4.2 Improving Gauss-Seidel - SOR" from the script (page 28).
	*/
	void SORCycle(
			GridFunction& p,
			const GridFunction& rhs,
			const Point& delta,
			const Ranges& inner_range,
			const Real& omega);

	/**
	* Compute one SOR iterator on every 'Red' or 'Black' cell of the pressure field.
	*/
	void SORCycleRedBlack(
			Domain& domain,
			const Real& omega,
			const Color& color);

}

#endif
