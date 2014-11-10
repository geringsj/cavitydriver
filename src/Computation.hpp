//! This class implements the computation
/*!
 * @author ...
 * @date 2014
 */

#ifndef Computation_HPP_
#define Computation_HPP_

// #include "Structs.hpp"
#include "Domain.hpp"

namespace Computation
{
	Real computeTimestep(
		Domain& domain,
		//const Real uMax, const Real vMax, const Point h, 
		const Real tau, const Real Re);

	//void setBoundaryX(GridFunction& x);

	void computeMomentumEquationsFGH(
		Domain& domain,
		//GridFunction& f, GridFunction& g,
		//GridFunction& u, GridFunction& v, 
		//GridFunction& gx, GridFunction& gy, const Point delta, 
		const Real deltaT, const Real Re);

	void computeRighthandSide(
		Domain& domain,
		//GridFunction& rhs, GridFunction& f, GridFunction& g, const Point delta, 
		const Real deltaT);

	/* solve for p in SOR here */

	void computeNewVelocities(
		Domain& domain,
		//GridFunction& u, GridFunction& v,
		//GridFunction& f, GridFunction& g,
		//GridFunction& p, const Point delta, 
		const Real deltaT);
};

#endif
