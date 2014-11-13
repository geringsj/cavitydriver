//! This namespace implements the computation
/*!
 * @author becherml, friesfn, geringsj
 * @date 2014
 */

#ifndef Computation_hpp
#define Computation_hpp

#include "Domain.hpp"

/* TODO: document this !! */

namespace Computation
{
	Real computeTimestep(
			Domain& domain,
			const Real tau, const Real Re);
	
	void computeMomentumEquationsFGH(
			Domain& domain,
			const Real deltaT, const Real Re);
	
	void computeRighthandSide(
			Domain& domain,
			const Real deltaT);
	
	void computeNewVelocities(
			Domain& domain,
			const Real deltaT);
};

#endif
