
#ifndef Computation_hpp
#define Computation_hpp

#include "Domain.hpp"

/** 
 * Implements functions for computing several steps of the simulation algorithm.
 *
 * @author becherml, friesfn, geringsj
 * @date 11/2014
 *
 * Since the Computation only operates on input values and has no need 
 * for private members, it is a namespace.
 *
 * All functions expect a Domain to operate on and some values 
 * from the simulation parameters. 
 *
 * The needed Grids to operate/compute on are retrieved from the domain 
 * dimension-wise using an easy interface, which allows the implementation 
 * of the Computation functions to be short and simple.
 *
 * Internally, for computation of derivatives, functions from the 
 * Derivatives namespace are used.
 */
namespace Computation
{

	/** 
	 * Compute the timestep deltaT that is to be used in further computations.
	 *
	 * Implements chapter "3.2.5 Time Step - Stability Conditions" formulas from 
	 * the script (page 25). 
	 *
	 * Needed grid width dx, height dy and uMax/vMax are retrieved from Domain.
	 */
	Real computeTimestep(
			Domain& domain, 
			/**< Domain where to get dy, dx, uMax, vMax from. */
			const Real tau, 
			/**< Some safety factor from simulation parameters SimParams. */
			const Real Re 
			/**< The Reynolds number from simulation parameters SimParams. */
			); 


	/** 
	 * Compute the "preliminary velocities" F, G (, H), 
	 * which result from the Momentum Equations given U and V.
	 *
	 * Implements chapter "3.2.1 Momentum Equations" formulas from the script (page 19f).
	 * U and V (and W) are retrieved from the given Domain.
	 * F and G ( and H) are stored in the given Domain.
	 */
	void computePreliminaryVelocities(
			Domain& domain, 
			/**< Domain where to get U, V, (and W) and where to store F and G (and H). */
			const Real deltaT, 
			/**< Time step as computed by Computation::computeTimestep. */
			const Real Re, 
			/**< The Reynolds number from simulation parameters SimParams. */
			const Real alpha 
			/**< Factor for weighted mean of 'original central difference 
			 * and the Donor-Cell Scheme'. See script chapter "3.2.4", page 22ff. 
			 * This value is stored in SimParams. */
			);


	/** 
	 * Computes the right-hand side of the discrete Poisson equation for p^(n+1).
	 * 
	 * Implements Equation (3.3) on page 21 of the script 
	 * (chapter "3.2.2 Continuity Equation".
	 * The right-hand side RHS is computed as sum of backward differences of F and G,
	 * divided by deltaT.
	 * RHS is stored in the given Domain.
	 */
	void computeRighthandSide(
			Domain& domain,
			/**< Where F and G are read from, where right-hand side RHS is stored. */
			const Real deltaT
			/**< Time step as computed by Computation::computeTimestep. */
			);
	

	/** 
	 * Computes new velocities U^(n+1) and V^(n+1) (and W^(n+1)).
	 *
	 * Implements equation (3.1) from the script, chapter 
	 * "3.2.1 Momentum Equations" (page 20).
	 *
	 * New velocities are computed from old velocities U^(n), V^(n) (, W^(n)) 
	 * and preliminary velocities F^(n), G^(n) (, H^(n))
	 * using pressure p^(n+1) and timestep deltaT.
	 */
	void computeNewVelocities(
			Domain& domain,
			/**< Domain where to get F, G, (H) and pressure P, 
			 * also where to store new velocities U, V (, W). */
			const Real deltaT
			/**< Time step as computed by Computation::computeTimestep. */
			);
};

#endif
