//! This class implements the computation
/*!
 * @author ...
 * @date 2014
 */

#include "typedef.hpp"

#ifndef Computation_HPP_
#define Computation_HPP_



class Computation
{
public:
	/** ctor */
	Computation();
	/** dtor */
	~Computation();
	RealType computeTimestep(RealType uMax, RealType vMax, RealType tau, 
		const PointType& h, RealType Re);
	void computeNewVelocities(GridFunction* u, GridFunction* v,
		GridFunctionType& f, GridFunctionType& g,
		GridFunctionType& p, const PointType& delta, RealType deltaT);
	void computeMomentumEquations(GridFunctionType& f, GridFunctionType& g,
		GridFunction* u, GridFunction* v, GridFunctionType& gx,
		GridFunctionType& gy, const PointType& delta, RealType deltaT);
	void setBoundaryX(GridFunction x);
	void computeRighthandSide(GridFunction* rhs, GridFunctionType& f,
		GridFunctionType& g, const PointType& delta, RealType deltaT);
};

#endif