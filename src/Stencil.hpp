//! This is still kind of a "Stencil", but now it is on steroids
/*!
 * @author becherml, friesfn, geringsj
 * @date 2014
 */

#ifndef Derivatives_hpp
#define Derivatives_hpp 

#include "GridFunction.hpp"
#include <functional>

namespace Derivatives
{
	/* those enums are expected when demanding a derivative operator, 
	 * so the functions know what to do (to get the dimensions right and stuff) */
	enum class Function {
		U = 1,
		V = 2,
		W = 3,
		P = 4
	};
	typedef Function F;

	enum class Direction {
		xb=-1,
		yb=-2,
		zb=-3,

		xf=1, 
		yf=2, 
		zf=3, 

		xx=4, 
		yy=5, 
		zz=6,
	};
	typedef Direction d;

	/* functions to get first derivatives (_f_orward and _b_ackward) 
	 * and second derivatives */
	std::function<Real(int, int)> getDerivative(GridFunction& gf, Delta d, Direction df);
	std::function<Real(int, int)> getDerivative(GridFunction& gf, Delta d, int df);

	/* first derivative of functions product 
	 * first function f1 says where to evaluate: exactly at grid point on f1(i,j)
	 * second function, f2, determines direction of derivative
	 * example: f2==V  => derive in y direction */
	std::function<Real(int, int)> getProductFirstDerivative(
			GridFunction& gf1, GridFunction& gf2, Delta d, 
			Function f1, Function f2);
	std::function<Real(int, int)> getProductFirstDerivative(
			GridFunction& gf1, GridFunction& gf2, Delta d, 
			int f1, int f2);

	/* donor cell scheme functions */
	/* TODO: understand DCS and implement functions */
	std::function<Real(int, int)> getProductFirstDerivativeDCS(
			GridFunction& gf1, GridFunction& gf2, Delta d, 
			Function f1, Function f2);
	std::function<Real(int, int)> getProductFirstDerivativeDCS(
			GridFunction& gf1, GridFunction& gf2, Delta d, 
			int f1, int f2);
};


#endif

