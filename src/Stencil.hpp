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
	/* here we have a wide variety of directions we 
	 * can take the derivative with respect to 
	 * think of it like in the formula: 
	 * (d F) / (d r) 
	 * where F is a function from Derivatives::Function
	 * and r is a direction from Derivatives::Direction
	 *
	 * rf / rb are the (first) forward / backward derivative, respectively 
	 * rr is the central second derivative 
	 * _r is the first derivative of a product of functions F and G: d(F * G)/dr
	 *
	 * the getDerivative() functions return the demanded derivatives as functions 
	 * that evaluate the given GridFunction using the given gridspacing delta 
	 * in an artbitrary gridpoint (i,j) and return the value of the 
	 * derivative in that point*/
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

		_x=7, 
		_y=8, 
		_z=9
	};

	/* TODO: document our supported derivative formats, if they are not 
	 * easy to derive from the definitions of the enums */

	std::function<Real(int, int)> getDerivative(GridFunction& gf, Delta d, Direction df);
	std::function<Real(int, int)> getDerivative(GridFunction& gf, Delta d, int df);

	std::function<Real(int, int)> getDerivative(
			GridFunction& gf1, GridFunction& gf2, Delta d, 
			Function f1, Function f2, Direction df);
	std::function<Real(int, int)> getDerivative(
			GridFunction& gf1, GridFunction& gf2, Delta d, 
			int f1, int f2, int df);
};


#endif

