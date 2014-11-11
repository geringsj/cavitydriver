//! This is still kind of a "Stencil", but now it is on steroids
/*!
 * @author becherml, friesfn, geringsj
 * @date 2014
 */

#ifndef Derivatives_hpp
#define Derivatives_hpp 

#include "Structs.hpp"
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

	/* attention: we don't handle weird Function-Direction cases here, 
	 * only the ones we will use later
	 */

	std::function<Real(int,int)> genFd(
			GridFunction& gf, Delta d,
			Direction df)
	{
		Real delta = 0.0;
		int dimX=0, dimY=0, dimZ=0;
		int add1=0, add2=0;
		switch(df) /* set right dimension */
		{
			case Direction::xf:
			case Direction::xb:
				delta = d.x;
				dimX = 1;
				dimY = 0;
				dimZ = 0;
				break;
			case Direction::yf:
			case Direction::yb:
				delta = d.y;
				dimX = 0;
				dimY = 1;
				dimZ = 0;
				break;
			case Direction::zf:
			case Direction::zb:
				delta = d.z;
				dimX = 0;
				dimY = 0;
				dimZ = 1;
				break;
			default:
				break;
		}
		switch(df) /* set right derivative */
		{
			case Direction::xf:
			case Direction::yf:
			case Direction::zf:
				/* forward */
				add1 = +1; /* index on the left */
					/* minus */
				add2 = 0; /* index on the right */
				break;
			case Direction::xb:
			case Direction::yb:
			case Direction::zb:
				/* backward */
				add1 = 0; /* index on the left */
					/* minus */
				add2 = -1; /* index the right */
				break;
			default:
				return [](int i, int j) { return i*j*0.0; };
				break;
		}
		return [&gf, delta, add1, add2, dimX, dimY, dimZ](int i, int j)
			{
				return (gf( i+ add1*dimX , j + add1*dimY) - gf( i + add2*dimX, j + add2*dimY)) / delta;
			};
	}

	std::function<Real(int,int)> genFdd(
			GridFunction& gf, Delta d,
			Direction df)
	{
		Real delta = 0.0;
		int dimX=0, dimY=0, dimZ=0;
		int add1=+1, add2=0, add3=-1;
		switch(df) /* set right dimension */
		{
			case Direction::xx:
				delta = d.x;
				dimX = 1;
				dimY = 0;
				dimZ = 0;
				break;
			case Direction::yy:
				delta = d.y;
				dimX = 0;
				dimY = 1;
				dimZ = 0;
				break;
			case Direction::zz:
				delta = d.z;
				dimX = 0;
				dimY = 0;
				dimZ = 1;
				break;
			default:
				break;
		}
		return [&gf, delta, add1, add2, add3, dimX, dimY, dimZ](int i, int j)
			{
				return 
					(gf( i+ add1*dimX , j + add1*dimY) 
					 - 2.0*gf( i + add2*dimX, j + add2*dimY) 
					 + gf( i + add3*dimX, j + add3*dimY)) 
					/ (delta*delta);
			};
	}

	std::function<Real(int, int)> genFG_d(
			GridFunction& gf1, GridFunction& gf2, 
			Delta d, 
			Function f1, Function f2, Direction df)
	{
		Real delta = 0.0;
		int f1dimX=0, f1dimY=0, f1dimZ=0;
		int f2dimX=0, f2dimY=0, f2dimZ=0;
		switch(df) /* set right dimension */
		{
			case Direction::_x:
				delta = d.x;
				break;
			case Direction::_y:
				delta = d.y;
				break;
			case Direction::_z:
				delta = d.z;
				break;
			default:
				break;
		}
		switch(f1)
		{
			case Function::U:
				f1dimX = 1;
				f1dimY = 0;
				f1dimZ = 0;
				break;
			case Function::V:
				f1dimX = 0;
				f1dimY = 1;
				f1dimZ = 0;
				break;
			case Function::W:
				f1dimX = 0;
				f1dimY = 0;
				f1dimZ = 1;
				break;
			default:
				break;
		}
		switch(f2)
		{
			case Function::U:
				f2dimX = 1;
				f2dimY = 0;
				f2dimZ = 0;
				break;
			case Function::V:
				f2dimX = 0;
				f2dimY = 1;
				f2dimZ = 0;
				break;
			case Function::W:
				f2dimX = 0;
				f2dimY = 0;
				f2dimZ = 1;
				break;
			default:
				break;
		}
		int add=+1, sub=-1;
		return [&gf1, &gf2, delta, add, sub, 
				 f1dimX, f1dimY, f1dimZ,
				 f2dimX, f2dimY, f2dimZ
			](int i, int j)
			{
				Real fgp = 
					(gf1( i , j ) + gf1( i + f2dimX*add , j + f2dimY*add)) 
					* 
					(gf2( i , j ) + gf2( i + f1dimX*add, j + f1dimY*add));
				Real fgm = 
					(gf1( i + f2dimX*sub, j + f2dimY*sub) + gf1( i , j )) 
					* 
					(gf2( i + f2dimX*sub, j + f2dimY*sub) 
					 + 
					 gf2( i + f2dimX*sub + f1dimX*add, j + f2dimY*sub + f1dimY*add));
				return (fgp - fgm) / (4.0*delta);
			};
	}


	/* call those from outside */

	std::function<Real(int, int)> getDerivative(
			GridFunction& gf, 
			Delta d, 
			Direction df)
	{
		switch(df)
		{
			case Direction::xf:
			case Direction::yf:
			case Direction::zf:
			case Direction::xb:
			case Direction::yb:
			case Direction::zb:
				return Derivatives::genFd(gf, d, df);
				break;
			case Direction::xx:
			case Direction::yy:
			case Direction::zz:
				return genFdd(gf, d, df);
				break;
			default:
				return [](int i, int j){ return i*j*0.0; };
				break;
		};
	}
	std::function<Real(int, int)> getDerivative(
			GridFunction& gf, 
			Delta d, 
			int df)
	{
		Direction dff = static_cast<Derivatives::Direction>(df);
		return Derivatives::getDerivative(gf, d, dff);
	}

	std::function<Real(int, int)> getDerivative(
			GridFunction& gf1, GridFunction& gf2, 
			Delta d, 
			Function f1, Function f2, Direction df)
	{
		switch(df)
		{
			case Direction::_x:
			case Direction::_y:
			case Direction::_z:
				return genFG_d(gf1, gf2, d, f1, f2, df);
				break;
			default:
				return [](int i, int j){ return i*j*0.0; };
				break;
		}
	}
	std::function<Real(int, int)> getDerivative(
			GridFunction& gf1, GridFunction& gf2, 
			Delta d, 
			int f1, int f2, int df)
	{
		Function ff1 = static_cast<Derivatives::Function>(f1);
		Function ff2 = static_cast<Derivatives::Function>(f2);
		Direction dff = static_cast<Derivatives::Direction>(df);
		return getDerivative(gf1, gf2, d, ff1, ff2, dff);
	}

};


#endif

