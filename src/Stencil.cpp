#include "Stencil.hpp"

#include "Debug.hpp"
#include <cmath>


namespace Derivatives
{
	namespace {
		std::function<Real(int,int)> genFd(
				GridFunction& gf, Delta d,
				Direction df)
		{
			Real delta = 0.0;
			int isX=0, isY=0, isZ=0;
			switch(df) /* set right dimension */
			{
				case Direction::xf:
				case Direction::xb:
					delta = d.x;
					isX = 1;
					isY = 0;
					isZ = 0;
					break;
				case Direction::yf:
				case Direction::yb:
					delta = d.y;
					isX = 0;
					isY = 1;
					isZ = 0;
					break;
				case Direction::zf:
				case Direction::zb:
					delta = d.z;
					isX = 0;
					isY = 0;
					isZ = 1;
					break;
				default:
					log_err("dimension could not be assigned");
					break;
			}
			int f = 0, b = 0;
			switch(df) /* set right derivative */
			{
				case Direction::xf:
				case Direction::yf:
				case Direction::zf:
					f = +1;
					b = 0;
					break;
				case Direction::xb:
				case Direction::yb:
				case Direction::zb:
					f = 0;
					b = -1;
					break;
				default:
					log_err("returning ZERO derivative operator");
					return [](int i, int j) { return i*j*0.0; };
					break;
			}
			return [&gf, delta, f, b, isX, isY, isZ](int i, int j)
				{
					return 
						(
						 gf( i + f*isX, j + f*isY)
							-
						 gf( i + b*isX, j + b*isY)
						) / delta;
				};
		}

		std::function<Real(int,int)> genFdd(
				GridFunction& gf, Delta d,
				Direction df)
		{
			Real delta = 0.0;
			int isX=0, isY=0, isZ=0;
			switch(df) /* set right dimension */
			{
				case Direction::xx:
					delta = d.x;
					isX = 1;
					isY = 0;
					isZ = 0;
					break;
				case Direction::yy:
					delta = d.y;
					isX = 0;
					isY = 1;
					isZ = 0;
					break;
				case Direction::zz:
					delta = d.z;
					isX = 0;
					isY = 0;
					isZ = 1;
					break;
				default:
					log_err("dimension could not be assigned");
					break;
			}
			return [&gf, delta, isX, isY, isZ](int i, int j)
				{
					return
						(
						 gf( i + isX, j + isY)
							+
						 gf( i - isX, j - isY)
							-
						 2.0*gf( i , j )
						)
						/ (delta*delta);
				};
		}

		std::function<Real(int, int)> genFG_g(
				GridFunction& gf, GridFunction& gg, Delta d, 
				Function f, Function g)
		{
			Real delta = 0.0;
			int FinX=0, FinY=0, FinZ=0;
			int GinX=0, GinY=0, GinZ=0;
			switch(f)
			{
				case Function::U:
					FinX = 1;
					FinY = 0;
					FinZ = 0;
					break;
				case Function::V:
					FinX = 0;
					FinY = 1;
					FinZ = 0;
					break;
				case Function::W:
					FinX = 0;
					FinY = 0;
					FinZ = 1;
					break;
				default:
					log_err("function dimension could not be assigned");
					break;
			}
			switch(g)
			{ /* g determines the direction of the derivative */
				case Function::U:
					delta = d.x;
					GinX = 1;
					GinY = 0;
					GinZ = 0;
					break;
				case Function::V:
					delta = d.y;
					GinX = 0;
					GinY = 1;
					GinZ = 0;
					break;
				case Function::W:
					delta = d.z;
					GinX = 0;
					GinY = 0;
					GinZ = 1;
					break;
				default:
					log_err("function dimension could not be assigned");
					break;
			}
			return [&gf, &gg, delta, 
					 FinX, FinY, FinZ,
					 GinX, GinY, GinZ
				](int i, int j)
				{
					/* just ask me to explain it to you. 
					 * right now i'm too tired to write a clear explanation. */

					Real fgp = /* this part goes in plus gg-direction */
						(
						 gf( i , j ) 
							+ 
						 gf( i + GinX , j  + GinY ) 
						)
						*
						(
						 gg( i , j ) 
							+ 
						 gg( i + FinX , j + FinY )
						);
					Real fgm = /* this part goes in minus gg-direction */
						(
						 gf( i , j ) 
							+ 
						 gf( i - GinX , j - GinY ) 
						)
						*
						(
						 gg( i - GinX , j - GinY ) 
						 + 
						 gg( i + FinX - GinX , j + FinY - GinY )
						);

					return (fgp - fgm) / (4.0*delta);
				};
		}

		/* DONOR CELL SCHEME */
		std::function<Real(int, int)> genFG_gdc(
				GridFunction& gf, GridFunction& gg, Delta d, 
				Function f, Function g)
		{
			Real delta = 0.0;
			int FinX=0, FinY=0, FinZ=0;
			int GinX=0, GinY=0, GinZ=0;
			switch(f)
			{
				case Function::U:
					FinX = 1;
					FinY = 0;
					FinZ = 0;
					break;
				case Function::V:
					FinX = 0;
					FinY = 1;
					FinZ = 0;
					break;
				case Function::W:
					FinX = 0;
					FinY = 0;
					FinZ = 1;
					break;
				default:
					log_err("function dimension could not be assigned");
					break;
			}
			switch(g)
			{ /* g determines the direction of the derivative */
				case Function::U:
					delta = d.x;
					GinX = 1;
					GinY = 0;
					GinZ = 0;
					break;
				case Function::V:
					delta = d.y;
					GinX = 0;
					GinY = 1;
					GinZ = 0;
					break;
				case Function::W:
					delta = d.z;
					GinX = 0;
					GinY = 0;
					GinZ = 1;
					break;
				default:
					log_err("function dimension could not be assigned");
					break;
			}
			return [&gf, &gg, delta, 
					 FinX, FinY, FinZ,
					 GinX, GinY, GinZ
				](int i, int j)
				{
					/* just ask me to explain it to you. 
					 * right now i'm too tired to write a clear explanation. */

					Real fgp = /* this part goes in plus gg-direction */
						(
						 gf( i , j ) 
							+ 
						 gf( i + GinX , j  + GinY ) 
						)
						*
					std::fabs(
						 gg( i , j ) 
							- 
						 gg( i + FinX , j + FinY )
						);
					Real fgm = /* this part goes in minus gg-direction */
						(
						 gf( i , j ) 
							+ 
						 gf( i - GinX , j - GinY ) 
						)
						*
					std::fabs(
						 gg( i - GinX , j - GinY ) 
							- 
						 gg( i + FinX - GinX , j + FinY - GinY )
						);

					return (fgp - fgm) / (4.0*delta);
				};
		}
	}

	/* call those from outside */


	/* the following two functions only take a direction enum */
	std::function<Real(int, int)> getDerivative(
			GridFunction& gf, Delta d, 
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
				log_err("returning ZERO derivative operator");
				return [](int i, int j){ return i*j*0.0; };
				break;
		};
	}
	std::function<Real(int, int)> getDerivative(
			GridFunction& gf, Delta d, 
			int df)
	{
		Direction dff = static_cast<Derivatives::Direction>(df);
		return Derivatives::getDerivative(gf, d, dff);
	}

	std::function<Real(int, int)> getProductFirstDerivative(
			GridFunction& gf1, GridFunction& gf2, Delta d, 
			Function f1, Function f2)
	{
		return genFG_g(gf1, gf2, d, f1, f2);
	}
	std::function<Real(int, int)> getProductFirstDerivative(
			GridFunction& gf1, GridFunction& gf2, Delta d, 
			int f1, int f2)
	{
		Function ff1 = static_cast<Derivatives::Function>(f1);
		Function ff2 = static_cast<Derivatives::Function>(f2);
		return getProductFirstDerivative(gf1, gf2, d, ff1, ff2);
	}

	std::function<Real(int, int)> getProductFirstDerivativeDCS(
			GridFunction& gf1, GridFunction& gf2, Delta d, 
			Function f1, Function f2)
	{
		return [](int i, int j){ return i*j*0.0; };
		return genFG_gdc(gf1, gf2, d, f1, f2);
	}
	std::function<Real(int, int)> getProductFirstDerivativeDCS(
			GridFunction& gf1, GridFunction& gf2, Delta d, 
			int f1, int f2)
	{
		Function ff1 = static_cast<Derivatives::Function>(f1);
		Function ff2 = static_cast<Derivatives::Function>(f2);
		return getProductFirstDerivative(gf1, gf2, d, ff1, ff2);
	}

};

