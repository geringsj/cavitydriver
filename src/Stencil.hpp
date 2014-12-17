
#ifndef Derivatives_hpp
#define Derivatives_hpp 

#include "GridFunction.hpp"
#include <functional>

/*!
 */

/**
 * This is (still) (kind of) a "Stencil", but (now) it is on steroids.
 * 
 * @author becherml, friesfn, geringsj
 * @date 11/2014
 *
 * The Derivatives namespace generates derivative operators that are demanded 
 * by the user, and returns those operators as functions D(i,j) to the user.
 * This makes things easier for the implementation and usage of derivatives.
 *
 * The following derivative types are supported (see Derivatives::Direction):
 * - first derivative, backward, in x, y (and z) direction
 * - first derivative, forward, in x, y (and z) direction
 * - secons derivative, central, in x, y (and z) direction
 * 
 * For the Momentum Equations we provide first derivatives of a product of functions F and G:
 * - ( F * G ) evaluated at the grid position of F(i,j) and derived in 
 *   the direction that is associated with G.
 * So the order of functions matters: (F*G) =/= (G*F).
 * F and G can be different, or the same function. 
 */
namespace Derivatives
{
	/** Supported directions for first and second derivatives of a function F. */
	enum class Direction {
		xb=-1, 
		/**< First backward derivative in x direction. */
		yb=-2,
		/**< First backward derivative in y direction. */
		zb=-3,
		/**< First backward derivative in z direction. */

		xf=1, 
		/**< First forward derivative in x direction. */
		yf=2, 
		/**< First forward derivative in y direction. */
		zf=3, 
		/**< First forward derivative in z direction. */

		xx=4, 
		/**< Second central derivative in x direction. */
		yy=5, 
		/**< Second central derivative in y direction. */
		zz=6
		/**< Second central derivative in z direction. */
	};
	typedef Direction d;

	/** Supported functions for derivatives of products (F*G). */
	enum class Function {
		U = 1, 
		/**< Corresponds to direction x. */
		V = 2,
		/**< Corresponds to direction y. */
		W = 3
		/**< Corresponds to direction z. */
	};
	typedef Function F;

	/** Function to get first or second derivatives.
	 * @return Returns a function 'Real D(int i, int j)' that evaluates as 
	 * the demanded derivative at (i,j).
	 */
	std::function<Real(int, int)> getDerivative(
			GridFunction& gf, 
			/**< GridFunction to take the derivative of */
			Delta d, 
			/**< Grid spacing delta, needed for finite difference. */
			Direction df
			/**< The direction in which the returned derivative operator should derive. */
			);

	/** Function to get first or second derivatives.
	 * @return Returns a function 'Real D(int i, int j)' that evaluates as 
	 * the demanded derivative at (i,j).
	 *
	 * This one takes an int as derivative direction. 
	 * See the Derivatives::Direction enum for accepted values of df.
	 */
	std::function<Real(int, int)> getDerivative(
			GridFunction& gf, 
			/**< GridFunction to take the derivative of */
			Delta d, 
			/**< Grid spacing delta, needed for finite difference. */
			int df
			/**< Integer that has a value corresponding to the demanded 
			 * derivative operator (see Derivatives::Direction enum). */
			);

	/** Function to get first derivative of product of functions (F*G).
	 * @return Returns a function 'Real D(int i, int j)' that evaluates as 
	 * the demanded derivative at (i,j).
	 *
	 * The semantics here are a little tricky: D(i,j) will evaluate at the
	 * grid position of F(i,j), and the derivative will be taken in the direction
	 * that is associated with G. 
	 * Example: G==V => derive in y direction.
	 * The order of functions matters: (F*G) =/= (G*F).
	 * F and G can be different, or the same (grid-)function. 
	 */
	std::function<Real(int, int)> getProductFirstDerivative(
			GridFunction& gfF, 
			/**< GridFunction for first function, F, of the product. */
			GridFunction& gfG, 
			/**< GridFunction for second function, G, of the product. */
			Delta d, 
			/**< Grid spacing delta, needed for finite difference. */
			Function F, 
			/**< Semantics of first function F: is it U or V (or W) ? */
			Function G
			/**< Semantics of second function G: is it U or V (or W) ? */
			);

	/** Function to get first derivative of product of functions (F*G).
	 * @return Returns a function 'Real D(int i, int j)' that evaluates as 
	 * the demanded derivative at (i,j).
	 *
	 * Just like the other getProductFirstDerivative, but this one takes 
	 * integer values corresponding to the enums for F and G (see Derivatives::Function).
	 */
	std::function<Real(int, int)> getProductFirstDerivative(
			GridFunction& gfF, 
			GridFunction& gfG, 
			Delta d, 
			int F, 
			int G
			);

	/** Returns first derivative of (F*G) with 'Donor-Cell Scheme' applied.
	 * Semantics and usage same as in Derivatives::getProductFirstDerivative.
	 */
	std::function<Real(int, int)> getProductFirstDerivativeDCS(
			GridFunction& gfF, 
			GridFunction& gfG, 
			Delta d, 
			Function F, 
			Function G
			);

	/** Returns first derivative of (F*G) with 'Donor-Cell Scheme' applied.
	 * Semantics and usage same as in Derivatives::getProductFirstDerivative.
	 */
	std::function<Real(int, int)> getProductFirstDerivativeDCS(
			GridFunction& gfF, 
			GridFunction& gfG, 
			Delta d, 
			int F, 
			int G
			);
};










#include <cmath>
/* =========================================================================== */

struct DerivativesgetDerivative {
	private:
	const GridFunction& gf;
	const Real delta; 
	int const f; 
	int const b;
	int const isX; 
	int const isY; 
	int const isZ;
	Derivatives::Direction const df;

	public:
	DerivativesgetDerivative(
			const GridFunction& gf, 
			const Delta d, 
			const Derivatives::Direction df)
		:
			gf(gf),
			delta(
					(
						df==Derivatives::Direction::xf || 
						df==Derivatives::Direction::xb || 
						df==Derivatives::Direction::xx 
						)?(d.x):
					(
					 (
						df==Derivatives::Direction::yf || 
						df==Derivatives::Direction::yb || 
						df==Derivatives::Direction::yy 
						)?(d.y):
					 (
					  (
						df==Derivatives::Direction::zf || 
						df==Derivatives::Direction::zb || 
						df==Derivatives::Direction::zz 
						)?(d.z):(0.0) 
					  ) 
					 ) 
					),
			f(
					(
					 df==Derivatives::Direction::xf ||
					 df==Derivatives::Direction::yf ||
					 df==Derivatives::Direction::zf
					 )?(1):
						(0)
					),
			b(
					(
					 df==Derivatives::Direction::xf ||
					 df==Derivatives::Direction::yf ||
					 df==Derivatives::Direction::zf
					 )?(0):
						(-1)
					),
			isX( 
					(
						df==Derivatives::Direction::xf || 
						df==Derivatives::Direction::xb || 
						df==Derivatives::Direction::xx 
						)?(1):(0) 
					),
			isY(
					(
						df==Derivatives::Direction::yf || 
						df==Derivatives::Direction::yb || 
						df==Derivatives::Direction::yy 
						)?(1):(0) 
					),
			isZ(
					(
						df==Derivatives::Direction::zf || 
						df==Derivatives::Direction::zb || 
						df==Derivatives::Direction::zz 
						)?(1):(0) 
					),
			df(df)
	{
	}
	DerivativesgetDerivative(
			GridFunction& gf, Delta d, 
			int df)
		: DerivativesgetDerivative(gf, d, static_cast<Derivatives::Direction>(df))
	{}

	Real operator()(const int i, const int j) const
	{
		switch(df) /* set right derivative */
		{
			case Derivatives::Direction::xf:
			case Derivatives::Direction::yf:
			case Derivatives::Direction::zf:
			case Derivatives::Direction::xb:
			case Derivatives::Direction::yb:
			case Derivatives::Direction::zb:
				return 
					( 
					 gf( i+ f*isX , j + f*isY) 
						- 
					 gf( i + b*isX, j + b*isY) 
					) / delta;

				break;
			case Derivatives::Direction::xx:
			case Derivatives::Direction::yy:
			case Derivatives::Direction::zz:
				return 
					(
					 gf( i + isX , j + isY) 
						+ 
					 gf( i - isX, j - isY)
						- 
					 2.0*gf( i , j ) 
					)
					/ (delta*delta);
				break;
			default:
				return 0.0*i*j;
				break;
		}
	}
};

/* =========================================================================== */

struct DerivativesgetProductFirstDerivative {
	private:
	const GridFunction& gf;
	const GridFunction& gg;
	const Real delta;
					 
	const int FinX;
	const int FinY;
	const int FinZ;
					 
	const int GinX;
	const int GinY;
	const int GinZ;

	public:
	DerivativesgetProductFirstDerivative(
			const GridFunction& gf1, const GridFunction& gf2, const Delta d, 
			const Derivatives::Function f1, const Derivatives::Function f2) :
		gf(gf1),
		gg(gf2),
		delta(
				(
					Derivatives::Function::U == f2
					)?(d.x):
				(
				 (
					Derivatives::Function::V == f2
					)?(d.y):
				 (
				  (
					Derivatives::Function::W == f2
					)?(d.z):(0.0) 
				  ) 
				 ) 
				),
		FinX( ( f1 == Derivatives::Function::U )?(1):(0) ),
		FinY( ( f1 == Derivatives::Function::V )?(1):(0) ),
		FinZ( ( f1 == Derivatives::Function::W )?(1):(0) ),
		GinX( ( f2 == Derivatives::Function::U )?(1):(0) ),
		GinY( ( f2 == Derivatives::Function::V )?(1):(0) ),
		GinZ( ( f2 == Derivatives::Function::W )?(1):(0) )
	{}
	DerivativesgetProductFirstDerivative(
			const GridFunction& gf1, 
			const GridFunction& gf2, 
			const Delta d, 
			const int f1, 
			const int f2) :
		DerivativesgetProductFirstDerivative(gf1, gf2, d, 
				static_cast<Derivatives::Function>(f1), 
				static_cast<Derivatives::Function>(f2))
	{}

	Real operator()(const int i, const int j) const
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
	}
};

/* =========================================================================== */

struct DerivativesgetProductFirstDerivativeDCS {
	private:
	const GridFunction& gf;
	const GridFunction& gg;
	const Real delta;
					 
	const int FinX;
	const int FinY;
	const int FinZ;
					 
	const int GinX;
	const int GinY;
	const int GinZ;

	public:
	DerivativesgetProductFirstDerivativeDCS (
			const GridFunction& gf1, const GridFunction& gf2, const Delta d, 
			const Derivatives::Function f1, const Derivatives::Function f2) :
		gf(gf1),
		gg(gf2),
		delta(
				(
					Derivatives::Function::U == f2
					)?(d.x):
				(
				 (
					Derivatives::Function::V == f2
					)?(d.y):
				 (
				  (
					Derivatives::Function::W == f2
					)?(d.z):(0.0) 
				  ) 
				 ) 
				),
		FinX( ( f1 == Derivatives::Function::U )?(1):(0) ),
		FinY( ( f1 == Derivatives::Function::V )?(1):(0) ),
		FinZ( ( f1 == Derivatives::Function::W )?(1):(0) ),
		GinX( ( f2 == Derivatives::Function::U )?(1):(0) ),
		GinY( ( f2 == Derivatives::Function::V )?(1):(0) ),
		GinZ( ( f2 == Derivatives::Function::W )?(1):(0) )
	{}
	DerivativesgetProductFirstDerivativeDCS (
			const GridFunction& gf1, 
			const GridFunction& gf2, 
			const Delta d, 
			const int f1, 
			const int f2) :
		DerivativesgetProductFirstDerivativeDCS(gf1, gf2, d, 
				static_cast<Derivatives::Function>(f1), 
				static_cast<Derivatives::Function>(f2))
	{}

	Real operator()(const int i, const int j) const
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
	}
};





#endif

