
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


#endif

