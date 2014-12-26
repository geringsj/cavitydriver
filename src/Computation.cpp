#include "Computation.hpp"
#include "Debug.hpp"
#include "Stencil.hpp"

#include <cmath>

namespace Computation
{

	Real computeTimestep(Domain& domain, const Real tau, const Real Re)
	{
		Real uMax = domain.getVelocity().m_u.getMaxValueGridFunction();
		Real vMax = domain.getVelocity().m_v.getMaxValueGridFunction();
		Delta maxVels(uMax, vMax);

		return computeTimestepFromMaxVelocities(maxVels, domain.getDelta(), tau, Re);
	}

	Real computeTimestepFromMaxVelocities(
			Delta maxVelocities, Delta cellsDelta, const Real tau, const Real Re)
	{
		Real uMax = maxVelocities.x;
		Real vMax = maxVelocities.y;

		Real dxx = pow(cellsDelta.x, 2.0);
		Real dyy = pow(cellsDelta.y, 2.0);

		Real result = std::fmin(
				(dxx * dyy * Re) / (2.0*(dxx+dyy)), std::fmin(
					cellsDelta.x / std::fabs(uMax),
					cellsDelta.y / std::fabs(vMax) ));

		return tau * result; /* tau is some safety factor in (0,1] */
	}


	void computePreliminaryVelocities(
			Domain& domain, const Real deltaT, const Real Re, const Real alpha)
	{
		for(uint D=0; D<DIMENSIONS; D++)
		{
			/* get second derivatives for momentum equation for dimension D */
			auto Fxx =
				Derivatives::getDerivative(domain.getVelocity()[D], domain.getDelta(),
						Derivatives::Direction::xx);
			auto Fyy = 
				Derivatives::getDerivative(domain.getVelocity()[D], domain.getDelta(),
						Derivatives::Direction::yy);

			/* derivative of product of functions */
			auto FFf = Derivatives::getProductFirstDerivative(
					domain.getVelocity()[D],domain.getVelocity()[D],domain.getDelta(),
					1+D,1+D);
			/* donor cell scheme correction operator */
			auto FFfdc = Derivatives::getProductFirstDerivativeDCS(
					domain.getVelocity()[D],domain.getVelocity()[D],domain.getDelta(),
					1+D,1+D);

			/* get product derivatives / donor cell scheme operator for next dimension
			 * with respect to D */
			int Gdim = (D+1)%DIMENSIONS;
			auto FGg = Derivatives::getProductFirstDerivative(
					domain.getVelocity()[D],domain.getVelocity()[Gdim],domain.getDelta(),
					1+D,1+Gdim);

			auto FGgdc = Derivatives::getProductFirstDerivativeDCS(
					domain.getVelocity()[D],domain.getVelocity()[Gdim],domain.getDelta(),
					1+D,1+Gdim);

			/* the formula: */
			for_vecrange(i,j,domain.getInnerRanges()[D])
			{
				domain.getPreliminaryVelocity()[D]( i,j ) =
					domain.getVelocity()[D]( i,j ) 
					+ deltaT*
					( 
					 ( Fxx(i,j) + Fyy(i,j) )/Re
					 - ( FFf(i,j) + alpha*FFfdc(i,j) ) - ( FGg(i,j) + alpha*FGgdc(i,j) )
					 + domain.g(D)
					);
			}
		}
	}



	void computeRighthandSide(
			Domain& domain,
			const Real deltaT)
	{
		auto Fxb = Derivatives::getDerivative(domain.F(),domain.getDelta(),
				Derivatives::Direction::xb);

		auto Gyb = Derivatives::getDerivative(domain.G(),domain.getDelta(),
				Derivatives::Direction::yb);

		for_vecrange(i,j,domain.getInnerRangeP())
		{
			domain.rhs()(i,j) = ( Fxb(i,j) + Gyb(i,j) )/deltaT;
		}
	}


	void computeNewVelocities(
			Domain& domain,
			const Real deltaT)
	{
		for(uint D=0; D<DIMENSIONS; D++)
		{
			auto Pdf = Derivatives::getDerivative(domain.p(),domain.getDelta(),
					(1+D)); /* D+1 ~ first forward derivative in (D+1) direction (1 => x) */

			for_vecrange(i,j,domain.getInnerRanges()[D])
			{
				domain.getVelocity()[D]( i,j ) =
					domain.getPreliminaryVelocity()[D]( i,j ) - deltaT*Pdf(i,j);
			}
		}
	}


};


/* TODO: sorry, but the following functions are full of indexing errors. */
// /**
//  * TODO:
//  * Please move to IO after we have one that is able to deal with MPI stuff.
//  */
// void ComputeVorticity(GridFunction& vorticity, Domain domain)
// {
// 	GridFunction& u = domain.u();
// 	GridFunction& v = domain.v();
// 	Real delta_x = domain.getDelta().x;
// 	Real delta_y = domain.getDelta().y;
// 	/* watch out: on all other grids starting from i=0/j=0 means taking with you
// 	 * the mostly useless left/lower boundary of that grid */
// 	for (int i = 0; i < domain.getDimension().i; i++)
// 	{
// 		for (int j = 0; j < domain.getDimension().j; j++)
// 		{
// 			vorticity(i, j) = (u(i, j + 1) - u(i, j)) / delta_y - (v(i + 1, j) - v(i, j)) / delta_x;
// 		}
// 	}
// }
// 
// void ComputeFlowFunction(GridFunction& flow, Domain domain)
// {
// 	GridFunction u = domain.u(); /* TODO: WROOONG */
// 	Real delta_y = domain.getDelta().y;
// 	/* Set the border values */
// 	for (int i = 0; i < domain.getDimension().i; i++)
// 	{
// 		flow(i, 0) = 0.0;
// 	}
// 	for (int i = 0; i < domain.getDimension().i; i++)
// 	{ /* TODO: why j=_1_ here ? */
// 		for (int j = 1; j < domain.getDimension().j; j++)
// 		{
// 			flow(i, j) = flow(i, j - 1) + u(i, j) * delta_y;
// 		}
// 	}
// }
