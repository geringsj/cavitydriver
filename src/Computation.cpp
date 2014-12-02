#include "Computation.hpp"
#include "Debug.hpp"
#include "Stencil.hpp"

#include <cmath>

namespace Computation
{


Real computeTimestep(
		Domain& domain,
		const Real tau, const Real Re)
{
	Real uMax = domain.getVelocity().m_u.getMaxValueGridFunction();
	Real vMax = domain.getVelocity().m_v.getMaxValueGridFunction();

	Real dxx = pow(domain.getDelta().x, 2.0);
	Real dyy = pow(domain.getDelta().y, 2.0);

	Real result = std::fmin(
			(dxx * dyy * Re) / (2.0*(dxx+dyy)), std::fmin(
			domain.getDelta().x / std::fabs(uMax),
			domain.getDelta().y / std::fabs(vMax) ));

	return tau * result; /* tau is some safety factor in (0,1] */
}


void computePreliminaryVelocities(
		Domain& domain,
		const Real deltaT, const Real Re, const Real alpha)
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
		forall(i,j,domain.getBeginInnerDomains(), domain.getEndInnerDomain()[D])
		{
			/* this point is used for the outer force g_D in dimension D, 
			 * it tells the outer force function where in the domain it should evaluate */
			Point gpoint = Point(
					/* here we are off by 1 or 2 because of the boundary */
						Real(i)/domain.getEndInnerDomain()[D][0], 
						Real(j)/domain.getEndInnerDomain()[D][1]);

			/* the formula (for real): */
			domain.getPreliminaryVelocity()[D]( i,j ) =
				domain.getVelocity()[D]( i,j ) 
				+ deltaT*
				( 
				 ( Fxx(i,j) + Fyy(i,j) )/Re
				  - ( FFf(i,j) + alpha*FFfdc(i,j) ) - ( FGg(i,j) + alpha*FGgdc(i,j) )
				 + domain.g(D,gpoint)
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
	
	forall(i,j,domain.getBeginInnerDomains(), domain.getEndInnerDomain()[3])
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

		forall(i,j,domain.getBeginInnerDomains(), domain.getEndInnerDomain()[D])
		{
			domain.getVelocity()[D]( i,j ) =
				domain.getPreliminaryVelocity()[D]( i,j ) - deltaT*Pdf(i,j);
		}
	}
}


};



/**
 * TODO:
 * Please move to IO after we have one that is able to deal with MPI stuff.
 */
void ComputeVorticity(GridFunction& vorticity, Domain domain)
{
	GridFunction u = domain.u();
	GridFunction v = domain.v();
	Real delta_x = domain.getDelta().x;
	Real delta_y = domain.getDelta().y;
	for (int i = 0; i < domain.getDimension().i; i++)
	{
		for (int j = 0; j < domain.getDimension().j; j++)
		{
			vorticity(i, j) = (u(i, j + 1) - u(i, j)) / delta_y - (v(i + 1, j) - v(i, j)) / delta_x;
		}
	}
}

void ComputeFlowFunction(GridFunction& flow, Domain domain)
{
	GridFunction u = domain.u();
	Real delta_y = domain.getDelta().y;
	/* Set the border values */
	for (int i = 0; i < domain.getDimension().i; i++)
	{
		flow(i, 0) = 0.0;
	}
	for (int i = 0; i < domain.getDimension().i; i++)
	{
		for (int j = 1; j < domain.getDimension().j; j++)
		{
			flow(i, j) = flow(i, j - 1) + u(i, j) * delta_y;
		}
	}
}