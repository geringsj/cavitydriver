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
	for(uint d=0; d<DIMENSIONS; d++)
	{
		auto Fxx = 
			Derivatives::getDerivative(domain.getVelocity()[d], domain.getDelta(),
					Derivatives::Direction::xx);
		auto Fyy = 
			Derivatives::getDerivative(domain.getVelocity()[d], domain.getDelta(),
					Derivatives::Direction::yy);
		
		/* 
		 * TODO: 
		 * 	drink a fine smoky whisky when putting the following code into service
		 *
		 * auto FFf = Derivatives::getProductFirstDerivative(
		 * 	domain.getVeolcity()[d],domain.getVeolcity()[d],domain.getDelta(),
		 * 	1+d,1+d);
		 *
		 * auto FFfdc = Derivatives::getProductFirstDerivativeDCS(
		 * 	domain.getVeolcity()[d],domain.getVeolcity()[d],domain.getDelta(),
		 * 	1+d,1+d);
		 *
		 * int Gdim = (d+1)%DIMENSIONS;
		 * auto FGg = Derivatives::getProductFirstDerivative(
		 * 	domain.getVeolcity()[d],domain.getVeolcity()[Gdim],domain.getDelta(),
		 * 	1+d,1+Gdim);
		 *
		 * auto FGgdc = Derivatives::getProductFirstDerivativeDCS(
		 * 	domain.getVeolcity()[d],domain.getVeolcity()[Gdim],domain.getDelta(),
		 * 	1+d,1+Gdim);
		 */

		/* the formula: */
		forall(i,j,domain.getBeginInnerDomains(), domain.getEndInnerDomain()[d])
		{
			Point gpoint = Point(
					/* here we are off by 1 or 2 because of the boundary */
						Real(i)/domain.getEndInnerDomain()[d][0], 
						Real(j)/domain.getEndInnerDomain()[d][1]);


			/* TODO: later, delete from here: */
/* ===========================================================================*/
			Real fph, fmh, gph, gmh, dd;
			if(d==0)
			{
				/* U*U / x */
				fph = gph = domain.u()(i,j) + domain.u()(i+1,j) ;
				fmh = gmh = domain.u()(i,j) + domain.u()(i-1,j) ;
				dd = domain.getDelta().x;
			}
			else{
				/* V*V / y */
				fph = gph = domain.v()(i,j) + domain.v()(i,j+1) ;
				fmh = gmh = domain.v()(i,j) + domain.v()(i,j-1) ;
				dd = domain.getDelta().y;
			}
			Real FF_df = (fph*gph - fmh*gmh) / (4.0 * dd);

			//if(d==0)
			//{
			//	/* V*U / x */
			//	fph = domain.u()(i,j) + domain.u()(i+1,j) ;
			//	gph = domain.v()(i+1,j) + domain.v()(i+1,j-1) ;
			//	fmh = domain.u()(i,j) + domain.u()(i-1,j) ;
			//	gmh = domain.v()(i,j) + domain.v()(i,j-1) ;
			//	dd = domain.getDelta().x;
			//}
			//else{
			//	/* V*U / y */
			//	fph = domain.u()(i,j) + domain.u()(i,j+1) ;
			//	gph = domain.v()(i,j) + domain.v()(i+1,j) ;
			//	fmh = domain.u()(i,j) + domain.u()(i,j-1) ;
			//	gmh = domain.v()(i,j-1) + domain.v()(i+1,j-1) ;
			//	dd = domain.getDelta().y;
			//}
			//Real FG_dg = (fph*gph - fmh*gmh) / (4.0 * dd);
			
			Real FG_dg;
			Real u_c =domain.u()(i,j);
			Real u_n =domain.u()(i,j + 1);
			Real u_s =domain.u()(i,j - 1);
			Real u_w = domain.u()(i - 1,j);
			Real u_nw = domain.u()(i - 1,j + 1);
			Real v_c = domain.v()(i,j);
			Real v_e = domain.v()(i + 1,+ j);
			Real v_s = domain.v()(i,j - 1);
			Real v_w = domain.v()(i - 1,j);
			Real v_se = domain.v()(i + 1,j - 1);
			Point h = domain.getDelta();
			if(d==1)
			{
				/* V*U / x */

				FG_dg =
					(1.0 / h[0]) * ((((u_c + u_n)*(v_c + v_e)) / 4.0) - (((u_w + u_nw)*(v_w + v_c)) / 4.0)) +
					(alpha / h[0]) *(((std::abs(u_c + u_n)*(v_c - v_e)) / 4.0) - ((std::abs(u_w + u_nw)*(v_w - v_c)) / 4.0));
			}
			else{
				/* V*U / y */

				FG_dg =
					(1.0 / h[1]) * ((((v_c + v_e)*(u_c + u_n)) / 4.0) - (((v_s + v_se)*(u_s + u_c)) / 4.0)) +
					(alpha / h[1]) *(((std::abs(v_c + v_e)*(u_c - u_n)) / 4.0) - ((std::abs(v_s + v_se)*(u_s - u_c)) / 4.0));
			}
/* ===========================================================================*/
			/* TODO: later, delete until here */

			domain.getPreliminaryVelocity()[d]( i,j ) =
				domain.getVelocity()[d]( i,j ) 
				+ deltaT*
				( 
				 ( Fxx(i,j) + Fyy(i,j) )/Re

				 - FF_df - FG_dg 
				 /* TODO: replace the line above by:
				  *
				  * 	- ( FFf(i,j) + alpha*FFfdc(i,j) ) - ( FGg(i,j) + alpha*FGgdc(i,j) )
				  *
				  * and comment in the function calls for the operators at the 
				  * start of this function 
				  *
				  * (don't forget the whisky) */

				 + domain.g(d,gpoint)
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
	for(uint d=0; d<DIMENSIONS; d++)
	{
		auto Pdf = Derivatives::getDerivative(domain.p(),domain.getDelta(),
				(1+d)); /* d+1 ~ first forward derivative in (d+1) direction (1 => x) */

		forall(i,j,domain.getBeginInnerDomains(), domain.getEndInnerDomain()[d])
		{
			domain.getVelocity()[d]( i,j ) =
				domain.getPreliminaryVelocity()[d]( i,j ) - deltaT*Pdf(i,j);
		}
	}
}


};
