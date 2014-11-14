#include "Computation.hpp"
#include "Debug.hpp"

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
		const Real deltaT, const Real Re)
{
	for(uint d=0; d<DIMENSIONS; d++)
	{
		// auto Fxx = 
		// 	Derivatives::getDerivative(
		// 			domain.getVeolcity()[d],
		// 			domain.getDelta(),
		// 			Derivatives::Direction::xx);
		// auto Fyy = 
		// 	Derivatives::getDerivative(
		// 			domain.getVeolcity()[d],
		// 			domain.getDelta(),
		// 			Derivatives::Direction::yy);
		
		/* TODO:
		 * are those derivatives right? 
		 * might explain the weird output 
		 * dimension is tricky here */
		// auto FF_df = 
		// 	Derivatives::getDerivative(
		// 			domain.getVeolcity()[d],
		// 			domain.getVeolcity()[d],
		// 			domain.getDelta(),
		// 			1+d, 
		// 			1+d,
		// 			1+6+d); /* 7==Derivatives::Direction::_x */
		// int Gdim = (d+1)%DIMENSIONS;
		// auto FG_dg = 
		// 	Derivatives::getDerivative(
		// 			domain.getVeolcity()[d],
		// 			domain.getVeolcity()[Gdim],
		// 			domain.getDelta(),
		// 			1+d,
		// 			1+Gdim,
		// 			1+6+Gdim);

		/* the formula: */
		forall(i,j,domain.getBeginInnerDomains(), domain.getEndInnerDomain()[d])
		{
			Point gpoint = Point(
					/* here we are off by 1 or 2 because of the boundary */
						Real(i)/domain.getEndInnerDomain()[d][0], 
						Real(j)/domain.getEndInnerDomain()[d][1]);

			Real fph, fmh, gph, gmh, dd;
			Real Fxx, Fyy;

			if(d==0)
			{
				Fxx = (domain.u()(i+1,j) - 2.0*domain.u()(i,j) + domain.u()(i-1,j)) / (pow(domain.getDelta().x,2.0));
				Fyy = (domain.u()(i,j+1) - 2.0*domain.u()(i,j) + domain.u()(i,j-1)) / (pow(domain.getDelta().y,2.0));
			}
			else
			{
				Fxx = (domain.v()(i+1,j) - 2.0*domain.v()(i,j) + domain.v()(i-1,j)) / (pow(domain.getDelta().x,2.0));
				Fyy = (domain.v()(i,j+1) - 2.0*domain.v()(i,j) + domain.v()(i,j-1)) / (pow(domain.getDelta().y,2.0));
			}

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

			if(d==0)
			{
				/* V*U / x */
				fph = domain.u()(i,j) + domain.u()(i+1,j) ;
				gph = domain.v()(i+1,j) + domain.v()(i+1,j-1) ;
				fmh = domain.u()(i,j) + domain.u()(i-1,j) ;
				gmh = domain.v()(i,j) + domain.v()(i,j-1) ;
				dd = domain.getDelta().x;
			}
			else{
				/* V*U / y */
				fph = domain.u()(i,j) + domain.u()(i,j+1) ;
				gph = domain.v()(i,j) + domain.v()(i+1,j) ;
				fmh = domain.u()(i,j) + domain.u()(i,j-1) ;
				gmh = domain.v()(i,j-1) + domain.v()(i+1,j-1) ;
				dd = domain.getDelta().y;
			}
			Real FG_dg = (fph*gph - fmh*gmh) / (4.0 * dd);

			domain.getPreliminaryVelocity()[d](i, j) =
				domain.getVelocity()[d](i, j) 
				+ deltaT*
				( 
				 ((Fxx/*(i,j)*/ + Fyy/*(i,j)*/ )/Re) 
				 - FF_df/*(i,j)*/ - FG_dg/*(i,j)*/ + domain.g(d,gpoint)
				);
		}
	}
}



void computeRighthandSide(
		Domain& domain,
		const Real deltaT)
{
	//auto Fxb = Derivatives::getDerivative(
	//		domain.F(),
	//		domain.getDelta(),
	//		Derivatives::Direction::xb);
	//auto Gyb = Derivatives::getDerivative(
	//		domain.G(),
	//		domain.getDelta(),
	//		Derivatives::Direction::yb);
	
	forall(i,j,domain.getBeginInnerDomains(), domain.getEndInnerDomain()[3])
	{
		Real Fxb = (domain.F()(i,j) - domain.F()(i-1,j)) / (domain.getDelta().x);
		Real Gyb = (domain.G()(i,j) - domain.G()(i,j-1)) / (domain.getDelta().y);

		domain.rhs()(i,j) = ( Fxb/*(i,j)*/ + Gyb/*(i,j)*/ )/deltaT;
	}
}


void computeNewVelocities(
		Domain& domain,
		const Real deltaT)
{
	for(uint d=0; d<DIMENSIONS; d++)
	{
		//auto Pdb = Derivatives::getDerivative(
		//		domain.p(),
		//		domain.getDelta(),
		//		(int)(1+d));

		forall(i,j,domain.getBeginInnerDomains(), domain.getEndInnerDomain()[d])
		{
			Real Pdf;
			
			if(d==0)
				Pdf = (domain.p()(i+1,j) - domain.p()(i,j))/(domain.getDelta().x);
			else
				Pdf = (domain.p()(i,j+1) - domain.p()(i,j))/(domain.getDelta().y);

			domain.getVelocity()[d]( i,j ) =
				domain.getPreliminaryVelocity()[d]( i,j ) - deltaT*Pdf/*(i,j)*/;
		}
	}
}


};
