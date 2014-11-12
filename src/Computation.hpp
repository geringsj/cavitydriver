//! This namespace implements the computation
/*!
 * @author becherml, friesfn, geringsj
 * @date 2014
 */

#ifndef Computation_hpp
#define Computation_hpp

#include "Domain.hpp"
#include "Stencil.hpp"
#include "Debug.hpp"

#include <cmath>

namespace Computation
{
Real computeTimestep(
		Domain& domain,
		const Real tau, const Real Re)
{
	Real result;
	Real uMax, vMax;

	uMax = domain.getVeolcity().m_u.getMaxValueGridFunction();//domain.getBeginInnerDomains(), domain.getEndInnerDomainU());
	vMax = domain.getVeolcity().m_v.getMaxValueGridFunction();//domain.getBeginInnerDomains(), domain.getEndInnerDomainV());
	result = std::fmin(
			domain.getDelta().x*domain.getDelta().x *
			domain.getDelta().y*domain.getDelta().y *Re
			/
			(2.0*(domain.getDelta().x*domain.getDelta().x+domain.getDelta().y*domain.getDelta().y)),

			std::fmin(
				domain.getDelta().x/std::abs(uMax),
				domain.getDelta().y/std::abs(vMax)
				));

	return tau * result; /* tau is some safety factor in (0,1] */
}

void computeMomentumEquationsFGH(
		Domain& domain,
		const Real deltaT, const Real Re)
{
	for(uint d=0; d<DIMENSIONS; d++)
	{
		auto Fxx = 
			Derivatives::getDerivative(
					domain.getVeolcity()[d],
					domain.getDelta(),
					Derivatives::Direction::xx);
		auto Fyy = 
			Derivatives::getDerivative(
					domain.getVeolcity()[d],
					domain.getDelta(),
					Derivatives::Direction::yy);
		/* TODO:
		 * are those derivatives right? 
		 * might explain the weird output 
		 * dimension is tricky here */
		auto FF_df = 
			Derivatives::getDerivative(
					domain.getVeolcity()[d],
					domain.getVeolcity()[d],
					domain.getDelta(),
					1+d, 
					1+d,
					1+6+d); /* 7==Derivatives::Direction::_x */
		int Gdim = (d+1)%DIMENSIONS;
		auto FG_dg = 
			Derivatives::getDerivative(
					domain.getVeolcity()[d],
					domain.getVeolcity()[Gdim],
					domain.getDelta(),
					1+d,
					1+Gdim,
					1+6+Gdim);

		/* the formula: */
		for(int i=domain.getBeginInnerDomains()[0]; 
				i<= domain.getEndInnerDomain()[d][0]; i++)
			for(int j=domain.getBeginInnerDomains()[1]; 
					j<= domain.getEndInnerDomain()[d][1]; j++)
		{
			Point gpoint = Point(
					/* here we are off by 1 or 2 because of the boundary */
						Real(i)/domain.getEndInnerDomain()[d][0], 
						Real(j)/domain.getEndInnerDomain()[d][1]);

			domain.getPreliminaryVeolcity()[d](i, j) =
				domain.getVeolcity()[d](i, j) + 
				(deltaT/Re)*( Fxx(i,j) + Fyy(i,j) )
				- FF_df(i,j) - FG_dg(i,j)
				+ domain.g(d,gpoint);
		}
	}
}

void computeRighthandSide(
		Domain& domain,
		const Real deltaT)
{
	auto Fxb = Derivatives::getDerivative(
			domain.F(),
			domain.getDelta(),
			Derivatives::Direction::xb);
	auto Gyb = Derivatives::getDerivative(
			domain.G(),
			domain.getDelta(),
			Derivatives::Direction::yb);

	for(int i=domain.getBeginInnerDomains()[0]; 
			i<= domain.getEndInnerDomain()[3][0]; i++)
		for(int j=domain.getBeginInnerDomains()[1]; 
				j<= domain.getEndInnerDomain()[3][1]; j++)
	{
		domain.p()(i,j) = ( Fxb(i,j) + Gyb(i,j) )/deltaT;
	}
}

void computeNewVelocities(
		Domain& domain,
		const Real deltaT)
{

	for(uint d=0; d<DIMENSIONS; d++)
	{
		auto Pdb = Derivatives::getDerivative(
				domain.p(),
				domain.getDelta(),
				-(int)(1+d));

		for(int i=domain.getBeginInnerDomains()[0]; 
				i<= domain.getEndInnerDomain()[d][0]; i++)
			for(int j=domain.getBeginInnerDomains()[1]; 
					j<= domain.getEndInnerDomain()[d][1]; j++)
		{
			domain.getVeolcity()[d]( i,j ) =
				domain.getPreliminaryVeolcity()[d]( i,j ) - deltaT*Pdb(i,j);
		}
	}
}
};

#endif
