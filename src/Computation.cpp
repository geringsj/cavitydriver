#include "Computation.hpp"
#include "Stencil.hpp"

#include <cmath>



Real Computation::computeTimestep(
		Domain& domain,
		//const Real uMax, const Real vMax, const Point h, 
		const Real tau, const Real Re)
{
	Real result;
	Real uMax, vMax;

	uMax = domain.getVeolcity().m_u.getMaxValueGridFunction();
	vMax = domain.getVeolcity().m_v.getMaxValueGridFunction();

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

//void Computation::setBoundaryX(GridFunction& x){}

void Computation::computeMomentumEquationsFGH(
		Domain& domain,
		//GridFunction& f, GridFunction& g,
		//GridFunction& u, GridFunction& v, 
		//GridFunction& gx, GridFunction& gy, const Point delta, 
		const Real deltaT, const Real Re)
{
	for(uint d=0; d<DIMENSIONS; d++)
	{
		auto Fxx = 
			Derivatives::getDerivative(
					domain.getVeolcity()[d],
					domain.getDelta(),
					4);
		auto Fyy = 
			Derivatives::getDerivative(
					domain.getVeolcity()[d],
					domain.getDelta(),
					5);
		auto FF_df = 
			Derivatives::getDerivative(
					domain.getVeolcity()[d],
					domain.getVeolcity()[d],
					domain.getDelta(),
					1+d, 
					1+d,
					1+d);
		auto FG_dg = 
			Derivatives::getDerivative(
					domain.getVeolcity()[d],
					domain.getVeolcity()[d],
					domain.getDelta(),
					1+d,
					1+((d+2)%3), /* TODO ! */
					1+((d+2)%3));

		/* reset inner (boundary not needed here) */
		domain.getPreliminaryVeolcity()[d].setGridFunction
			(domain.getBeginInnerDomains(),domain.getEndInnerDomain()[d],0.0);

		/* the formula: */
		for(int i=domain.getBeginInnerDomains()[0]; 
				i<= domain.getEndInnerDomain()[d][0]; i++)
			for(int j=domain.getBeginInnerDomains()[1]; 
					j<= domain.getEndInnerDomain()[d][1]; j++)
		{
			domain.getPreliminaryVeolcity()[d](i, j) =
				domain.getVeolcity()[d](i, j) + 
				(deltaT/Re)*( Fxx(i,j) + Fyy(i,j) )
				- FF_df(i,j) - FG_dg(i,j)
				+ domain.gx(
						Point(
						Real(i)/domain.getEndInnerDomain()[d][0]+2, 
						Real(j)/domain.getEndInnerDomain()[d][1]+2));
		}
	}
}

void Computation::computeRighthandSide(
		Domain& domain,
		//GridFunction& rhs, GridFunction& f, GridFunction& g, const Point delta, 
		const Real deltaT)
{
	auto Fxf = Derivatives::getDerivative(
			domain.F(),
			domain.getDelta(),
			Derivatives::Direction::xf);
	auto Gyf = Derivatives::getDerivative(
			domain.G(),
			domain.getDelta(),
			Derivatives::Direction::yf);

	for(int i=domain.getBeginInnerDomains()[0]; 
			i<= domain.getEndInnerDomain()[3][0]; i++)
		for(int j=domain.getBeginInnerDomains()[1]; 
				j<= domain.getEndInnerDomain()[3][1]; j++)
	{
		domain.p()(i,j) = ( Fxf(i,j) + Gyf(i,j) )/deltaT;
	}
}

void Computation::computeNewVelocities(
		Domain& domain,
		//GridFunction& u, GridFunction& v,
		//GridFunction& f, GridFunction& g,
		//GridFunction& p, const Point delta, 
		const Real deltaT)
{

	for(uint d=0; d<DIMENSIONS; d++)
	{
		auto Pdb = Derivatives::getDerivative(
				domain.p(),
				domain.getDelta(),
				1+d);

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

