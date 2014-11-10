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
		const Real deltaT)
{
	for(uint d=0; d<DIMENSIONS; d++)
	{
		/* reset inner (boundary not needed here) */
		domain.getPreliminaryVeolcity()[d].setGridFunction
			(domain.getBeginInnerDomains(),domain.getEndInnerDomain()[d],0.0);

		/* add laplace */
		domain.getPreliminaryVeolcity()[d];
			domain.getVeolcity()[d]; 
		/* times 1/Re */
		/* minus mixed terms */
		/* plus external force */
		/* times deltaT */
		/* plus old direction value */
	}
}

void Computation::computeRighthandSide(
		Domain& domain,
		//GridFunction& rhs, GridFunction& f, GridFunction& g, const Point delta, 
		const Real deltaT)
{
}

void Computation::computeNewVelocities(
		Domain& domain,
		//GridFunction& u, GridFunction& v,
		//GridFunction& f, GridFunction& g,
		//GridFunction& p, const Point delta, 
		const Real deltaT)
{
}

