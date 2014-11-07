#include "Computation.hpp"

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
			domain.getDelta()[0]*domain.getDelta()[0] *
			domain.getDelta()[1]*domain.getDelta()[1]*Re
			/
			(2.0*(domain.getDelta()[0]*domain.getDelta()[0]+domain.getDelta()[1]*domain.getDelta()[1])),

			std::fmin(
				domain.getDelta()[0]/std::abs(uMax),
				domain.getDelta()[1]/std::abs(vMax)
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

