#include "src/IO.hpp"
#include "src/Domain.hpp"
#include "src/Computation.hpp"
#include "src/Solver.hpp"
#include "src/Debug.hpp"

#define _USE_MATH_DEFINES
#include <math.h>

Real u(Index index, GridFunction& gf, Dimension dim, Simparam& simparam)
{
	Real value = simparam.ui;
	if (index.i == 0) value = 0.0;
	if (index.i == dim.i) value = 0.0;
	if (index.j == 0) value = -gf(index.i, index.j + 1);
	if (index.j == dim.j) value = 2.0 - gf(index.i, index.j - 1);
	return value;
}

Real v(Index index, GridFunction& gf, Dimension dim, Simparam& simparam)
{
	Real value = simparam.vi;
	if (index.i == 0) value = -gf(index.i+1, index.j);
	if (index.i == dim.i) value = -gf(index.i-1, index.j);
	if (index.j == 0) value = 0.0;
	if (index.j == dim.j) value = 0.0;
	return value;
}

int main(int argc, char** argv)
{
	IO io(argc, argv);
	Simparam simparam = io.readInputfile();
	simparam.writeSimParamToSTDOUT();

	Dimension dim;
	dim.i = simparam.iMax;
	dim.j = simparam.jMax;
	Delta delta;
	delta.x = simparam.xLength;
	delta.y = simparam.yLength;
	Domain domain(dim, delta,
		std::bind(u, std::placeholders::_1, std::placeholders::_2, 
			std::placeholders::_3, std::ref(simparam)),
		std::bind(v, std::placeholders::_1, std::placeholders::_2, 
			std::placeholders::_3, std::ref(simparam)),
		[](Index i, GridFunction& gf, Dimension dim)
			{return 0.0*i.i*gf.getGridDimension().i*dim.i; }, 
		[&simparam](Index i, GridFunction& gf, Dimension dim)
			{return simparam.pi + 0.0*(i.i*gf.getGridDimension().i*dim.i); });

	/* main loop */
	/* TODO:
	 * Calculate omega for grid dimension
	 */
	Real t = 0.0, dt = simparam.deltaT, res, h;
	h = std::fmin(simparam.xLength, simparam.yLength);
	if (simparam.xLength == simparam.yLength) simparam.omg = 2.0 /(1.0 + sin(0.9*M_PI/(h)));
	debug("omega: %f", simparam.omg);
	int it, step=0;
	while (t < simparam.tEnd)
	{
		dt = Computation::computeTimestep(domain, simparam.tau, simparam.re);
		t += dt;
		debug("dt: %f t/tmx: %f", dt, t / simparam.tEnd);
		domain.setVelocitiesBoundaries();
		Computation::computeMomentumEquationsFGH(domain, dt, simparam.re);

		domain.setPreliminaryVelocitiesBoundaries();
		domain.setPressureBoundaries();

		Computation::computeRighthandSide(domain, dt);

		it = 0;
		do
		{
			res = Solver::computeResidual(
					domain.p(), domain.rhs(), delta, 
					domain.getBeginInnerDomains(), domain.getEndInnerDomainP());
			Solver::SORCycle(
					domain.p(), domain.rhs(), delta, 
					domain.getBeginInnerDomains(), domain.getEndInnerDomainP(), simparam.omg);
			it++;
		} while (it < simparam.iterMax && res > simparam.eps);

		Computation::computeNewVelocities(domain, dt);
		domain.setVelocitiesBoundaries();
		io.writeVTKFile(
				domain.getDimension(), domain.u(), domain.v(), domain.p(), delta, step);

		step++;
	}

	return 0;
}
