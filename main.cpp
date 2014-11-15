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
	delta.x = simparam.xLength / simparam.iMax;
	delta.y = simparam.yLength / simparam.jMax;

	Domain domain(dim, delta,
		/* U */std::bind(u, std::placeholders::_1, std::placeholders::_2, 
			std::placeholders::_3, std::ref(simparam)),
		/* V */std::bind(v, std::placeholders::_1, std::placeholders::_2, 
			std::placeholders::_3, std::ref(simparam)),
		/* W==0 */[](Index i, GridFunction& gf, Dimension dim)
			{return 0.0*i.i*gf.getGridDimension().i*dim.i; }, 
		/* P==0 */[&simparam](Index i, GridFunction& gf, Dimension dim)
			{return simparam.pi + 0.0*(i.i*gf.getGridDimension().i*dim.i); });

	/* main loop */
	Real t = 0.0, dt = simparam.deltaT, res;
	//http://userpages.umbc.edu/~gobbert/papers/YangGobbertAML2007.pdf
	//Real h = 1.0 / simparam.iMax;// std::fmin(simparam.xLength, simparam.yLength);
	//if (simparam.xLength == simparam.yLength) simparam.omg = 2.0 /(1.0 + sin(M_PI*(h)));
	//debug("omega: %f", simparam.omg);

	int it = 0, step=0;
	io.writeVTKFile(domain.getDimension(), domain.u(), domain.v(), domain.p(), delta, step);
	step++;

	while (t < simparam.tEnd)
	{
		log_info("- Round %i", step);
		if (m_log_info) log_info("= = = = = = = = = = = = = = =");
		if (m_log_info) log_info("Grids with borders: ");
		if (m_log_info) { log_info("U:"); domain.u().printSTDOUT(); }
		if (m_log_info) { log_info("F:"); domain.F().printSTDOUT(); }
		if (m_log_info) { log_info("V:"); domain.v().printSTDOUT();	}
		if (m_log_info) { log_info("G:"); domain.G().printSTDOUT();	}
		if (m_log_info) { log_info("P:"); domain.p().printSTDOUT();	}
		if (m_log_info) { log_info("RHS:"); domain.rhs().printSTDOUT(); }
		if (m_log_info) log_info("= = = = = = = = = = = = = = =");

		dt = Computation::computeTimestep(domain, simparam.tau, simparam.re); 
		t += dt;
		log_info("-- dt=%f | t/tmax=%f", dt, t / simparam.tEnd);

		domain.setVelocitiesBoundaries();
		domain.setPreliminaryVelocitiesBoundaries();

		Computation::computePreliminaryVelocities(domain, dt, simparam.re, simparam.alpha);
		Computation::computeRighthandSide(domain, dt);

		do
		{
			domain.setPressureBoundaries();

			Solver::SORCycle(
					domain.p(), domain.rhs(), delta, 
					domain.getBeginInnerDomains(), domain.getEndInnerDomainP(), simparam.omg);

			res = Solver::computeResidual(
					domain.p(), domain.rhs(), delta, 
					domain.getBeginInnerDomains(), domain.getEndInnerDomainP());
			it++;
		} while (it < simparam.iterMax && res > simparam.eps); 
		//if (m_log_info) 
			log_info("-- Solver done: it=%i (max:%i)| res=%f (max:%f)",
					it, simparam.iterMax, res, simparam.eps);
		it = 0;

		Computation::computeNewVelocities(domain, dt);
		//domain.setVelocitiesBoundaries();
		//domain.setPressureBoundaries();
		io.writeVTKFile(
				domain.getDimension(), domain.u(), domain.v(), domain.p(), delta, step);
		step++;
	}
	return 0;
}
