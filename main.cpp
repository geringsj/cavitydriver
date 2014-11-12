#include "src/IO.hpp"
#include "src/Domain.hpp"
#include "src/Computation.hpp"
#include "src/Solver.hpp"
#include "src/Debug.hpp"

#define FORALL(I,J,SI,EI) 
void printGrid(GridFunction& gf, Index SI, Index EI)
{
	for(int I=EI[1]+1; I>=SI[1]-1; I--)
	{
		for(int J=SI[0]-1; J<=EI[0]+1; J++)
		{
			std::cout << gf(J,I) << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

Real u(Index index, GridFunction& gf, Dimension dim, IO& io)
{
	Real value = io.getUi();
	if (index.i == 0) value = 0.0;
	if (index.i == dim.i) value = 0.0;
	if (index.j == 0) value = -gf(index.i, index.j + 1);
	if (index.j == dim.j) value = 2.0 - gf(index.i, index.j - 1);
	return value;
}

Real v(Index index, GridFunction& gf, Dimension dim, IO& io)
{
	Real value = io.getVi();
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
	io.writeSimParamToSTDOUT();

	Dimension dim;
	dim.i = simparam.iMax;
	dim.j = simparam.jMax;
	Delta delta;
	delta.x = simparam.xLength;
	delta.y = simparam.yLength;
	Domain domain(dim, delta,
		std::bind(u, std::placeholders::_1, std::placeholders::_2, 
			std::placeholders::_3, std::ref(io)),
		std::bind(v, std::placeholders::_1, std::placeholders::_2, 
			std::placeholders::_3, std::ref(io)),
		[](Index i, GridFunction& gf, Dimension dim)
			{return 0.0*i.i*gf.getGridDimension().i*dim.i; }, 
		[&simparam](Index i, GridFunction& gf, Dimension dim)
			{return simparam.pi + 0.0*(i.i*gf.getGridDimension().i*dim.i); });

	/* main loop */
	/* TODO:
	 * Calculate omega for grid dimension
	 */
	Real t = 0.0, dt = simparam.deltaT, res;
	int it, step=0;
	io.writeVTKFile(
			domain.getDimension(), domain.u(), domain.v(), domain.p(), delta, step);
	step++;
	while (t < simparam.tEnd)
	{
		log_info("Round %i, here are the Grids (with borders): ", step);
		log_info("U:");
		printGrid(domain.u(),
				domain.getBeginInnerDomains(),domain.getEndInnerDomainU());
		log_info("V:");
		printGrid(domain.v(),
				domain.getBeginInnerDomains(),domain.getEndInnerDomainV());
		log_info("P:");
		printGrid(domain.p(),
				domain.getBeginInnerDomains(),domain.getEndInnerDomainP());
		log_info("F:");
		printGrid(domain.F(),
				domain.getBeginInnerDomains(),domain.getEndInnerDomainU());
		log_info("G:");
		printGrid(domain.G(),
				domain.getBeginInnerDomains(),domain.getEndInnerDomainV());
		log_info("RHS:");
		printGrid(domain.rhs(),
				domain.getBeginInnerDomains(),domain.getEndInnerDomainP());
		std::cin.get();

		dt = Computation::computeTimestep(domain, simparam.tau, simparam.re);
		t += dt;
		std::cout << "dt: " << dt << " , t/tmax: " << t / simparam.tEnd << std::endl;

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
		domain.setPressureBoundaries();
		io.writeVTKFile(
				domain.getDimension(), domain.u(), domain.v(), domain.p(), delta, step);

		step++;
	}

	return 0;
}
