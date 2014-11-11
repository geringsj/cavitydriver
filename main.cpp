#include "src/IO.hpp"
#include "src/Domain.hpp"
#include "src/Computation.hpp"
#include "src/Solver.hpp"

// default values for the output directory and the settings file
const char *output = "out";
const char *settings = "inputvals";

/**
 * Print the help text in order to show what
 * arguments can be processed.
 */
void dieSynopsis() {
	std::cerr << "Synopsis: IO "
		<< "--out \"output directory\" --set \"path to the settings file\" ";
	exit(-1);
}

/**
 * Parse the input arguments. If there are no arguments
 * the default values are used.
 * @param the length of the input
 * @param the arguments
 */
void parseArguments(int argc, char** argv)
{
	--argc; ++argv; // We don't care about program name;
	while (argc > 0 && argv[0][0] == '-')
	{
		if (argv[0] == (std::string) "--help" || argv[0] == (std::string) "-h")
		{
			dieSynopsis();
		}
		else if (argv[0] == (std::string)"--out")
		{
			if (argv[1] != NULL && argv[1] != (std::string)"--set")
			{
				output = argv[1];
				--argc; ++argv;
				--argc; ++argv;
			}
			else --argc; ++argv;
		}
		else if (argv[0] == (std::string)"--set")
		{
			if (argv[1] != NULL)
			{
				settings = argv[1];
				--argc; ++argv;
				--argc; ++argv;
			}
			else --argc; ++argv;
		}
	}
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
	parseArguments(argc, argv);
	std::cout << settings << " " << output << std::endl;
	IO io(settings, output);
	//IO io((char*)"inputvals", (char*)"out");
	io.writeSimParamToSTDOUT();

	Dimension dim;
	dim.i = io.getIMax();
	dim.j = io.getJMax();
	Delta delta;
	delta.x = io.getXLength();
	delta.y = io.getYLength();
	Domain domain(dim, delta,
		std::bind(u, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::ref(io)),
		std::bind(v, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::ref(io)),
		[](Index i, GridFunction& gf, Dimension dim){return 0.0; }, 
		[&io](Index i, GridFunction& gf, Dimension dim){return io.getPi(); });

	/* main loop */
	/* todo:
	 * Calculate omega for grid dimension
	 * 
	 */
	Real t = 0.0;
	Real dt;
	Real res;
	int it = 0;
	while (t < io.getTEnd())
	{
		dt = Computation::computeTimestep(domain, io.getTau(), io.getRe());
		t += dt;
		std::cout << "dt: " << dt << " , t/tmax: " << t / io.getTEnd() << std::endl;

		domain.setVelocitiesBoundaries();
		Computation::computeMomentumEquationsFGH(domain, dt, io.getRe());
		domain.setPreliminaryVelocitiesBoundaries();
		domain.setPressureBoundaries();
		Computation::computeRighthandSide(domain, dt);

		it = 0;
		do
		{
			res = Solver::computeResidual(domain.p(), domain.rhs(), delta, domain.getBeginInnerDomains(), domain.getEndInnerDomainP());
			Solver::SORCycle(domain.p(), domain.rhs(), delta, domain.getBeginInnerDomains(), domain.getEndInnerDomainP(), io.getOmg());
			it++;
		} while (it < io.getIterMax() && res > io.getEps());

		Computation::computeNewVelocities(domain, dt);
		domain.setVelocitiesBoundaries();
		io.writeVTKFile(domain.getDimension(), domain.u(), domain.v(), domain.p(), delta, t);
	}

	return 0;
}
