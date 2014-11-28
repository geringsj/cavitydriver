#include "src/IO.hpp"
#include "src/Domain.hpp"
#include "src/Computation.hpp"
#include "src/Solver.hpp"
#include "src/Debug.hpp"
/* communication test */
#include "src/Communication.hpp"

#define _USE_MATH_DEFINES
#include <math.h>

#include <chrono>

Real u_bounds(Index index, GridFunction& gf, Dimension dim, SimParams& simparam)
{
	Real value = simparam.ui;
	if (index.i == 0) value = 0.0;
	if (index.i == dim.i) value = 0.0;
	if (index.j == 0) value = -gf(index.i, index.j + 1);
	if (index.j == dim.j) value = 2.0 - gf(index.i, index.j - 1);
	return value;
}

Real v_bounds(Index index, GridFunction& gf, Dimension dim, SimParams& simparam)
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
	/* time measurement variables */
	std::chrono::steady_clock::time_point t_frame_start, t_frame_end;
	std::chrono::steady_clock::time_point t_sor_start, t_sor_end;
	std::chrono::duration<double> time_span;
	double t_sor_avg = 0.0;
	double t_frame_avg = 0.0;

	/* init IO/parameters */
	int io_argc = argc;			
	char** io_argv = argv;
	/** 
	 * The IO parser "eats" the arguments leaving them empty after the constructor.
	 * We need these arguments for the MPI_Init call so we give the IO a copy of
	 * the arguments.
	 */
	IO io(io_argc, io_argv);
	SimParams simparam = io.readInputfile();
	simparam.writeSimParamsToSTDOUT();

	/* init problem dimensions and grid spacing delta */
	Dimension dim;
	dim.i = simparam.iMax;
	dim.j = simparam.jMax;
	Delta delta;
	delta.x = simparam.xLength / simparam.iMax;
	delta.y = simparam.yLength / simparam.jMax;

	/* communication test */
	Communication comm = Communication(dim, argc, argv);

	/* init domain, which holds all grids and knows about their dimensions */
	Domain domain(dim, delta,
		/* 
		 * you can ignore the next few lines. 
		 * they set the boundary conditions on the domain by passing boundary functions */
		/* U */std::bind(u_bounds, std::placeholders::_1, std::placeholders::_2, 
			std::placeholders::_3, std::ref(simparam)),
		/* V */std::bind(v_bounds, std::placeholders::_1, std::placeholders::_2, 
			std::placeholders::_3, std::ref(simparam)),
		/* W==0 */[](Index i, GridFunction& gf, Dimension dim)
			{return 0.0*i.i*gf.getGridDimension().i*dim.i; }, 
		/* P==0 */[&simparam](Index i, GridFunction& gf, Dimension dim)
			{return simparam.pi + 0.0*(i.i*gf.getGridDimension().i*dim.i); });
	log_info("My Domain starts at Color: %s", 
			(domain.getDomainFirstCellColor()==Color::Red)?("Red"):("Black"));

	/* next: omega and time parameters */
	Real h = 1.0 / simparam.iMax;// std::fmin(simparam.xLength, simparam.yLength);
	// concerning h, see: http://userpages.umbc.edu/~gobbert/papers/YangGobbertAML2007.pdf
	if (simparam.xLength == simparam.yLength) 
		simparam.omg = 2.0 /(1.0 + sin(M_PI*(h))); //debug("omega: %f", simparam.omg);

	Real t = 0.0, dt = simparam.deltaT, res;
	int it = 0, step=0;

	/* write initial state of velocities and pressure */
	io.writeVTKFile(domain.getDimension(), domain.u(), domain.v(), domain.p(), delta, step);
	step++;
	Real nextWrite = 0.0;

	/* main loop */
	while (t < simparam.tEnd)
	{
		t_frame_start = std::chrono::steady_clock::now();

		/* the following block is for debugging purposes */
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


		/* the magic starts here */
		dt = Computation::computeTimestep(domain, simparam.tau, simparam.re); 
		t += dt;
		log_info("-- dt=%f | t/tmax=%f", dt, t / simparam.tEnd);

		domain.setVelocitiesBoundaries();
		Computation::computePreliminaryVelocities(domain, dt, simparam.re, simparam.alpha);
		domain.setPreliminaryVelocitiesBoundaries();

		Computation::computeRighthandSide(domain, dt);

		t_sor_start = std::chrono::steady_clock::now();
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
		t_sor_end = std::chrono::steady_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t_sor_end-t_sor_start);
		t_sor_avg += time_span.count();
		//log_info("SOR solver time: %f seconds",time_span.count());

		if (m_log_info) 
			log_info("-- Solver done: it=%i (max:%i)| res=%f (max:%f)",
					it, simparam.iterMax, res, simparam.eps);
		it = 0;

		Computation::computeNewVelocities(domain, dt);

		nextWrite += dt;
		if(nextWrite > simparam.deltaVec)
		{
			io.writeVTKFile(
				domain.getDimension(), domain.u(), domain.v(), domain.p(), delta, step);
			nextWrite = 0;
		}
		step++;

		t_frame_end = std::chrono::steady_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t_frame_end-t_frame_start);
		t_frame_avg += time_span.count();
		//log_info("Overall frame time: %f seconds",time_span.count());
	}

	/* output average time per frame and pressure computation per frame */
	log_info("Average frame time: %f seconds",t_frame_avg / (double)step-1);
	log_info("Average SOR time: %f seconds",t_sor_avg / (double)step-1);

	/* end of magic */
	return 0;
}
