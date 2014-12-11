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

void mainLoop()
{
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
	IO io(argc, argv);
	SimParams simparam = io.readInputfile();
	//simparam.writeSimParamsToSTDOUT();

	/* init problem dimensions and grid spacing delta */
	Dimension global_dim;
	global_dim.i = simparam.iMax;
	global_dim.j = simparam.jMax;
	Delta delta;
	delta.x = simparam.xLength / simparam.iMax;
	delta.y = simparam.yLength / simparam.jMax;

	/* This is just for testing, delete it if we are done with testing. */
	Communication communication = Communication(global_dim);
	//if(!communication.getProcessValid())
	//	log_info("I am a useless process. ");
	//else
	//	log_info("I am a good process. Yaay.");

	Dimension local_dim = communication.getLocalDimensions();

	/* init domain, which holds all grids and knows about their dimensions */
	Domain domain(local_dim, delta,
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
			{return simparam.pi + 0.0*(i.i*gf.getGridDimension().i*dim.i); },
		/* outer forces */
			simparam.gx, simparam.gy, 0.0,
		/* boundaries and color pattern */
			communication.getBoundaryCompetence(), communication.getFirstCellColor());

	log_info("process %i has competence over boundaries: %s %s %s %s",
			communication.getRank(),
			((communication.getBoundaryCompetence().Up) ?
			("Up ") : ("")),
			((communication.getBoundaryCompetence().Right) ?
			("Right ") : ("")),
			((communication.getBoundaryCompetence().Down) ?
			("Down ") : ("")),
			((communication.getBoundaryCompetence().Left) ?
			("Left") : (""))
			);
	log_info("process %i has end indices: p=(%i,%i), u=(%i,%i), v=(%i,%i), firstColor=%s",
			communication.getRank(),
			domain.getEndInnerDomainP().i, domain.getEndInnerDomainP().j,
			domain.getEndInnerDomainU().i, domain.getEndInnerDomainU().j,
			domain.getEndInnerDomainV().i, domain.getEndInnerDomainV().j,
		(domain.getDomainFirstCellColor() == Color::Red) ? ("Red") : ("Black"));

	/* next: omega and time parameters */
	Real h = 1.0 / simparam.iMax;// std::fmin(simparam.xLength, simparam.yLength);
	// concerning h, see: http://userpages.umbc.edu/~gobbert/papers/YangGobbertAML2007.pdf
	if (simparam.xLength == simparam.yLength) 
		simparam.omg = 2.0 /(1.0 + sin(M_PI*(h))); 

	Real t = 0.0, dt = simparam.deltaT, res;
	int it = 0, step=0;

	/* write initial state of velocities and pressure */
	io.writeVTKSlaveFile(domain,step,communication,simparam);
	if(communication.getRank() == 0)
		io.writeVTKMasterFile(global_dim,step,communication);
	step++;
	Real nextWrite = 0.0;

	std::chrono::steady_clock::time_point t_start = std::chrono::steady_clock::now();

	/* main loop */
	while (t < simparam.tEnd)
	{
		t_frame_start = std::chrono::steady_clock::now();

		if(communication.getRank() == 0)
			log_info("- Round %i", step);

		/* the magic starts here */
		//dt = Computation::computeTimestep(domain, simparam.tau, simparam.re); 
		Delta maxVelocities = 
			Delta(domain.u().getMaxValueGridFunction(), domain.v().getMaxValueGridFunction());
		maxVelocities = communication.getGlobalMaxVelocities(maxVelocities);

		dt = Computation::computeTimestepFromMaxVelocities
			(maxVelocities, domain.getDelta(), simparam.tau, simparam.re);
		t += dt;
		if(communication.getRank() == 0)
			log_info("-- dt=%f | t/tmax=%f", dt, t / simparam.tEnd);

		domain.setVelocitiesBoundaries();
		communication.exchangeGridBoundaryValues(domain, Communication::Handle::Velocities);

		Computation::computePreliminaryVelocities(domain, dt, simparam.re, simparam.alpha);

		domain.setPreliminaryVelocitiesBoundaries();
		communication.exchangeGridBoundaryValues(domain, Communication::Handle::PreliminaryVelocities);

		Computation::computeRighthandSide(domain, dt);

		t_sor_start = std::chrono::steady_clock::now();
		do
		{
			domain.setPressureBoundaries();

			//Solver::SORCycle(
			//		domain.p(), domain.rhs(), delta, 
			//		domain.getBeginInnerDomains(), domain.getEndInnerDomainP(), simparam.omg);

			Solver::SORCycleRedBlack(domain, simparam.omg, Color::Red);
			communication.exchangeGridBoundaryValues(domain, Communication::Handle::Pressure);

			Solver::SORCycleRedBlack(domain, simparam.omg, Color::Black);
			communication.exchangeGridBoundaryValues(domain, Communication::Handle::Pressure);

			//res = Solver::computeResidual(
			//		domain.p(), domain.rhs(), delta, 
			//		domain.getBeginInnerDomains(), domain.getEndInnerDomainP(), global_dim);
			res = Solver::computeSquaredResidual(
					domain.p(), domain.rhs(), delta, 
					domain.getBeginInnerDomains(), domain.getEndInnerDomainP(), global_dim);
			res = communication.getGlobalResidual(res);
			it++;
		} while (communication.checkGlobalFinishSOR
				(it < simparam.iterMax && res > simparam.eps));
		t_sor_end = std::chrono::steady_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t_sor_end-t_sor_start);
		t_sor_avg += time_span.count();
		//log_info("SOR solver time: %f seconds",time_span.count());

		if(communication.getRank() == 0)
			log_info("-- Solver done: it=%i (max:%i)| res=%f (max:%f)",
						it, simparam.iterMax, res, simparam.eps);
		it = 0;

		Computation::computeNewVelocities(domain, dt);

		nextWrite += dt;
		if(nextWrite > simparam.deltaVec)
		{
			nextWrite = 0;

			//io.writeVTKFile(
			//	domain, step);

			// write all slave files
			io.writeVTKSlaveFile(domain,step,communication,simparam);

			// only rank 0 writes a master file
			if(communication.getRank() == 0)
				io.writeVTKMasterFile(global_dim,step,communication);
		}
		step++;

		t_frame_end = std::chrono::steady_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t_frame_end-t_frame_start);
		t_frame_avg += time_span.count();
		//log_info("Overall frame time: %f seconds",time_span.count());
	}

	std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t_end-t_start);
	log_info("Overall time: %f seconds",time_span);

	/* output average time per frame and pressure computation per frame */
	log_info("Average frame time: %f seconds",t_frame_avg / (double)step-1);
	log_info("Average SOR time: %f seconds",t_sor_avg / (double)step-1);

	/* end of magic */
	return 0;
}
