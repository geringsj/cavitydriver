
#include "src/VTKOutput.hpp"
#include "src/SimulationParameters.hpp"
#include "src/Domain.hpp"
#include "src/Computation.hpp"
#include "src/Solver.hpp"
#include "src/Debug.hpp"
#include "src/Communication.hpp"
//#include "src/CavityPainter.hpp"

#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>

#include <fstream>

#include <chrono>


int main(int argc, char** argv)
{
	argc = argc*(1+0*(**argv)); /* just to get rid of some warnings */

	/* time measurement variables */
	std::chrono::steady_clock::time_point t_frame_start, t_frame_end;
	std::chrono::steady_clock::time_point t_sor_start, t_sor_end;
	std::chrono::duration<double> time_span;
	double t_sor_avg = 0.0;
	double t_frame_avg = 0.0;

	/* init simulation parameters */
	std::string inputvals("inputvals");
	if(argc == 2)
		inputvals = std::string(argv[1]);
	SimulationParameters simparam( inputvals );
	//simparam.writeToSTDOUT();

	/* init problem dimensions and grid spacing delta */
	Dimension global_dim;
	global_dim.i = simparam.iMax;
	global_dim.j = simparam.jMax;
	Delta delta;
	delta.x = simparam.xLength / simparam.iMax;
	delta.y = simparam.yLength / simparam.jMax;

	/* init communication of processes */
	Communication communication = Communication(global_dim);

	if(communication.getRank() == 0)
	log_info("[P%i] loading inputvals from file \"%s\"",
		communication.getRank(), inputvals.c_str());
	if(communication.getRank() == 0)
	log_info("[P%i] tasks count: %i | doing %ix%i procs-grid | global inner: [(%i,%i),(%i,%i)]", 
		communication.getRank(),
		communication.getProcsCount(),
		communication.getProcsGridDim().i,
		communication.getProcsGridDim().j,
		communication.getGlobalInnerRange().begin.i,
		communication.getGlobalInnerRange().begin.j,
		communication.getGlobalInnerRange().end.i,
		communication.getGlobalInnerRange().end.j );

	log_info("[P%i] procs-grid position: (%i,%i) | local inner: [(%i,%i),(%i,%i)] | competences : %s%s%s%s",
		communication.getRank(),
		communication.getProcsGridPosition().i,
		communication.getProcsGridPosition().j,
		communication.getLocalInnerRange().begin.i,
		communication.getLocalInnerRange().begin.j,
		communication.getLocalInnerRange().end.i,
		communication.getLocalInnerRange().end.j,
		((communication.getBoundaryCompetence().Up) ? ("Up ") : ("")),
		((communication.getBoundaryCompetence().Right) ? ("Right ") : ("")),
		((communication.getBoundaryCompetence().Down) ? ("Down ") : ("")),
		((communication.getBoundaryCompetence().Left) ? ("Left") : ("")) );

	Dimension local_dim = communication.getLocalDimension();

	/* init domain, which holds all grids and knows about their dimensions */
	Domain domain(local_dim, delta,
		/* init boundary, for this we need the inner range of the local process
		 * w.r.t. the range of the global domain. and the competences. */
		Boundary(
			communication.getLocalInnerRange(), 
			communication.getBoundaryCompetence(), 
			simparam.boundary_conditions),
		/* outer forces */
		simparam.gx, simparam.gy,
		/* initial grid values */
		simparam.ui, simparam.vi, simparam.pi,
		/* color pattern */
		communication.getFirstCellColor());

	Real global_fluidCells = 
		communication.getGlobalFluidCellsCount(domain.getFluidCellsCount());
	//debug("FluidCells: %f", global_fluidCells);

	//simparam.writeSettingsFile("inputvals_"+simparam.name);

	log_info("[P%i] range p=(%i,%i), firstColor=%s, subRangesCount: p=%lu, u=%lu, v=%lu",
		communication.getRank(),
		domain.getWholeInnerRange().end.i, domain.getWholeInnerRange().end.j,
		(domain.getDomainFirstCellColor() == Color::Red) ? ("Red") : ("Black"),
		domain.getInnerRangeP().size(),
		domain.getInnerRangeU().size(),
		domain.getInnerRangeV().size() );

	debug("ranges P:");
	for(auto& r : domain.getInnerRangeP())	
		debug("[(%i,%i)(%i,%i)]",
		r.begin.i,
		r.begin.j,
		r.end.i,
		r.end.j );

	debug("ranges U:");
	for(auto& r : domain.getInnerRangeU())	
		debug("[(%i,%i)(%i,%i)]",
		r.begin.i,
		r.begin.j,
		r.end.i,
		r.end.j );

	debug("ranges V:");
	for(auto& r : domain.getInnerRangeV())	
		debug("[(%i,%i)(%i,%i)]",
		r.begin.i,
		r.begin.j,
		r.end.i,
		r.end.j );

	/* next: omega and time parameters */
	Real h = 1.0 / std::min(simparam.iMax, simparam.jMax); // be careful
	// concerning h, see: http://userpages.umbc.edu/~gobbert/papers/YangGobbertAML2007.pdf
	simparam.omg = 2.0 /(1.0 + sin(M_PI*(h))); 

	Real t = 0.0, dt, res;
	int it = 0, step=0;

	/* write initial state of velocities and pressure */
	VTKOutput vtkoutput(domain, "out", communication);
	vtkoutput.writeVTKFile(0.0);
	Real nextVTKWrite = 0.0;

	std::chrono::steady_clock::time_point t_start = std::chrono::steady_clock::now();

	/* main loop */
	while(t < simparam.tEnd)
	{
		t_frame_start = std::chrono::steady_clock::now();

		/* the magic starts here */
		//dt = Computation::computeTimestep(domain, simparam.tau, simparam.re); 
		Delta maxVelocities = 
			Delta(domain.u().getMaxValueGridFunction(), domain.v().getMaxValueGridFunction());
		maxVelocities = communication.getGlobalMaxVelocities(maxVelocities);

		dt = Computation::computeTimestepFromMaxVelocities
			(maxVelocities, domain.getDelta(), simparam.tau, simparam.re);
		t += dt;

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
#ifdef WITHMPI
			Solver::SORCycleRedBlack(domain, simparam.omg, Color::Red);
			communication.exchangeGridBoundaryValues(domain,Communication::Handle::Pressure);
			Solver::SORCycleRedBlack(domain, simparam.omg, Color::Black);
			communication.exchangeGridBoundaryValues(domain,Communication::Handle::Pressure);
#else
			Solver::SORCycle(domain.p(), domain.rhs(), delta,
					domain.getInnerRangeP(), simparam.omg);
#endif

			res = Solver::computeSquaredResidual(
				domain.p(), domain.rhs(), delta, domain.getInnerRangeP(), global_fluidCells);
			res = communication.getGlobalResidual(res);
			it++;
		} while (communication.checkGlobalFinishSOR
			(it < simparam.iterMax && res > simparam.eps));
		t_sor_end = std::chrono::steady_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t_sor_end-t_sor_start);
		t_sor_avg += time_span.count();
		//log_info("SOR solver time: %f seconds",time_span.count());

		/* write new line of current statistics to stdout */
		if(communication.getRank() == 0)
		{
			printf("[INFO] - Round %i: %.2f%% | dt=%f | Solver: it=%i / res=%f%s\r",
				step, (t / simparam.tEnd)*100., dt, it, res, 
				((t < simparam.tEnd)) ? ("") : ("\n"));
			fflush( stdout );
		}
		it = 0;

		Computation::computeNewVelocities(domain, dt);

		nextVTKWrite += dt;
		if(nextVTKWrite > simparam.deltaVec)
		{
			nextVTKWrite = 0.0;
			vtkoutput.writeVTKFile(dt);
		}
		step++;

		t_frame_end = std::chrono::steady_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t_frame_end-t_frame_start);
		t_frame_avg += time_span.count();
		//log_info("Overall frame time: %f seconds",time_span.count());
	}

	std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t_end-t_start);

	/* output time: overall, per frame and pressure computation per frame */
	log_info("[P%i] Overall time: %fs | avg. frame time: %fs | avg. SOR time: %fs",
		communication.getRank(),
		time_span.count(),
		t_frame_avg / (double)(step-1),
		t_sor_avg / (double)(step-1) );

	/* end of magic */
	return 0;
}
