
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

#include <chrono>


int main(int argc, char** argv)
{
	argc = argc*(**argv); // just to get rid of some warnings */

	/* time measurement variables */
	std::chrono::steady_clock::time_point t_frame_start, t_frame_end;
	std::chrono::steady_clock::time_point t_sor_start, t_sor_end;
	std::chrono::duration<double> time_span;
	double t_sor_avg = 0.0;
	double t_frame_avg = 0.0;

	/* init simulation parameters */
	SimulationParameters simparam("inputvals");

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
	log_info("[P%i] number of tasks: %i - going for a %ix%i processors grid - global inner: [(%i,%i),(%i,%i)]", 
		communication.getRank(),
		communication.getProcsCount(),
		communication.getProcsGridDim().i,
		communication.getProcsGridDim().j,
		communication.getGlobalInnerRange().begin.i,
		communication.getGlobalInnerRange().begin.j,
		communication.getGlobalInnerRange().end.i,
		communication.getGlobalInnerRange().end.j );

	log_info("[P%i] procs-grid position: (%i,%i) - local inner: [(%i,%i),(%i,%i)] - competences : %s%s%s%s",
		communication.getRank(),
		communication.getProcsGridPosition().i,
		communication.getProcsGridPosition().j,
		communication.getLocalInnerRange().begin.i,
		communication.getLocalInnerRange().begin.j,
		communication.getLocalInnerRange().end.i,
		communication.getLocalInnerRange().end.j,
		((communication.getBoundaryCompetence().Up) ?
		("Up ") : ("")),
		((communication.getBoundaryCompetence().Right) ?
		("Right ") : ("")),
		((communication.getBoundaryCompetence().Down) ?
		("Down ") : ("")),
		((communication.getBoundaryCompetence().Left) ?
		("Left") : (""))
		);

	Dimension local_dim = communication.getLocalDimension();

	/* init domain, which holds all grids and knows about their dimensions */
	Domain domain(local_dim, delta,
		/* init boundary, for this we need the inner range of the local process
		 * w.r.t. the range of the global domain. and the competences. */
		Boundary(communication.getLocalInnerRange(),
			communication.getBoundaryCompetence()),
		/* outer forces */
		simparam.gx, simparam.gy,
		/* initial grid values */
		simparam.ui, simparam.vi, simparam.pi,
		/* color pattern */
		communication.getFirstCellColor());

	log_info("[P%i] range p=(%i,%i), firstColor=%s, subRangesCount: p=%lu, u=%lu, v=%lu",
		communication.getRank(),
		domain.getWholeInnerRange().end.i, domain.getWholeInnerRange().end.j,
		(domain.getDomainFirstCellColor() == Color::Red) ? ("Red") : ("Black"),
		domain.getInnerRangeP().size(),
		domain.getInnerRangeU().size(),
		domain.getInnerRangeV().size() );

	/* next: omega and time parameters */
	Real h = 1.0 / simparam.iMax;// std::fmin(simparam.xLength, simparam.yLength);
	// concerning h, see: http://userpages.umbc.edu/~gobbert/papers/YangGobbertAML2007.pdf
	if (simparam.xLength == simparam.yLength) 
		simparam.omg = 2.0 /(1.0 + sin(M_PI*(h))); 

	Real t = 0.0, dt, res;
	int it = 0, step=0;

	/* write initial state of velocities and pressure */
	VTKOutput vtkoutput(domain, "out", communication);
	vtkoutput.writeVTKFile();
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

			//Solver::SORCycle(domain.p(), domain.rhs(), delta, 
			//		domain.getInnerRangeP(), simparam.omg);
			Solver::SORCycleRedBlack(domain, simparam.omg, Color::Red);
			communication.exchangeGridBoundaryValues(domain,Communication::Handle::Pressure);
			Solver::SORCycleRedBlack(domain, simparam.omg, Color::Black);
			communication.exchangeGridBoundaryValues(domain,Communication::Handle::Pressure);

			//res = Solver::computeResidual(
			//		domain.p(), domain.rhs(), delta, 
			//		domain.getInnerRangeP(), global_dim);
			res = Solver::computeSquaredResidual(
					domain.p(), domain.rhs(), delta, domain.getInnerRangeP(), global_dim);
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
			printf("\r[INFO] - Round %i: t/tmax=%f | dt=%f | Solver: it=%i (max:%i) | res=%f (max:%f) %s",
				step, t / simparam.tEnd, dt, it, simparam.iterMax, res, simparam.eps,
				((t < simparam.tEnd)) ? ("") : ("\n"));
			fflush( stdout );
		}
		it = 0;

		Computation::computeNewVelocities(domain, dt);

		nextVTKWrite += dt;
		if(nextVTKWrite > simparam.deltaVec)
		{
			nextVTKWrite = 0.0;
			vtkoutput.writeVTKFile();
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
