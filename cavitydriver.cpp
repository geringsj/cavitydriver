
#include "src/VTKOutput.hpp"
#include "src/SimulationParameters.hpp"
#include "src/Domain.hpp"
#include "src/Computation.hpp"
#include "src/Solver.hpp"
#include "src/Debug.hpp"
#include "src/Communication.hpp"

#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>

#include <fstream>

#include <chrono>

#ifdef WITHUNCER
#if defined(__linux)
#include "sys/stat.h"
#endif
#if defined(_WIN64)
#include <Windows.h>
#endif
#if defined(_WIN32)
#include <Windows.h>
#endif
void checkOutputPath()
{
	std::string filename;
	filename.append("./uncertainty/");

	// Test if the directory exits, if not create it.
	std::filebuf test_dir;
	test_dir.open(const_cast < char *>(filename.c_str()), std::ios::out);
	if (!test_dir.is_open())
	{
		// Directory doesn't exist.
#if defined(_WIN64)
		CreateDirectory(filename.c_str(), NULL);
#elif defined(_WIN32)
		CreateDirectory(filename.c_str(), NULL);
#elif defined(__linux)
		mkdir(filename.c_str(), 0700);
#endif
	}
}

void printUncertainty(Domain& domain, Real time, std::string name)
{
	std::string filename;
	filename.append("./uncertainty/");
	filename.append(name);
	filename.append(".txt");

	std::string lines;
	std::string line;
	std::ifstream read_uncer(filename);
	if (read_uncer.is_open())
	{
		while (getline(read_uncer, line))
		{
			lines.append(line);
			lines.append("\n");
		}
		read_uncer.close();
	}

	Real u_val_120_5 = domain.u()(120, 5);
	Real u_val_64_64 = domain.u()(64, 64);
	Real u_val_5_120 = domain.u()(5, 120);
	Real v_val_120_5 = domain.v()(120, 5);
	Real v_val_64_64 = domain.v()(64, 64);
	Real v_val_5_120 = domain.v()(5, 120);

	std::string values = "";
	values.append(std::to_string(u_val_120_5));
	values.append(",");
	values.append(std::to_string(v_val_120_5));
	values.append(" | ");
	values.append(std::to_string(u_val_64_64));
	values.append(",");
	values.append(std::to_string(v_val_64_64));
	values.append(" | ");
	values.append(std::to_string(u_val_5_120));
	values.append(",");
	values.append(std::to_string(v_val_5_120));
	values.append(" | ");
	values.append(std::to_string(time));

	std::ofstream write_uncer(filename);
	if (write_uncer.is_open())
	{
		write_uncer << lines;
		write_uncer << values;
		write_uncer.close();
	}
}
#endif

int main(int argc, char** argv)
{
	argc = argc*(1+0*(**argv)); /* just to get rid of some warnings */

	/* time measurement variables */
	std::chrono::steady_clock::time_point t_frame_start, t_frame_end;
	std::chrono::steady_clock::time_point t_sor_start, t_sor_end;
	std::chrono::duration<double> time_span;
	double t_sor_avg = 0.0, t_frame_avg = 0.0;

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

	/* TODO: move statistics somewhere else */
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

	Real global_fluidCellsCount = static_cast<Real>(
		communication.getGlobalFluidCellsCount(domain.getFluidCellsCount()) );

	log_info("[P%i] range p=(%i,%i), firstColor=%s, subRangesCount: p=%lu, u=%lu, v=%lu",
		communication.getRank(),
		domain.getWholeInnerRange().end.i, domain.getWholeInnerRange().end.j,
		(domain.getDomainFirstCellColor() == Color::Red) ? ("Red") : ("Black"),
		domain.getInnerRangeP().size(),
		domain.getInnerRangeU().size(),
		domain.getInnerRangeV().size() );

	/* next: omega and time parameters */
	Real h = 1.0 / std::min(simparam.iMax, simparam.jMax); // be careful
	// concerning h, see: http://userpages.umbc.edu/~gobbert/papers/YangGobbertAML2007.pdf
	simparam.omg = 2.0 /(1.0 + sin(M_PI*(h))); 

	Real t = 0.0, dt = 0.0, res = 42.;
	int it = 0, step=0;

	/* write initial state of velocities and pressure */
	VTKOutput vtkoutput(domain, "out", communication);
	vtkoutput.writeVTKFile(0.0); /* first vtk frame: all zero */
	Real nextVTKWrite = 0.0;

#ifdef WITHUNCER
	checkOutputPath();
	printUncertainty(domain, t, simparam.name);
#endif

	std::chrono::steady_clock::time_point t_start = std::chrono::steady_clock::now();

	/* main loop */
	while(t < simparam.tEnd)
	{
		t_frame_start = std::chrono::steady_clock::now();

		/* the magic starts here */

		/* to get consistent maxvalues from new velocities, 
		 * first set boundaries, then compute max values for timestep */
		domain.setVelocitiesBoundaries();
		communication.exchangeGridBoundaryValues
			(domain, Communication::Handle::Velocities);

		/* maybe write vtk */
		
		if ((nextVTKWrite += dt) > simparam.deltaVec)
		{
			vtkoutput.writeVTKFile(dt); nextVTKWrite = 0.0;

#ifdef WITHUNCER
			printUncertainty(domain, t, simparam.name);
#endif
		}

		//dt = Computation::computeTimestep(domain, simparam.tau, simparam.re);
		Delta maxVelocities(domain.u().getMaxValue(), domain.v().getMaxValue());
		maxVelocities = communication.getGlobalMaxVelocities(maxVelocities);

		dt = Computation::computeTimestepFromMaxVelocities
			(maxVelocities, domain.getDelta(), simparam.tau, simparam.re);
		t += dt; /* for status output */

		Computation::computePreliminaryVelocities(domain, dt, simparam.re, simparam.alpha);
		domain.setPreliminaryVelocitiesBoundaries();
		communication.exchangeGridBoundaryValues
			(domain, Communication::Handle::PreliminaryVelocities);

		Computation::computeRighthandSide(domain, dt);

		t_sor_start = std::chrono::steady_clock::now();
		do
		{
			domain.setPressureBoundaries();
#ifdef WITHMPI
			Solver::SORCycleRedBlack(domain, simparam.omg, Color::Red);
			communication.exchangeGridBoundaryValues
				(domain,Communication::Handle::Pressure);

			Solver::SORCycleRedBlack(domain, simparam.omg, Color::Black);
			communication.exchangeGridBoundaryValues
				(domain,Communication::Handle::Pressure);
#else
			Solver::SORCycle(domain.p(), domain.rhs(), delta,
					domain.getInnerRangeP(), simparam.omg);
#endif

			res = Solver::computeSquaredResidual(
				domain.p(), domain.rhs(), delta,
				domain.getInnerRangeP(), global_fluidCellsCount);
			res = communication.getGlobalResidual(res);
		} while (communication.checkGlobalFinishSOR
			(++it < simparam.iterMax && res > simparam.eps));

		t_sor_end = std::chrono::steady_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>
			(t_sor_end-t_sor_start); t_sor_avg += time_span.count();
		//log_info("SOR solver time: %f seconds",time_span.count());

		/* write new line of current statistics to stdout */
		if(communication.getRank() == 0)
		{
			printf("[INFO] - Round %i: %.2f%% | dt=%f | Solver: it=%i / res=%f%s\r",
			step,(t/simparam.tEnd)*100.,dt,it,res,(t<simparam.tEnd)?(""):("\n"));
			fflush( stdout );
		}
		it = 0; step++;

		domain.setPressureBoundaries();
		communication.exchangeGridBoundaryValues
			(domain,Communication::Handle::Pressure);
		Computation::computeNewVelocities(domain, dt);

		t_frame_end = std::chrono::steady_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>
			(t_frame_end-t_frame_start); t_frame_avg += time_span.count();
		//log_info("Overall frame time: %f seconds",time_span.count());
	}

	std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t_end-t_start);

	/* output time: overall, per frame and pressure computation per frame */
	log_info("[P%i] Overall time: %fs | avg. frame time: %fs | avg. SOR time: %fs",
		communication.getRank(), time_span.count(),
		t_frame_avg/(double)(step-1), t_sor_avg/(double)(step-1) );

	/* end of magic */
	return 0;
}
