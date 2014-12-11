/**
 * Working with Windows, so we need to disable the warnings in order to uses
 * sprintf.
 */
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "IO.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

IO::IO (const char *input, const char *output)
{
  //Read the file with the simulations parameters
  settings = input;
  this->output = (char*)output;
}

IO::IO(int argc, char** argv)
{
	parseArguments(argc, argv);
	std::cout << "out dir: " << output << " settings file: " << settings << std::endl;
}

IO::~IO ()
{

}

void IO::dieSynopsis() {
	std::cerr << "Synopsis: IO "
		<< "--out \"output directory\" --set \"path to the settings file\" ";
	exit(-1);
}

void IO::parseArguments(int argc, char** argv)
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
		else if (argv[0] == (std::string)"--d")
		{
			m_debug = true;
			--argc; ++argv;
		}
		else if (argv[0] == (std::string)"--li")
		{
			m_log_info = true;
			--argc; ++argv;
		}
		else if (argv[0] == (std::string)"--le")
		{
			m_log_err = true;
			--argc; ++argv;
		}
		else if (argv[0] == (std::string)"--lw")
		{
			m_log_warn = true;
			--argc; ++argv;
		}
		else
		{
			--argc; ++argv;
		}
	}
}

SimParams IO::readInputfile()
{
  //Store the input parameters.
	std::ifstream file (settings, std::ios::in);// | std::ios::binary);

	if(! file.is_open() ) 
		std::cerr << "ERROR: could not open file \"" << settings << "\"" << "\"" << std::endl;

	std::stringstream buffer;
	buffer << file.rdbuf();
	std::size_t foundat;


	SimParams simparam;

	/* pretty much look for every option we support, one by one.
	 * this way we also ignore substring that have nothing to do with 
	 * configuration parameters at all, e.g. comments */

	foundat = buffer.str().find("xLength=");
	if(foundat != std::string::npos)
		simparam.xLength = ::strtod(&buffer.str()[foundat+8], 0);
	else
		std::cerr << "ERROR: could not find parameter 'xLength' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("yLength=");
	if(foundat != std::string::npos)
		simparam.yLength = ::strtod(&buffer.str()[foundat+8], 0);
	else
		std::cerr << "ERROR: could not find parameter 'yLength' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("iMax=");
	if(foundat != std::string::npos)
		simparam.iMax = ::strtol(&buffer.str()[foundat+5], 0, 10);
	else
		std::cerr << "ERROR: could not find parameter 'iMax' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("jMax=");
	if(foundat != std::string::npos)
		simparam.jMax = ::strtol(&buffer.str()[foundat+5], 0, 10);
	else
		std::cerr << "ERROR: could not find parameter 'jMax' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("tEnd=");
	if(foundat != std::string::npos)
		simparam.tEnd = ::strtod(&buffer.str()[foundat+5], 0);
	else
		std::cerr << "ERROR: could not find parameter 'tEnd' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("deltaT=");
	if(foundat != std::string::npos)
		simparam.deltaT = ::strtod(&buffer.str()[foundat+7], 0);
	else
		std::cerr << "ERROR: could not find parameter 'deltaT' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("tau=");
	if(foundat != std::string::npos)
		simparam.tau = ::strtod(&buffer.str()[foundat+4], 0);
	else
		std::cerr << "ERROR: could not find parameter 'tau' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("deltaVec=");
	if(foundat != std::string::npos)
		simparam.deltaVec = ::strtod(&buffer.str()[foundat+9], 0);
	else
		std::cerr << "ERROR: could not find parameter 'deltaVec' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("iterMax=");
	if(foundat != std::string::npos)
		simparam.iterMax = ::strtol(&buffer.str()[foundat+8], 0, 10);
	else
		std::cerr << "ERROR: could not find parameter 'iterMax' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("eps=");
	if(foundat != std::string::npos)
		simparam.eps = ::strtod(&buffer.str()[foundat+4], 0);
	else
		std::cerr << "ERROR: could not find parameter 'tau' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("omg=");
	if(foundat != std::string::npos)
		simparam.omg = ::strtod(&buffer.str()[foundat+4], 0);
	else
		std::cerr << "ERROR: could not find parameter 'omg' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("alpha=");
	if(foundat != std::string::npos)
		simparam.alpha = ::strtod(&buffer.str()[foundat+6], 0);
	else
		std::cerr << "ERROR: could not find parameter 'alpha' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("re=");
	if(foundat != std::string::npos)
		simparam.re = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 're' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("gx=");
	if(foundat != std::string::npos)
		simparam.gx = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 'gx' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("gy=");
	if(foundat != std::string::npos)
		simparam.gy = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 'gy' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("ui=");
	if(foundat != std::string::npos)
		simparam.ui = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 'ui' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("vi=");
	if(foundat != std::string::npos)
		simparam.vi = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 'vi' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("pi=");
	if(foundat != std::string::npos)
		simparam.pi = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 'pi' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("xProcs=");
	if (foundat != std::string::npos)
		simparam.xProcs = ::strtol(&buffer.str()[foundat + 7], 0, 10);
	else
		std::cerr << "ERROR: could not find parameter 'xProcs' in file \"" << settings << "\"" << std::endl;

	foundat = buffer.str().find("yProcs=");
	if (foundat != std::string::npos)
		simparam.yProcs = ::strtol(&buffer.str()[foundat + 7], 0, 10);
	else
		std::cerr << "ERROR: could not find parameter 'yProcs' in file \"" << settings << "\"" << std::endl;
	
	return simparam;
}

#if defined(__linux)
	#include "sys/stat.h"
#endif
void IO::checkOutputDir()
{
	std::string filename;
	filename.append("./");
	filename.append(output);
	filename.append("/");

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

void IO::writeVTKMasterFile(Communication& comm, int step)
{
	std::string filename;
	filename.append("./");
	filename.append(output);
	filename.append("/");
	filename.append("field_");
	filename.append(std::to_string(step));
	filename.append(".pvtr");

	std::filebuf fb;
	fb.open(const_cast < char *>(filename.c_str()), std::ios::out);
	std::ostream os(&fb);

	/* all these values are inclusive. like in 'for(int i=min; i <= max; i++)' */
	/* also these values should start and end at the INNER of the domain.
	 * when/if needed, the offset for boundaries will be added as -1/+1 */
	int xGlobMin = 1; 
	int xGlobMax = comm.getGlobalDimensions().i; 
	int yGlobMin = 1; 
	int yGlobMax = comm.getGlobalDimensions().j; 

	os << "<?xml version=\"1.0\"?>" << std::endl;
	os << "<VTKFile type=\"PRectilinearGrid\">" << std::endl;
	/* extends */
	os << "<PRectilinearGrid WholeExtent=\""
		<< xGlobMin <<" "<< xGlobMax << " " 
		<< yGlobMin <<" "<< yGlobMax 
		<< " 0 0\" GhostLevel=\"" << 0 << "\">"<< std::endl;
	os << "<PCoordinates>" << std::endl;
	os << "<PDataArray type=\"Float64\"/>" << std::endl;
	os << "<PDataArray type=\"Float64\"/>" << std::endl;
	os << "<PDataArray type=\"Float64\"/>" << std::endl;
	os << "</PCoordinates>" << std::endl;

	/* iterate the local domains and write out their extents */
	/* the following Locl variables are analogous to the above "Glob" ones, 
	 * but for the local domain of the process "localProcRank" */
	int xLoclbMin, xLoclbMax, yLoclbMin, yLoclbMax; 
	int localProcRank;

	/* loop through all processes */
	for (int x = 0; x < comm.getProcsGridDim().i; x++)
	{
		for (int y = 0; y < comm.getProcsGridDim().j; y++)
		{
			/* compute 'local' process rank and its dimensions from position in 
			 * grid of processes */
			localProcRank = x * comm.getProcsGridDim().j + y; 

			Index expctdCellsPerDim(
				 (comm.getGlobalDimensions().i / comm.getProcsGridDim().i) ,
				 (comm.getGlobalDimensions().j / comm.getProcsGridDim().j) );

			Dimension localOffsetToGlobalDomain(
				 expctdCellsPerDim.i * x ,
				 expctdCellsPerDim.j * y );

			Dimension localDomainDim(
			/* last process in x/y dimensions gets the rest of cells in that direction,
			 * that is why the following computation looks ugly */
				 expctdCellsPerDim.i 
				 	+ ( (x < comm.getProcsGridDim().i-1)
					?(0):(comm.getGlobalDimensions().i % comm.getProcsGridDim().i)) ,
				 expctdCellsPerDim.j 
				 	+ ( (y < comm.getProcsGridDim().j-1)
					?(0):(comm.getGlobalDimensions().j % comm.getProcsGridDim().j)) );

			/* here, write extends of sub-domain of every process Omega_{i,j} */
			/* for a consistent visualization ParaView wants 
			 * the subdomains to overlapp, this is why we extend the inner by 
			 * one in all directions. so we later need to write out the boundary too.
			 * because we will later write everything with respect to the pressure
			 * (which sits in the center of each cell), this will work out. 
			 * (we will need to interpolate velocities accordingly) */
			xLoclbMin = localOffsetToGlobalDomain.i + xGlobMin;
			xLoclbMax = xLoclbMin + localDomainDim.i -1; 
			/* 'localDomainDim.i-1'
			 * because we want to go from first cell of domain to last, but not too far */
			yLoclbMin = localOffsetToGlobalDomain.j + yGlobMin;
			yLoclbMax = yLoclbMin + localDomainDim.j -1;

			/* walk _on_ the boundary of domain of current process */
			xLoclbMin += -1;
			xLoclbMax += +1; 
			yLoclbMin += -1;
			yLoclbMax += +1;

			/* write extent of Domain of current process, 
			 * with respect to global Domain */
			os << "<Piece Extent=\"" 
				<< xLoclbMin << " " << xLoclbMax << " " 
				<< yLoclbMin << " " << yLoclbMax << " 0 0\" Source=\"field_" 
				<< std::to_string(step) << "_processor_"
				<< localProcRank << ".vtr\"/>" << std::endl;
		}
	}

	/* standard fields: */
	os << "<PPointData Vectors=\"field\" Scalars=\"p\">"<< std::endl; // , vorticity, stream
	os << "<PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"field\" format=\"ascii\"/>"<< std::endl;
	os << "<PDataArray type=\"Float64\" Name=\"p\" format=\"ascii\"/>"<< std::endl;

	/* new fields: vorticity and stream */
	//os << "<PDataArray type=\"Float64\" Name=\"vorticity\" format=\"ascii\"/>"<< std::endl;
	//os << "<PDataArray type=\"Float64\" Name=\"stream\" format=\"ascii\"/>"<< std::endl;

	os << "</PPointData>" << std::endl;
	os << "</PRectilinearGrid>" << std::endl;
	os << "</VTKFile>" << std::endl;
}

void IO::writeVTKSlaveFile(Domain& domain, Communication& comm, int step)
{

	Real cellDeltaX = domain.getDelta().x;
	Real cellDeltaY = domain.getDelta().y;
	int localProcRank = comm.getRank();

	std::string filename;
	filename.append("./");
	filename.append(output);
	filename.append("/");

	filename.append("field_");
	filename.append(std::to_string(step));
	filename.append("_processor_");
	filename.append(std::to_string(localProcRank));
	filename.append(".vtr");

	std::filebuf fb;
	fb.open(const_cast < char *>(filename.c_str()), std::ios::out);
	std::ostream os(&fb);

	int xGlobMin = 1; 
	//int xGlobMax = comm.getGlobalDimensions().i; 
	int yGlobMin = 1; 
	//int yGlobMax = comm.getGlobalDimensions().j; 

	int xLoclbMin, xLoclbMax, yLoclbMin, yLoclbMax; 
	/* set the four indices like in WriteVTKMasterFile */
	xLoclbMin = comm.getOwnOffsetToGlobalDomain().i + xGlobMin;
	xLoclbMax = xLoclbMin + domain.getDimension().i -1; 
	yLoclbMin = comm.getOwnOffsetToGlobalDomain().j + yGlobMin;
	yLoclbMax = yLoclbMin + domain.getDimension().j -1;
	/* walk _on_ the boundary of domain of current process */
	xLoclbMin += -1;
	xLoclbMax += +1; 
	yLoclbMin += -1;
	yLoclbMax += +1;

	os << "<?xml version=\"1.0\"?>" << std::endl;
	os << "<VTKFile type=\"RectilinearGrid\">" << std::endl;
	// whole extent = whole extent in master
	// os << "<RectilinearGrid WholeExtent=\"0 " << s.iMax - 1 << " 0 "
	// 		<< s.jMax - 1 << " 0 0\" GhostLevel=\"" << stencilwidth << "\">"
	// 		<< std::endl;

	// whole extent = piece extent in slave
	os << "<RectilinearGrid WholeExtent=\"" 
		<< xLoclbMin << " " << xLoclbMax << " " 
		<< yLoclbMin << " " << yLoclbMax
		<< " 0 0\" GhostLevel=\"" << 1 << "\">" << std::endl;

	os << "<Piece Extent=\"" 
		<< xLoclbMin << " " << xLoclbMax << " " 
		<< yLoclbMin << " " << yLoclbMax
		<< " 0 0\">" << std::endl;

	os << "<Coordinates>" << std::endl;
	os << "<DataArray type=\"Float64\" format=\"ascii\">" << std::endl;

	/* now we write the coordinates of the grid points in x and y direction */
	for (int i = xLoclbMin-1; i <= xLoclbMax+1; i++)
		os << std::scientific << i * cellDeltaX << " ";

	os << std::endl;
	os << "</DataArray>" << std::endl;

	os << "<DataArray type=\"Float64\" format=\"ascii\">" << std::endl;
	for (int j = yLoclbMin-1; j <= yLoclbMax+1; j++)
	{
		os << std::scientific << j * cellDeltaY << " ";
	}
	os << std::endl;
	os << "</DataArray>" << std::endl;

	os << "<DataArray type=\"Float64\" format=\"ascii\">0 0</DataArray>" << std::endl;
	os << "</Coordinates>" << std::endl;

	os << "<PointData Vectors=\"field\" Scalars=\"p\">" << std::endl; // , vorticity, stream
	os << "<DataArray Name=\"field\" NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">" << std::endl;

	/* write grid values with boundaries.
	 * interpolate u and v to be at p position of cell.
	 * we have to make up values on the boundary just to satisfy overlapping
	 * conditions of sub domains in the VTK file. */
	Real u_inter, v_inter;
	for (int j = 1 -1; j <= domain.getEndInnerDomainP().i +1; j++)
	{
		for (int i = 1 -1; i <= domain.getEndInnerDomainP().j +1; i++)
		{
			if(i == 0)
				u_inter = (domain.u()(i, j) + domain.u()(i , j)) / 2.0;
			else 
			if(i == domain.getEndInnerDomainP().i+1)
				u_inter = (domain.u()(i - 1, j) + domain.u()(i - 1, j)) / 2.0;
			else
				u_inter = (domain.u()(i, j) + domain.u()(i - 1, j)) / 2.0;


			if(j == 0 )
				v_inter = (domain.v()(i, j) + domain.v()(i, j)) / 2.0;
			else
			if(j == domain.getEndInnerDomainP().j +1)
				v_inter = (domain.v()(i, j - 1) + domain.v()(i, j - 1)) / 2.0;
			else
				v_inter = (domain.v()(i, j) + domain.v()(i, j - 1)) / 2.0;

			os << std::scientific <<
				u_inter << " " <<
				v_inter << " " <<
				0.0 << std::endl;
		}
		os << std::endl;
	}
	os << "</DataArray>" << std::endl;

	os << "<DataArray type=\"Float64\" Name=\"p\" format=\"ascii\">"
		<< std::endl;
	
	/* write out p as it is */
	for (int j = 1 -1; j <= domain.getEndInnerDomainP().i +1; j++)
	{
		for (int i = 1 -1; i <= domain.getEndInnerDomainP().j +1; i++)
		{
			os << std::scientific << domain.p()(i, j) << " ";
		}
		os << std::endl;
	}

	os << "</DataArray>" << std::endl;

//	os << "<DataArray type=\"Float64\" Name=\"vorticity\" format=\"ascii\">"
//		<< std::endl;
//
//	/* TODO TODO */
//	for (int i = 0; i < iMax + 2; i++)
//	for (int j = 1 -1; j <= domain.getEndInnerDomainP().i +1; j++)
//	{
//		os << std::scientific << 0.0 << " ";
//	}
//	os << std::endl;
//	for (int j = 1; j < jMax+1; j++)
//	{
//		os << std::scientific << 0.0 << " ";
//		for (int i = 1; i < iMax+1; i++)
//		{
//			Real zeta = (domain.u()(i, j + 1) - domain.u()(i, j)) / deltaY
//				- (domain.v()(i + 1, j) - domain.v()(i, j)) / deltaX;
//			os << std::scientific << zeta << " ";
//		}
//		os << std::scientific << 0.0 << std::endl;
//	}
//	// we need one additional line of vorticity...
//	for (int i = 0; i < iMax + 2; i++)
//	{
//		os << std::scientific << 0.0 << " ";
//	}
//	os << std::endl;
//
//	os << "</DataArray>" << std::endl;
//
//	os << "<DataArray type=\"Float64\" Name=\"stream\" format=\"ascii\">"
//		<< std::endl;
//
//	Dimension stream_dim(iMax + 2, jMax + 2);
//	GridFunction stream(stream_dim);
//
//	// we need one initial line of stream...
//	for (int i = 0; i < jMax + 2; i++)
//	{
//		stream(i, 0) = 0.0;
//	}
//	for (int j = 1; j < jMax + 2; j++)
//	{
//		for (int i = 0; i < iMax + 2; i++)
//		{
//			stream(i, j) = stream(i, j - 1)
//				+ domain.u()(i, j) * deltaY;
//		}
//	}
//	os << std::endl;
//	for (int j = 0; j < jMax + 2; j++)
//	{
//		//os << std::scientific << 0.0 << " ";
//		for (int i = 0; i < iMax + 2; i++)
//		{
//			os << std::scientific << stream(i, j) << " ";
//		}
//		os << std::endl;
//	}
//	
//	os << "</DataArray>" << std::endl;


	os << "</PointData>" << std::endl;
	os << "</Piece>" << std::endl;

	os << "</RectilinearGrid>" << std::endl;
	os << "</VTKFile>" << std::endl;
}

void
IO::writeVTKFile (Domain& domain, int step)
{
  Real deltaX = domain.getDelta()[0];
  Real deltaY = domain.getDelta()[1];

  int iMax = domain.getDimension()[0];// w.r.t. inner of P
  int jMax = domain.getDimension()[1];// w.r.t. inner of P

  std::string filename;
  filename.append("./");
  filename.append (output);
  filename.append("/");

  // Test if the directory exits, if not create it.
  std::filebuf test_dir;
  test_dir.open(const_cast < char *>(filename.c_str()), std::ios::out);

  if (!test_dir.is_open())
  {
	  // Directory doesn't exist.
      #if defined(_WIN64)
			CreateDirectory(filename.c_str(),NULL);
      #elif defined(_WIN32)
			CreateDirectory(filename.c_str(),NULL);
      #elif defined(__linux)
			mkdir(filename.c_str(),0700);
	  #endif
  }
  
  filename.append ("field_");
  filename.append (std::to_string(step));
  filename.append (".vts");

  std::filebuf fb;
  fb.open (const_cast < char *>(filename.c_str ()), std::ios::out);
  std::ostream os (&fb);

  os << "<?xml version=\"1.0\"?>" << std::endl
    << "<VTKFile type=\"StructuredGrid\">" << std::endl
    << "<StructuredGrid WholeExtent=\""
    << "0" << " " << (iMax-1) << " "
    << "0" << " " << (jMax-1) << " "
    << "0" << " " << "0" << " "
    << "\" GhostLevel=\"" << "1" << "\">" << std::endl
    << "<Piece Extent=\""
    << "0" << " " << (iMax-1) << " "
    << "0" << " " << (jMax-1) << " "
    << "0" << " " << "0" << " "
    << "\">" << std::endl
    << "<Points>" << std::endl
    <<
    "<DataArray type=\"Float64\" format=\"ascii\" NumberOfComponents=\"3\"> "
    << std::endl;
  for (int i = 1; i <= iMax; ++i)
  {
	  for (int j = 1; j <= jMax; ++j)
	  {
		  os << std::scientific << i * deltaX << " " << j *
			  deltaY << " " << 0.0 << std::endl;
	  }
  }
  os << "</DataArray>" << std::endl
    << "</Points>" << std::endl
    << "<PointData Vectors=\"field\"  Scalars=\"P\">"
    << std::endl <<
    "<DataArray Name=\"field\" NumberOfComponents=\"3\" type=\"Float64\" >" <<
    std::endl;
  for (int i = 1; i <= iMax; ++i)
  {
	  for (int j = 1; j <= jMax; ++j)
	  {
		  os << std::scientific << 
			  (domain.u()(i, j) + domain.u()(i - 1, j)) / 2.0
			  /*interpolateVelocityU(x, y, u, delta)*/ << " " <<
			  (domain.v()(i, j) + domain.v()(i, j - 1)) / 2.0
			  /*interpolateVelocityV(x, y, v, delta)*/ << " " << 
			  0.0 << std::endl;
	  }
  }
  os << "</DataArray>" << std::endl
    << "<DataArray type=\"Float64\" Name=\"P\" format=\"ascii\">" <<
    std::endl;
  for (int i = 1; i <= iMax; ++i)
  {
	  for (int j = 1; j <= jMax; ++j)
	  {
		  os << std::scientific << domain.p()(i, j) << " ";
	  }
	  os << std::endl;
  }

  os << "</DataArray>" << std::endl
    << "</PointData>" << std::endl
    << "</Piece>" << std::endl
    << "</StructuredGrid>" << std::endl << "</VTKFile>" << std::endl;
  fb.close ();
}
