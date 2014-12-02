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


#define Element(field,ic) ((field)((ic)[0],(ic)[1]))

Real
  IO::interpolateVelocityU (Real x, Real y, GridFunction & u,
			    const Point & delta)
{

  Real deltaX = delta[0];
  Real deltaY = delta[1];

  Index index;

  // Computation of u(x,y)
  index[0] = ((int) (x / deltaX)) + 1;
  index[1] = ((int) ((y + (deltaY / 2)) / deltaY)) + 1;

  // The coordinates of the cell corners

  Real x1 = (index[0] - 1) * deltaX;
  Real x2 = index[0] * deltaX;
  Real y1 = ((index[1] - 1) - 0.5) * deltaY;
  Real y2 = (index[1] - 0.5) * deltaY;

  Index offset;

  offset[0] = index[0] - 1;
  offset[1] = index[1] - 1;

  Real u1 = Element(u, offset);	// datafields->u->getField ()[i - 1][j - 1];

  offset[0] = index[0];
  offset[1] = index[1] - 1;

  Real u2 = Element(u, offset);	//datafields->u->getField ()[i][j - 1];

  offset[0] = index[0] - 1;
  offset[1] = index[1];

  Real u3 = Element(u, offset);	//datafields->u->getField ()[i - 1][j];
  Real u4 = Element(u, index);

  Real
    uInterploated =
    (1.0 / (deltaX * deltaY)) * ((x2 - x) * (y2 - y) *
				 u1 + (x - x1) * (y2 -
						  y) *
				 u2 + (x2 - x) * (y -
						  y1) *
				 u3 + (x - x1) * (y - y1) * u4);

  return uInterploated;
  //return 0.0;
}


Real
  IO::interpolateVelocityV (Real x, Real y, GridFunction & v,
			    const Point & delta)
{
  Real deltaX = delta[0];
  Real deltaY = delta[1];

  // Computation of v(x,y)
  Index index;
  index[0] = ((int) ((x + (deltaX / 2)) / deltaX)) + 1;
  index[1] = ((int) (y / deltaY)) + 1;

  // The coordinates of the cell corners

  Real x1 = ((index[0] - 1) - 0.5) * deltaX;
  Real x2 = (index[0] - 0.5) * deltaX;
  Real y1 = (index[1] - 1) * deltaY;
  Real y2 = index[1] * deltaY;

  Index offset;

  offset[0] = index[0] - 1;
  offset[1] = index[1] - 1;

  Real v1 = Element (v, offset);	//datafields->v->getField ()[i - 1][j - 1];

  offset[0] = index[0];
  offset[1] = index[1] - 1;

  Real v2 = Element (v, offset);	//datafields->v->getField ()[i][j - 1];

  offset[0] = index[0] - 1;
  offset[1] = index[1];

  Real v3 = Element (v, offset);	//datafields->v->getField ()[i - 1][j];


  Real v4 = Element (v, index);	//datafields->v->getField ()[i][j];

  Real
    vInterpolated =
    (1.0 / (deltaX * deltaY)) * ((x2 - x) * (y2 - y) *
				 v1 + (x - x1) * (y2 -
						  y) *
				 v2 + (x2 - x) * (y -
						  y1) *
				 v3 + (x - x1) * (y - y1) * v4);
  return vInterpolated;
  //return 0.0;
}

#if defined(__linux)
	#include "sys/stat.h"
#endif

void IO::writeVTKMasterFile(const Index & griddimension, 
	int step, Communication &comm)
{
	int stencilwidth = 0;

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

	filename.append("field_");
	filename.append(std::to_string(step));
	filename.append(".pvtr");

	std::filebuf fb;
	fb.open(const_cast < char *>(filename.c_str()), std::ios::out);
	std::ostream os(&fb);

	// iLocalMax & jLocalMax are the size of the area..
	int iMax = griddimension.i; 
	int jMax = griddimension.j; 
	os << "<?xml version=\"1.0\"?>" << std::endl;
	os << "<VTKFile type=\"PRectilinearGrid\">" << std::endl;
	os << "<PRectilinearGrid WholeExtent=\"0 " << (iMax - 1) << " 0 "
		<< (jMax - 1) << " 0 0\" GhostLevel=\"" << stencilwidth << "\">"
		<< std::endl;
	os << "<PCoordinates>" << std::endl;
	os << "<PDataArray type=\"Float64\"/>" << std::endl;
	os << "<PDataArray type=\"Float64\"/>" << std::endl;
	os << "<PDataArray type=\"Float64\"/>" << std::endl;
	os << "</PCoordinates>" << std::endl;

	// iterate the dimensions, calculating the rank of the specific coords and the respective variables...

	iMax = comm.getLocalDimensions()[0]-1; // s.iLocalMax - 1; // w.r.t. inner of P
	jMax = comm.getLocalDimensions()[1]-1; // s.jLocalMax - 1; // w.r.t. inner of P
	for (int x = 0; x < griddimension.i; x++)
	{
		for (int y = 0; y < griddimension.j; y++)
		{
			int curRank;
			int coords[2] = { x, y };
			comm.getRankByCoords(coords, curRank);
			int x1 = x * iMax - stencilwidth;
			int x2 = (x + 1) * iMax + stencilwidth - 1;
			int x3 = y * jMax - stencilwidth;
			int x4 = (y + 1) * jMax + stencilwidth - 1;

			os << "<Piece Extent=\"" << x1 << " " << x2 << " " << x3 << " "
				<< x4 << " 0 0\" Source=\"field_" << std::to_string(step) << "_processor_"
				<< curRank << ".vtr\"/>" << std::endl;
		}
	}

	os << "<PPointData Vectors=\"field\" Scalars=\"p, vorticity, stream\">"
		<< std::endl;
	os
		<< "<PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"field\" format=\"ascii\"/>"
		<< std::endl;
	os << "<PDataArray type=\"Float64\" Name=\"p\" format=\"ascii\"/>"
		<< std::endl;
	os << "<PDataArray type=\"Float64\" Name=\"vorticity\" format=\"ascii\"/>"
		<< std::endl;
	os << "<PDataArray type=\"Float64\" Name=\"stream\" format=\"ascii\"/>"
		<< std::endl;
	os << "</PPointData>" << std::endl;

	os << "</PRectilinearGrid>" << std::endl;
	os << "</VTKFile>" << std::endl;
}

void IO::writeVTKSlaveFile(Domain& domain, int step,
	Communication &comm, SimParams sim_params)
{
	int stencilwidth = 0;
	int iMax = comm.getLocalDimensions()[0] - 1; // s.iLocalMax - 1;
	int jMax = comm.getLocalDimensions()[1] - 1; // s.jLocalMax - 1;

	Real deltaX = domain.getDelta()[0];
	Real deltaY = domain.getDelta()[1];

	int rank = comm.getRank();

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

	filename.append("field_");
	filename.append(std::to_string(step));
	filename.append("_processor_");
	filename.append(std::to_string(rank));
	filename.append(".vtr");

	std::filebuf fb;
	fb.open(const_cast < char *>(filename.c_str()), std::ios::out);
	std::ostream os(&fb);

	//int coords[2];
	//comm.getOwnCoords(coords); //TODO obtain this information via getProcsGridPosition and us an Index here
	Index coords = comm.getProcsGridPosition();
	int x = coords[0];
	int y = coords[1];

	int x1 = x * iMax - stencilwidth;
	int x2 = (x + 1) * iMax + stencilwidth - 1;
	int x3 = y * jMax - stencilwidth;
	int x4 = (y + 1) * jMax + stencilwidth - 1;

	os << "<?xml version=\"1.0\"?>" << std::endl;
	os << "<VTKFile type=\"RectilinearGrid\">" << std::endl;
	// whole extent = whole extent in master
	// os << "<RectilinearGrid WholeExtent=\"0 " << s.iMax - 1 << " 0 "
	// 		<< s.jMax - 1 << " 0 0\" GhostLevel=\"" << stencilwidth << "\">"
	// 		<< std::endl;

	// whole extent = piece extent in slave
	os << "<RectilinearGrid WholeExtent=\"" << x1 << " " << x2 << " " << x3 << " "
		<< x4 << " 0 0\" GhostLevel=\"" << stencilwidth << "\">"
		<< std::endl;


	os << "<Piece Extent=\"" << x1 << " " << x2 << " " << x3 << " " << x4
		<< " 0 0\">" << std::endl;

	os << "<Coordinates>" << std::endl;

	os << "<DataArray type=\"Float64\" format=\"ascii\">" << std::endl;

	Real iOffset, jOffset;

	int iStart = x1 + stencilwidth;
	int jStart = x3 + stencilwidth;
	if (iStart == 0)
	{
		iOffset = 0.;
	}
	else
	{
		iOffset = (sim_params.xLength / sim_params.iMax) * iStart;
	}

	if (jStart == 0)
	{
		jOffset = 0.;
	}
	{
		jOffset = (sim_params.yLength / sim_params.jMax) * jStart;
	}

	for (int i = 0; i <= iMax + 1; i++)
	{
		os << std::scientific << i * deltaX + iOffset << " ";
	}
	os << std::endl;
	os << "</DataArray>" << std::endl;

	os << "<DataArray type=\"Float64\" format=\"ascii\">" << std::endl;
	for (int j = 0; j <= jMax + 1; j++)
	{
		os << std::scientific << j * deltaY + jOffset << " ";
	}
	os << std::endl;
	os << "</DataArray>" << std::endl;

	os << "<DataArray type=\"Float64\" format=\"ascii\">0 0</DataArray>"
		<< std::endl;

	os << "</Coordinates>" << std::endl;

	os << "<PointData Vectors=\"field\" Scalars=\"p, vorticity, stream\">"
		<< std::endl;
	os
		<< "<DataArray Name=\"field\" NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">"
		<< std::endl;
	for (int j = 0; j <= jMax + 1; j++)
	{
		for (int i = 0; i <= iMax + 1; ++i)
		{
			os << std::scientific << domain.u()(i,j) << " " << domain.v()(i,j) << " " << 0.
				<< " ";
		}
		os << std::endl;
	}
	os << "</DataArray>" << std::endl;

	os << "<DataArray type=\"Float64\" Name=\"p\" format=\"ascii\">"
		<< std::endl;
	for (int j = 1; j <= jMax + 1; j++)
	{
		for (int i = 1; i <= iMax + 2; i++)
		{
			os << std::scientific << (domain.p()(i,j) + domain.p()(i,j + 1)) / 2. << " ";
		}
		os << std::endl;
	}

	// print last line of pressure:
	for (int i = 1; i <= iMax + 2; i++)
	{
		os << std::scientific << domain.p()(i,jMax + 2);
	}
	os << "</DataArray>" << std::endl;

	os << "<DataArray type=\"Float64\" Name=\"vorticity\" format=\"ascii\">"
		<< std::endl;
	for (int i = 0; i <= iMax + 1; i++)
	{
		os << std::scientific << 0.0 << " ";
	}
	os << std::endl;

	for (int j = 1; j <= jMax; j++)
	{
		os << std::scientific << 0.0 << " ";
		for (int i = 1; i <= iMax; i++)
		{
			Real zeta = (domain.u()(i, j + 1) - domain.u()(i, j)) / deltaY
				- (domain.v()(i + 1, j) - domain.v()(i, j)) / deltaX;
			os << std::scientific << zeta << " ";
		}
		os << std::scientific << 0.0 << std::endl;
	}
	// we need one additional line of vorticity...
	for (int i = 0; i <= jMax + 1; i++)
	{
		os << std::scientific << 0.0 << " ";
	}
	os << std::endl;
	os << "</DataArray>" << std::endl;

	os << "<DataArray type=\"Float64\" Name=\"stream\" format=\"ascii\">"
		<< std::endl;
	Dimension stream_dim(iMax + 2, jMax + 2);
	GridFunction stream(stream_dim);

	for (int i = 1; i <= iMax + 1; i++)
	{
		for (int j = 1; j <= jMax + 1; j++)
		{
			stream(i,j) = stream(i,j - 1)
				+ domain.u()(i, j) * deltaY;
		}
	}

	// we need one initial line of stream...
	for (int i = 0; i <= jMax + 1; i++)
	{
		os << std::scientific << 0.0 << " ";
	}
	os << std::endl;
	for (int j = 1; j <= jMax + 1; j++)
	{
		os << std::scientific << 0.0 << " ";
		for (int i = 1; i <= iMax + 1; i++)
		{
			os << std::scientific << stream(i,j) << " ";
		}
		os << std::endl;
	}

	os << "</DataArray>" << std::endl;
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
