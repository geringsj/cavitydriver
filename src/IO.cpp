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
  readInputfile();
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
		else
		{
			--argc; ++argv;
		}
	}
}

Simparam IO::readInputfile()
{
  //Store the input parameters.
	std::ifstream file (settings, std::ios::in);// | std::ios::binary);

	if(! file.is_open() ) 
		std::cerr << "ERROR: could not open file \"" << settings << "\"" << "\"" << std::endl;

	std::stringstream buffer;
	buffer << file.rdbuf();
	std::size_t foundat;

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

void 
IO::writeSimParamToSTDOUT()
{
	std::cout << "SimParam: " << std::endl <<
		"xLength=" <<
		simparam.xLength << std::endl <<
		"yLength=" <<
		simparam.yLength << std::endl <<
		"iMax=" <<
		simparam.iMax << std::endl <<
		"jMax=" <<
		simparam.jMax << std::endl <<
		"tEnd=" <<
		simparam.tEnd << std::endl <<
		"deltaT=" <<
		simparam.deltaT << std::endl <<
		"tau=" <<
		simparam.tau << std::endl <<
		"deltaVec=" <<
		simparam.deltaVec << std::endl <<
		"iterMax=" <<
		simparam.iterMax << std::endl <<
		"eps=" <<
		simparam.eps << std::endl <<
		"omg=" <<
		simparam.omg << std::endl <<
		"alpha=" <<
		simparam.alpha << std::endl <<
		"re=" <<
		simparam.re << std::endl <<
		"gx=" <<
		simparam.gx << std::endl <<
		"gy=" <<
		simparam.gy << std::endl <<
		"ui=" <<
		simparam.ui << std::endl <<
		"vi=" <<
		simparam.vi << std::endl <<
		"pi=" <<
		simparam.pi << std::endl ;
		/* TODO remove comment when implementing MPI stuff */
		// <<
		// "xProcs=" <<
		// simparam.xProcs << std::endl <<
		// "yProcs=" <<
		// simparam.yProcs << std::endl;
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
void
IO::writeVTKFile (const Index & griddimension, GridFunction & u,
		  GridFunction & v, GridFunction & p,
		  const Point & delta, int step)
{
  Real deltaX = delta[0];
  Real deltaY = delta[1];

  int iMax = griddimension[0];// w.r.t. inner of P
  int jMax = griddimension[1];// w.r.t. inner of P

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
  for (int i = 0; i < iMax; ++i)
  {
	  for (int j = 0; j < jMax; ++j)
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
			  (u(i,j)+u(i-1,j)) / 2.0
			  /*interpolateVelocityU(x, y, u, delta)*/ << " " <<
			  (v(i,j)+v(i,j-1)) / 2.0
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
		  os << std::scientific << p(i,j) << " ";
	  }
	  os << std::endl;
  }

  os << "</DataArray>" << std::endl
    << "</PointData>" << std::endl
    << "</Piece>" << std::endl
    << "</StructuredGrid>" << std::endl << "</VTKFile>" << std::endl;
  fb.close ();
}
