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

IO::IO (char *input, char *output)
{
  //Read the file with the simulations parameters
  readInputfile(input);
  this->output = output;
}

IO::~IO ()
{

}

void
IO::readInputfile (char *filename)
{
  //Store the input parameters.
	std::ifstream file (filename, std::ios::in);// | std::ios::binary);

	if(! file.is_open() ) 
		std::cerr << "ERROR: could not open file \"" << filename << "\"" << "\"" << std::endl;

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
		std::cerr << "ERROR: could not find parameter 'xLength' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("yLength=");
	if(foundat != std::string::npos)
		simparam.yLength = ::strtod(&buffer.str()[foundat+8], 0);
	else
		std::cerr << "ERROR: could not find parameter 'yLength' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("iMax=");
	if(foundat != std::string::npos)
		simparam.iMax = ::strtol(&buffer.str()[foundat+5], 0, 10);
	else
		std::cerr << "ERROR: could not find parameter 'iMax' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("jMax=");
	if(foundat != std::string::npos)
		simparam.jMax = ::strtol(&buffer.str()[foundat+5], 0, 10);
	else
		std::cerr << "ERROR: could not find parameter 'jMax' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("tEnd=");
	if(foundat != std::string::npos)
		simparam.tEnd = ::strtod(&buffer.str()[foundat+5], 0);
	else
		std::cerr << "ERROR: could not find parameter 'tEnd' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("deltaT=");
	if(foundat != std::string::npos)
		simparam.deltaT = ::strtod(&buffer.str()[foundat+7], 0);
	else
		std::cerr << "ERROR: could not find parameter 'deltaT' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("tau=");
	if(foundat != std::string::npos)
		simparam.tau = ::strtod(&buffer.str()[foundat+4], 0);
	else
		std::cerr << "ERROR: could not find parameter 'tau' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("deltaVec=");
	if(foundat != std::string::npos)
		simparam.deltaVec = ::strtod(&buffer.str()[foundat+9], 0);
	else
		std::cerr << "ERROR: could not find parameter 'deltaVec' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("iterMax=");
	if(foundat != std::string::npos)
		simparam.iterMax = ::strtol(&buffer.str()[foundat+8], 0, 10);
	else
		std::cerr << "ERROR: could not find parameter 'iterMax' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("eps=");
	if(foundat != std::string::npos)
		simparam.eps = ::strtod(&buffer.str()[foundat+4], 0);
	else
		std::cerr << "ERROR: could not find parameter 'tau' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("omg=");
	if(foundat != std::string::npos)
		simparam.omg = ::strtod(&buffer.str()[foundat+4], 0);
	else
		std::cerr << "ERROR: could not find parameter 'omg' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("alpha=");
	if(foundat != std::string::npos)
		simparam.alpha = ::strtod(&buffer.str()[foundat+6], 0);
	else
		std::cerr << "ERROR: could not find parameter 'alpha' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("re=");
	if(foundat != std::string::npos)
		simparam.re = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 're' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("gx=");
	if(foundat != std::string::npos)
		simparam.gx = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 'gx' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("gy=");
	if(foundat != std::string::npos)
		simparam.gy = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 'gy' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("ui=");
	if(foundat != std::string::npos)
		simparam.ui = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 'ui' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("vi=");
	if(foundat != std::string::npos)
		simparam.vi = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 'vi' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("pi=");
	if(foundat != std::string::npos)
		simparam.pi = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 'pi' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("xProcs=");
	if (foundat != std::string::npos)
		simparam.xProcs = ::strtol(&buffer.str()[foundat + 7], 0, 10);
	else
		std::cerr << "ERROR: could not find parameter 'xProcs' in file \"" << filename << "\"" << std::endl;

	foundat = buffer.str().find("yProcs=");
	if (foundat != std::string::npos)
		simparam.yProcs = ::strtol(&buffer.str()[foundat + 7], 0, 10);
	else
		std::cerr << "ERROR: could not find parameter 'yProcs' in file \"" << filename << "\"" << std::endl;
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
		simparam.pi << std::endl <<
		"xProcs=" <<
		simparam.xProcs << std::endl <<
		"yProcs=" <<
		simparam.yProcs << std::endl;
}


#define Element(field,ic) ((field)[(ic)[0]][(ic)[1]])

RealType
  IO::interpolateVelocityU (RealType x, RealType y, GridFunctionType & u,
			    const PointType & delta)
{

  RealType deltaX = delta[0];
  RealType deltaY = delta[1];

  MultiIndexType index;

  // Computation of u(x,y)
  index[0] = ((int) (x / deltaX)) + 1;
  index[1] = ((int) ((y + (deltaY / 2)) / deltaY)) + 1;

  // The coordinates of the cell corners

  RealType x1 = (index[0] - 1) * deltaX;
  RealType x2 = index[0] * deltaX;
  RealType y1 = ((index[1] - 1) - 0.5) * deltaY;
  RealType y2 = (index[1] - 0.5) * deltaY;

  MultiIndexType offset;

  offset[0] = index[0] - 1;
  offset[1] = index[1] - 1;

  RealType u1 = Element (u, offset);	// datafields->u->getField ()[i - 1][j - 1];

  offset[0] = index[0];
  offset[1] = index[1] - 1;

  RealType u2 = Element (u, offset);	//datafields->u->getField ()[i][j - 1];

  offset[0] = index[0] - 1;
  offset[1] = index[1];

  RealType u3 = Element (u, offset);	//datafields->u->getField ()[i - 1][j];
  RealType u4 = Element (u, index);

  RealType
    uInterploated =
    (1.0 / (deltaX * deltaY)) * ((x2 - x) * (y2 - y) *
				 u1 + (x - x1) * (y2 -
						  y) *
				 u2 + (x2 - x) * (y -
						  y1) *
				 u3 + (x - x1) * (y - y1) * u4);

  return uInterploated;
}


RealType
  IO::interpolateVelocityV (RealType x, RealType y, GridFunctionType & v,
			    const PointType & delta)
{
  RealType deltaX = delta[0];
  RealType deltaY = delta[1];

  // Computation of v(x,y)
  MultiIndexType index;
  index[0] = ((int) ((x + (deltaX / 2)) / deltaX)) + 1;
  index[1] = ((int) (y / deltaY)) + 1;

  // The coordinates of the cell corners

  RealType x1 = ((index[0] - 1) - 0.5) * deltaX;
  RealType x2 = (index[0] - 0.5) * deltaX;
  RealType y1 = (index[1] - 1) * deltaY;
  RealType y2 = index[1] * deltaY;

  MultiIndexType offset;

  offset[0] = index[0] - 1;
  offset[1] = index[1] - 1;

  RealType v1 = Element (v, offset);	//datafields->v->getField ()[i - 1][j - 1];

  offset[0] = index[0];
  offset[1] = index[1] - 1;

  RealType v2 = Element (v, offset);	//datafields->v->getField ()[i][j - 1];

  offset[0] = index[0] - 1;
  offset[1] = index[1];

  RealType v3 = Element (v, offset);	//datafields->v->getField ()[i - 1][j];


  RealType v4 = Element (v, index);	//datafields->v->getField ()[i][j];

  RealType
    vInterpolated =
    (1.0 / (deltaX * deltaY)) * ((x2 - x) * (y2 - y) *
				 v1 + (x - x1) * (y2 -
						  y) *
				 v2 + (x2 - x) * (y -
						  y1) *
				 v3 + (x - x1) * (y - y1) * v4);
  return vInterpolated;
}

void
IO::writeVTKFile (const MultiIndexType & griddimension, GridFunctionType & u,
		  GridFunctionType & v, GridFunctionType & p,
		  const PointType & delta, int step)
{
  RealType deltaX = delta[0];
  RealType deltaY = delta[1];

  IndexType iMax = griddimension[0] - 1;
  IndexType jMax = griddimension[1] - 1;

  char numstr[21];
  sprintf (numstr, "%d", step);
  std::string filename;
  filename.append ("./");
  filename.append ("field_");
  filename.append (numstr);
  filename.append (".vts");

  std::filebuf fb;
  fb.open (const_cast < char *>(filename.c_str ()), std::ios::out);
  std::ostream os (&fb);

  os << "<?xml version=\"1.0\"?>" << std::endl
    << "<VTKFile type=\"StructuredGrid\">" << std::endl
    << "<StructuredGrid WholeExtent=\""
    << "0" << " " << (iMax - 1) << " "
    << "0" << " " << (jMax - 1) << " "
    << "0" << " " << "0" << " "
    << "\" GhostLevel=\"" << "1" << "\">" << std::endl
    << "<Piece Extent=\""
    << "0" << " " << (iMax - 1) << " "
    << "0" << " " << (jMax - 1) << " "
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
  for (int i = 0; i < iMax; ++i)
    {
      RealType x = i * deltaX;

      for (int j = 0; j < jMax; ++j)
	{
	  RealType y = j * deltaY;

	  os << std::scientific << interpolateVelocityU (x, y, u,
							 delta) << " " <<
	    interpolateVelocityV (x, y, v, delta) << " " << 0. << std::endl;
	}

    }
  os << "</DataArray>" << std::endl
    << "<DataArray type=\"Float64\" Name=\"P\" format=\"ascii\">" <<
    std::endl;
  for (int i = 0; i <= iMax; ++i)
    {
      for (int j = 0; j <= jMax; ++j)
	{
	  os << std::scientific << p[i][j] << " ";

	}
      os << std::endl;

    }

  os << "</DataArray>" << std::endl
    << "</PointData>" << std::endl
    << "</Piece>" << std::endl
    << "</StructuredGrid>" << std::endl << "</VTKFile>" << std::endl;
  fb.close ();
}
