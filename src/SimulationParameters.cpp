
#include "SimulationParameters.hpp"
#include "Debug.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

#define _USE_MATH_DEFINES
#include <math.h>

void SimulationParameters::init()
{
	/* just some default values */

	this->xLength = 1.0;
	this->yLength = 1.0;

	this->iMax = 64;
	this->jMax = 64;

	this->tau = 0.5;
	this->tEnd = 16.5;
	this->deltaT = 0.1;
	this->deltaVec = 0.2;

	this->eps = 0.001;
	this->omg = 1.7;;
	this->alpha = 0.9;
	this->iterMax = 100;

	this->re = 1000;
	this->gx = 0.0;
	this->gy = 0.0;

	this->ui = 0.0;
	this->vi = 0.0;
	this->pi = 0.0;

	/* extended simulation parameters */
	this->name = std::string("");
	this->KarmanAngle = 2.0 * M_PI / 8.0;
	this->KarmanObjectWidth = 1.5 * this->yLength / this->jMax;

	this->useComplexGeometry = 0;

	this->boundary_conditions.clear();
}

SimulationParameters::SimulationParameters(std::string settingsfile)
{
	this->init();
	this->readInputfile(settingsfile);
}
SimulationParameters::SimulationParameters() { this->init(); }

void SimulationParameters::writeToSTDOUT()
{
	std::cout << "SimulationParameters: " << std::endl <<
		"xLength=" <<
		xLength << std::endl <<
		"yLength=" <<
		yLength << std::endl <<
		"iMax=" <<
		iMax << std::endl <<
		"jMax=" <<
		jMax << std::endl <<
		"tEnd=" <<
		tEnd << std::endl <<
		"deltaT=" <<
		deltaT << std::endl <<
		"tau=" <<
		tau << std::endl <<
		"deltaVec=" <<
		deltaVec << std::endl <<
		"iterMax=" <<
		iterMax << std::endl <<
		"eps=" <<
		eps << std::endl <<
		"omg=" <<
		omg << std::endl <<
		"alpha=" <<
		alpha << std::endl <<
		"re=" <<
		re << std::endl <<
		"gx=" <<
		gx << std::endl <<
		"gy=" <<
		gy << std::endl <<
		"ui=" <<
		ui << std::endl <<
		"vi=" <<
		vi << std::endl <<
		"pi=" <<
		pi << std::endl <<
		"name=" << "\"" << name << "\"" << std::endl <<
		"KarmanAngle=" <<
		KarmanAngle << std::endl <<
		"KarmanObjectWidth=" <<
		KarmanObjectWidth << std::endl <<
		"useComplexGeometry=" <<
		useComplexGeometry << std::endl;
}

void SimulationParameters::writeSettingsFile(
		std::string settingsfile)
{
	std::ofstream fs(settingsfile); 
	
	if(!fs)
	{
		std::cerr << "ERROR: could not open file \""
			<< settingsfile << "\"" << std::endl;
	return;
	}

	fs << "-- Extended Inputvals for SimulationParameters.\n";
	fs << "-- When a parameter is given more than once,\n";
	fs << "-- the first occurrence of the parameter in the file will be used.\n\n";

	fs << "name=\"" << this->name << "\"\n\n";
	fs << "useComplexGeometry=" << std::scientific << this->useComplexGeometry << " <- set to 0 for driven cavity\n\n";

	fs << "-- Domain parameters\n";
	fs << "xLength=" << std::scientific << this->xLength << "\n";
	fs << "yLength=" << std::scientific << this->yLength << "\n";

	fs << "iMax=" << std::scientific << this->iMax << "\n";
	fs << "jMax=" << std::scientific << this->jMax << "\n\n";

	fs << "-- Time and Simulation parameters\n";
	fs << "tEnd=" << std::scientific << this->tEnd << "\n";
	fs << "tau=" << std::scientific << this->tau << "\n";
	fs << "deltaVec=" << std::scientific << this->deltaVec << "\n";
	fs << "deltaT=" << std::scientific << this->deltaT << "\n";
	fs << "iterMax=" << std::scientific << this->iterMax << "\n";
	fs << "eps=" << std::scientific << this->eps << "\n";
	fs << "omg=" << std::scientific << this->omg << "\n";
	fs << "alpha=" << std::scientific << this->alpha << "\n";
	fs << "re=" << std::scientific << this->re << "\n\n";


	fs << "-- Outside forces on the system\n";
	fs << "gx=" << std::scientific << this->gx << "\n";
	fs << "gy=" << std::scientific << this->gy << "\n\n";

	fs << "-- Initial field values\n";
	fs << "ui=" << std::scientific << this->ui << "\n";
	fs << "vi=" << std::scientific << this->vi << "\n";
	fs << "pi=" << std::scientific << this->pi << "\n\n";

	fs << "-- Karman Parameters\n";
	fs << "KarmanAngle=" << std::scientific << this->KarmanAngle << "\n";
	fs << "KarmanObjectWidth=" << std::scientific << this->KarmanObjectWidth << "\n\n";

	fs << "-- do NOT edit this by hand! \n-- this is going to be parsed as information for boundary conditions when handling complex geometry\n";
	fs << "BoundaryPieces=\"";
	std::size_t curr = 0;
	while(curr < this->boundary_conditions.size())
	{
		if(curr!=0)
			fs << "|";

		fs << std::scientific 
			<< static_cast<int>(this->boundary_conditions[curr].direction) << ",";
		fs << std::scientific 
			<< static_cast<int>(this->boundary_conditions[curr].condition) << ",";
		fs << std::scientific 
			<< static_cast<int>(this->boundary_conditions[curr].gridtype) << ",";
		fs << std::scientific << this->boundary_conditions[curr].condition_value << ",";
		fs << std::scientific << this->boundary_conditions[curr].range.begin.i << ",";
		fs << std::scientific << this->boundary_conditions[curr].range.begin.j << ",";
		fs << std::scientific << this->boundary_conditions[curr].range.end.i << ",";
		fs << std::scientific << this->boundary_conditions[curr].range.end.j;
		curr++;
	}
	fs << "\"\n";

	fs.close();
}

static bool findReal(Real& ret, const std::string& str, const std::string& buffer, const std::string& settingsfile)
{
	std::size_t foundat;
	foundat = buffer.find( str.c_str() );
	if(foundat != std::string::npos)
	{
		ret = (Real)strtod(&buffer[foundat+str.size()], 0);
		return true;
	}
	else
		std::cerr << "ERROR: could not find parameter \"" 
			<< str << "\" in file \"" << settingsfile << "\"" << std::endl;
	return false;
}

static Boundary::BoundaryPiece readSingleBP(std::string& bps)
{
	Real tmp[8];
	std::size_t tpos = 0;
	int tmppos = 0;
	std::string commadel(",");
	while ((tpos = bps.find(commadel)) != std::string::npos) 
	{
		auto nextnumber = bps.substr(0, tpos);
		tmp[tmppos++] = strtod( &nextnumber[0] , 0);
		bps.erase(0, tpos + commadel.length());
	}
	tmp[tmppos++] = strtod( &bps[0] , 0);
	Boundary::Direction dir(static_cast<Boundary::Direction>((int)tmp[0]));
	Boundary::Condition cond(static_cast<Boundary::Condition>((int)tmp[1]));
	Boundary::Grid grid(static_cast<Boundary::Grid>((int)tmp[2]));
	Real condv(tmp[3]);
	Range range(
		Index(static_cast<int>(tmp[4]), static_cast<int>(tmp[5])),
		Index(static_cast<int>(tmp[6]), static_cast<int>(tmp[7])) );

	return Boundary::BoundaryPiece(dir,cond,grid,condv,range);
}

static void readBoundaryPieces(
	std::vector<Boundary::BoundaryPiece>& bconds,
	std::string& boundps
	)
{ /* read in: (Direction=int, Condition=int, Grid=int, CondValue=Real, Range=4int) */
	// Boundary::BoundaryPiece bp;
	if(boundps.size() < 8)
		return;

	std::string delimiter = "|";

	/* taken from: 
	 * http://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c */
	std::size_t pos = 0;
	std::string token;
	while ((pos = boundps.find(delimiter)) != std::string::npos) {
		token = boundps.substr(0, pos);
		/* process token */
		//std::cout << token << std::endl;
		bconds.push_back( readSingleBP(token) );
		/* erase token from string */
		boundps.erase(0, pos + delimiter.length());
	}
	bconds.push_back( readSingleBP(boundps) );
	//std::cout << "[INFO] using given boundaries information in inputvals\n";
}

bool SimulationParameters::readInputfile(std::string settingsfile)
{
	std::ifstream file (settingsfile, std::ios::in);

	if(! file.is_open() ) 
	{
		std::cerr << "[ERROR] could not open inputvals file \"" << settingsfile << "\"" << "\"" << std::endl;
		std::exit(EXIT_FAILURE); // educate the user
	}
	//std::cout << "[INFO] loading inputvals from file \"" << settingsfile << "\"\n";

	std::stringstream ssbuffer;
	ssbuffer << file.rdbuf();
	std::string buffer = ssbuffer.str();
	Real tmp;

	/* pretty much look for every option we support, one by one.
	 * this way we also ignore substring that have nothing to do with 
	 * configuration parameters at all, e.g. comments */

	 findReal(this->xLength, std::string("xLength="), buffer, settingsfile);
	 findReal(this->yLength, std::string("yLength="), buffer, settingsfile);

	 if(!findReal(tmp, std::string("iMax="), buffer, settingsfile))
	 findReal(tmp, std::string("xCells="), buffer, settingsfile);
	 this->iMax = static_cast<int>(tmp);
	 if(!findReal(tmp, std::string("jMax="), buffer, settingsfile))
	 findReal(tmp, std::string("yCells="), buffer, settingsfile);
	 this->jMax = static_cast<int>(tmp);

	 findReal(this->tEnd, std::string("tEnd="), buffer, settingsfile);
	 findReal(this->deltaT, std::string("deltaT="), buffer, settingsfile);
	 findReal(this->tau, std::string("tau="), buffer, settingsfile);

	 if(!findReal(this->deltaVec, std::string("deltaVec="), buffer, settingsfile))
	 findReal(this->tDeltaWrite, std::string("tDeltaWrite="), buffer, settingsfile);

	 findReal(tmp, std::string("iterMax="), buffer, settingsfile);
	 this->iterMax = static_cast<int>(tmp);

	 findReal(this->eps, std::string("eps="), buffer, settingsfile);
	 findReal(this->omg, std::string("omg="), buffer, settingsfile);
	 findReal(this->alpha, std::string("alpha="), buffer, settingsfile);

	 findReal(this->re, std::string("re="), buffer, settingsfile);
	 findReal(this->gx, std::string("gx="), buffer, settingsfile);
	 findReal(this->gy, std::string("gy="), buffer, settingsfile);
	 findReal(this->ui, std::string("ui="), buffer, settingsfile);
	 findReal(this->vi, std::string("vi="), buffer, settingsfile);
	 findReal(this->pi, std::string("pi="), buffer, settingsfile);

	 findReal(this->KarmanAngle, std::string("KarmanAngle="), buffer, settingsfile);
	 findReal(this->KarmanObjectWidth, 
		std::string("KarmanObjectWidth="), buffer, settingsfile);

	/* get name */
	std::size_t foundat, foundat2;
	std::string namestr("name=\"");
	foundat = buffer.find( namestr );
	if(foundat != std::string::npos)
	{
		foundat = foundat+namestr.size();
		foundat2 = buffer.find("\"",foundat);
		if(foundat2 != std::string::npos)
			this->name = buffer.substr(foundat, foundat2-foundat);
	}
	else
		std::cerr << "ERROR: could not find parameter \"" 
			<< namestr << "\" in file \"" << settingsfile << "\"" << std::endl;

	/* read in and handle conditionals on 'useGeometryFile' */
	 if(findReal(tmp, std::string("useComplexGeometry="), buffer, settingsfile))
	 	this->useComplexGeometry= static_cast<int>(tmp);

	/* handle read-in of possible boundary values for more complex setting */
	 if(this->useComplexGeometry)
	 {
		std::string boundps("BoundaryPieces=\"");
		foundat = buffer.find( boundps );
		if(foundat != std::string::npos)
		{
			foundat = foundat+boundps.size();
			boundps = std::string("");
			foundat2 = buffer.find("\"",foundat);
			if(foundat2 != std::string::npos)
			{
				boundps = buffer.substr(foundat, foundat2-foundat);
			}
			else
			{
				std::cerr << "ERROR: could not read Boundary Pieces in file \"" 
					<< settingsfile << "\"" << std::endl;
			}
		}
		else
		{ /* error reading boundary pieces! */
			boundps = std::string("");
			std::cerr << "ERROR: could not read Boundary Pieces in file \"" 
				<< settingsfile << "\"" << std::endl;
		}
		readBoundaryPieces(this->boundary_conditions, boundps);
	 }
	
	return true;
}
