
#include "SimulationParameters.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

SimulationParameters::SimulationParameters(std::string settingsfile)
{
	this->readInputfile(settingsfile);
}

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
		pi << std::endl;
}


bool SimulationParameters::readInputfile(std::string settingsfile)
{
	std::ifstream file (settingsfile, std::ios::in);

	if(! file.is_open() ) 
	{
		std::cerr << "ERROR: could not open file \"" << settingsfile << "\"" << "\"" << std::endl;
		return false;
	}

	std::stringstream buffer;
	buffer << file.rdbuf();
	std::size_t foundat;

	/* pretty much look for every option we support, one by one.
	 * this way we also ignore substring that have nothing to do with 
	 * configuration parameters at all, e.g. comments */

	foundat = buffer.str().find("xLength=");
	if(foundat != std::string::npos)
		this->xLength = ::strtod(&buffer.str()[foundat+8], 0);
	else
		std::cerr << "ERROR: could not find parameter 'xLength' in file \"" << settingsfile << "\"" << std::endl;

	foundat = buffer.str().find("yLength=");
	if(foundat != std::string::npos)
		this->yLength = ::strtod(&buffer.str()[foundat+8], 0);
	else
		std::cerr << "ERROR: could not find parameter 'yLength' in file \"" << settingsfile << "\"" << std::endl;

	foundat = buffer.str().find("iMax=");
	if(foundat != std::string::npos)
		this->iMax = ::strtol(&buffer.str()[foundat+5], 0, 10);
	else
		std::cerr << "ERROR: could not find parameter 'iMax' in file \"" << settingsfile << "\"" << std::endl;

	foundat = buffer.str().find("jMax=");
	if(foundat != std::string::npos)
		this->jMax = ::strtol(&buffer.str()[foundat+5], 0, 10);
	else
		std::cerr << "ERROR: could not find parameter 'jMax' in file \"" << settingsfile << "\"" << std::endl;

	foundat = buffer.str().find("tEnd=");
	if(foundat != std::string::npos)
		this->tEnd = ::strtod(&buffer.str()[foundat+5], 0);
	else
		std::cerr << "ERROR: could not find parameter 'tEnd' in file \"" << settingsfile << "\"" << std::endl;

	foundat = buffer.str().find("deltaT=");
	if(foundat != std::string::npos)
		this->deltaT = ::strtod(&buffer.str()[foundat+7], 0);
	else
		std::cerr << "ERROR: could not find parameter 'deltaT' in file \"" << settingsfile << "\"" << std::endl;

	foundat = buffer.str().find("tau=");
	if(foundat != std::string::npos)
		this->tau = ::strtod(&buffer.str()[foundat+4], 0);
	else
		std::cerr << "ERROR: could not find parameter 'tau' in file \"" << settingsfile << "\"" << std::endl;

	foundat = buffer.str().find("deltaVec=");
	if(foundat != std::string::npos)
		this->deltaVec = ::strtod(&buffer.str()[foundat+9], 0);
	else
		std::cerr << "ERROR: could not find parameter 'deltaVec' in file \"" << settingsfile << "\"" << std::endl;

	foundat = buffer.str().find("iterMax=");
	if(foundat != std::string::npos)
		this->iterMax = ::strtol(&buffer.str()[foundat+8], 0, 10);
	else
		std::cerr << "ERROR: could not find parameter 'iterMax' in file \"" << settingsfile << "\"" << std::endl;

	foundat = buffer.str().find("eps=");
	if(foundat != std::string::npos)
		this->eps = ::strtod(&buffer.str()[foundat+4], 0);
	else
		std::cerr << "ERROR: could not find parameter 'tau' in file \"" << settingsfile << "\"" << std::endl;

	foundat = buffer.str().find("omg=");
	if(foundat != std::string::npos)
		this->omg = ::strtod(&buffer.str()[foundat+4], 0);
	else
		std::cerr << "ERROR: could not find parameter 'omg' in file \"" << settingsfile << "\"" << std::endl;

	foundat = buffer.str().find("alpha=");
	if(foundat != std::string::npos)
		this->alpha = ::strtod(&buffer.str()[foundat+6], 0);
	else
		std::cerr << "ERROR: could not find parameter 'alpha' in file \"" << settingsfile << "\"" << std::endl;

	foundat = buffer.str().find("re=");
	if(foundat != std::string::npos)
		this->re = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 're' in file \"" << settingsfile << "\"" << std::endl;

	foundat = buffer.str().find("gx=");
	if(foundat != std::string::npos)
		this->gx = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 'gx' in file \"" << settingsfile << "\"" << std::endl;

	foundat = buffer.str().find("gy=");
	if(foundat != std::string::npos)
		this->gy = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 'gy' in file \"" << settingsfile << "\"" << std::endl;

	foundat = buffer.str().find("ui=");
	if(foundat != std::string::npos)
		this->ui = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 'ui' in file \"" << settingsfile << "\"" << std::endl;

	foundat = buffer.str().find("vi=");
	if(foundat != std::string::npos)
		this->vi = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 'vi' in file \"" << settingsfile << "\"" << std::endl;

	foundat = buffer.str().find("pi=");
	if(foundat != std::string::npos)
		this->pi = ::strtod(&buffer.str()[foundat+3], 0);
	else
		std::cerr << "ERROR: could not find parameter 'pi' in file \"" << settingsfile << "\"" << std::endl;
	
	return true;
}
