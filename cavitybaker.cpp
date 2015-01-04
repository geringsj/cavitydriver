#include "src/CavityRenderer.hpp"
#include "src/SimulationParameters.hpp"
#include "src/Debug.hpp"

#include "src/Bakery.hpp"

#include <string>

#if defined(__linux)
	#include "external/include/docopt.h"
#endif

void getCLIDouble(Real& real, std::map<std::string, docopt::value>& args, std::string string)
{
	real = (args[string].isString()) ? std::strtod(args[string].asString().c_str(),0) : (real);
}
void getCLIInt(int& _int, std::map<std::string, docopt::value> args, std::string string)
{
	_int = (args[string].isString()) ? std::stoi(args[string].asString().c_str(),0) : (_int);
}

void initCLI(int& setting, bool& gui, SimulationParameters& simparam, int argc, char** argv)
{
#if defined(__linux)
	static const char USAGE[] = R"(
Usage: cavitybaker SETTING [options]

SETTING:
 A name identifying the setting to generate.
 Can be one of the following (case sensitive):
   DrivenCavity (alternatives: DC , 0)
   ChannelFlow (CF , 1)
   ChannelFlowUpperHalf (CFUH , 2)
   StepFlow (SF , 3)
   ObstacleChannelFlow (OCF , 4)

Options:
  -h --help         Show this message.
  --gui             Opens the OpenGL user interface for interactive
                    parameters input. Only available when program was
                    compiled with OpenGL support.
  --name=STRING     Name of the file in which the configuration parameters
                    will be saved. If an empty string is given, the 
                    parameters file will be written to stdout. [default: ""]
  --xLength=FLOAT   Length of domain in x-direction. [default: 1.0]
  --yLength=FLOAT   Length of domain in y-direction. [default: 1.0]
  --iMax=INT        Number of cells the discrete grid will have 
                    in x-direction. [default: 64]
  --jMax=INT        Number of cells the discrete grid will have 
                    in y-direction. [default: 64]
  --tau=FLOAT       Scaling factor for timestep computations 
                    in (0,1]. [default: 0.5]
  --tEnd=FLOAT      Length of simulation in seconds. [default: 16.5]
  --deltaT=FLOAT    Timestep size. Will be ignored.[default: 0.1]
  --deltaVec=FLOAT  Time between two VTK Outputs of current 
                    state of simulation. [default: 0.2]
  --eps=FLOAT       Upper bound for residual for SOR solver. [default: 0.001]
  --omg=FLOAT       Omega factor for SOR calculation. [default: 1.7]
  --alpha=FLOAT     Mixing factor for Donor-Cell scheme 
                    for next velocities. [default: 0.9]
  --iterMax=INT     Maximum number of iterations for SOR solver. [default: 100]
  --re=INT          Reynolds number to be used. [default: 1000]
  --gx=FLOAT        External force in x direction. [default: 0.0]
  --gy=FLOAT        External force in y direction. [default: 0.0]
  --ui=FLOAT        Constant initial value for U-velocity field. [default: 0.0]
  --vi=FLOAT        Constant initial value for V-velocity field. [default: 0.0]
  --pi=FLOAT        Constant initial value for Pressure field. [default: 0.0]
  --ka=FLOAT        Angle of object for Karman Vortex Street SETTING. 
                    Will be ignored if other SETTING than 'ObstacleChannelFlow' 
                    is used. Specified in radiant. [default: 0.7854]
  --kow=FLOAT       Width of object for Karman Vortex Street SETTING. 
                    Will be ignored if other SETTING than 'ObstacleChannelFlow' 
                    is used. [default: 1.5*yLength/jMax]
)";

	std::map<std::string, docopt::value> args = 
		docopt::docopt(USAGE, { argv + 1, argv + argc });

	/* move cli domain params into a start simparams */
	//for(auto const& arg : args)
	//	std::cout << arg.first << " is " <<  arg.second << "- isString:"<<arg.second.isString()<<std::endl;

	if(
		args["SETTING"] == docopt::value(std::string("DrivenCavity")) ||
		args["SETTING"] == docopt::value(std::string("DC")) ||
		args["SETTING"] == docopt::value(std::string("0"))
		)
		setting = 0;
	else if(
		args["SETTING"] == docopt::value(std::string("ChannelFlow")) ||
		args["SETTING"] == docopt::value(std::string("CF")) ||
		args["SETTING"] == docopt::value(std::string("1"))
		)
		setting = 1;
	else if(
		args["SETTING"] == docopt::value(std::string("ChannelFlowUpperHalf")) ||
		args["SETTING"] == docopt::value(std::string("CFUH")) ||
		args["SETTING"] == docopt::value(std::string("2"))
		)
		setting = 2;
	else if(
		args["SETTING"] == docopt::value(std::string("StepFlow")) ||
		args["SETTING"] == docopt::value(std::string("SF")) ||
		args["SETTING"] == docopt::value(std::string("3"))
		)
		setting = 3;
	else if(
		args["SETTING"] == docopt::value(std::string("ObstacleChannelFlow")) ||
		args["SETTING"] == docopt::value(std::string("OCF")) ||
		args["SETTING"] == docopt::value(std::string("4"))
		)
		setting = 4;
	else setting = -1;

	gui = args["--gui"].asBool();

	if( args["--name"] != docopt::value(std::string("\"\"")) )
		simparam.name = args["--name"].asString();

	getCLIDouble(simparam.yLength, args, "--yLength");
	getCLIDouble(simparam.xLength, args, "--xLength");

	getCLIInt(simparam.iMax, args, "--iMax");
	getCLIInt(simparam.jMax, args, "--jMax");

	getCLIDouble(simparam.KarmanAngle, args, "--ka");
	getCLIDouble(simparam.KarmanObjectWidth, args, "--kow");

	/* TODO: rest of parameters */

#else
	/* on windows, default to GUI */
	gui = true;
#endif
}

void overwritePostBakeryParams(SimulationParameters& newsp, SimulationParameters& clisp)
{
	newsp.name = clisp.name;

	/* TODO: rest of parameters */
}

int main(int argc, char** argv)
{
	argc = argc*(1+0*(**argv)); /* just to get rid of some warnings */

	SimulationParameters simparam;
	bool gui = false;
	int setting = 0;

	/* parse cli values: use gui? inputvals-filename? what other values to set? */
	initCLI(setting, gui, simparam, argc, argv);

	/* let bakery create boundaries 'n stuff using start simparams */
	if(setting < 0 || setting > 4)
	{
		std::cout << "Invalid SETTING parameter.\n";
		exit(EXIT_FAILURE);
	}

	SimulationParameters newparam = Bakery::get(static_cast<Bakery::Setting>(setting), simparam);

	/* replace non-domain parameters in new simparams, to overwrite init-defaults */
	overwritePostBakeryParams(newparam, simparam);

	/* maybe give simparams to gui - gui feedback? */

	if(gui)
	{
		/**
		 * This is just a test but the range should
		 * match the inner range of p.
		 */
		Index begin = Index(1, 1);
		Index end = Index(newparam.iMax, newparam.jMax);
		Range range = Range(begin, end);

		CavityRenderer cavity_renderer;
		if (cavity_renderer.initBakeryVis(640, 480, newparam))
		{
			cavity_renderer.createGrid(range);
			cavity_renderer.paint();
		}
		/* TODO: get/verwrite newparams from gui */
	}

	if(newparam.name == "")
	{
		newparam.writeToSTDOUT();
	}
	else
	{
		newparam.writeSettingsFile(newparam.name);
		log_info("Config file written to \"%s\".", newparam.name.c_str());
	}

	return 0;
}
