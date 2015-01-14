#include "src/SimulationParameters.hpp"
#include "src/Debug.hpp"
#include "src/Bakery.hpp"
#include "src/CavityRenderer.hpp"

#include <future>
#include <string>
#include <thread>
#include <iostream>

#include "optionparser.h"
static option::ArgStatus NonEmpty(const option::Option& option, bool msg)
{
  if (option.arg != 0 && option.arg[0] != 0)
    return option::ARG_OK;

  if (msg) log_err("Option \"%s\" requires a non-empty argument",option.name);
  return option::ARG_ILLEGAL;
}
enum optionIndex
{ 
	UNKNOWN,
	HELP,
	GUI,
	NAME,
	XLENGTH,
	YLENGTH,
	IMAX,
	JMAX,
	TAU,
	TEND,
	DELTAT,
	DELTAVEC,
	EPS,
	OMG,
	ALPHA,
	ITERMAX,
	RE,
	GX,
	GY,
	UI,
	VI,
	PI,
	KA,
	KOW,
	INFLOW
};
const option::Descriptor usage[] =
{
	{UNKNOWN, 0, "", "", option::Arg::None, 
R"(
Usage: cavitybaker SETTING [options]

Examples:
  Generate a driven cavity with RE=500 in a 3x1.5 box:
    cavitybaker DrivenCavity --re=500 -xLength==3 --yLength=1.5

  Generate a Karman street with object width beeing object 
  height(=yLength/2) with 1000 max iterations for SOR:
    cavitybaker 4 --yLength=1 --xLength=5 --kow=0.5 --iterMax=1000

SETTING:
 A name identifying the setting to generate.
 Can be one of the following (case sensitive):
   DrivenCavity (alternatives: DC , 0)
   ChannelFlow (CF , 1)
   ChannelFlowUpperHalf (CFUH , 2)
   StepFlow (SF , 3)
   ObstacleChannelFlow (OCF , KarmanStreet , KS , 4)

Options:)" },
	{HELP, 0, "h", "help", option::Arg::None,
"  -h --help  	Show this message." },
 {GUI, 0, "g", "gui", option::Arg::None,
"  -g --gui  \tOpens the OpenGL user interface for interactive parameters input. Only available when program was compiled with OpenGL support." },
 {NAME, 0, "", "name", NonEmpty, 
"  --name=STRING  \tName of the file in which the configuration parameters will be saved. If no name is given, the parameters file will be written to stdout. [default: no name]" },

 {XLENGTH, 0, "", "xLength", NonEmpty,
"  --xLength=FLOAT  \tLength of domain in x-direction. [default: 1.0]" },
 {YLENGTH, 0, "", "yLength", NonEmpty,
"  --yLength=FLOAT  \tLength of domain in y-direction. [default: 1.0 in DC / at least 5.0*xLength in flows]" },
 {IMAX, 0, "", "iMax", NonEmpty,
"  --iMax=INT  \tNumber of cells the discrete grid will have in x-direction. [default: 64]" },
 {JMAX, 0, "", "jMax", NonEmpty,
"  --jMax=INT  \tNumber of cells the discrete grid will have in y-direction. [default: 64]" },
 {TAU, 0, "", "tau", NonEmpty,
"  --tau=FLOAT  \tScaling factor for timestep computations in (0,1]. [default: 0,5]" },
 {TEND, 0, "", "tEnd", NonEmpty,
"  --tEnd=FLOAT  \tLength of simulation in seconds. [default: 16.5]" },
 {DELTAT, 0, "", "deltaT", NonEmpty,
"  --deltaT=FLOAT  \tTimestep size. Will be ignored." },
 {DELTAVEC, 0, "", "deltaVec", NonEmpty,
"  --deltaVec=FLOAT  \tTime in seconds between two VTK Outputs of current state of simulation. [default: 0.2]" },
 {EPS, 0, "", "", NonEmpty,
"  --eps=FLOAT  \tUpper bound for residual for SOR solver. [default: 0.001]" },
 {OMG, 0,  "", "omg", NonEmpty,
"  --omg=FLOAT  \tOmega factor for SOR calculation. [default: 1.7]" },
 {ALPHA, 0, "", "alpha", NonEmpty,
"  --alpha=FLOAT  \tMixing factor for Donor-Cell scheme for next velocities. [default: 0.9]" },
 {ITERMAX, 0, "", "iterMax", NonEmpty,
"  --iterMax=INT  \tMaximum number of iterations for SOR solver. [default: 100]" },
 {RE, 0, "", "re", NonEmpty,
"  --re=FLOAT\tReynolds number to be used. [default: 1000.0]" },
 {GX, 0, "", "gx", NonEmpty,
"  --gx=FLOAT  \tExternal force in x direction. [default: 0.0]" },
 {GY, 0, "", "gy", NonEmpty,
"  --gy=FLOAT  \tExternal force in y direction. [default: 0.0]" },
 {UI, 0, "", "ui", NonEmpty,
"  --ui=FLOAT  \tConstant initial value for U-velocity field. [default: 0.0]" },
 {VI, 0, "", "vi", NonEmpty,
"  --vi=FLOAT  \tConstant initial value for V-velocity field. [default: 0.0]" },
 {PI, 0, "", "pi", NonEmpty,
"  --pi=FLOAT  \tConstant initial value for Pressure field. [default: 0.0]" },
 {KA, 0, "", "ka", NonEmpty,
"  --ka=FLOAT  \tAngle of object for Karman Vortex Street SETTING specified in radiant. Will be ignored if other SETTING than 'ObstacleChannelFlow is used. [default: 45° in radiant]" },
 {KOW, 0, "", "kow", NonEmpty,
"  --kow=FLOAT  \tWidth of object for Karman Vortex Street SETTING. Will be ignored if other SETTING than 'ObstacleChannelFlow is used. [default: 2.5*5*deltaX]" },
 {INFLOW, 0, "", "inflow", NonEmpty,
"  --inflow=FLOAT  \tValue of INFLOW boundary condition in all cases. For DrivenCavity this is the U velocity at the upper boundary. For all other settings (the flows), this is the pressure INFLOW at the left boundary, resulting in beeing the pressure difference between left and right boundary. [default: 1.0 in DC / 0.5 in flows]" },
 {UNKNOWN, 0, 0, 0, 0, 0 }
};

void getMaybeCLIDouble(Real& real, option::Option* options, optionIndex ind)
{
	if(options[ind])
	{
		real = std::strtod(options[ind].arg,0);
		std::cout << "[INFO] found " << options[ind].name+2 << std::endl;
	}
}
void getMaybeCLIInt(int& _int, option::Option* options, optionIndex ind)
{
	if(options[ind])
	{
		_int = std::stoi(options[ind].arg,0);
		std::cout << "[INFO] found " << options[ind].name+2 << std::endl;
	}
}


int initCLI(int& setting, bool& gui, Real& inflowVal, SimulationParameters& simparam, int argc, char** argv)
{
	int doexit = 0;
	argc-=(argc>0); argv+=(argc>0); // skip program name argv if present
	option::Stats  stats(usage, argc, argv);
	//debug("options_max: %i, buffer_max: %i", stats.options_max, stats.buffer_max);
	//debug("argc: %i", argc);
	option::Option* options = new option::Option[stats.options_max];
	stats.buffer_max = argc;
	option::Option* buffer  = new option::Option[stats.buffer_max]; /* !! */
	option::Parser parse(true, usage, argc, argv, options, buffer);

	if(parse.error())
	{
		std::cout << "An options-parser error occurred.\n";
		doexit = EXIT_FAILURE+1;
	}
	if(options[HELP])
	{
		option::printUsage(std::cout, usage);
		doexit = EXIT_SUCCESS+1;
	}
	if(argc == 0 || !parse.nonOptionsCount())
	{
		std::cout << "Missing SETTING parameter. Call with '--help' for options.\n";
		doexit = EXIT_FAILURE+1;
	}
	if(parse.nonOptionsCount() > 1)
	{
		std::cout << "Too many non-option or unkown parameters, i am confused.\n";
		for(int i=0; i<parse.nonOptionsCount(); i++)
			std::cout << "  '" << parse.nonOption(i) << "'\n";
		doexit = EXIT_FAILURE+1;
	}
	if(doexit)
	{
		delete[] options; delete[] buffer;
		return doexit;
	}


	std::string settingstr = std::string(parse.nonOption(0));
	setting = -1;
	if(
			"0" == settingstr ||
			"DrivenCavity" == settingstr ||
			"DC" == settingstr
		)
 	setting = 0;
	if(
			"1" == settingstr ||
			"ChannelFlow" == settingstr ||
			"CF" == settingstr
		)
 	setting = 1;
	if(
			"2" == settingstr ||
			"ChannelFlowUpperHalf" == settingstr ||
			"CFUH" == settingstr
		)
 	setting = 2;
	if(
			"3" == settingstr ||
			"StepFlow" == settingstr ||
			"SF" == settingstr
		)
 	setting = 3;
	if(
			"4" == settingstr ||
			"ObstacleChannelFlow" == settingstr ||
			"OCF" == settingstr ||
			"KarmanStreet" == settingstr ||
			"KS" == settingstr
		)
 	setting = 4;

	gui = static_cast<bool>(options[GUI]);
	if(!gui && (setting < 0 || setting > 4))
	{
		std::cout << "Invalid SETTING parameter '"<< std::string(settingstr) 
			<<"'. Call with '--help' for options.\n";
		doexit = EXIT_FAILURE+1;
	}
	if(doexit)
	{
		delete[] options; delete[] buffer;
		return doexit;
	}

	if(setting >= 0)
		std::cout << "[INFO] found SETTING=" << setting << std::endl;
	if(gui)
		std::cout << "[INFO] found GUI switch\n";

	if(options[NAME])
	{
		simparam.name = std::string(options[NAME].arg);
		std::cout << "[INFO] found name=" << simparam.name << std::endl;
	}

	/* overwrite those on pre-bakery processing, because
	 * bakery uses them */
	getMaybeCLIDouble(simparam.yLength, options, YLENGTH);
	getMaybeCLIDouble(simparam.xLength, options, XLENGTH);
	getMaybeCLIInt(simparam.iMax, options, IMAX);
	getMaybeCLIInt(simparam.jMax, options, JMAX);
	getMaybeCLIDouble(simparam.KarmanAngle, options, KA);
	getMaybeCLIDouble(simparam.KarmanObjectWidth, options, KOW);

	/* overwrite those on post-bakery processing */
	getMaybeCLIDouble(simparam.tau, options, TAU);
	getMaybeCLIDouble(simparam.tEnd, options, TEND);
	getMaybeCLIDouble(simparam.deltaT, options, DELTAT);
	getMaybeCLIDouble(simparam.deltaVec, options, DELTAVEC);
	getMaybeCLIDouble(simparam.eps, options, EPS);
	getMaybeCLIDouble(simparam.omg, options, OMG);
	getMaybeCLIDouble(simparam.alpha, options, ALPHA);
	getMaybeCLIDouble(simparam.re, options, RE);
	getMaybeCLIDouble(simparam.gx, options, GX);
	getMaybeCLIDouble(simparam.gy, options, GY);
	getMaybeCLIDouble(simparam.ui, options, UI);
	getMaybeCLIDouble(simparam.vi, options, VI);
	getMaybeCLIDouble(simparam.pi, options, PI);
	getMaybeCLIInt(simparam.iterMax, options, ITERMAX);

	switch(setting)
	{
		case 0: /* Driven Cavity */
			inflowVal = 1.0;
		case 1: /* ChannelFlow */
		case 2: /* ChannelFlow Upper Half */
		case 3: /* Step Flow */
		case 4: /* Karman Object Street */
			inflowVal = 0.5;
			break;
		default: /* all other: fail */
			break;
	}
	getMaybeCLIDouble(inflowVal, options, INFLOW);

	delete[] options; delete[] buffer;
	return doexit;
}

void overwritePostBakeryParams(SimulationParameters& newsp, SimulationParameters& clisp)
{
	newsp.name = clisp.name;

	/* if user specifies another re than given in exercise sheet, use users re */
	newsp.re = clisp.re;
}

/**
 * Helper function that encloses complete OpenGL execution.
 * If done like this, we don't have to clean up OpenGL objects explicitly.
 */
void runVisualization(MTQueue<SimulationParameters>& inbox, MTQueue<SimulationParameters>& outbox,
	int window_width, int window_height, SimulationParameters& initial_params)
{
	CavityRenderer cavity_renderer(inbox,outbox);
	cavity_renderer.initBakeryVis(window_width,window_height,initial_params);
	cavity_renderer.paint();
}

int main(int argc, char** argv)
{
	argc = argc*(1+0*(**argv)); /* just to get rid of some warnings */

	SimulationParameters simparam;
	bool gui = true;
	int setting = -1;
	Real inflowVal = 0.0;

	/* parse cli values: use gui? inputvals-filename? what other values to set? */
	if(int doexit = initCLI(setting, gui, inflowVal, simparam, argc, argv) )
		exit(doexit-1);

	/* let bakery create boundaries 'n stuff using start simparams */
	SimulationParameters newparam = Bakery::get(static_cast<Bakery::Setting>(setting), inflowVal, simparam);
	newparam.useComplexGeometry = setting;

	/* replace non-domain parameters in new simparams, to overwrite init-defaults */
	overwritePostBakeryParams(newparam, simparam);

	/* give simparams to gui - gui feedback? */
	if(gui)
	{
		MTQueue<SimulationParameters> outbox;
		MTQueue<SimulationParameters> inbox;

		// Create renderer. Obviously the renderer's inbox is the outbox on this side.
		auto render_execution = std::async(std::launch::async, &runVisualization, std::ref(outbox), std::ref(inbox), 640, 480, std::ref(newparam));

		/* TODO: get/overwrite newparams from gui */
		while(true)
		{
			SimulationParameters received_params;
			if(inbox.tryPop(received_params,std::chrono::milliseconds(500)))
			{
				std::cout<<"I got the goods!"<<std::endl;

				newparam = Bakery::get(static_cast<Bakery::Setting>(received_params.useComplexGeometry), inflowVal, received_params);

				outbox.push(newparam);
			}

			auto status = render_execution.wait_for(std::chrono::seconds(0));
			if(status == std::future_status::ready)
				break;
		}

	}

	if(newparam.name == "")
	{
		newparam.writeToSTDOUT();
	}
	else
	{
		newparam.writeSettingsFile(newparam.name);
		log_info(">> config file written to \"%s\"", newparam.name.c_str());
	}

	return 0;
}
