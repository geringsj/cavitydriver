#include "src/CavityRenderer.hpp"
#include "src/SimulationParameters.hpp"
#include "src/Debug.hpp"

#include "src/Bakery.hpp"

#include <string>

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
	KOW
};
const option::Descriptor usage[] =
{
	{UNKNOWN, 0, "", "", option::Arg::None, 
R"(
Usage: cavitybaker SETTING [options]

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
"  --name=STRING  \tName of the file in which the configuration parameters will be saved. If no name is given, the parameters file will be written to stdout." },

 {XLENGTH, 0, "", "xLength", NonEmpty,
"  --xLength=FLOAT  \tLength of domain in x-direction." },
 {YLENGTH, 0, "", "yLength", NonEmpty,
"  --yLength=FLOAT  \tLength of domain in y-direction." },
 {IMAX, 0, "", "iMax", NonEmpty,
"  --iMax=INT  \tNumber of cells the discrete grid will have in x-direction." },
 {JMAX, 0, "", "jMax", NonEmpty,
"  --jMax=INT  \tNumber of cells the discrete grid will have in y-direction." },
 {TAU, 0, "", "tau", NonEmpty,
"  --tau=FLOAT  \tScaling factor for timestep computations in (0,1]." },
 {TEND, 0, "", "tEnd", NonEmpty,
"  --tEnd=FLOAT  \tLength of simulation in seconds." },
 {DELTAT, 0, "", "deltaT", NonEmpty,
"  --deltaT=FLOAT  \tTimestep size. Will be ignored." },
 {DELTAVEC, 0, "", "deltaVec", NonEmpty,
"  --deltaVec=FLOAT  \tTime between two VTK Outputs of current state of simulation." },
 {EPS, 0, "", "", NonEmpty,
"  --eps=FLOAT  \tUpper bound for residual for SOR solver." },
 {OMG, 0,  "", "omg", NonEmpty,
"  --omg=FLOAT  \tOmega factor for SOR calculation." },
 {ALPHA, 0, "", "alpha", NonEmpty,
"  --alpha=FLOAT  \tMixing factor for Donor-Cell scheme for next velocities." },
 {ITERMAX, 0, "", "iterMax", NonEmpty,
"  --iterMax=INT  \tMaximum number of iterations for SOR solver." },
 {RE, 0, "", "re", NonEmpty,
"  --re=INT  \tReynolds number to be used." },
 {GX, 0, "", "gx", NonEmpty,
"  --gx=FLOAT  \tExternal force in x direction." },
 {GY, 0, "", "gy", NonEmpty,
"  --gy=FLOAT  \tExternal force in y direction." },
 {UI, 0, "", "ui", NonEmpty,
"  --ui=FLOAT  \tConstant initial value for U-velocity field." },
 {VI, 0, "", "vi", NonEmpty,
"  --vi=FLOAT  \tConstant initial value for V-velocity field." },
 {PI, 0, "", "pi", NonEmpty,
"  --pi=FLOAT  \tConstant initial value for Pressure field." },
 {KA, 0, "", "ka", NonEmpty,
"  --ka=FLOAT  \tAngle of object for Karman Vortex Street SETTING specified in radiant. Will be ignored if other SETTING than 'ObstacleChannelFlow is used." },
 {KOW, 0, "", "kow", NonEmpty,
"  --kow=FLOAT  \tWidth of object for Karman Vortex Street SETTING. Will be ignored if other SETTING than 'ObstacleChannelFlow is used." }
};

void getMaybeCLIDouble(Real& real, option::Option* options, optionIndex ind)
{
	if(options[ind])
		real = std::strtod(options[ind].arg,0);
}
void getMaybeCLIInt(int& _int, option::Option* options, optionIndex ind)
{
	if(options[ind])
		_int = std::stoi(options[ind].arg,0);
}

void initCLI(int& setting, bool& gui, SimulationParameters& simparam, int argc, char** argv)
{
  argc-=(argc>0); argv+=(argc>0); // skip program name argv if present
  option::Stats  stats(usage, argc, argv);
  option::Option* options = new option::Option[stats.options_max];
  option::Option* buffer  = new option::Option[stats.buffer_max*100]; /* !! */
  option::Parser parse(true, usage, argc, argv, options, buffer);

	if(parse.error())
	{
		std::cout << "An options-parser error occurred.\n";
		exit(EXIT_FAILURE);
	}


	if(options[HELP])
	{
		option::printUsage(std::cout, usage);
		exit(EXIT_SUCCESS);
	}
	if(argc == 0 || !parse.nonOptionsCount())
	{
		std::cout << "Missing SETTING parameter. Call with '--help' for options.\n";
		exit(EXIT_FAILURE);
	}

	std::string settingstr = std::string(parse.nonOption(0));
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

	gui = static_cast<bool>(options[GUI].count());
	if(!gui && (setting < 0 || setting > 4))
	{
		std::cout << "Invalid SETTING parameter. Call with '--help' for options.\n";
		exit(EXIT_FAILURE);
	}

	if(options[NAME])
		simparam.name = std::string(options[NAME].arg);

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

	delete[] options;
	delete[] buffer;
}

void overwritePostBakeryParams(SimulationParameters& newsp, SimulationParameters& clisp)
{
	newsp.name = clisp.name;

	/* if user specifies another re than given in exercise sheet, use users re */
	newsp.re = clisp.re;
}

int main(int argc, char** argv)
{
	argc = argc*(1+0*(**argv)); /* just to get rid of some warnings */

	SimulationParameters simparam;
	bool gui = false;
	int setting = -1;

	/* parse cli values: use gui? inputvals-filename? what other values to set? */
	initCLI(setting, gui, simparam, argc, argv);

	/* let bakery create boundaries 'n stuff using start simparams */

	/**
	 * Karman:
	 * Re: 1000, max iter: 5000
	 * Step:
	 * Re: 1000, max iter:  100
	 */
	SimulationParameters newparam = Bakery::get(static_cast<Bakery::Setting>(setting), simparam);
	newparam.useComplexGeometry = setting;

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
