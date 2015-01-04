#include "src/CavityRenderer.hpp"
#include "src/SimulationParameters.hpp"
#include "src/Debug.hpp"

#if defined(__linux)
	#include "external/include/docopt.h"
#endif

int main(int argc, char** argv)
{
	argc = argc*(1+0*(**argv)); /* just to get rid of some warnings */

	SimulationParameters simparam;
	bool gui = false;
	int setting = 0;

	/* parse cli values: use gui? inputvals-filename? what other values to set? */
#if defined(__linux)
	static const char USAGE[] = R"(
Usage: cavitybaker SETTING [options]

SETTING:
 A name identifying the setting to generate.
 Can be one of the following (case insensitive):
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
  --name=STRING     Nam of the file in which the configuration parameters
                    will be saved. If an empty string is given, the 
                    parameters file will be writen to stdout. [default: ""]
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
                    is used. Specified in radiant. [default: 0.7854 == PI/4.0]
  --kow=FLOAT       Width of object for Karman Vortex Street SETTING. 
                    Will be ignored if other SETTING than 'ObstacleChannelFlow' 
                    is used. [default: 1.5*yLength/jMax]
)";

	std::map<std::string, docopt::value> args = 
		docopt::docopt(USAGE, { argv + 1, argv + argc });

	/* move cli domain params into a start simparams */
	for(auto const& arg : args)
		std::cout << arg.first << " is " <<  arg.second << std::endl;
#else
#endif

	/* let bakery create boundaries 'n stuff using start simparams */

	/* replace non-domain parameters in simparams, to overwrite init-defaults */

	/* maybe give simparams to gui - gui feedback? */

	if(gui)
	{
	/**
	 * This is just a test but the range should
	 * match the inner range of p.
	 */
	Index begin = Index(1, 1);
	Index end = Index(simparam.iMax, simparam.jMax);
	Range range = Range(begin, end);

		CavityRenderer cavity_renderer;
		if (cavity_renderer.initBakeryVis(640, 480, simparam))
		{
			cavity_renderer.createGrid(range);
			cavity_renderer.paint();
		}
	}

	return 0;
}
