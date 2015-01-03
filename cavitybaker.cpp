#include "src/CavityRenderer.hpp"
#include "src/SimulationParameters.hpp"
#include "src/Debug.hpp"

#if defined(__linux)
	#include "external/include/docopt.h"
#endif

int main(int argc, char** argv)
{
	argc = argc*(1+0*(**argv)); /* just to get rid of some warnings */

	/* parse cli values: use gui? inputvals-filename? what other values to set? */
#if defined(__linux)
	static const char USAGE[]=
		"
		Usage: cavitybaker (SETTING) ... // TODO
		";
	std::map<std::string, docopt::value> args = 
		docopt::docopt(USAGE, { argv + 1, argv + argc });
#else
#endif

	/* move cli domain params into a start simparams */

	/* let bakery create boundaries 'n stuff using start simparams */

	/* replace non-domain parameters in simparams, to overwrite init-defaults */

	/* maybe give simparams to gui - gui feedback? */

	/**
	 * This is just a test but the range should
	 * match the inner range of p.
	 */
	SimulationParameters simparam("inputvals");
	Index begin = Index(1, 1);
	Index end = Index(simparam.iMax, simparam.jMax);
	Range range = Range(begin, end);
	simparam.writeToSTDOUT();
	/**
	 * This is just a test but the range should
	 * match the inner range of p.
	 */

	CavityRenderer cavity_renderer;
	if (cavity_renderer.initBakeryVis(640, 480, simparam))
	{
		cavity_renderer.createGrid(range);
		cavity_renderer.paint();
	}

	return 0;
}
