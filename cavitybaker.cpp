#include "src/CavityRenderer.hpp"
#include "src/SimulationParameters.hpp"
#include "src/Debug.hpp"

int main(int argc, char** argv)
{
	argc = argc*(1+0*(**argv)); // just to get rid of some warnings */

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
