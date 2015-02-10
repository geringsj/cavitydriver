#include "src/CavityRenderer.hpp"
#include "src/SimulationParameters.hpp"
#include "src/Debug.hpp"

int main(int argc, char** argv)
{
	argc = argc*(1+0*(**argv)); /* just to get rid of some warnings */

	if(argc != 2)
	{
		log_err("provide a config file to the program!");
		exit(EXIT_FAILURE);
	}
	SimulationParameters simparam(argv[1]);


	// Don't really need these, but current implementation of Renderer is quite chatty
	MTQueue<SimulationParameters> outbox;
	MTQueue<SimulationParameters> inbox;

	CavityRenderer cavity_renderer(inbox,outbox);
	if( !cavity_renderer.initPainterVis(800,450,simparam,"out/field_") )
	{
		log_err("could not start painter");
		exit(EXIT_FAILURE);
	}
	cavity_renderer.paint();

	return 0;
}
