#include "src/IO.hpp"
#include "src/Domain.hpp"

// default values for the output directory and the settings file
const char *output = "out";
const char *settings = "inputvals";

/**
 * Print the help text in order to show what
 * arguments can be processed.
 */
void dieSynopsis() {
	std::cerr << "Synopsis: IO "
		<< "--out \"output directory\" --set \"path to the settings file\" ";
	exit(-1);
}

/**
 * Parse the input arguments. If there are no arguments
 * the default values are used.
 * @param the length of the input
 * @param the arguments
 */
void parseArguments(int argc, char** argv)
{
	--argc; ++argv; // We don't care about program name;
	while (argc > 0 && argv[0][0] == '-')
	{
		if (argv[0] == (std::string) "--help" || argv[0] == (std::string) "-h")
		{
			dieSynopsis();
		}
		else if (argv[0] == (std::string)"--out")
		{
			if (argv[1] != NULL && argv[1] != (std::string)"--set")
			{
				output = argv[1];
				--argc; ++argv;
				--argc; ++argv;
			}
			else --argc; ++argv;
		}
		else if (argv[0] == (std::string)"--set")
		{
			if (argv[1] != NULL)
			{
				settings = argv[1];
				--argc; ++argv;
				--argc; ++argv;
			}
			else --argc; ++argv;
		}
	}
}

int main(int argc, char** argv)
{
	parseArguments(argc, argv);
	std::cout << settings << " " << output << std::endl;
	IO io(settings, output);
	//IO io((char*)"inputvals", (char*)"out");
	io.writeSimParamToSTDOUT();

	//Dimension dim;
	//Delta delta;
	//Domain domain(dim,delta);

	return 0;
}
