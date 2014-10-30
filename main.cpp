#include "src/IO.hpp"

int main(int argc, char** argv)
{
	IO io((char*)"inputvals", (char*)"out");
	io.writeSimParamToSTDOUT();

	return 0;
}
