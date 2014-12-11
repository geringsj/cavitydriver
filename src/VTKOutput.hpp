#ifndef VTKWriter_hpp
#define VTKWriter_hpp

#include "Domain.hpp"
#include "Communication.hpp"

#include <string>

class VTKOutput{
private:

	unsigned int framestep;

	const std::string output_path;

	Domain& domain;

	const Communication* const communication;
	
	void checkOutputPath();
	
	void writeVTKMasterFile();
	void writeVTKSlaveFile();

public:

	VTKOutput(Domain& domain, const std::string outputpath);
	VTKOutput(Domain& domain, const std::string outputpath, const Communication& comm);
	~VTKOutput();

	void writeVTKFile();

	void writeVTKFileParallel();
};

#endif
