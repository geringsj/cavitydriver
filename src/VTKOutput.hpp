#ifndef VTKWriter_hpp
#define VTKWriter_hpp

#define NOMINMAX // this is for the windows folks

#include "Domain.hpp"
#include "Communication.hpp"

#include <algorithm>
#include <string>

class VTKOutput{
private:

	unsigned int framestep;

	const std::string output_path;

	Domain& domain;

	const Communication* const communication;

	std::vector<std::vector<Point>> m_particles;
	
	void checkOutputPath();
	
	void writeVTKSingleFile();

	void writeVTKMasterFile();
	void writeVTKSlaveFile();

	void writeParticleVTPFile();
	void writeStreamlineVTPFile();

	void interpolateUV(Point p, Real& u, Real& v);

public:

	VTKOutput(Domain& domain, const std::string outputpath);
	VTKOutput(Domain& domain, const std::string outputpath, const Communication& comm);
	~VTKOutput();

	void writeVTKFile();
};

#endif
