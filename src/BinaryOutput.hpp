#ifndef BinaryOutput_hpp
#define BinaryOutput_hpp

#include "Domain.hpp"

#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

class BinaryOutput
{
public:
	BinaryOutput(Domain& domain, const std::string outputpath);
	~BinaryOutput();
	void write();
private:
	Domain& m_domain;
	const std::string m_output_path;

	unsigned int m_framestep;

	void writePfm();
};

#endif