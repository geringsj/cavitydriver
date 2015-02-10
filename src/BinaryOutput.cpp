#include "BinaryOutput.hpp"


BinaryOutput::BinaryOutput(Domain& domain, const std::string outputpath) : m_domain(domain), m_output_path(outputpath), m_framestep(0)
{
}

BinaryOutput::~BinaryOutput()
{

}

void BinaryOutput::write()
{
	writePfm();
	m_framestep++;
}

void BinaryOutput::writePfm()
{
	/* generate filename with current framestep */
	std::string filename;
	filename.append("./");
	filename.append (this->m_output_path);
	filename.append("/");
	filename.append ("field_");
	filename.append (std::to_string(this->m_framestep));
	filename.append (".pfm");

	int xMin = this->m_domain.getWholeInnerRange().begin.i;
	int xMax = this->m_domain.getWholeInnerRange().end.i;
	int yMin = this->m_domain.getWholeInnerRange().begin.j;
	int yMax = this->m_domain.getWholeInnerRange().end.j;

	std::ofstream  out;
	out.open(const_cast < char *>(filename.c_str ()));

	out << "PF" << std::endl;
	//out << (xMax-xMin)+1 << " " << (yMax-yMin)+1 << std::endl;
	out << (xMax-xMin)+3 << " " << (yMax-yMin)+3 << std::endl;
	out << -1 << std::endl; // TODO assuming little endian here... use better (actually check) in future

	out.close();
	
	//unsigned int data_size = ((xMax-xMin)+1) * ((yMax-yMin)+1) * 3;
	unsigned int cell_entries = 4;
	unsigned int data_size = ((xMax-xMin)+3) * ((yMax-yMin)+3) * cell_entries;
	float* data = new float[data_size];

	unsigned int data_index = 0;

	//for(int j = yMin; j <= yMax; ++j)
	//{
	//	for(int i = xMin; i <= xMax; ++i)
	//	{
	for(int j = yMin-1; j <= yMax+1; ++j)
	{
		for(int i = xMin-1; i <= xMax+1; ++i)
		{
			if(i == xMin-1)
				data[data_index++] = (float)m_domain.u()(i, j);
			else if(i == xMax+1)
				data[data_index++] = (float)m_domain.u()(i-1, j);
			else
				data[data_index++] = (float)(m_domain.u()(i, j) + m_domain.u()(i - 1, j)) / 2.0;

			if(j==yMin-1)
				data[data_index++] = (float)m_domain.v()(i, j);
			else if(j==yMax+1)
				data[data_index++] = (float)m_domain.v()(i, j-1);
			else
				data[data_index++] = (float)(m_domain.v()(i, j) + m_domain.v()(i, j - 1)) / 2.0;

			data[data_index++] = (float)m_domain.p()(i,j);

			data[data_index++] = (float)m_domain.t()(i,j);
		}
	}

	out.open (const_cast < char *>(filename.c_str ()), std::ios::app | std::ios::binary);
	out.write(reinterpret_cast<const char*>(data), sizeof(float)*data_size);
	out.close();

	delete[] data;
}