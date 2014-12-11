
#include "VTKOutput.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

VTKOutput::VTKOutput(
		Domain& domain, 
		const std::string outputpath) 
	:framestep(0), output_path(outputpath), domain(domain), communication(NULL)
{
	this->checkOutputPath();
}

VTKOutput::VTKOutput(
		Domain& domain, 
		const std::string outputpath, 
		const Communication& comm)
	:framestep(0), output_path(outputpath), domain(domain), communication(&comm)
{
	this->checkOutputPath();
}

VTKOutput::~VTKOutput()
{
}

#if defined(__linux)
	#include "sys/stat.h"
#endif
void VTKOutput::checkOutputPath()
{
	std::string filename;
	filename.append("./");
	filename.append(this->output_path);
	filename.append("/");

	// Test if the directory exits, if not create it.
	std::filebuf test_dir;
	test_dir.open(const_cast < char *>(filename.c_str()), std::ios::out);
	if (!test_dir.is_open())
	{
		// Directory doesn't exist.
		#if defined(_WIN64)
			CreateDirectory(filename.c_str(), NULL);
		#elif defined(_WIN32)
			CreateDirectory(filename.c_str(), NULL);
		#elif defined(__linux)
			mkdir(filename.c_str(), 0700);
		#endif
	}
}

void VTKOutput::writeVTKFile()
{
	/* generate filename with current framestep */
	std::string filename;
  	filename.append("./");
  	filename.append (this->output_path);
  	filename.append("/");
  	filename.append ("field_");
  	filename.append (std::to_string(this->framestep));
  	filename.append (".vts");

  	/* for next VTK write, increase framestep */
  	this->framestep++;

  	/* open stream for VTK XML data */
  	std::filebuf fb;
  	fb.open (const_cast < char *>(filename.c_str ()), std::ios::out);
  	std::ostream os (&fb);

  	/* write extent of grid, the Min/Max values are _inclusive_, 
  	 * like in for(int i=Min; i<=Max; ...) */
  	int xMin = this->domain.getBeginInnerDomains().i;
  	int xMax = this->domain.getEndInnerDomainP().i;
  	int yMin = this->domain.getBeginInnerDomains().j;
  	int yMax = this->domain.getEndInnerDomainP().j;
	os 
	<< "<?xml version=\"1.0\"?>" << std::endl
	<< "<VTKFile type=\"StructuredGrid\">" << std::endl
	<< "<StructuredGrid WholeExtent=\""
	<< xMin << " " << xMax << " "
	<< yMin << " " << yMax << " "
	<< "0" << " " << "0" << " "
	<< "\" GhostLevel=\"" << "0" << "\">" << std::endl
	<< "<Piece Extent=\""
	<< xMin << " " << xMax << " "
	<< yMin << " " << yMax << " "
	<< "0" << " " << "0" << " "
	<< "\">" << std::endl
	
	/* write floating point coordinates of grid cell points */
	<< "<Points>" << std::endl 
	<< "<DataArray type=\"Float64\" format=\"ascii\" NumberOfComponents=\"3\"> " << std::endl;
	Real deltaX = this->domain.getDelta().x;
  	Real deltaY = this->domain.getDelta().y;
	for(int i = xMin; i <= xMax; i++)
	{
		for(int j = yMin; j <= yMax; j++)
		{
			os << std::scientific 
				<< i * deltaX << " " 
				<< j * deltaY << " " 
				<< 0.0 << std::endl;
		}
	}
	os 
	<< "</DataArray>" << std::endl
	<< "</Points>" << std::endl

	/* write inner of vector field U/V 
	 * here we interpolate U and V to be at the grid-position of P,
	 * so we 'implicitly' touch the boundary values of U and V
	 * when looping using the ranges of P */
	<< "<PointData Vectors=\"Field\"  Scalars=\"P\">" << std::endl 
	<< "<DataArray Name=\"Field\" NumberOfComponents=\"3\" type=\"Float64\" >" << std::endl;
	for(int i = xMin; i <= xMax; ++i)
	{
		for(int j = yMin; j <= yMax; ++j)
		{
			os << std::scientific 
				<< (domain.u()(i, j) + domain.u()(i - 1, j)) / 2.0 << " " 
				<< (domain.v()(i, j) + domain.v()(i, j - 1)) / 2.0 << " " 
				<< 0.0 << std::endl;
		}
	}

	/* write inner of P, no weird stuff happening here */
	os
	<< "</DataArray>" << std::endl
	<< "<DataArray type=\"Float64\" Name=\"P\" format=\"ascii\">" << std::endl;
	for(int i = xMin; i <= xMax; ++i)
	{
		for(int j = yMin; j <= yMax; ++j)
		{
			os << std::scientific << domain.p()(i, j) << " ";
		}
		os << std::endl;
	}
	os 
	<< "</DataArray>" << std::endl
	<< "</PointData>" << std::endl
	
	<< "</Piece>" << std::endl
	<< "</StructuredGrid>" << std::endl 
	<< "</VTKFile>" << std::endl;
	fb.close ();
}


void VTKOutput::writeVTKFileParallel()
{
	if(this->communication->getRank() == 0)
		this->writeVTKMasterFile();
	this->writeVTKSlaveFile();
	this->framestep++;
}


void VTKOutput::writeVTKMasterFile()
{
	std::string filename;
	filename.append("./");
	filename.append(this->output_path);
	filename.append("/");
	filename.append("field_");
	filename.append(std::to_string(this->framestep));
	filename.append(".pvtr");

	std::filebuf fb;
	fb.open(const_cast < char *>(filename.c_str()), std::ios::out);
	std::ostream os(&fb);

	/* all these values are inclusive. like in 'for(int i=min; i <= max; i++)' */
	/* also these values should start and end at the INNER of the domain.
	 * when/if needed, the offset for boundaries will be added as -1/+1 */
	int xGlobMin = 1; 
	int xGlobMax = this->communication->getGlobalDimensions().i; 
	int yGlobMin = 1; 
	int yGlobMax = this->communication->getGlobalDimensions().j; 

	os << "<?xml version=\"1.0\"?>" << std::endl;
	os << "<VTKFile type=\"PRectilinearGrid\">" << std::endl;
	/* write extends */
	os << "<PRectilinearGrid WholeExtent=\""
		<< xGlobMin <<" "<< xGlobMax << " " 
		<< yGlobMin <<" "<< yGlobMax 
		<< " 0 0\" GhostLevel=\"" << 0 << "\">"<< std::endl;
	os << "<PCoordinates>" << std::endl;
	os << "<PDataArray type=\"Float64\"/>" << std::endl;
	os << "<PDataArray type=\"Float64\"/>" << std::endl;
	os << "<PDataArray type=\"Float64\"/>" << std::endl;
	os << "</PCoordinates>" << std::endl;

	/* iterate the local domains and write out their extents */
	/* the following Locl variables are analogous to the above "Glob" ones, 
	 * but for the local domain of the process "localProcRank" */
	int xLoclbMin, xLoclbMax, yLoclbMin, yLoclbMax; 
	int localProcRank;
	/* loop through all processes */
	for (int x = 0; x < this->communication->getProcsGridDim().i; x++)
	{
		for (int y = 0; y < this->communication->getProcsGridDim().j; y++)
		{
			/* compute 'local' process rank and its dimensions from position in 
			 * grid of processes */
			localProcRank = x * this->communication->getProcsGridDim().j + y; 

			Index expctdCellsPerDim(
					(this->communication->getGlobalDimensions().i 
					 / this->communication->getProcsGridDim().i) ,
					(this->communication->getGlobalDimensions().j 
					 / this->communication->getProcsGridDim().j) );

			Dimension localOffsetToGlobalDomain(
					expctdCellsPerDim.i * x ,
					expctdCellsPerDim.j * y );

			Dimension localDomainDim(
					/* last process in x/y dimensions gets the rest of cells in that direction,
					 * that is why the following computation looks ugly */
					expctdCellsPerDim.i 
					+ ( (x < this->communication->getProcsGridDim().i-1)
						?(0):(this->communication->getGlobalDimensions().i 
							% this->communication->getProcsGridDim().i)) ,
					expctdCellsPerDim.j 
					+ ( (y < this->communication->getProcsGridDim().j-1)
						?(0):(this->communication->getGlobalDimensions().j 
							% this->communication->getProcsGridDim().j)) );

			/* here, write extends of sub-domain of every process Omega_{i,j} */
			/* for a consistent visualization ParaView wants 
			 * the subdomains to overlapp, this is why we extend the inner by 
			 * one in all directions. so we later need to write out the boundary too.
			 * because we will later write everything with respect to the pressure
			 * (which sits in the center of each cell), this will work out. 
			 * (we will need to interpolate velocities accordingly) */
			xLoclbMin = localOffsetToGlobalDomain.i + xGlobMin;
			xLoclbMax = xLoclbMin + localDomainDim.i -1; 
			/* 'localDomainDim.i-1'
			 * because we want to go from first cell of domain to last, but not too far */
			yLoclbMin = localOffsetToGlobalDomain.j + yGlobMin;
			yLoclbMax = yLoclbMin + localDomainDim.j -1;

			/* walk _on_ the boundary of domain of current process */
			xLoclbMin += -1;
			xLoclbMax += +1; 
			yLoclbMin += -1;
			yLoclbMax += +1;

			/* write extent of Domain of current process, 
			 * with respect to global Domain */
			os << "<Piece Extent=\"" 
				<< xLoclbMin << " " << xLoclbMax << " " 
				<< yLoclbMin << " " << yLoclbMax << " 0 0\" Source=\"field_" 
				<< std::to_string(this->framestep) << "_processor_"
				<< localProcRank << ".vtr\"/>" << std::endl;
		}
	}

	/* standard fields: */
	os << "<PPointData Vectors=\"Field\" Scalars=\"P\">"<< std::endl; // , Vorticity, Stream
	os << "<PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Field\" format=\"ascii\"/>"<< std::endl;
	os << "<PDataArray type=\"Float64\" Name=\"P\" format=\"ascii\"/>"<< std::endl;

	/* new fields: vorticity and stream */
	//os << "<PDataArray type=\"Float64\" Name=\"Vorticity\" format=\"ascii\"/>"<< std::endl;
	//os << "<PDataArray type=\"Float64\" Name=\"Stream\" format=\"ascii\"/>"<< std::endl;

	os << "</PPointData>" << std::endl;
	os << "</PRectilinearGrid>" << std::endl;
	os << "</VTKFile>" << std::endl;
}



void VTKOutput::writeVTKSlaveFile()
{
	int localProcRank = this->communication->getRank();
	std::string filename;
	filename.append("./");
	filename.append(this->output_path);
	filename.append("/");
	filename.append("field_");
	filename.append(std::to_string(this->framestep));
	filename.append("_processor_");
	filename.append(std::to_string(localProcRank));
	filename.append(".vtr");

	std::filebuf fb;
	fb.open(const_cast < char *>(filename.c_str()), std::ios::out);
	std::ostream os(&fb);

	int xGlobMin = 1; 
	//int xGlobMax = this->communication->getGlobalDimensions().i; 
	int yGlobMin = 1; 
	//int yGlobMax = this->communication->getGlobalDimensions().j; 

	int xLoclbMin, xLoclbMax, yLoclbMin, yLoclbMax; 
	/* set the four indices like in WriteVTKMasterFile */
	xLoclbMin = this->communication->getOwnOffsetToGlobalDomain().i + xGlobMin;
	xLoclbMax = xLoclbMin + domain.getDimension().i -1; 
	yLoclbMin = this->communication->getOwnOffsetToGlobalDomain().j + yGlobMin;
	yLoclbMax = yLoclbMin + domain.getDimension().j -1;
	/* walk _on_ the boundary of domain of current process */
	xLoclbMin += -1;
	xLoclbMax += +1; 
	yLoclbMin += -1;
	yLoclbMax += +1;

	os << "<?xml version=\"1.0\"?>" << std::endl;
	os << "<VTKFile type=\"RectilinearGrid\">" << std::endl;
	// whole extent = whole extent in master
	// os << "<RectilinearGrid WholeExtent=\"0 " << s.iMax - 1 << " 0 "
	// 		<< s.jMax - 1 << " 0 0\" GhostLevel=\"" << stencilwidth << "\">"
	// 		<< std::endl;

	// whole extent = piece extent in slave
	os << "<RectilinearGrid WholeExtent=\"" 
		<< xLoclbMin << " " << xLoclbMax << " " 
		<< yLoclbMin << " " << yLoclbMax
		<< " 0 0\" GhostLevel=\"" << 1 << "\">" << std::endl;

	os << "<Piece Extent=\"" 
		<< xLoclbMin << " " << xLoclbMax << " " 
		<< yLoclbMin << " " << yLoclbMax
		<< " 0 0\">" << std::endl;

	/* now we write the coordinates of the grid points in x and y direction */
	os << "<Coordinates>" << std::endl;
	os << "<DataArray type=\"Float64\" format=\"ascii\">" << std::endl;
	Real cellDeltaX = domain.getDelta().x;
	Real cellDeltaY = domain.getDelta().y;
	for (int i = xLoclbMin; i <= xLoclbMax; i++)
		os << std::scientific << i * cellDeltaX << " ";
	os << std::endl;
	os << "</DataArray>" << std::endl;
	os << "<DataArray type=\"Float64\" format=\"ascii\">" << std::endl;
	for (int j = yLoclbMin; j <= yLoclbMax; j++)
		os << std::scientific << j * cellDeltaY << " ";
	os << std::endl;
	os << "</DataArray>" << std::endl;
	os << "<DataArray type=\"Float64\" format=\"ascii\">0 0</DataArray>" << std::endl;
	os << "</Coordinates>" << std::endl;

	os << "<PointData Vectors=\"Field\" Scalars=\"P\">" << std::endl; // , Vorticity, Stream
	os << "<DataArray Name=\"Field\" NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">" << std::endl;

	/* write grid values with boundaries.
	 * interpolate u and v to be at p position of cell.
	 * we have to make up values on the boundary just to satisfy overlapping
	 * conditions of sub domains in the VTK file. */
	/* the {x,y}Locl{Min,Max} are w.r.t. the offset of the global domain, 
	 * so we can not use them as loop boundaries here */
	Real u_inter, v_inter;
	for (int j = 1 -1; j <= domain.getEndInnerDomainP().i +1; j++)
	{
		for (int i = 1 -1; i <= domain.getEndInnerDomainP().j +1; i++)
		{
			if(i == 0)
				u_inter = (domain.u()(i, j) + domain.u()(i , j)) / 2.0;
			else 
				if(i == domain.getEndInnerDomainP().i+1)
					u_inter = (domain.u()(i - 1, j) + domain.u()(i - 1, j)) / 2.0;
				else
					u_inter = (domain.u()(i, j) + domain.u()(i - 1, j)) / 2.0;

			if(j == 0 )
				v_inter = (domain.v()(i, j) + domain.v()(i, j)) / 2.0;
			else
				if(j == domain.getEndInnerDomainP().j +1)
					v_inter = (domain.v()(i, j - 1) + domain.v()(i, j - 1)) / 2.0;
				else
					v_inter = (domain.v()(i, j) + domain.v()(i, j - 1)) / 2.0;

			os << std::scientific <<
				u_inter << " " <<
				v_inter << " " <<
				0.0 << std::endl;
		}
		os << std::endl;
	}
	os << "</DataArray>" << std::endl;

	/* write out p as it is */
	os << "<DataArray type=\"Float64\" Name=\"P\" format=\"ascii\">"<< std::endl;
	for (int j = 1 -1; j <= domain.getEndInnerDomainP().i +1; j++)
	{
		for (int i = 1 -1; i <= domain.getEndInnerDomainP().j +1; i++)
		{
			os << std::scientific << domain.p()(i, j) << " ";
		}
		os << std::endl;
	}

	os << "</DataArray>" << std::endl;

	//	os << "<DataArray type=\"Float64\" Name=\"vorticity\" format=\"ascii\">"
	//		<< std::endl;
	//
	//	/* TODO TODO */
	//	for (int i = 0; i < iMax + 2; i++)
	//	for (int j = 1 -1; j <= domain.getEndInnerDomainP().i +1; j++)
	//	{
	//		os << std::scientific << 0.0 << " ";
	//	}
	//	os << std::endl;
	//	for (int j = 1; j < jMax+1; j++)
	//	{
	//		os << std::scientific << 0.0 << " ";
	//		for (int i = 1; i < iMax+1; i++)
	//		{
	//			Real zeta = (domain.u()(i, j + 1) - domain.u()(i, j)) / deltaY
	//				- (domain.v()(i + 1, j) - domain.v()(i, j)) / deltaX;
	//			os << std::scientific << zeta << " ";
	//		}
	//		os << std::scientific << 0.0 << std::endl;
	//	}
	//	// we need one additional line of vorticity...
	//	for (int i = 0; i < iMax + 2; i++)
	//	{
	//		os << std::scientific << 0.0 << " ";
	//	}
	//	os << std::endl;
	//
	//	os << "</DataArray>" << std::endl;
	//
	//	os << "<DataArray type=\"Float64\" Name=\"stream\" format=\"ascii\">"
	//		<< std::endl;
	//
	//	Dimension stream_dim(iMax + 2, jMax + 2);
	//	GridFunction stream(stream_dim);
	//
	//	// we need one initial line of stream...
	//	for (int i = 0; i < jMax + 2; i++)
	//	{
	//		stream(i, 0) = 0.0;
	//	}
	//	for (int j = 1; j < jMax + 2; j++)
	//	{
	//		for (int i = 0; i < iMax + 2; i++)
	//		{
	//			stream(i, j) = stream(i, j - 1)
	//				+ domain.u()(i, j) * deltaY;
	//		}
	//	}
	//	os << std::endl;
	//	for (int j = 0; j < jMax + 2; j++)
	//	{
	//		//os << std::scientific << 0.0 << " ";
	//		for (int i = 0; i < iMax + 2; i++)
	//		{
	//			os << std::scientific << stream(i, j) << " ";
	//		}
	//		os << std::endl;
	//	}
	//	
	//	os << "</DataArray>" << std::endl;


	os << "</PointData>" << std::endl;
	os << "</Piece>" << std::endl;

	os << "</RectilinearGrid>" << std::endl;
	os << "</VTKFile>" << std::endl;
}

