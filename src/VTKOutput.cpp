
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

	// add empty particle groups based on domain size
	//	Real delta_y = this->domain.getDelta().y;
	int dim_y = this->domain.getDimension().j;

	for(int i=0; i<dim_y/2; i++)
		m_particles.push_back(std::vector<Point>());

}

VTKOutput::VTKOutput(
		Domain& domain, 
		const std::string outputpath, 
		const Communication& comm)
	:framestep(0), output_path(outputpath), domain(domain), communication(&comm)
{
	this->checkOutputPath();

	// add empty particle groups based on domain size
	//	Real delta_y = this->domain.getDelta().y;
	int dim_y = this->domain.getDimension().j;

	for(int i=0; i < (dim_y/2); i++)
		m_particles.push_back(std::vector<Point>());
}

VTKOutput::~VTKOutput()
{
}

#if defined(__linux)
	#include "sys/stat.h"
#endif
#if defined(_WIN64)
	#include <Windows.h>
#endif
#if defined(_WIN32)
	#include <Windows.h>
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

void VTKOutput::writeVTKSingleFile()
{
	/* generate filename with current framestep */
	std::string filename;
	filename.append("./");
	filename.append (this->output_path);
	filename.append("/");
	filename.append ("field_");
	filename.append (std::to_string(this->framestep));
	filename.append (".vts");

	/* open stream for VTK XML data */
	std::filebuf fb;
	fb.open (const_cast < char *>(filename.c_str ()), std::ios::out);
	std::ostream os (&fb);

	/* write extent of grid, the Min/Max values are _inclusive_, 
	 * like in for(int i=Min; i<=Max; ...) */
	int xMin = this->domain.getWholeInnerRange().begin.i;
	int xMax = this->domain.getWholeInnerRange().end.i;
	int yMin = this->domain.getWholeInnerRange().begin.j;
	int yMax = this->domain.getWholeInnerRange().end.j;
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
	for(int j = yMin; j <= yMax; j++)
	{
		for(int i = xMin; i <= xMax; i++)
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
	for(int j = yMin; j <= yMax; ++j)
	{
		for(int i = xMin; i <= xMax; ++i)
		{
			os << std::scientific 
				<< (domain.u()(i, j) + domain.u()(i - 1, j)) / 2.0 << " " 
				<< std::scientific 
				<< (domain.v()(i, j) + domain.v()(i, j - 1)) / 2.0 << " " 
				<< 0.0 << std::endl;
		}
	}

	/* write inner of P, no weird stuff happening here */
	os
	<< "</DataArray>" << std::endl
	<< "<DataArray type=\"Float64\" Name=\"P\" format=\"ascii\">" << std::endl;
	for(int j = yMin; j <= yMax; ++j)
	{
		for(int i = xMin; i <= xMax; ++i)
		{
			os << std::scientific << domain.p()(i, j) << " ";
		}
		os << std::endl;
	}

	/* write the vorticity */
	os 
	<< "</DataArray>" << std::endl;


	os << "<DataArray type=\"Float64\" Name=\"Vorticity\" format=\"ascii\">" << std::endl;
	for (int j = yMin; j <= yMax; ++j)
	{
		for (int i = xMin; i <= xMax; ++i)
		{
			os
			<< std::scientific <<
			(domain.u()(i, j) - domain.u()(i, j - 1)) / domain.getDelta().y -
			(domain.v()(i, j) - domain.v()(i - 1, j)) / domain.getDelta().x
			<< " ";
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


void VTKOutput::writeVTKFile(Real timestep)
{
#ifdef WITHMPI
	if(this->communication->getRank() == 0)
		this->writeVTKMasterFile();
	this->writeVTKSlaveFile();
#else
	this->writeVTKSingleFile();
	this->writeParticleVTPFile(timestep);
	this->writeStreamlineVTPFile();
#endif
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
	int xGlobMin = communication->getGlobalInnerRange().begin.i;
	int xGlobMax = communication->getGlobalInnerRange().end.i;
	int yGlobMin = communication->getGlobalInnerRange().begin.j;
	int yGlobMax = communication->getGlobalInnerRange().end.j;

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
	/* loop through all processes, here we assume that processes have ranks in
	 * [0 .. Communication.getProcsCount() ) */
	for(int rank = 0; rank < this->communication->getProcsCount(); rank++)
	{
		/* get sub-range of process from communication, and the rank of process */
		int localRank = rank;
		Range localSubRange = 
			this->communication->getProcLocalInnerRange(localRank);

		/* here, write extends of sub-domain of every process Omega_{i,j} */
		/* for a consistent visualization ParaView wants 
		 * the subdomains to overlapp, this is why we extend the inner by 
		 * one in all directions. so we later need to write out the boundary too.
		 * because we will later write everything with respect to the pressure
		 * (which sits in the center of each cell), this will work out. 
		 * (we will need to interpolate velocities accordingly) */
		xLoclbMin = localSubRange.begin.i ;
		xLoclbMax = localSubRange.end.i ;
		yLoclbMin = localSubRange.begin.j;
		yLoclbMax = localSubRange.end.j;
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
			<< localRank << ".vtr\"/>" << std::endl;
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

	int xLoclbMin, xLoclbMax, yLoclbMin, yLoclbMax; 
	/* set the four indices like in WriteVTKMasterFile */
	/* get sub-range of process from communication */
	Range localSubRange = this->communication->getLocalInnerRange();
	xLoclbMin = localSubRange.begin.i ;
	xLoclbMax = localSubRange.end.i ;
	yLoclbMin = localSubRange.begin.j;
	yLoclbMax = localSubRange.end.j;
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
		os << std::scientific << i * cellDeltaX + cellDeltaX/2.0 << " ";
	os << std::endl;
	os << "</DataArray>" << std::endl;
	os << "<DataArray type=\"Float64\" format=\"ascii\">" << std::endl;
	for (int j = yLoclbMin; j <= yLoclbMax; j++)
		os << std::scientific << j * cellDeltaY + cellDeltaY/2.0 << " ";
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
  	int xBegin = this->domain.getWholeInnerRange().begin.i -1;
  	int xEnd = this->domain.getWholeInnerRange().end.i +1;
  	int yBegin = this->domain.getWholeInnerRange().begin.j -1;
  	int yEnd = this->domain.getWholeInnerRange().end.j +1;
	/* note: the loops first write in x direction, then in y direction */
	for (int j = yBegin; j <= yEnd; j++)
	{
		for (int i = xBegin; i <= xEnd; i++)
		{
			if(i == xBegin)
				u_inter = (domain.u()(i, j) + domain.u()(i , j)) / 2.0;
			else 
				if(i == xEnd)
					u_inter = (domain.u()(i - 1, j) + domain.u()(i - 1, j)) / 2.0;
				else
					u_inter = (domain.u()(i, j) + domain.u()(i - 1, j)) / 2.0;

			if(j == yBegin)
				v_inter = (domain.v()(i, j) + domain.v()(i, j)) / 2.0;
			else
				if(j == yEnd)
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
	/* note: the loops first write in x direction, then in y direction */
	for (int j = yBegin; j <= yEnd; j++)
	{
		for (int i = xBegin; i <= xEnd; i++)
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


void VTKOutput::writeParticleVTPFile(Real timestep)
{
	// compute new positions for existing particles
	for(auto& particle_group : m_particles)
	{
		for(auto& particle : particle_group)
		{
			Real u,v;
			interpolateUV(particle,u,v);

			particle.x += timestep * u;
			particle.y += timestep * v;
		}
	}

	// add new particles at start position
	Real start_x = this->domain.getDelta().x;
	Real start_y = 0.0;
	for(auto& particle_group : m_particles)
	{
		start_y += 2.0* (this->domain.getDelta().y);

		particle_group.push_back(Point(start_x,start_y,0.0));
	}


	/* generate filename with current framestep */
	std::string filename;
  	filename.append("./");
  	filename.append (this->output_path);
  	filename.append("/");
  	filename.append ("particles_");
  	filename.append (std::to_string(this->framestep));
  	filename.append (".vtp");

  	/* open stream for VTK XML data */
  	std::filebuf fb;
  	fb.open (const_cast < char *>(filename.c_str ()), std::ios::out);
  	std::ostream os (&fb);

	os 
	<< "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl
	<< "<PolyData>" << std::endl;
	
	for(auto& particle_group : m_particles)
	{
		os
		<< "<Piece NumberOfPoints=\""<< particle_group.size() <<"\" "
		<< "NumberOfVerts=\""<< particle_group.size() <<"\" "
		<< "NumberOfLines=\""<< particle_group.size()-1<<"\""
		<<">"
		<<std::endl;

		os
		<< "<Points>" <<std::endl
		<< "<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" <<std::endl;

		for(auto& particle : particle_group)
		{
			os
			<< particle.x <<" "<< particle.y <<" "<< particle.z <<std::endl;
		}

		os
		<< "</DataArray>" <<std::endl
		<< "</Points>" <<std::endl;

		os
		<< "<Verts>" <<std::endl
		<< "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" <<std::endl;
		for(long i=0; i<(long)particle_group.size(); i++)
			os<< i <<std::endl;
		os
		<< "</DataArray>" <<std::endl
		<< "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" <<std::endl;
		for(long i=1; i<=(long)particle_group.size(); i++)
			os<< i <<std::endl;
		os
		<< "</DataArray>" <<std::endl
		<< "</Verts>" <<std::endl;

		os
		<< "<Lines>" <<std::endl
		<< "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" <<std::endl;
		for(long i=0; i<(long)particle_group.size(); i++)
			os<< i <<" "<< i+1 <<std::endl;
		os
		<< "</DataArray>" <<std::endl
		<< "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" <<std::endl;
		for(long i=2; i<=(long)(particle_group.size()*2); i=i+2)
			os<< i <<std::endl;
		os
		<< "</DataArray>" <<std::endl
		<< "</Lines>" <<std::endl;

		os
		<< "</Piece>" <<std::endl;
	}

	os
	<< "</PolyData>" << std::endl
	<< "</VTKFile>" << std::endl;
	fb.close ();
}

void VTKOutput::writeStreamlineVTPFile()
{
	int dim_x = this->domain.getDimension().i;
	int dim_y = this->domain.getDimension().j;
	Real dx = this->domain.getDelta().x;
	Real dy = this->domain.getDelta().y;
	// verrrryyy average cell size
	Real cell_size = std::sqrt(dx*dx + dy*dy);
	Real domain_size_x = (Real)dim_x * dx;
	Real domain_size_y = (Real)dim_y * dy;

	std::vector<std::vector<Point>> streamlines;

	for(int i=0; i<dim_y/2; i++)
		streamlines.push_back(std::vector<Point>());

	// build streamlines
	Real start_x = dx;
	Real start_y = 0.0;
	for(auto& streamline : streamlines)
	{
		start_y += 2.0*dy;
		streamline.push_back(Point(start_x,start_y,0.0));

		while( 0.0 < streamline.back().x && streamline.back().x < domain_size_x &&
				0.0 < streamline.back().y && streamline.back().y < domain_size_y )
		{
			Point p = streamline.back();
			Real u,v;
			interpolateUV(p,u,v);

			if(!((u+v)>0.0))
				break;

			Real scaling = cell_size / std::sqrt(u*u + v*v);

			Point p_new(p.x + scaling*u, p.y + scaling*v,0.0);
			streamline.push_back(p_new);

			if((int)streamline.size() > std::max(dim_x,dim_y))
				break;
		}
	}

	/* generate filename with current framestep */
	std::string filename;
  	filename.append("./");
  	filename.append (this->output_path);
  	filename.append("/");
  	filename.append ("streamlines_");
  	filename.append (std::to_string(this->framestep));
  	filename.append (".vtp");

  	/* open stream for VTK XML data */
  	std::filebuf fb;
  	fb.open (const_cast < char *>(filename.c_str ()), std::ios::out);
  	std::ostream os (&fb);

	os 
	<< "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl
	<< "<PolyData>" << std::endl;
	
	for(auto& streamline : streamlines)
	{
		os
		<< "<Piece NumberOfPoints=\""<< streamline.size() <<"\" "
		<< "NumberOfLines=\""<< streamline.size()-1<<"\""
		<<">"
		<<std::endl;

		os
		<< "<Points>" <<std::endl
		<< "<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" <<std::endl;

		for(auto& vertex : streamline)
		{
			os
			<< vertex.x <<" "<< vertex.y <<" "<< vertex.z <<std::endl;
		}

		os
		<< "</DataArray>" <<std::endl
		<< "</Points>" <<std::endl;

		os
		<< "<Lines>" <<std::endl
		<< "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" <<std::endl;
		for(long i=0; i<(long)streamline.size(); i++)
			os<< i <<" "<< i+1 <<std::endl;
		os
		<< "</DataArray>" <<std::endl
		<< "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" <<std::endl;
		for(long i=2; i<=(long)(streamline.size()*2); i=i+2)
			os<< i <<std::endl;
		os
		<< "</DataArray>" <<std::endl
		<< "</Lines>" <<std::endl;

		os
		<< "</Piece>" <<std::endl;
	}

	os
	<< "</PolyData>" << std::endl
	<< "</VTKFile>" << std::endl;
	fb.close ();
}

void VTKOutput::interpolateUV(Point p, Real& u, Real& v)
{
	Real dx = this->domain.getDelta().x;
	Real dy = this->domain.getDelta().y;

	// step1: find grid cell
	int i = (int)(p.x/dx)+1;
	int j = (int)((p.y+(dy/2.0))/dy)+1;

	// step2: find neighbour cells
	Real x_1 = (Real)(i-1) * dx;
	Real x_2 = (Real)i * dx;
	Real y_1 = (Real)(j-1)*dy-(dy/2.0);
	Real y_2 = (Real)(j)*dy-(dy/2.0);

	// step3: bilinear interpolation

	u = 1.0/(dx*dy) * ( (x_2 - p.x)*(y_2 - p.y) * domain.u()(i-1,j-1) +
						(p.x - x_1)*(y_2 - p.y) * domain.u()(i,j-1) +
						(x_2 - p.x)*(p.y - y_1) * domain.u()(i-1,j) +
						(p.x - x_1)*(p.y - y_1) * domain.u()(i,j) );

	v = 1.0/(dx*dy) * ( (x_2 - p.x)*(y_2 - p.y) * domain.v()(i-1,j-1) +
						(p.x - x_1)*(y_2 - p.y) * domain.v()(i,j-1) +
						(x_2 - p.x)*(p.y - y_1) * domain.v()(i-1,j) +
						(p.x - x_1)*(p.y - y_1) * domain.v()(i,j) );
}
