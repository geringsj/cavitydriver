#ifndef Domain_h
#define Domain_h

#include "Structs.hpp"
#include "GridFunction.hpp"


struct Grid3D
{
	Grid3D(Dimension dim) : 
		m_u(dim.i+1,dim.j+2),
		m_v(dim.i+2,dim.j+1),
		m_w(0,0) {}
	//Grid3D(Dimension dimU, Dimension dimV, Dimension dimW) : 
	//	m_u(dimU.i+1 ,dimU.j+2 ,dimU.k+2),
	//	m_v(dimV.i+2 ,dimV.j+1 ,dimV.k+2),
	//	m_w(dimW.i+2 ,dimW.j+2 ,dimW.k+1) {}

	GridFunction m_u;
	GridFunction m_v;
	GridFunction m_w;

	GridFunction& operator[](uint index)
	{
		switch (index)
		{
		case 1:
			return m_u;
		case 2:
			return m_v;
		case 3:
			return m_w;
		default:
			return m_u;
			//return GridFunction(); 
			//exit(-1); //not to sure about this
		}
	}
};

/**
 * The class represents the complete simulation domain.
 * Therefore it stores and manages all grid and their properties.
 */
class Domain
{
public:
	Domain(Dimension dimension, Point delta);
	Domain(Dimension dimension, Point delta, Point extForces);
	~Domain();

	Dimension getDimension() { return m_dimension; }
	Point getDelta() { return m_delta; }

	Grid3D& getVeolcity() { return m_velocity; }

	GridFunction& p() { return m_p; }
	GridFunction& u() { return m_velocity.m_u; }
	GridFunction& v() { return m_velocity.m_v; }
	GridFunction& w() { return m_velocity.m_w; }
	GridFunction& F() { return m_preliminary_velocities_FGH.m_u ; }
	GridFunction& G() { return m_preliminary_velocities_FGH.m_v ; }
	GridFunction& H() { return m_preliminary_velocities_FGH.m_w ; }
	GridFunction& rhs() { return m_rhs; }
	Real gx() { return m_gx; }
	Real gy() { return m_gy; }
	Real gz() { return m_gz; }

	void setBoundary(); 
	/* TODO use suer-given boundary functions for each side ? 
	 * => dynamic class member function for boundary */

private:
	/**
	 * Size of the domain
	 */
	Dimension m_dimension;
	/**
	 * Marks the inner of the domain 
	 */
	Dimension m_inner_begin, m_inner_end;

	/**
	 * Grid spacing
	 */
	Point m_delta;

	/**
	 * Pressure Grid
	 */
	GridFunction m_p;
	/**
	 * Velocity grid
	 */
	Grid3D m_velocity;
	/**
	 * Other grids
	 */
	Grid3D m_preliminary_velocities_FGH;
	GridFunction m_rhs;
	/**
	 * External forces 
	 */
	Real m_gx, m_gy, m_gz; /* TODO are these supposed to be GridFunctions? */
};

#endif
