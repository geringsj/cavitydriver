#ifndef Domain_h
#define Domain_h

#include "Structs.hpp"
#include "typedef.hpp"
#include "GridFunction.hpp"

/**
 * The class represents the complete simulation domain.
 * Therefore it stores and manages all grid and their properties.
 */
class Domain
{
public:
	Domain(Dimension3D dimension, Point3D delta);
	~Domain();

	Dimension3D getDimension() { return m_dimension; }
	Point3D getDelta() { return m_delta; }

	Grid3D& getVeolcity() { return m_velocity; }

	GridFunction& p() { return m_p; }
	GridFunction& u() { return m_velocity.m_u; }
	GridFunction& v() { return m_velocity.m_v; }
	GridFunction& w() { return m_velocity.m_w; }
	GridFunction& F() { return m_F; }
	GridFunction& G() { return m_G; }
	GridFunction& H() { return m_H; }
	GridFunction& rhs() { return m_rhs; }
	GridFunction& gx() { return m_gx; }
	GridFunction& gy() { return m_gy; }
	GridFunction& gz() { return m_gz; }

private:
	/**
	 * Size of the domain
	 */
	Dimension3D m_dimension;
	/**
	 * Grid spacing
	 */
	Point3D m_delta;

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
	GridFunction m_F, m_G, m_H, m_rhs, m_gx, m_gy, m_gz;
};

#endif