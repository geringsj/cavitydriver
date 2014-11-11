#ifndef Domain_h
#define Domain_h

#include "Structs.hpp"
#include "GridFunction.hpp"

#include <functional>

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
 * Therefore it stores and manages all grids and their properties.
 */
class Domain
{
public:
	Domain(Dimension dimension, Point delta,
		std::function<Real(Index,GridFunction&,Dimension)> in_u,
		std::function<Real(Index,GridFunction&,Dimension)> in_v,
		std::function<Real(Index,GridFunction&,Dimension)> in_w,
		std::function<Real(Point)> in_gx = [](Point coord)->Real{ return coord.x*0.0; },
		std::function<Real(Point)> in_gy = [](Point coord)->Real{ return coord.x*0.0; },
		std::function<Real(Point)> in_gz = [](Point coord)->Real{ return coord.x*0.0; }
		);

	~Domain();

	Dimension getDimension() { return m_dimension; }
	Dimension getBeginInnerDomains() {return m_inner_begin; }
	Dimension* getEndInnerDomain() {return m_inner_end; }
	Dimension getEndInnerDomainU() {return m_inner_end[0]; }
	Dimension getEndInnerDomainV() {return m_inner_end[1]; }
	Dimension getEndInnerDomainW() {return m_inner_end[2]; }
	Dimension getEndInnerDomainP() {return m_inner_end[3]; }
	Point getDelta() { return m_delta; }

	Grid3D& getVeolcity() { return m_velocities; }
	Grid3D& getPreliminaryVeolcity() { return m_preliminary_velocities_FGH; }

	GridFunction& p() { return m_p; }
	GridFunction& u() { return m_velocities.m_u; }
	GridFunction& v() { return m_velocities.m_v; }
	GridFunction& w() { return m_velocities.m_w; }
	GridFunction& F() { return m_preliminary_velocities_FGH.m_u ; }
	GridFunction& G() { return m_preliminary_velocities_FGH.m_v ; }
	GridFunction& H() { return m_preliminary_velocities_FGH.m_w ; }
	GridFunction& rhs() { return m_p_rhs; }

	Real gx(Point coord) { return m_infunc_gx(coord); }
	Real gy(Point coord) { return m_infunc_gy(coord); }
	Real gz(Point coord) { return m_infunc_gz(coord); }

	void setVelocitiesBoundaries()
	{
		for(uint d=0; d<DIMENSIONS; d++)
		{
			for(int i=m_inner_begin[0]-1; i<=m_inner_end[d][0]+1; i+=(m_inner_end[d][0]+1-(m_inner_begin[0]-1)))
				for(int j=m_inner_begin[1]-1; j<=m_inner_end[d][1]+1; j+=(m_inner_end[d][1]+1-(m_inner_begin[1]-1)))
					for(int k=m_inner_begin[2]-1; k<=m_inner_end[d][2]+1; k+=(m_inner_end[d][2]+1-(m_inner_begin[2]-1)))
			{
				Dimension current(i,j,k);
				if(d==0)
					u()(i,j,k) = m_infunc_u(current);
				if(d==1)
					v()(i,j,k) = m_infunc_v(current);
				if(d==3)
					w()(i,j,k) = m_infunc_w(current);
			}
		}
	}

private:
	/**
	 * Size of the domain
	 */
	Dimension m_dimension;
	/**
	 * Marks the inner of the domain 
	 */
	Dimension m_inner_begin, m_inner_end[4];

	/**
	 * Grid spacing
	 */
	Point m_delta;

	/**
	 * Pressure Grid
	 */
	GridFunction m_p;
	GridFunction m_p_rhs;

	/**
	 * Velocities grid
	 */
	Grid3D m_velocities;
	Grid3D m_preliminary_velocities_FGH;

	/**
	 * Input functions used to (re)set boundaries and starting conditions
	 */
	std::function<Real(Index)> m_infunc_u;
	std::function<Real(Index)> m_infunc_v;
	std::function<Real(Index)> m_infunc_w;

	/**
	 * External forces
	 */
	std::function<Real(Point)> m_infunc_gx;
	std::function<Real(Point)> m_infunc_gy;
	std::function<Real(Point)> m_infunc_gz;
};

#endif
