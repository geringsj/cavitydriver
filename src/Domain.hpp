//! This file implements the domain we will be working in 
/*!
 * @author becherml, friesfn, geringsj
 * @date 2014
 */

#ifndef Domain_hpp
#define Domain_hpp

#include "Structs.hpp"
#include "GridFunction.hpp"

#include <functional>

/* TODO: move INTO Domain class or mvoe to Structs.hpp ? 
 * nobody uses 3D grids, except for the Domain */
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
		case 0:
			return m_u;
		case 1:
			return m_v;
		case 2:
			return m_w;
		default:
			return m_u; /* yes, this is wrong. print warning? */
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
		/* the functions for setting boundary and initial grid values */
		std::function<Real(Index,GridFunction&,Dimension)> in_u,
		std::function<Real(Index,GridFunction&,Dimension)> in_v,
		std::function<Real(Index,GridFunction&,Dimension)> in_w,
		std::function<Real(Index, GridFunction&, Dimension)> in_p,
		/* if no outer forces are given, we assume they are zero. 
		 * TODO: can this be done in .cpp with the preset values still beeing used? 
		 * or maybe split into two constructors ? 
		 * => constructor delegation would be nice to use*/
		std::function<Real(Point)> in_gx = 
			[](Point coord)->Real{ return coord.x*0.0; }, 
		std::function<Real(Point)> in_gy = 
			[](Point coord)->Real{ return coord.x*0.0; },
		std::function<Real(Point)> in_gz = 
			[](Point coord)->Real{ return coord.x*0.0; }
		);

	~Domain();

	/* TODO: document for doxygen and move implementation into .cpp ?
	 * on the other hand: explaining multiple simple similiar functions is stupid, 
	 * better have the user see what is done? */
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

	Real g(int dim, Point coord) 
	{ 
		if(dim == 0) return gx(coord);
		if(dim == 1) return gy(coord);
		return gz(coord);
	}


	Real gx(Point coord) { return m_forcefunc_gx(coord); }
	Real gy(Point coord) { return m_forcefunc_gy(coord); }
	Real gz(Point coord) { return m_forcefunc_gz(coord); }

	/* TODO: definitely move THIS and boundary functions to .cpp */
#define LEFT(start,end)	for(int i = m_inner_begin[0]-1, j = start; j <= end; j++)
#define RIGHT(start,end,d) for(int i = m_inner_end[d][0] + 1, j = start; j <= end; j++)
#define TOP(start,end,d) for(int j = m_inner_end[d][1] + 1, i = start; i <= end; i++)
#define BOTTOM(start,end) for(int j = m_inner_begin[1]-1, i = start; i <= end; i++)

#define VELOCITIESBOUNDARIES(i,j,current,d) \
	current = Dimension(i, j);\
	if (d == 0)\
		u()(i, j)/*,k)*/ = m_borderfunc_u(current);\
	if (d == 1)\
		v()(i, j)/*,k)*/ = m_borderfunc_v(current);\
	if (d == 3)\
		w()(i, j)/*,k)*/ = m_borderfunc_w(current);

#define PRELIMINARYVELOCITIESBOUNDARIES(i,j,d) \
	if (d == 0)\
		F()(i, j) = u()(i, j);\
	if (d == 1)\
		G()(i, j) = v()(i, j);\
	if (d == 3)\
		H()(i, j) = w()(i, j);

#define PRESSUREBOUNDARIES(i,j) \
	if (i == 0) p()(i, j) = p()(i + 1, j);\
	if (i == m_dimension[0] + 1) p()(i, j) = p()(i - 1, j);\
	if (j == m_dimension[1] + 1) p()(i, j) = p()(i, j - 1);\
	if (j == 0) p()(i, j) = p()(i, j + 1);


	void setVelocitiesBoundaries()
	{
		Dimension current;
		for(uint d=0; d<DIMENSIONS; d++)
		{
			LEFT(m_inner_begin[1] - 1, m_inner_end[d][1] + 1)
			{
				VELOCITIESBOUNDARIES(i, j, current, d);
			}
			RIGHT(m_inner_begin[1] - 1, m_inner_end[d][1] + 1, d) 
			{
				VELOCITIESBOUNDARIES(i, j, current, d);
			}
			TOP(m_inner_begin[0] - 1, m_inner_end[d][0] + 1, d)
			{
				VELOCITIESBOUNDARIES(i, j, current, d);
			}
			BOTTOM(m_inner_begin[0] - 1, m_inner_end[d][0] + 1)
			{
				VELOCITIESBOUNDARIES(i, j, current, d);
			}
		}
	}

	void setPreliminaryVelocitiesBoundaries()
	{
		for (uint d = 0; d<DIMENSIONS; d++)
		{
			LEFT(m_inner_begin[1] - 1, m_inner_end[d][1] + 1)
			{
				PRELIMINARYVELOCITIESBOUNDARIES(i, j, d);
			}
			RIGHT(m_inner_begin[1] - 1, m_inner_end[d][1] + 1, d)
			{
				PRELIMINARYVELOCITIESBOUNDARIES(i, j, d);
			}
			TOP(m_inner_begin[0] - 1, m_inner_end[d][0] + 1, d)
			{
				PRELIMINARYVELOCITIESBOUNDARIES(i, j, d);
			}
			BOTTOM(m_inner_begin[0] - 1, m_inner_end[d][0] + 1)
			{
				PRELIMINARYVELOCITIESBOUNDARIES(i, j, d);
			}
		}
	}

	void setPressureBoundaries()
	{
		LEFT(m_inner_begin[1] - 1, m_dimension[1] + 1)
		{
			PRESSUREBOUNDARIES(i, j);
		}
		RIGHT(m_inner_begin[1] - 1, m_dimension[1] + 1, 3)
		{
			PRESSUREBOUNDARIES(i, j);
		}
		TOP(m_inner_begin[0] - 1, m_dimension[0] + 1, 3)
		{
			PRESSUREBOUNDARIES(i, j);
		}
		BOTTOM(m_inner_begin[0] - 1, m_dimension[0] + 1)
		{
			PRESSUREBOUNDARIES(i, j);
		}
	}

private:
	/* TODO document in a non-stupid manner. 
	 * also mention: we are prepared to do 3D, but don't do it yet */
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
	 * Pressure Grid and and its RHS that is needed for computation.
	 */
	GridFunction m_p;
	GridFunction m_p_rhs;

	/**
	 * Velocities grids
	 */
	Grid3D m_velocities;
	Grid3D m_preliminary_velocities_FGH;

	/**
	 * Input functions used to (re)set boundaries and starting conditions
	 */
	std::function<Real(Index)> m_borderfunc_u;
	std::function<Real(Index)> m_borderfunc_v;
	std::function<Real(Index)> m_borderfunc_w;

	/**
	 * External forces
	 */
	std::function<Real(Point)> m_forcefunc_gx;
	std::function<Real(Point)> m_forcefunc_gy;
	std::function<Real(Point)> m_forcefunc_gz;
};

#endif
