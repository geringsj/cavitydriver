
#ifndef Domain_hpp
#define Domain_hpp

#include "Structs.hpp"
#include "GridFunction.hpp"
#include "Boundary.hpp"

#include <vector>

typedef std::vector<Range> Ranges;
#define for_vecrange(F,S,R) for(auto&M:R)for_range(F,S,M)


/** 
 * Implements the Domain we will be working on, managing dimensions, grids and properties.
 * 
 * @author becherml, friesfn, geringsj
 * @date 02/2015
 *
 * The Domain class represents the domain the simulation is computed on.
 * Therefore it stores and manages all grids and their properties.
 *
 * Sorry for not documenting every method in detail. 
 * We just have a bunch of getter-methods to get the grids and their properties.
 * Also the Domain knows how to (re)set the boundaries of all grids ("itself").
 */
class Domain
{
private: 
	/** 
	 * We need this struct for things to be easier.
	 *
	 * Because U, V (, W) and F, G (, H) have special needs for their dimensions, 
	 * this struct is used to encapsulate the initialisation and storage of those grids.
	 */
	struct Grid3D
	{
		Grid3D(Dimension dim, Boundary::Competence bndrycomp) : 
			m_u(Dimension(dim.i+2-(int)bndrycomp.Right,dim.j+2)),
			m_v(Dimension(dim.i+2,dim.j+2-(int)bndrycomp.Up)),
			m_w(Dimension(0,0)) {}
	
		GridFunction m_u;
		GridFunction m_v;
		GridFunction m_w;
	
		GridFunction& operator[](uint index)
		{
			switch (index)
			{
			case 0:
				return m_u;
				break;
			case 1:
				return m_v;
				break;
			case 2:
				return m_w;
				break;
			default:
				return m_u; /* yes, this is wrong. print warning? */
				break;
			}
		}
	};

public:
	Domain(
		Dimension dimension, Delta delta,

		Boundary boundary,

		/* field forces */
		Real in_gx = 0.0,
		Real in_gy = 0.0,
		//Real in_gz = 0.0,

		/* initial grid values */
		Real in_uinit = 0.0,
		Real in_vinit = 0.0,
		//Real in_winit = 0.0,
		Real in_pinit = 0.0,
		Real in_tinit = 0.0,

		/* color for SOR Red/Black pattern */
		Color firstCellColor = Color::Red
		);

	~Domain();

	/* TODO: document for doxygen and move implementation into .cpp ?
	 * on the other hand: explaining multiple simple similiar functions is stupid, 
	 * better have the user see what is done? */
	Dimension getDimension() const { return m_dimension; }

	const Range getWholeInnerRange() const { return m_whole_inner_range; }
	const Ranges* getInnerRanges() const { return m_inner_ranges; }
	const Ranges& getInnerRangeU() const { return m_inner_ranges[0]; }
	const Ranges& getInnerRangeV() const { return m_inner_ranges[1]; }
	//Rang&es getInnerRangeW() const { return m_inner_ranges[2]; };
	const Ranges& getInnerRangeP() const { return m_inner_ranges[3]; }
	const Ranges& getInnerRangeT() const { return m_inner_ranges[4]; }

	Point getDelta() const { return m_delta; }

	Grid3D& getVelocity() { return m_velocities; }
	Grid3D& getPreliminaryVelocity() { return m_preliminary_velocities_FGH; }

	GridFunction& p() { return m_p; }
	GridFunction& u() { return m_velocities.m_u; }
	GridFunction& v() { return m_velocities.m_v; }
	GridFunction& w() { return m_velocities.m_w; }
	GridFunction& F() { return m_preliminary_velocities_FGH.m_u ; }
	GridFunction& G() { return m_preliminary_velocities_FGH.m_v ; }
	GridFunction& H() { return m_preliminary_velocities_FGH.m_w ; }
	GridFunction& rhs() { return m_p_rhs; }
	GridFunction& t() { return m_temperature; }
	GridFunction& tcache() { return m_temperature_computationcache; }

	Real g(int dim) const
	{ if(dim == 0) return gx(); if(dim == 1) return gy(); return gy();/*TODO later do gz*/ }

	Real gx() const { return m_force_gx; }
	Real gy() const { return m_force_gy; }
	//Real gz() const { return m_force_gz; }

	Color getDomainFirstCellColor(){ return m_FirstCellColor; };

	int getFluidCellsCount() const
	{
		int sum = 0;
		for_vecrange(i,j,m_inner_ranges[3])
			sum += 1;
		return sum;
	}

	void setPressureBoundaries();
	void setPreliminaryVelocitiesBoundaries();
	void setVelocitiesBoundaries();
	void setTemperatureBoundaries();

private:
	/**
	 * Size of the domain.
	 */
	Dimension m_dimension;

	Boundary m_boundary;

	/**
	 * Gives the ranges of the inner of fields.
	 * Array entries correspond to:
	 * [0] -> inner Range of U
	 * [1] -> inner Range of V
	 * [2] -> inner Range of W (maye we will need it ...)
	 * [3] -> inner Range of pressure P
	 * [4] -> inner Range of temperature T
	 */
	Ranges m_inner_ranges[5];
	Range m_whole_inner_range;

	/**
	 * Color of first inner ((1,1)) Cell of domain.
	 */
	Color m_FirstCellColor;

	/**
	 * Grid spacing.
	 */
	Point m_delta;

	/**
	 * Pressure Grid and and its RHS that is needed for computation.
	 */
	GridFunction m_p;
	GridFunction m_p_rhs;

	/**
	 * Velocities grids.
	 */
	Grid3D m_velocities;
	Grid3D m_preliminary_velocities_FGH;

	/**
	 * External forces
	 */
	Real m_force_gx;
	Real m_force_gy;
	//Real m_force_gz;

	/**
	 * Temperature
	 */
	GridFunction m_temperature;
	GridFunction m_temperature_computationcache;
};

#endif
