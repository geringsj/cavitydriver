
#ifndef Domain_hpp
#define Domain_hpp

#include "Structs.hpp"
#include "GridFunction.hpp"

#include <functional>


/** 
 * Implements the Domain we will be working on, managing dimensions, grids and properties.
 * 
 * @author becherml, friesfn, geringsj
 * @date 11/2014
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
public:
	struct Boundary {
		/* to be used as booleans ! */
		int Up;
		int Down;
		int Left;
		int Right;

		Boundary() : Up(1), Down(1), Left(1), Right(1) {};
		Boundary(int u, int d, int l, int r) 
			: Up(u), Down(d), Left(l), Right(r) {};

		bool operator[](int i)
		{
			switch(i)
			{
				case 0: /* U */
					return Right;
					break;
				case 1: /* V */
					return Up;
					break;
				case 2: /* W */
					return 0; /* TODO Front/Back */
					break;
				default:
					return 0;
					break;
			}
		}
	};

private: 
	/** 
	 * We need this struct for things to be easier.
	 *
	 * Because U, V (, W) and F, G (, H) have special needs for their dimensions, 
	 * this struct is used to encapsulate the initialisation and storage of those grids.
	 */
	struct Grid3D
	{
		Grid3D(Dimension dim, Boundary bndry=Boundary()) : 
			m_u(Dimension(dim.i+2-(int)bndry.Right,dim.j+2)),
			m_v(Dimension(dim.i+2,dim.j+2-(int)bndry.Up)),
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
	Domain(Dimension dimension, Delta delta,

		/* the functions for setting boundary and initial grid values */
		std::function<Real(Index,GridFunction&,Dimension)> in_u,
		std::function<Real(Index,GridFunction&,Dimension)> in_v,
		std::function<Real(Index,GridFunction&,Dimension)> in_w,
		std::function<Real(Index, GridFunction&, Dimension)> in_p,
		Real in_gx = 0.0,
		Real in_gy = 0.0,
		Real in_gz = 0.0,
		Domain::Boundary bndry = Boundary(),
		Color firstCellColor = Color::Red
		);

	~Domain();

	/* TODO: document for doxygen and move implementation into .cpp ?
	 * on the other hand: explaining multiple simple similiar functions is stupid, 
	 * better have the user see what is done? */
	Dimension getDimension() const { return m_dimension; }

	Range* getInnerRanges() { return m_inner_ranges; };
	Range getInnerRangeU() const { return m_inner_ranges[0]; };
	Range getInnerRangeV() const { return m_inner_ranges[1]; };
	Range getInnerRangeW() const { return m_inner_ranges[2]; };
	Range getInnerRangeP() const { return m_inner_ranges[3]; };

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

	Real g(int dim) const
	{ if(dim == 0) return gx(); if(dim == 1) return gy(); return gz(); }

	Real gx() const { return m_force_gx; }
	Real gy() const { return m_force_gy; }
	Real gz() const { return m_force_gz; }

	Color getDomainFirstCellColor(){ return m_FirstCellColor; };

	void setPressureBoundaries();
	void setPreliminaryVelocitiesBoundaries();
	void setVelocitiesBoundaries();

private:
	/**
	 * Size of the domain.
	 */
	Dimension m_dimension;
	/**
	 * Marks the inner of the domain.
	 */
	Dimension m_inner_begin, m_inner_end[4]; /* TODO: get rid of these */

	/**
	 * Gives the ranges of the inner of fields.
	 * Array entries correspond to:
	 * [0] -> inner Range of U
	 * [1] -> inner Range of V
	 * [2] -> inner Range of W (maye we will need it ...)
	 * [3] -> inner Range of pressure P
	 */
	Range m_inner_ranges[4];

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
	 * Input functions used to (re)set boundaries and starting conditions
	 */
	std::function<Real(Index)> m_borderfunc_u;
	std::function<Real(Index)> m_borderfunc_v;
	std::function<Real(Index)> m_borderfunc_w;

	Boundary m_boundary;

	/**
	 * External forces
	 */
	Real m_force_gx;
	Real m_force_gy;
	Real m_force_gz;
};

#endif
