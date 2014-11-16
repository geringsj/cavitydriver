
#ifndef GridFunction_hpp
#define GridFunction_hpp


#include "Structs.hpp"
#include <iomanip>


/**
 * Better use this makro to walk on the grid.
 */
#define forall(F,S,B,E) for(int F=B[0];F<=E[0];F++)for(int S=B[1];S<=E[1];S++)


/** 
 * Provides indexed access to a 2D (3D) grid structure of Real values.
 * 
 * @author becherml, friesfn, geringsj
 * @date 11/2014
 *
 * The GridFunction maintains a grid of dimensions X*Y(*Z) (for now, we set Z=1).
 * The first grid element is at Index position (0,0,0).
 * One shall NOT pass around GridFunctions per value, 
 * because no sane copy constructors are implemented.
 * Again, do NOT move or pass GridFunctions per value. 
 * Just don't.
 */
class GridFunction
{
private:
	Dimension dimension; /**< Dimensions of the grid. */
	Real* grid; /**< Raw memory of the grid as array of Real values. */

public:
	/**
	 * Initialize Grid with certain dimensions.
	 */
	GridFunction(
			const Dimension griddimension
			/**< Dimensions of the grid. */
			);
	virtual ~GridFunction();

	/** 
	 * Provides write access to the grid via the [] operator.
	 * @return Value at Index position ind, per reference.
	 */
	Real& operator[](Index ind);
	/** 
	 * Provides read access to the grid via the [] operator.
	 * @return Value of grid at Index position ind, per value;
	 */
	Real operator[](Index ind) const;

	/** 
	 * Provides 2D write access to the grid via the () operator.
	 * @return Value at position (i,j), per reference.
	 */
	Real& operator()(int i, int j);
	/** 
	 * Provides 2D read access to the grid via the () operator.
	 * @return Value at position (i,j), per value.
	 */
	Real operator()(int i, int j) const;

	/** 
	 * @return Dimensions the grid was initialized with.
	 */
	Dimension getGridDimension() const;

	/** 
	 * Determines the maximum value out of all grid entries.
	 * @return Maximum grid value.
	 */
	Real getMaxValueGridFunction();

	/** Prints 2D grid to stdout in a nicely formatted manner, for debugging purposes. */
	void printSTDOUT();
};

#endif
