
#ifndef GridFunction_hpp
#define GridFunction_hpp


#include "Structs.hpp"
#include <iomanip>


/**
 * Better use this makro to walk on the grid.
 */
#define forall(F,S,B,E) for(int S=B[1];S<=E[1];S++)for(int F=B[0];F<=E[0];F++)

/* F beeing i (in x direction) and S beeing j (in y direction) */
#define for_range(F,S,R) for(int S=R.begin.j;S<=R.end.j;S++)for(int F=R.begin.i;F<=R.end.i;F++)


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
	Range wholeRange; /**< Range in which access to grid points is defined. All bounds are inclusive, like in for(int i=begin; i <= end; ...)  */
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
	Real& operator()(const Index& d);
	//Real& operator[](const Index& ind);
	/** 
	 * Provides read access to the grid via the [] operator.
	 * @return Value of grid at Index position ind, per value;
	 */
	Real operator()(const Index& d) const;
	//Real operator[](const Index& ind) const;

	/** 
	 * Provides 2D write access to the grid via the () operator.
	 * @return Value at position (i,j), per reference.
	 */
	Real& operator()(const int& i, const int& j);
	/** 
	 * Provides 2D read access to the grid via the () operator.
	 * @return Value at position (i,j), per value.
	 */
	Real operator()(const int& i, const int& j) const;

	/** 
	 * @return Dimensions the grid was initialized with.
	 */
	Dimension getGridDimension() const;

	/** 
	 * Determines the maximum of all absolute values of grid entries.
	 * @return Maximum of absolute grid values.
	 */
	Real getMaxAbsValue();

	/** Prints 2D grid to stdout in a nicely formatted manner, for debugging purposes. */
	void printSTDOUT();
};

#endif
