//! The GridFunction ("Grid") is our bitch of choice in handling raw-memory access
/*!
 * @author becherml, friesfn, geringsj
 * @date 2014
 */

#ifndef GridFunction_hpp
#define GridFunction_hpp

#include "Structs.hpp"
#include <iomanip>

/* the gridfunction just does what it is told to do: maintain a grid 
 * of size X * Y  ( * Z - for now we set Z==1 )
 * basically all this class needs to do is allocate the memory
 * and give us indexed access to it in 2 or 3 dimensions
 * be it using operator()(int i, int j, int k) or some sort of
 * operator[](Index index) */

/**
 * Better use this makro to walk on the grid (including beginning B and ending E
 */
#define forall(F,S,B,E) for(int F=B[0];F<=E[0];F++)for(int S=B[1];S<=E[1];S++)

class GridFunction
{
private:
	Dimension dimension;
	Real* grid;

public:
	GridFunction(const Index griddimension);
	virtual ~GridFunction();

	Real& operator[](Index d);
	Real operator[](Index d) const;
	Real& operator()(int i, int j);
	Real operator()(int i, int j) const;

	Index getGridDimension() const;
	Real getMaxValueGridFunction();
	void printSTDOUT();
};

#endif
