//! The GridFunction ("Grid") is our bitch of choice in handling raw-memory access
/*!
 * @author becherml, friesfn, geringsj
 * @date 2014
 */

#ifndef GridFunction_hpp
#define GridFunction_hpp

#include "Structs.hpp"
#include "Debug.hpp"


/* the gridfunction just does what it is told to do: maintain a grid 
 * of size X * Y * Z - for now we set Z==1
 * basically all this class needs to do is allocate the memory
 * and give us indexed access to it in 2 or 3 dimensions
 * be it using operator()(int i, int j, int k) or some sort of
 * operator[](Index index) */
/* TODO: when i think of it, almost all of the following {set,get}GridFunction()
 * methods are pretty useless, we really just care for indexed memory access,
 * getting the MAX value of the grid and maybe the dimensions
 */

class GridFunction
{
private:
	Dimension dimension;
	Real* grid;

	void init(const uint dimX, const uint dimY, const uint dimZ);

public:
	GridFunction(const uint dimX, const uint dimY, const uint dimZ = 1);
	GridFunction(const Index griddimension);
	virtual ~GridFunction();

	inline Real& operator[](Dimension d){
		return this->operator()(d.i, d.j, d.k);
	}
	inline Real operator[](Dimension d) const {
		return this->operator()(d.i, d.j, d.k);
	}
	inline Real& operator[](int i){
		return this->grid[i];
	}
	inline Real operator[](int i) const {
		return this->grid[i];
	}
	inline Real& operator()(int i){
		return this->grid[i];
	}
	inline Real operator()(int i) const {
		return this->grid[i];
	}
	inline Real& operator()(int i, int j){
		return this->grid[i * dimension.j + j];
	}
	inline Real operator()(int i, int j) const { 
		return this->grid[i * dimension.j + j];
	}
	inline Real& operator()(int i, int j, int k){ /* TODO !? */
		return this->grid[i * dimension.j * dimension.k + j * dimension.k + k];
	}
	inline Real operator()(int i, int j, int k) const {
		return this->grid[i * dimension.j * dimension.k + j * dimension.k + k];
	}

	Index getGridDimension() const;
	void setGridFunction(
		const Index begin, const Index end, 
		const Real value);
	void scaleGridFunction(
		const Index begin, const Index end, 
		const Real scale);
	void divideGridFunction(
		const Index begin, const Index end, 
		const Real factor);
	void setGridFunction(
		const Index begin, const Index end, 
		const Real factor, 
		const GridFunction sourcegridfunction);
	void setGridFunction(
		const Index begin, const Index end, 
		const Real factor, 
		const GridFunction sourcegridfunction, 
		const Index offset);
	void setGridFunction(
		const Index begin, const Index end, 
		const Real factor, 
		const GridFunction sourcegridfunction, 
		const Index offset, 
		const Real constant);
	void addToGridFunction (
		const Index begin, const Index end, 
		const Real factor, 
		const GridFunction sourcegridfunction);
	Real getMaxValueGridFunction(
		const Index begin, const Index end);
	Real getMaxValueGridFunction();
};

#endif
