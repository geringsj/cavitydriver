
#ifndef GridFunction_hpp
#define GridFunction_hpp

#include "Structs.hpp"

class GridFunction
{
private:
	Dimension dimension;
	Real* grid;

	void init(const uint dimX, const uint dimY, const uint dimZ);

public:
	GridFunction();
	GridFunction(const uint dimX, const uint dimY, const uint dimZ = 1);
	GridFunction(const Index griddimension);
	virtual ~GridFunction();

	inline Real& operator[](int i){ /* TODO multi dimension? */
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
