
#ifndef GridFunction_hpp
#define GridFunction_hpp

#include "typedef.hpp"

class GridFunction
{
private:
	GridFunctionType gridfunction;
	MultiIndexType griddimension;
	RealType* rawmemory;

public:
	GridFunction(const uint dimX, const uint dimY);
	GridFunction(const MultiIndexType griddimension);
	virtual ~GridFunction();

	GridFunctionType getGridFunction();	
	MultiIndexType getGridDimension();
	void setGridFunction(
		const MultiIndexType begin, 
		const MultiIndexType end, 
		const RealType value);
	void scaleGridFunction(
		const MultiIndexType begin, 
		const MultiIndexType end, 
		const RealType scale);
	void setGridFunction(
		const MultiIndexType begin, 
		const MultiIndexType end, 
		const RealType factor, 
		const GridFunctionType sourcegridfunction);
	void setGridFunction(
		const MultiIndexType begin, 
		const MultiIndexType end, 
		const RealType factor, 
		const GridFunctionType sourcegridfunction, 
		const MultiIndexType offset);
	void setGridFunction(
		const MultiIndexType begin, 
		const MultiIndexType end, 
		const RealType factor, 
		const GridFunctionType sourcegridfunction, 
		const MultiIndexType offset, 
		const RealType constant);
	void addToGridFunction (
		const MultiIndexType begin, 
		const MultiIndexType end, 
		const RealType factor, 
		const GridFunctionType sourcegridfunction);
	RealType getMaxValueGridFunction(
		const MultiIndexType begin, 
		const MultiIndexType end);
	RealType getMaxValueGridFunction();
};

#endif
