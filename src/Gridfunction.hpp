#include "typedef.hpp"

class Gridfunction
{
private:
	GridFunctionType gridfunction;
	MultiIndexType griddimension;
	RealType* rawmemory;

public:
	Gridfunction(const uint dimX, const uint dimY);
	Gridfunction(const MultiIndexType griddimension);
	virtual ~Gridfunction();

	GridFunctionType getGridfunction();	
	MultiIndexType getGridDimension();
	void setGridfunction(
		const MultiIndexType begin, 
		const MultiIndexType end, 
		const RealType value);
	void scaleGridfunction(
		const MultiIndexType begin, 
		const MultiIndexType end, 
		const RealType scale);
	void setGridfunction(
		const MultiIndexType begin, 
		const MultiIndexType end, 
		const RealType factor, 
		const GridFunctionType sourcegridfunction);
	void setGridfunction(
		const MultiIndexType begin, 
		const MultiIndexType end, 
		const RealType factor, 
		const GridFunctionType sourcegridfunction, 
		const MultiIndexType offset);
	void setGridfunction(
		const MultiIndexType begin, 
		const MultiIndexType end, 
		const RealType factor, 
		const GridFunctionType sourcegridfunction, 
		const MultiIndexType offset, 
		const RealType constant);
	void addToGridfunction (
		const MultiIndexType begin, 
		const MultiIndexType end, 
		const RealType factor, 
		const GridFunctionType sourcegridfunction);
	RealType getMaxValueGridfunction(
		const MultiIndexType begin, 
		const MultiIndexType end);
	RealType getMaxValueGridfunction();
};
