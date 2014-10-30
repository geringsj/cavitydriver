#include "GridFunction.hpp"

GridFunction::GridFunction(const uint dimX, const uint dimY)
{
	griddimension[0] = dimX;
	griddimension[1] = dimY;
	/* TODO dim == 0 ? => warning */

	this->rawmemory = new RealType[dimX*dimY];
	this->gridfunction = new RealType*[dimY];
	for(uint i=0; i<dimY; i++)
		this->gridfunction[i] = rawmemory + i*dimX;
}

GridFunction::GridFunction(const MultiIndexType griddimension) 
	: GridFunction(griddimension[0], griddimension[1])
{
}

GridFunction::~GridFunction()
{
	for(int i=0; i<this->griddimension[1]; i++)
		this->gridfunction[i] = NULL;
	delete[] this->gridfunction;
	this->gridfunction = NULL;
	delete[] this->rawmemory;
	this->rawmemory = NULL;
}

GridFunctionType GridFunction::getGridFunction()
{
	return this->gridfunction;
}

MultiIndexType GridFunction::getGridDimension()
{
	return this->griddimension;
}

#define forall(F,S,B,E) for(int F=B[0];F<=E[0];F++)for(int S=B[1];S<=E[1];S++)

void GridFunction::setGridFunction(
	const MultiIndexType begin, 
	const MultiIndexType end, 
	const RealType value)
{
	forall(i,j,begin,end)
	{
		this->gridfunction[i][j] = value;
	}
}

void GridFunction::scaleGridFunction(
	const MultiIndexType begin, 
	const MultiIndexType end, 
	const RealType scale)
{
	forall(i,j,begin,end)
	{
		this->gridfunction[i][j] *= scale;
	}
}

void GridFunction::setGridFunction(
	const MultiIndexType begin, 
	const MultiIndexType end, 
	const RealType factor, 
	const GridFunctionType sourcegridfunction)
{
	forall(i,j,begin,end)
	{
		this->gridfunction[i][j] = factor*sourcegridfunction[i][j];
	}
}

void GridFunction::setGridFunction(
	const MultiIndexType begin, 
	const MultiIndexType end, 
	const RealType factor, 
	const GridFunctionType sourcegridfunction, 
	const MultiIndexType offset)
{
	forall(i,j,begin,end)
	{
		this->gridfunction[i][j] = 
			factor*sourcegridfunction[offset[0]+i][offset[1]+j];
	}
}

void GridFunction::setGridFunction(
	const MultiIndexType begin, 
	const MultiIndexType end, 
	const RealType factor, 
	const GridFunctionType sourcegridfunction, 
	const MultiIndexType offset, 
	const RealType constant)
{
	forall(i,j,begin,end)
	{
		this->gridfunction[i][j] = 
			factor*sourcegridfunction[offset[0]+i][offset[1]+j]+constant;
	}
}

void GridFunction::addToGridFunction(
	const MultiIndexType begin, 
	const MultiIndexType end, 
	const RealType factor, 
	const GridFunctionType sourcegridfunction)
{
	forall(i,j,begin,end)
	{
		this->gridfunction[i][j] += factor*sourcegridfunction[i][j];
	}
}

RealType GridFunction::getMaxValueGridFunction(
	const MultiIndexType begin, 
	const MultiIndexType end)
{
	RealType max = this->gridfunction[0][0];

	forall(i,j,begin,end)
	{
		if(max < this->gridfunction[i][j])
			max = this->gridfunction[i][j];
	}
	return max;
}

RealType GridFunction::getMaxValueGridFunction()
{
	MultiIndexType begin; begin[0]=0; begin[1]=0;
	MultiIndexType end; 
	end[0]=this->griddimension[0]; end[1]=griddimension[1];
	/* TODO off by one? */
	return getMaxValueGridFunction(begin,end);
}

