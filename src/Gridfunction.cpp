#include "Gridfunction.hpp"

Gridfunction::Gridfunction(const uint dimX, const uint dimY)
{
	griddimension[0] = dimX;
	griddimension[1] = dimY;
	/* TODO dim == 0 ? => warning */

	this->rawmemory = new RealType[dimX*dimY];
	this->gridfunction = new RealType*[dimY];
	for(uint i=0; i<dimY; i++)
		this->gridfunction[i] = rawmemory + i*dimX;
}

Gridfunction::Gridfunction(const MultiIndexType griddimension) 
	: Gridfunction(griddimension[0], griddimension[1])
{
}

Gridfunction::~Gridfunction()
{
	for(int i=0; i<this->griddimension[1]; i++)
		this->gridfunction[i] = NULL;
	delete[] this->gridfunction;
	this->gridfunction = NULL;
	delete[] this->rawmemory;
	this->rawmemory = NULL;
}

GridFunctionType Gridfunction::getGridfunction()
{
	return this->gridfunction;
}

MultiIndexType Gridfunction::getGridDimension()
{
	return this->griddimension;
}

#define forall(F,S,B,E) for(int F=B[0];F<=E[0];F++)for(int S=B[1];S<=E[1];S++)

void Gridfunction::setGridfunction(
	const MultiIndexType begin, 
	const MultiIndexType end, 
	const RealType value)
{
	forall(i,j,begin,end)
	{
		this->gridfunction[i][j] = value;
	}
}

void Gridfunction::scaleGridfunction(
	const MultiIndexType begin, 
	const MultiIndexType end, 
	const RealType scale)
{
	forall(i,j,begin,end)
	{
		this->gridfunction[i][j] *= scale;
	}
}

void Gridfunction::setGridfunction(
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

void Gridfunction::setGridfunction(
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

void Gridfunction::setGridfunction(
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

void Gridfunction::addToGridfunction(
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

RealType Gridfunction::getMaxValueGridfunction(
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

RealType Gridfunction::getMaxValueGridfunction()
{
	MultiIndexType begin; begin[0]=0; begin[1]=0;
	MultiIndexType end; 
	end[0]=this->griddimension[0]; end[1]=griddimension[1];
	/* TODO off by one? */
	return getMaxValueGridfunction(begin,end);
}

