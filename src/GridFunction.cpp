#include "GridFunction.hpp"

void GridFunction::init(const uint dimX, const uint dimY, const uint dimZ = 1)
{
	dimension[0] = dimX;
	dimension[1] = dimY;
	dimension[2] = dimZ;

	this->grid = new Real[dimX*dimY*dimZ];
}

GridFunction::GridFunction(const uint dimX, const uint dimY, const uint dimZ)
{
	this->init(dimX, dimY, dimZ);
}

GridFunction::GridFunction()
{
}

GridFunction::GridFunction(Index griddimension) 
{
	this->init(griddimension[0], griddimension[1], griddimension[2]);
}

GridFunction::~GridFunction()
{
	delete[] this->grid;
	this->grid = NULL;
}

Index GridFunction::getGridDimension() const
{
	return this->dimension;
}

#define forall(F,S,B,E) for(int F=B[0];F<=E[0];F++)for(int S=B[1];S<=E[1];S++)

void GridFunction::setGridFunction(
	const Index begin, 
	const Index end, 
	const Real value)
{
	forall(i,j,begin,end)
	{
		this->operator()(i,j) = value;
	}
}

void GridFunction::scaleGridFunction(
	const Index begin, 
	const Index end, 
	const Real scale)
{
	forall(i,j,begin,end)
	{
		this->operator()(i,j) *= scale;
	}
}

void GridFunction::setGridFunction(
	const Index begin, 
	const Index end, 
	const Real factor, 
	const GridFunction sourcegridfunction)
{
	forall(i,j,begin,end)
	{
		this->operator()(i,j) = factor*sourcegridfunction(i,j);
	}
}

void GridFunction::setGridFunction(
	Index begin, 
	Index end, 
	const Real factor, 
	const GridFunction sourcegridfunction, 
	Index offset)
{
	forall(i,j,begin,end)
	{
		this->operator()(i,j) = 
			factor*sourcegridfunction(offset[0]+i,offset[1]+j);
	}
}

void GridFunction::setGridFunction(
	const Index begin, 
	const Index end, 
	const Real factor, 
	const GridFunction sourcegridfunction, 
	const Index offset, 
	const Real constant)
{
	forall(i,j,begin,end)
	{
		this->operator()(i,j) = 
			factor*sourcegridfunction(offset[0]+i,offset[1]+j)+constant;
	}
}

void GridFunction::addToGridFunction(
	const Index begin, 
	const Index end, 
	const Real factor, 
	const GridFunction sourcegridfunction)
{
	forall(i,j,begin,end)
	{
		this->operator()(i,j) += factor*sourcegridfunction(i,j);
	}
}

Real GridFunction::getMaxValueGridFunction(
	const Index begin, 
	const Index end)
{
	Real max = this->operator()(0,0);

	forall(i,j,begin,end)
	{
		if(max < this->operator()(i,j))
			max = this->operator()(i,j);
	}
	return max;
}

Real GridFunction::getMaxValueGridFunction()
{
	Index begin; begin[0]=1; begin[1]=1;
	Index end; 
	end[0]=this->dimension[0]-1; end[1]=this->dimension[1]-1;
	return getMaxValueGridFunction(begin,end);
}

