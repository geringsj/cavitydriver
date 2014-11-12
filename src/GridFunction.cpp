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

Real GridFunction::getMaxValueGridFunction(
	const Index begin, 
	const Index end)
{
	Real max = this->operator()(begin.i, begin.j);

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

