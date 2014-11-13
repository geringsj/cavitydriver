#include "GridFunction.hpp"

GridFunction::GridFunction(Index dims) 
{
	dimension[0] = dims.i;
	dimension[1] = dims.j;
	dimension[2] = 1;

	this->grid = new Real[dimension[0]*dimension[1]*dimension[2]];
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

Real GridFunction::getMaxValueGridFunction()
{
	Index begin(0,0);
	Real max = this->operator()(begin.i, begin.j);

	forall(i,j,begin,this->dimension)
	{
		if(max < this->operator()(i,j))
			max = this->operator()(i,j);
	}
	return max;
}

Real& GridFunction::operator[](Index d){
	return this->operator()(d.i, d.j);
}
Real GridFunction::operator[](Index d) const {
	return this->operator()(d.i, d.j);
}
Real& GridFunction::operator()(int i, int j){
	return this->grid[i * dimension.j + j];
}
Real GridFunction::operator()(int i, int j) const { 
	return this->grid[i * dimension.j + j];
}

//Real& operator()(int i, int j, int k){ /* TODO !? */
//	return this->grid[i * dimension.j * dimension.k + j * dimension.k + k];
//}
//Real operator()(int i, int j, int k) const {
//	return this->grid[i * dimension.j * dimension.k + j * dimension.k + k];
//}
