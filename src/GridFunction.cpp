#include "GridFunction.hpp"

GridFunction::GridFunction(Index dims) 
{
	dimension[0] = dims.i;
	dimension[1] = dims.j;
	dimension[2] = 1;

	this->wholeRange = Range(Index(0,0), Index(dims.i-1, dims.j-1));

	this->grid = new Real[dimension[0]*dimension[1]*dimension[2]];

	for(int i=0; i<dimension[0]*dimension[1]*dimension[2]; i++)
		this->grid[i] = 0.0;
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

Real GridFunction::getMaxValue()
{
	Real max = this->operator()(0,0);

	for_range(i,j,this->wholeRange)
	{
		if(max < this->operator()(i,j))
			max = this->operator()(i,j);
	}
	return max;
}

//Real& GridFunction::operator[](const Index& d){
//	return this->operator()(d.i, d.j);
//}
//Real GridFunction::operator[](const Index& d) const {
//	return this->operator()(d.i, d.j);
//}
Real& GridFunction::operator()(const Index& d){
	return this->operator()(d.i, d.j);
}
Real GridFunction::operator()(const Index& d) const {
	return this->operator()(d.i, d.j);
}
Real& GridFunction::operator()(const int& i, const int& j){
	return this->grid[j * dimension.i + i];
	//return this->grid[i * dimension.j + j];
}
Real GridFunction::operator()(const int& i, const int& j) const { 
	return this->grid[j * dimension.i + i];
	//return this->grid[i * dimension.j + j];
}

/* WARNING: here, the indices move as follows: first k, then j, i is slowest, 
 * so beware when extending to 3D */
//Real& operator()(int i, int j, int k){ /* TODO !? */
//	return this->grid[i * dimension.j * dimension.k + j * dimension.k + k];
//}
//Real operator()(int i, int j, int k) const {
//	return this->grid[i * dimension.j * dimension.k + j * dimension.k + k];
//}

#include <iostream>
void GridFunction::printSTDOUT()
{
	Index SI(0,0), EI(this->dimension.i-1,this->dimension.j-1);
	for(int J=EI[1]; J>=SI[1]; J--)
	{
		for(int I=SI[0]; I<=EI[0]; I++)
		{
			std::cout <<std::setprecision(5)<<std::fixed << this->operator()(I,J) << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}
