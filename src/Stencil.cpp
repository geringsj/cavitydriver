#include "Stencil.hpp"

Stencil::Stencil(int stencilwidth, const PointType h)
	: stencilwidth(stencilwidth) , h(h) 
{
}

Stencil::~Stencil()
{

}

void Stencil::ApplyStencilOperator(const MultiIndexType& gridbegin, 
								   const MultiIndexType& gridreadend,
								   const MultiIndexType& gridwritebegin,
								   const MultiIndexType& gridwriteend, 
								   const GridFunction sourcegridfunction,
								   GridFunction imagegridfunction)
{

}

void Stencil::setFxxStencil()
{
	stencil[-1][0] = 1.0 / (h[0] * h[0]);
	stencil[1][0] = 1.0 / (h[0] * h[0]);
	stencil[0][-1] = 0.0;
	stencil[0][1] = 0.0;
	stencil[0][0] = -2.0 / (h[0] * h[0]);
}

void Stencil::setFxyStencil()
{

}

void Stencil::setFyxStencil()
{

}

void Stencil::setFyyStencil()
{

}