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

void Stencil::setFxxStencil(const MultiIndexType& begin, 
							const MultiIndexType& end, RealType value)
{

}

void Stencil::setFxxStencil(const MultiIndexType& begin, 
							const MultiIndexType& end, RealType factor,
							MultiIndexType& offset)
{

}

//void Stencil::setFxxStencil(const MultiIndexType& begin, 
//							const MultiIndexType& end, RealType factor)
//{
//
//}

void Stencil::setFxxStencil(const MultiIndexType& begin, 
							const MultiIndexType& end, RealType factor,
							GridFunctionType& sourcegridfunction)
{

}

void Stencil::setFxxStencil(const MultiIndexType& begin, 
							const MultiIndexType& end, RealType factor,
							GridFunctionType& sourcegridfunction, 
							MultiIndexType& offset)
{

}

void Stencil::setFxxStencil(const MultiIndexType& begin,
							const MultiIndexType& end, RealType factor,
							GridFunctionType& sourcegridfunction,
							MultiIndexType& offset, RealType constant)
{

}

//void Stencil::setFxxStencil(const MultiIndexType& begin, 
//							const MultiIndexType& end, RealType factor, 
//							GridFunctionType& sourcegridfunction)
//{
//
//}

void Stencil::setFxxStencil(const MultiIndexType& begin,
							const MultiIndexType& end)
{

}
