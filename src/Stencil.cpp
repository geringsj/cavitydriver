#include "Stencil.hpp"

Stencil::Stencil(int stencilwidth, const PointType h)
	: stencilwidth(stencilwidth) , h(h) 
{
	/* ignore stencilwidth ... 
	 * we don't take kindly to freakish stencils around here! */
	stencil = new RealType*[3];
	stencil++;
	stencil[-1] = &this->rawmemory[1]; // shift by 1 for [-1] access
	stencil[0] = &this->rawmemory[4];
	stencil[1] = &this->rawmemory[7];
}

Stencil::~Stencil()
{
	stencil[-1] = NULL;
	stencil[0] = NULL;
	stencil[1] = NULL;
	stencil--;
	delete[] stencil;
}

void Stencil::ApplyStencilOperator(const MultiIndexType& gridbegin, 
								   const MultiIndexType& gridreadend,
								   const MultiIndexType& gridwritebegin,
								   const MultiIndexType& gridwriteend, 
								   const GridFunction sourcegridfunction,
								   GridFunction imagegridfunction)
{
}

/* stencil is set as stencil[x][y],
 * as in 'U_i,j' with i moving in x direction (and j in y)
 */
void Stencil::setFxStencil()
{
	stencil[-1][-1] = 0.0;
	stencil[-1][0]  = 0.0;
	stencil[-1][1]  = 0.0;
	stencil[0][-1]  = 0.0;
	stencil[0][0]   = 0.0;
	stencil[0][1]   = 0.0;
	stencil[1][-1]  = 0.0;
	stencil[1][0]   = 0.0;
	stencil[1][1]   = 0.0;
}

void Stencil::setFyStencil()
{
	stencil[-1][-1] = 0.0;
	stencil[-1][0]  = 0.0;
	stencil[-1][1]  = 0.0;
	stencil[0][-1]  = 0.0;
	stencil[0][0]   = 0.0;
	stencil[0][1]   = 0.0;
	stencil[1][-1]  = 0.0;
	stencil[1][0]   = 0.0;
	stencil[1][1]   = 0.0;
}

void Stencil::setFFyStencil()
{
	stencil[-1][-1] = 0.0;
	stencil[-1][0]  = 0.0;
	stencil[-1][1]  = 0.0;
	stencil[0][-1]  = 0.0;
	stencil[0][0]   = 0.0;
	stencil[0][1]   = 0.0;
	stencil[1][-1]  = 0.0;
	stencil[1][0]   = 0.0;
	stencil[1][1]   = 0.0;
}
void Stencil::setFFxStencil()
{
	stencil[-1][-1] = 0.0;
	stencil[-1][0]  = 0.0;
	stencil[-1][1]  = 0.0;
	stencil[0][-1]  = 0.0;
	stencil[0][0]   = 0.0;
	stencil[0][1]   = 0.0;
	stencil[1][-1]  = 0.0;
	stencil[1][0]   = 0.0;
	stencil[1][1]   = 0.0;
}

void Stencil::setFxxStencil()
{
	stencil[-1][-1] = 0.0;
	stencil[-1][0]  = 1.0 / (h[0] * h[0]);
	stencil[-1][1]  = 0.0;
	stencil[0][-1]  = 0.0;
	stencil[0][0]   = -2.0 / (h[0] * h[0]);
	stencil[0][1]   = 0.0;
	stencil[1][-1]  = 0.0;
	stencil[1][0]   = 1.0 / (h[0] * h[0]);
	stencil[1][1]   = 0.0;
}

void Stencil::setFyyStencil()
{
	stencil[-1][-1] = 0.0;
	stencil[-1][0]  = 0.0;
	stencil[-1][1]  = 0.0;
	stencil[0][-1]  = 1.0 / (h[1] * h[1]);
	stencil[0][0]   = -2.0 / (h[1] * h[1]);
	stencil[0][1]   = 1.0 / (h[1] * h[1]);
	stencil[1][-1]  = 0.0;
	stencil[1][0]   = 0.0;
	stencil[1][1]   = 0.0;
}

