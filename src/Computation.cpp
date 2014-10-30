#include "Computation.hpp"

Computation::Computation()
{

}

Computation::~Computation()
{

}

RealType Computation::computeTimestep(RealType uMax, RealType vMax,
									  RealType tau, const PointType& h,
									  RealType Re)
{

}

void Computation::computeNewVelocities(GridFunction* u, GridFunction* v,
									   GridFunctionType& f, 
									   GridFunctionType& g,
									   GridFunctionType& p, 
									   const PointType& delta, RealType deltaT)
{

}

void Computation::computeMomentumEquations(GridFunctionType& f, 
										   GridFunctionType& g, 
										   GridFunction* u, GridFunction* v,
										   GridFunctionType& gx, 
										   GridFunctionType& gy, 
										   const PointType& delta, 
										   RealType deltaT)
{

}

void Computation::setBoundaryX(GridFunction x)
{

}

void Computation::computeRighthandSide(GridFunction* rhs, GridFunctionType& f,
									   GridFunctionType& g, 
									   const PointType& delta, RealType deltaT)
{

}