#ifndef Solver_hpp
#define Solver_hpp

#include "GridFunction.hpp"
#include "Structs.hpp"

#include <math.h>

namespace Solver
{

	Real computeResidual(
			GridFunction& p,
			GridFunction& rhs,
			const Point& delta,
			Dimension inner_begin, Dimension inner_end)
	{
		Real res;

		Real numerator = 0.0;
		Real denominator = 
			Real((inner_end.i-inner_begin.i+1) * (inner_end.j-inner_begin.j+1));

		Real dxx_numerator;
		Real dxx_denominator = (delta.x * delta.x);

		Real dyy_numerator;
		Real dyy_denominator = (delta.y * delta.y);
		for (int i = inner_begin.i; i <= inner_end.i; i++)
		{
			for (int j = inner_begin.j; j <= inner_end.j; j++)
			{
				dxx_numerator = p(i+1,j) - 2.0 * p(i,j) + p(i-1,j);
				dyy_numerator = p(i,j+1) - 2.0 * p(i,j) + p(i,j-1);
				numerator = numerator + 
					pow(
						(dxx_numerator / dxx_denominator + 
						 dyy_numerator / dyy_denominator - rhs(i,j)
						 )
					,2.0);
			}
		}
		res = sqrt(numerator / denominator);
		return res;
	}
	
	void SORCycle(
			GridFunction& p,
			GridFunction& rhs,
			const Point& delta,
			Dimension inner_begin, Dimension inner_end,
			Real omega)
	{
		Real dxx_numerator;
		Real dyy_numerator;
		Real dxx = pow(delta.x, 2.0);
		Real dyy = pow(delta.y, 2.0);
		Real old_value, new_value;

		for (int i = inner_begin.i; i <= inner_end.i; i++)
		{
			for (int j = inner_begin.j; j <= inner_end.j; j++)
			{
				old_value = p(i,j);
				dxx_numerator = p(i - 1,j) + p(i + 1,j);
				dyy_numerator = p(i,j - 1) + p(i,j + 1);

				new_value = (1. - omega) * old_value + 
					omega * ((dxx * dyy)/(2.0 * (dxx + dyy))) * 
						( (dxx_numerator/dxx) + (dyy_numerator/dyy) - rhs(i,j) );

				p(i,j) = new_value;
			}
		}
	}
	
}

#endif
