#include "Solver.hpp"
#include "Debug.hpp"

#include <cmath>

namespace Solver
{


	Real computeResidual(
			GridFunction& p,
			GridFunction& rhs,
			const Point& delta,
			Dimension inner_begin, Dimension inner_end)
	{
		Real numerator = 0.0;
		Real denominator = 
			Real((inner_end.i-inner_begin.i+1) * (inner_end.j-inner_begin.j+1));

		Real dxx = pow(delta.x, 2.0);
		Real dyy = pow(delta.y, 2.0);

		forall(i,j,inner_begin,inner_end)
		{
			Real pxx = (p(i+1,j) - 2.0*p(i,j) + p(i-1,j)) / dxx;
			Real pyy = (p(i,j+1) - 2.0*p(i,j) + p(i,j-1)) / dyy;
			Real help = pxx + pyy - rhs(i,j);
			numerator += pow(help,2.0);
		}

		Real res = sqrt(numerator / denominator);
		return res;
	}
	

	void SORCycle(
			GridFunction& p,
			GridFunction& rhs,
			const Point& delta,
			Dimension inner_begin, Dimension inner_end,
			Real omega)
	{
		Real dxx = pow(delta.x, 2.0);
		Real dyy = pow(delta.y, 2.0);

		forall(i,j,inner_begin,inner_end)
		{
			Real old_value = p(i,j);
			Real pxx = (p(i - 1,j) + p(i + 1,j)) / dxx;
			Real pyy = (p(i,j - 1) + p(i,j + 1)) / dyy;

			Real new_value = 
				(1. - omega) * old_value 
				+ 
				omega 
				* ((dxx*dyy)/(2.0*(dxx+dyy))) 
				* ( pxx + pyy - rhs(i,j) );

			p(i,j) = new_value;
		}
	}
	

}
