#include "Solver.hpp"
#include "Debug.hpp"

#include <cmath>

namespace Solver
{
	Real computeResidual(
			const GridFunction& p,
			const GridFunction& rhs,
			const Point& delta,
			const Dimension inner_begin, 
			const Dimension inner_end,
			const Dimension globalDims)
	{
		return 
			sqrt( computeSquaredResidual(p, rhs, delta, inner_begin, inner_end, globalDims) );
	}

	Real computeSquaredResidual(
			const GridFunction& p,
			const GridFunction& rhs,
			const Point& delta,
			const Dimension inner_begin, 
			const Dimension inner_end,
			const Dimension globalDims)
	{
		Real numerator = 0.0;
		Real denominator = (globalDims.i * globalDims.j);

		Real dxx = pow(delta.x, 2.0);
		Real dyy = pow(delta.y, 2.0);

		forall(i,j,inner_begin,inner_end)
		{
			Real pxx = (p(i+1,j) - 2.0*p(i,j) + p(i-1,j)) / dxx;
			Real pyy = (p(i,j+1) - 2.0*p(i,j) + p(i,j-1)) / dyy;

			Real help = pxx + pyy - rhs(i,j);

			numerator += pow(help,2.0);
		}

		return (numerator / denominator);
	}


	namespace {
		inline void evaluateSOR(
				GridFunction& p, 
				const GridFunction& rhs,
				const int& i, 
				const int& j, 
				const Real& dxx, 
				const Real& dyy, 
				const Real& omega)
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
	};

	void SORCycle(
			GridFunction& p,
			const GridFunction& rhs,
			const Point& delta,
			const Dimension& inner_begin, 
			const Dimension& inner_end,
			const Real& omega)
	{
		Real dxx = pow(delta.x, 2.0);
		Real dyy = pow(delta.y, 2.0);

		forall(i,j,inner_begin,inner_end)
		{
			evaluateSOR(p, rhs, i, j, dxx, dyy, omega);
			//Real old_value = p(i,j);
			//Real pxx = (p(i - 1,j) + p(i + 1,j)) / dxx;
			//Real pyy = (p(i,j - 1) + p(i,j + 1)) / dyy;

			//Real new_value = 
			//	(1. - omega) * old_value 
			//	+ 
			//	omega 
			//	* ((dxx*dyy)/(2.0*(dxx+dyy))) 
			//	* ( pxx + pyy - rhs(i,j) );

			//p(i,j) = new_value;
		}
	}
	
	
	void SORCycleRedBlack(
			Domain& domain,
			const Real& omega,
			const Color& color)
	{
		GridFunction& p = domain.p();
		const GridFunction& rhs = domain.rhs();
		const Delta& delta = domain.getDelta();
		const Dimension& inner_begin = domain.getBeginInnerDomains();
		const Dimension& inner_end = domain.getEndInnerDomainP();

		Real dxx = pow(delta.x, 2.0);
		Real dyy = pow(delta.y, 2.0);

		// Decided offset (in y-direction) based on color of first cell.
		// If color matches the color of the first cell, no offset is required in the first column,
		// else start with offset +1. Subsequently, the offset will be flipped between 0 and 1 after each column
		// in order to achieve a check-board pattern (red-black scheme).
		int offset = (color == domain.getDomainFirstCellColor()) ? 0 : 1;

		for(int i=inner_begin[0]; i<=inner_end[0]; i++)
		{
			for(int j=inner_begin[1]+offset; j<=inner_end[1]; j=j+2)// skip every other cell in y dimension
			{
				evaluateSOR(p, rhs, i, j, dxx, dyy, omega);
				//Real old_value = p(i,j);
				//Real pxx = (p(i - 1,j) + p(i + 1,j)) / dxx;
				//Real pyy = (p(i,j - 1) + p(i,j + 1)) / dyy;

				//Real new_value = 
				//	(1. - omega) * old_value 
				//	+ 
				//	omega 
				//	* ((dxx*dyy)/(2.0*(dxx+dyy))) 
				//	* ( pxx + pyy - rhs(i,j) );

				//p(i,j) = new_value;
			}
			// 'Flip' offset for next column
			offset = (offset+1)%2;
		}

	}

}
