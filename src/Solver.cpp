#include "Solver.hpp"
#include "Debug.hpp"

#include <cmath>

namespace Solver
{
	Real computeResidual(
			const GridFunction& p,
			const GridFunction& rhs,
			const Point delta,
			const Range inner_range,
			const Dimension globalDims)
	{
		return 
			sqrt( computeSquaredResidual(p, rhs, delta, inner_range, globalDims) );
	}

	Real computeSquaredResidual(
			const GridFunction& p,
			const GridFunction& rhs,
			const Point& delta,
			const Range inner_range,
			const Dimension globalDims)
	{
		Real numerator = 0.0;
		Real denominator = (globalDims.i * globalDims.j);
		Real dxx = pow(delta.x, 2.0);
		Real dyy = pow(delta.y, 2.0);

		for_range(i,j,inner_range)
		{
			Real pxx = (p(i+1,j) - 2.0*p(i,j) + p(i-1,j)) / dxx;
			Real pyy = (p(i,j+1) - 2.0*p(i,j) + p(i,j-1)) / dyy;
			Real help = pxx + pyy - rhs(i,j);

			numerator += pow(help,2.0);
		}
		return (numerator / denominator);
	}


	namespace {
		inline Real evaluateSOR(
				const GridFunction& p, 
				const GridFunction& rhs,
				const int& i, 
				const int& j, 
				const Real& dxx, 
				const Real& dyy, 
				const Real& omegaMinus1,
				const Real& omegaTimesDxxDyy)
		{
			const Real old_value = p(i,j);
			const Real pxx = (p(i - 1,j) + p(i + 1,j)) / dxx;
			const Real pyy = (p(i,j - 1) + p(i,j + 1)) / dyy;

			const Real new_value = 
				omegaMinus1 * old_value 
				+ 
				omegaTimesDxxDyy 
				* ( pxx + pyy - rhs(i,j) );

			return new_value;
		}
	};

	void SORCycle(
			GridFunction& p,
			const GridFunction& rhs,
			const Delta & delta,
			const Range& inner_range,
			const Real& omega)
	{
		const Real dxx = pow(delta.x, 2.0);
		const Real dyy = pow(delta.y, 2.0);
		const Real omegaMinus1 = (1. - omega);
		const Real omegaTimesDxxDyy = omega * ((dxx*dyy)/(2.0*(dxx+dyy)));

		for_range(i,j,inner_range)
		{
			p(i,j) = evaluateSOR(p, rhs, i, j, dxx, dyy, omegaMinus1, omegaTimesDxxDyy);
		}
	}
	
	
	void SORCycleRedBlack(
			Domain& domain,
			const Real& omega,
			const Color& color)
	{
		GridFunction& p = domain.p();
		const GridFunction& rhs = domain.rhs();
		const Delta delta = domain.getDelta();
		const Range inner_range = domain.getInnerRangeP();

		const Real dxx = pow(delta.x, 2.0);
		const Real dyy = pow(delta.y, 2.0);
		const Real omegaMinus1 = (1. - omega);
		const Real omegaTimesDxxDyy = omega * ((dxx*dyy)/(2.0*(dxx+dyy)));

		// Decided offset (in y-direction) based on color of first cell.
		// If color matches the color of the first cell, no offset is required 
		// in the first column, else start with offset +1. 
		// Subsequently, the offset will be flipped between 0 and 1 after each column
		// in order to achieve a check-board pattern (red-black scheme).
		int offset = (color == domain.getDomainFirstCellColor()) ? 0 : 1;

		for(int i=inner_range.begin[0]; i<=inner_range.end[0]; i++)
		{
			for(int j=inner_range.begin[1]+offset; j<=inner_range.end[1]; j=j+2)
				// skip every other cell in y dimension
			{
				p(i,j) = evaluateSOR(p, rhs, i, j, dxx, dyy, omegaMinus1, omegaTimesDxxDyy);
			}
			// 'Flip' offset for next column
			offset = (offset+1)%2;
		}
	}

}
