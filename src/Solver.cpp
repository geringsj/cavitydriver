#include "Solver.hpp"
#include "Debug.hpp"

#include <cmath>

namespace Solver
{
	Real computeResidual(
			const GridFunction& p,
			const GridFunction& rhs,
			const Point delta,
			const Ranges& inner_ranges,
			const Real global_fluidCells)
	{
		return 
			sqrt( computeSquaredResidual(p, rhs, delta, inner_ranges, global_fluidCells) );
	}

	Real computeSquaredResidual(
			const GridFunction& p,
			const GridFunction& rhs,
			const Point& delta,
			const Ranges& inner_ranges,
			const Real global_fluidCells)
	{
		Real numerator = 0.0;
		Real dxx = pow(delta.x, 2.0);
		Real dyy = pow(delta.y, 2.0);

		for_vecrange(i,j,inner_ranges)
		{
			Real pxx = (p(i+1,j) - 2.0*p(i,j) + p(i-1,j)) / dxx;
			Real pyy = (p(i,j+1) - 2.0*p(i,j) + p(i,j-1)) / dyy;
			Real help = pxx + pyy - rhs(i,j);

			numerator += pow(help,2.0);
		}
		return (numerator / global_fluidCells);
	}


	namespace {
		inline Real evaluateSOR(
				GridFunction& p, 
				const GridFunction& rhs,
				const int i, 
				const int j, 
				const Real dxx, 
				const Real dyy, 
				const Real omegaMinus1,
				const Real omegaTimesDxxDyy)
		{
			const Real old_value = p(i,j);
			const Real pxx = (p(i - 1,j) + p(i + 1,j)) / dxx;
			const Real pyy = (p(i,j - 1) + p(i,j + 1)) / dyy;

			return 
				omegaMinus1 * old_value 
				+ 
				omegaTimesDxxDyy 
				* ( pxx + pyy - rhs(i,j) );
		}
	};

	void SORCycle(
			GridFunction& p,
			const GridFunction& rhs,
			const Delta & delta,
			const Ranges& inner_ranges,
			const Real omega)
	{
		const Real dxx = pow(delta.x, 2.0);
		const Real dyy = pow(delta.y, 2.0);
		const Real OneMinusOmega = (1. - omega);
		const Real omegaTimesDxxDyy = omega * ((dxx*dyy)/(2.0*(dxx+dyy)));

		for_vecrange(i,j,inner_ranges)
		{
			p(i,j) = evaluateSOR(p, rhs, i, j, dxx, dyy, OneMinusOmega, omegaTimesDxxDyy);
		}
	}
	
	
	void SORCycleRedBlack(
			Domain& domain,
			const Real omega,
			const Color color)
	{
		GridFunction& p = domain.p();
		const GridFunction& rhs = domain.rhs();
		const Delta delta = domain.getDelta();
		const Ranges& inner_ranges = domain.getInnerRangeP();

		const Real dxx = pow(delta.x, 2.0);
		const Real dyy = pow(delta.y, 2.0);
		const Real OneMinusOmega = (1. - omega);
		const Real omegaTimesDxxDyy = omega * ((dxx*dyy)/(2.0*(dxx+dyy)));

		// Decide offset (in y-direction) based on color of first cell.
		// If color matches the color of the first cell, no offset is required 
		// in the first column, else start with offset +1. 
		// Subsequently, the offset will be flipped between 0 and 1 after each column
		// in order to achieve a check-board pattern (red-black scheme).

		//for(auto& r : inner_ranges)
		Range r = inner_ranges;
		{/* we need to recompute the offset for each subrange :( */
			const int xBegin = r.begin.i;
			const int yBegin = r.begin.j;
			const int xEnd = r.end.i;
			const int yEnd = r.end.j;

			int offset = (color == domain.getDomainFirstCellColor()) ? 0 : 1;
			int subRangeOffset = ((xBegin-1)%2+(yBegin-1)%2)%2; // sorry. this line exploits that the inner of domains starts at (1,1)
			offset = (offset + subRangeOffset) % 2;

			for(int j=yBegin; j<=yEnd; j++)
			{
				for(int i=xBegin+offset; i<=xEnd; i=i+2)// skip every other cell in x dimension
				{
					p(i,j) = 
						evaluateSOR(p, rhs, i, j, dxx, dyy, OneMinusOmega, omegaTimesDxxDyy);
				}
				offset = (offset+1)%2;// 'Flip' offset for next row
			}
		}
	}

}
