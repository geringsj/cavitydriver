#include "Solver.hpp"

namespace Solver
{

	RealType computeResidual(GridFunctionType& sourcegridfunction,
											GridFunctionType& rhs,
											const PointType& h
											/* There is no way
											 * to get iMax and
											 * jMax from a 
											 * GridFunctionType
											 * variable. So I
											 * added two int values
											 */
											 ,int iMax, int jMax)
	{
		RealType res;

		RealType numerator;
		RealType denominator = iMax * jMax;

		RealType dxx_numerator;
		RealType dxx_denominator = (h[0] * h[0]);

		RealType dyy_numerator;
		RealType dyy_denominator = (h[1] * h[1]);
		for (int i = 0; i < iMax; i++)
		{
			for (int j = 0; j < jMax; j++)
			{
				dxx_numerator = sourcegridfunction[i+1][j] - 
					2.0 * sourcegridfunction[i][j] + 
					sourcegridfunction[i-1][j];
				dyy_numerator = sourcegridfunction[i][j+1] - 
					2.0 * sourcegridfunction[i][j] + 
					sourcegridfunction[i][j-1];
				numerator = numerator + 
					pow((dxx_numerator / dxx_denominator + 
					dyy_numerator / dyy_denominator - 
					rhs[i][j]),2.0);
			}
		}
		res = sqrt(numerator / denominator);
		return res;
	}
	
	void SORCycle(GridFunction* gridfunction,
							GridFunctionType& rhs,
							const PointType& delta,
							RealType omega)
	{
		int iMax = gridfunction->getGridDimension()[0];
		int jMax = gridfunction->getGridDimension()[1];
		MultiIndexType begin;
		MultiIndexType end;
		RealType value;
		RealType dxx_numerator;
		RealType dyy_numerator;
		RealType dxx = pow(delta[0], 2.0);
		RealType dyy = pow(delta[1], 2.0);
		RealType old_value;
		for (int i = 0; i < iMax; i++)
		{
			for (int j = 0; j < jMax; j++)
			{
				/**
				 * Get the new value.
				 * Old value = p_{i,j}
				 * dxx_numerator = 
				 * p_{i-1,j} + p_{i+1,j}
				 * dyy_numerator = 
				 * p_{i,j-1} + p_{i,j+1}
				 * dxx = dx²
				 * dyy = dy²
				 */
				old_value = gridfunction->
					getGridFunction()[i][j];
				dxx_numerator = gridfunction->
					getGridFunction()[i - 1][j] + 
					gridfunction->
					getGridFunction()[i + 1][j];
				dyy_numerator = gridfunction->
					getGridFunction()[i][j - 1] + 
					gridfunction->
					getGridFunction()[i][j + 1];
				value = old_value + (1 - omega) * old_value
					+ omega * (dxx * dyy) / (2.0 * (dxx + dyy))
					* (dxx_numerator / dxx + dyy_numerator/dyy 
					- rhs[i][j]);
				/**
				 * I'm not sure about how we set p_{i,j} so I 
				 * set the begin and the end to the same value.
				 * If this implementation is correct only 
				 * p_{i,j} should be updated. If this is not 
				 * correct please fix it, thanks.
				 */
				begin[0] = i;
				begin[1] = j;
				end = begin;
				gridfunction->setGridFunction(begin, end, value);
			}
		}
	}
	
}
