
#ifndef Stencil2_hpp
#define Stencil2_hpp 

#include "typedef.hpp"
#include "GridFunction.hpp"
#include <functional>

class Stencil
{
public:
	enum class StencilOperator;
private:
	const PointType h;

	inline RealType Fx_forward(int i,int j, GridFunction& gf)
	{
		GridFunctionType f = gf.getGridFunction();
		return (f[i+1][j] - f[i][j]) / (h[0]);
	}
	inline RealType Fy_forward(const int i, const int j, GridFunction& gf)
	{
		GridFunctionType f = gf.getGridFunction();
		return (f[i][j+1] - f[i][j]) / (h[1]);
	}
	inline RealType Fx_backward(const int i, const int j, GridFunction& gf)
	{
		GridFunctionType f = gf.getGridFunction();
		return (f[i][j] - f[i-1][j]) / (h[0]);
	}
	inline RealType Fy_backward(const int i, const int j, GridFunction& gf)
	{
		GridFunctionType f = gf.getGridFunction();
		return (f[i][j] - f[i][j-1]) / (h[1]);
	}


	inline RealType Fxx(const int i, const int j, GridFunction& gf)
	{
		GridFunctionType f = gf.getGridFunction();
		return
			(f[i+1][j] - 2.0*f[i][j] + f[i-1][j]) / (h[0]*h[0]);
	}
	inline RealType Fyy(const int i, const int j, GridFunction& gf)
	{
		GridFunctionType f = gf.getGridFunction();
		return
			(f[i][j+1] - 2.0*f[i][j] + f[i][j-1]) / (h[1]*h[1]);
	}

#define square(X) ((X)*(X))
	inline RealType FFx(const int i, const int j, GridFunction& gf)
	{
		GridFunctionType f = gf.getGridFunction();
		return 
			(square(f[i+1][j] + f[i][j]) - square(f[i][j] + f[i-1][j]))
			/ (4.0*h[0]);
	}
	inline RealType FFy(const int i, const int j, GridFunction& gf)
	{
		GridFunctionType f = gf.getGridFunction();
		return 
			(square(f[i][j+1] + f[i][j]) - square(f[i][j] + f[i][j-1]))
			/ (4.0*h[1]);
	}

	inline RealType FGx(const int i, const int j, 
			GridFunction& gf1,
			GridFunction& gf2)
	{
		GridFunctionType f = gf1.getGridFunction();
		GridFunctionType g = gf2.getGridFunction();
		return 
			/* d(uv)/dx , should be okay */
			((f[i][j+1] + f[i][j])*(g[i+1][j] + g[i][j]) 
			 - (f[i-1][j+1] + f[i-1][j])*(g[i][j] + g[i-1][j]))
			/ (4.0*h[0]);
	}
	inline RealType FGy(const int i, const int j, 
			GridFunction& gf1,
			GridFunction& gf2)
	{
		GridFunctionType f = gf1.getGridFunction();
		GridFunctionType g = gf2.getGridFunction();
		return 
			((g[i+1][j] + g[i][j])*(f[i][j+1] + f[i][j]) 
			 - (g[i+1][j-1] + g[i][j-1])*(f[i][j] + f[i][j-1]))
			/ (4.0*h[1]);
	}

	inline std::function<RealType(const int i, const int j)>
		getStencilOperatorFunction(const StencilOperator so,
				GridFunction gf)
		{
			switch(so)
			{
				case StencilOperator::Fx_forward:
					return std::bind
						(&Stencil::Fx_forward, this,
						 std::placeholders::_1,
						 std::placeholders::_2,
						 std::ref(gf));
					break;
				case StencilOperator::Fy_forward:
					return std::bind
						(&Stencil::Fy_forward, this,
						 std::placeholders::_1,
						 std::placeholders::_2,
						 std::ref(gf));
					break;
				case StencilOperator::Fx_backward:
					return std::bind
						(&Stencil::Fx_backward, this,
						 std::placeholders::_1,
						 std::placeholders::_2,
						 std::ref(gf));
					break;
				case StencilOperator::Fy_backward:
					return std::bind
						(&Stencil::Fy_backward, this,
						 std::placeholders::_1,
						 std::placeholders::_2,
						 std::ref(gf));
					break;
				case StencilOperator::Fxx:
					return std::bind
						(&Stencil::Fxx, this,
						 std::placeholders::_1,
						 std::placeholders::_2,
						 std::ref(gf));
					break;
				case StencilOperator::Fyy:
					return std::bind
						(&Stencil::Fyy, this,
						 std::placeholders::_1,
						 std::placeholders::_2,
						 std::ref(gf));
					break;
				case StencilOperator::FFx:
					return std::bind
						(&Stencil::FFx, this,
						 std::placeholders::_1,
						 std::placeholders::_2,
						 std::ref(gf));
					break;
				case StencilOperator::FFy:
					return std::bind
						(&Stencil::FFy, this,
						 std::placeholders::_1,
						 std::placeholders::_2,
						 std::ref(gf));
					break;
				default: /* ehm */
					break;
			};
			return 0;
		}

	inline std::function<RealType(const int i, const int j)>
		getStencilOperatorFunction(const StencilOperator so,
				GridFunction gf1,
				GridFunction gf2)
		{
			switch(so)
			{
				case StencilOperator::FGx:
					return std::bind
						(&Stencil::FGx, this,
						 std::placeholders::_1,
						 std::placeholders::_2,
						 std::ref(gf1), std::ref(gf2));
					break;
				case StencilOperator::FGy:
					return std::bind
						(&Stencil::FGy, this,
						 std::placeholders::_1,
						 std::placeholders::_2,
						 std::ref(gf1), std::ref(gf2));
					break;
				default: /* uuhh... */
					break;
			};
			return 0;
		}

public:
	Stencil(const PointType h) : h(h) {};
	virtual ~Stencil();

#define forall(F,S,B,E) for(int F=B[0];F<=E[0];F++)for(int S=B[1];S<=E[1];S++)

	inline void applyStencilOperator(
		const StencilOperator stencil,
		const MultiIndexType& gridreadbegin, 
		const MultiIndexType& gridreadend,
		const MultiIndexType& gridwritebegin,
		const MultiIndexType& gridwriteend, 
		GridFunction& sourcegridfunction,
		GridFunction& targetgridfunction)
	{
		auto stenciloperator = 
			getStencilOperatorFunction(stencil, sourcegridfunction);
		forall(i,j,gridreadbegin,gridreadend)
		{
			targetgridfunction.getGridFunction()
					[gridwritebegin[0]+i-gridreadbegin[0]]
					[gridwritebegin[1]+j-gridreadbegin[1]] 
				= stenciloperator(i,j);
		}
	}

	inline void applyStencilOperator(
		const StencilOperator stencil,
		const MultiIndexType& gridreadbegin, 
		const MultiIndexType& gridreadend,
		const MultiIndexType& gridwritebegin,
		const MultiIndexType& gridwriteend, 
		GridFunction& sourcegridfunction1,
		GridFunction& sourcegridfunction2,
		GridFunction& targetgridfunction)
	{
		auto stenciloperator = 
			getStencilOperatorFunction(stencil, 
					sourcegridfunction1,sourcegridfunction2);
		forall(i,j,gridreadbegin,gridreadend)
		{
			targetgridfunction.getGridFunction()
					[gridwritebegin[0]+i-gridreadbegin[0]]
					[gridwritebegin[1]+j-gridreadbegin[1]] 
				= stenciloperator(i,j);
		}
	}

	enum class StencilOperator {
		Fx_forward,
		Fy_forward,
		Fx_backward,
		Fy_backward,
		Fxx,
		Fyy,
		FFx,
		FFy,
		FGx,
		FGy
	};
};

#endif
