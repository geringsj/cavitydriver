
#ifndef Stencil2_hpp
#define Stencil2_hpp 

#include "Structs.hpp"
#include "GridFunction.hpp"
#include <functional>

class Stencil
{
public:
	enum class StencilOperator;
private:
	const Point h;

	inline Real Fx_forward(int i,int j, GridFunction& f)
	{
		return (f[i+1][j] - f[i][j]) / (h[0]);
	}
	inline Real Fy_forward(const int i, const int j, GridFunction& f)
	{
		return (f[i][j+1] - f[i][j]) / (h[1]);
	}
	inline Real Fx_backward(const int i, const int j, GridFunction& f)
	{
		return (f[i][j] - f[i-1][j]) / (h[0]);
	}
	inline Real Fy_backward(const int i, const int j, GridFunction& f)
	{
		return (f[i][j] - f[i][j-1]) / (h[1]);
	}


	inline Real Fxx(const int i, const int j, GridFunction& f)
	{
		return
			(f[i+1][j] - 2.0*f[i][j] + f[i-1][j]) / (h[0]*h[0]);
	}
	inline Real Fyy(const int i, const int j, GridFunction& f)
	{
		return
			(f[i][j+1] - 2.0*f[i][j] + f[i][j-1]) / (h[1]*h[1]);
	}

#define square(X) ((X)*(X))
	inline Real FFx(const int i, const int j, GridFunction& f)
	{
		return 
			(square(f[i+1][j] + f[i][j]) - square(f[i][j] + f[i-1][j]))
			/ (4.0*h[0]);
	}
	inline Real FFy(const int i, const int j, GridFunction& f)
	{
		GridFunction f = gf.getGridFunction();
		return 
			(square(f[i][j+1] + f[i][j]) - square(f[i][j] + f[i][j-1]))
			/ (4.0*h[1]);
	}

	inline Real FGx(const int i, const int j, 
			GridFunction& f,
			GridFunction& g)
	{
		return 
			/* d(uv)/dx , should be okay */
			((f[i][j+1] + f[i][j])*(g[i+1][j] + g[i][j]) 
			 - (f[i-1][j+1] + f[i-1][j])*(g[i][j] + g[i-1][j]))
			/ (4.0*h[0]);
	}
	inline Real FGy(const int i, const int j, 
			GridFunction& f,
			GridFunction& g)
	{
		return 
			((g[i+1][j] + g[i][j])*(f[i][j+1] + f[i][j]) 
			 - (g[i+1][j-1] + g[i][j-1])*(f[i][j] + f[i][j-1]))
			/ (4.0*h[1]);
	}

	inline 
		std::function<Real(Stencil*, const int i, const int j, GridFunction& gf)>
	getStencilOperatorFunctionSingle(const StencilOperator so)
		{
			switch(so)
			{
				case StencilOperator::Fx_forward:
					return &Stencil::Fx_forward;
					break;
				case StencilOperator::Fy_forward:
					return &Stencil::Fy_forward;
					break;
				case StencilOperator::Fx_backward:
					return &Stencil::Fx_backward;
					break;
				case StencilOperator::Fy_backward:
					return &Stencil::Fy_backward;
					break;
				case StencilOperator::Fxx:
					return &Stencil::Fxx;
					break;
				case StencilOperator::Fyy:
					return &Stencil::Fyy;
					break;
				case StencilOperator::FFx:
					return &Stencil::FFx;
					break;
				case StencilOperator::FFy:
					return &Stencil::FFy;
					break;
				default: /* uuhh... */
					break;
			};
			return 0;
		}
	inline std::function<Real(Stencil*, const int i, const int j, GridFunction& gf, GridFunction& gf2)>
	getStencilOperatorFunctionDouble(const StencilOperator so)
		{
			switch(so)
			{
				case StencilOperator::FGx:
					return &Stencil::FGx;
					break;
				case StencilOperator::FGy:
					return &Stencil::FGy;
					break;
				default: /* uuhh... */
					break;
			};
			return 0;
		}

public:
	Stencil(const Point h) : h(h) {};
	virtual ~Stencil();

#define forall(F,S,B,E) for(int F=B[0];F<=E[0];F++)for(int S=B[1];S<=E[1];S++)

	inline void applyStencilOperator(
		const StencilOperator stencil,
		const Index& gridreadbegin, 
		const Index& gridreadend,
		const Index& gridwritebegin,
		const Index& gridwriteend, 
		GridFunction& sourcegrid,
		GridFunction& targetgrid)
	{
		auto stenciloperator = 
			getStencilOperatorFunctionSingle(stencil);
		forall(i,j,gridreadbegin,gridreadend)
		{
			targetgrid(
					gridwritebegin[0]+i-gridreadbegin[0],
					gridwritebegin[1]+j-gridreadbegin[1])
				= stenciloperator(this,i,j,sourcegrid);
		}
	}

	inline void applyStencilOperator(
		const StencilOperator stencil,
		const Index& gridreadbegin, 
		const Index& gridreadend,
		const Index& gridwritebegin,
		const Index& gridwriteend, 
		GridFunction& sourcegrid,
		GridFunction& sourcegrid2,
		GridFunction& targetgrid)
	{
		auto stenciloperator = 
			getStencilOperatorFunctionDouble(stencil);
		forall(i,j,gridreadbegin,gridreadend)
		{
			targetgrid(
					gridwritebegin[0]+i-gridreadbegin[0],
					gridwritebegin[1]+j-gridreadbegin[1])
				= stenciloperator(this,i,j,sourcegrid,sourcegrid2);
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

