
#ifndef Stencil2_hpp
#define Stencil2_hpp 

#include "Structs.hpp"
#include "GridFunction.hpp"
#include <functional>

namespace Derivatives
{
	enum class Function {
		U,
		V,
		W,
		P
	};
	enum class Direction {
		xf, xb,
		yf, yb,
		zf, zb,
		xx, yy, zz,
		_x, _y, _z
	};

	/* attention: we don't handle weird Function-Direction cases here, 
	 * only the ones we will use later
	 */

	std::function<Real(int,int)> genFd(
			GridFunction& gf, Delta d,
			Direction df)
	{
		Real delta;
		int dimX=0, dimY=0, dimZ=0;
		int add1=0, add2=0;
		switch(df) /* set right dimension */
		{
			case Direction::xf:
			case Direction::xb:
				delta = d.x;
				dimX = 1;
				dimY = 0;
				dimZ = 0;
				break;
			case Direction::yf:
			case Direction::yb:
				delta = d.y;
				dimX = 0;
				dimY = 1;
				dimZ = 0;
				break;
			case Direction::zf:
			case Direction::zb:
				delta = d.z;
				dimX = 0;
				dimY = 0;
				dimZ = 1;
				break;
			default:
				break;
		}
		switch(df) /* set right derivative */
		{
			case Direction::xf:
			case Direction::yf:
			case Direction::zf:
				/* forward */
				add1 = +1; /* index on the left */
					/* minus */
				add2 = 0; /* index on the right */
				break;
			case Direction::xb:
			case Direction::yb:
			case Direction::zb:
				/* backward */
				add1 = 0; /* index on the left */
					/* minus */
				add2 = -1; /* index the right */
				break;
			default:
				return [](int i, int j) { return i*j*0.0; };
				break;
		}
		return [&gf, delta, add1, add2, dimX, dimY, dimZ](int i, int j)
			{
				return (gf( i+ add1*dimX , j + add1*dimY) - gf( i + add2*dimX, j + add2*dimY)) / delta;
			};
	}

	std::function<Real(int,int)> genFdd(
			GridFunction& gf, Delta d,
			Direction df)
	{
		Real delta = 0.0;
		int dimX=0, dimY=0, dimZ=0;
		int add1=+1, add2=0, add3=-1;
		switch(df) /* set right dimension */
		{
			case Direction::xx:
				delta = d.x;
				dimX = 1;
				dimY = 0;
				dimZ = 0;
				break;
			case Direction::yy:
				delta = d.y;
				dimX = 0;
				dimY = 1;
				dimZ = 0;
				break;
			case Direction::zz:
				delta = d.z;
				dimX = 0;
				dimY = 0;
				dimZ = 1;
				break;
			default:
				break;
		}
		return [&gf, delta, add1, add2, add3, dimX, dimY, dimZ](int i, int j)
			{
				return 
					(gf( i+ add1*dimX , j + add1*dimY) 
					 - 2.0*gf( i + add2*dimX, j + add2*dimY) 
					 + gf( i + add3*dimX, j + add3*dimY)) 
					/ (delta*delta);
			};
	}

	std::function<Real(int, int)> genFG_d(
			GridFunction& gf1, GridFunction& gf2, 
			Delta d, 
			Function f1, Function f2, Direction df)
	{
		Real delta = 0.0;
		int f1dimX=0, f1dimY=0, f1dimZ=0;
		int f2dimX=0, f2dimY=0, f2dimZ=0;
		switch(df) /* set right dimension */
		{
			case Direction::xx:
				delta = d.x;
				break;
			case Direction::yy:
				delta = d.y;
				break;
			case Direction::zz:
				delta = d.z;
				break;
			default:
				break;
		}
		switch(f1)
		{
			case Function::U:
				f1dimX = 1;
				f1dimY = 0;
				f1dimZ = 0;
				break;
			case Function::V:
				f1dimX = 0;
				f1dimY = 1;
				f1dimZ = 0;
				break;
			case Function::W:
				f1dimX = 0;
				f1dimY = 0;
				f1dimZ = 1;
				break;
			default:
				break;
		}
		switch(f2)
		{
			case Function::U:
				f2dimX = 1;
				f2dimY = 0;
				f2dimZ = 0;
				break;
			case Function::V:
				f2dimX = 0;
				f2dimY = 1;
				f2dimZ = 0;
				break;
			case Function::W:
				f2dimX = 0;
				f2dimY = 0;
				f2dimZ = 1;
				break;
			default:
				break;
		}
		int add=+1, sub=-1;
		return [&gf1, &gf2, delta, add, sub, 
				 f1dimX, f1dimY, f1dimZ,
				 f2dimX, f2dimY, f2dimZ
			](int i, int j)
			{
				Real fgp = 
					(gf1( i , j ) + gf1( i + f2dimX*add , j + f2dimY*add)) 
					* 
					(gf2( i , j ) + gf2( i + f1dimX*add, j + f1dimY*add));
				Real fgm = 
					(gf1( i + f2dimX*sub, j + f2dimY*sub) + gf1( i , j )) 
					* 
					(gf2( i + f2dimX*sub, j + f2dimY*sub) 
					 + 
					 gf2( i + f2dimX*sub + f1dimX*add, j + f2dimY*sub + f1dimY*add));
				return (fgp - fgm) / (4.0*delta);
			};
	}


	/* call those from outside */

	std::function<Real(int, int)> getDerivative(
			GridFunction& gf, 
			Delta d, 
			Direction df)
	{
		switch(df)
		{
			case Direction::xf:
			case Direction::yf:
			case Direction::zf:
			case Direction::xb:
			case Direction::yb:
			case Direction::zb:
				return Derivatives::genFd(gf, d, df);
				break;
			case Direction::xx:
			case Direction::yy:
			case Direction::zz:
				return genFdd(gf, d, df);
				break;
			default:
				return [](int i, int j){ return i*j*0.0; };
				break;
		};
	}
	std::function<Real(int, int)> getDerivative(
			GridFunction& gf, 
			Delta d, 
			int df)
	{
		Direction dff = Direction::xf;
		switch(df)
		{
			case -3:
				dff = Direction::zb;
				break;
			case -2:
				dff = Direction::yb;
				break;
			case -1:
				dff = Direction::xb;
				break;
			case 3:
				dff = Direction::zf;
				break;
			case 2:
				dff = Direction::yf;
				break;
			case 1:
				dff = Direction::xf;
				break;
			case 6:
				dff = Direction::zz;
				break;
			case 5:
				dff = Direction::yy;
				break;
			case 4:
				dff = Direction::xx;
				break;
			default:
				break;
		};
		return Derivatives::getDerivative(gf, d, dff);
	}

	std::function<Real(int, int)> getDerivative(
			GridFunction& gf1, GridFunction& gf2, 
			Delta d, 
			Function f1, Function f2, Direction df)
	{
		switch(df)
		{
			case Direction::_x:
			case Direction::_y:
			case Direction::_z:
				return genFG_d(gf1, gf2, d, f1, f2, df);
				break;
			default:
				return [](int i, int j){ return i*j*0.0; };
				break;
		}
	}
	std::function<Real(int, int)> getDerivative(
			GridFunction& gf1, GridFunction& gf2, 
			Delta d, 
			int f1, int f2, int df)
	{
		Function ff1 = Derivatives::Function::U;
		Function ff2 = Derivatives::Function::U;
		Direction dff = Derivatives::Direction::_x;

		switch(f1)
		{
			case 1:
				ff1 = Function::U;
				break;
			case 2:
				ff1 = Function::V;
				break;
			case 3:
				ff1 = Function::W;
				break;
		}
		switch(f2)
		{
			case 1:
				ff2 = Function::U;
				break;
			case 2:
				ff2 = Function::V;
				break;
			case 3:
				ff2 = Function::W;
				break;
		}
		switch(df)
		{
			case 1:
				dff = Direction::_x;
				break;
			case 2:
				dff = Direction::_y;
				break;
			case 3:
				dff = Direction::_z;
				break;
		}
		return getDerivative(gf1, gf2, d, ff1, ff2, dff);
	}

};












//  class Stencil
//  {
//  public:
//  	enum class StencilOperator;
//  private:
//  	const Point h;
//  
//  	inline Real Fx_forward(int i,int j, GridFunction& f)
//  	{
//  		return (f[i+1][j] - f[i][j]) / (h[0]);
//  	}
//  	inline Real Fy_forward(const int i, const int j, GridFunction& f)
//  	{
//  		return (f[i][j+1] - f[i][j]) / (h[1]);
//  	}
//  	inline Real Fx_backward(const int i, const int j, GridFunction& f)
//  	{
//  		return (f[i][j] - f[i-1][j]) / (h[0]);
//  	}
//  	inline Real Fy_backward(const int i, const int j, GridFunction& f)
//  	{
//  		return (f[i][j] - f[i][j-1]) / (h[1]);
//  	}
//  
//  
//  	inline Real Fxx(const int i, const int j, GridFunction& f)
//  	{
//  		return
//  			(f[i+1][j] - 2.0*f[i][j] + f[i-1][j]) / (h[0]*h[0]);
//  	}
//  	inline Real Fyy(const int i, const int j, GridFunction& f)
//  	{
//  		return
//  			(f[i][j+1] - 2.0*f[i][j] + f[i][j-1]) / (h[1]*h[1]);
//  	}
//  
//  #define square(X) ((X)*(X))
//  	inline Real FFx(const int i, const int j, GridFunction& f)
//  	{
//  		return 
//  			(square(f[i+1][j] + f[i][j]) - square(f[i][j] + f[i-1][j]))
//  			/ (4.0*h[0]);
//  	}
//  	inline Real FFy(const int i, const int j, GridFunction& f)
//  	{
//  		GridFunction f = gf.getGridFunction();
//  		return 
//  			(square(f[i][j+1] + f[i][j]) - square(f[i][j] + f[i][j-1]))
//  			/ (4.0*h[1]);
//  	}
//  
//  	inline Real FGx(const int i, const int j, 
//  			GridFunction& f,
//  			GridFunction& g)
//  	{
//  		return 
//  			/* d(uv)/dx , should be okay */
//  			((f[i][j+1] + f[i][j])*(g[i+1][j] + g[i][j]) 
//  			 - (f[i-1][j+1] + f[i-1][j])*(g[i][j] + g[i-1][j]))
//  			/ (4.0*h[0]);
//  	}
//  	inline Real FGy(const int i, const int j, 
//  			GridFunction& f,
//  			GridFunction& g)
//  	{
//  		return 
//  			((g[i+1][j] + g[i][j])*(f[i][j+1] + f[i][j]) 
//  			 - (g[i+1][j-1] + g[i][j-1])*(f[i][j] + f[i][j-1]))
//  			/ (4.0*h[1]);
//  	}
//  
//  	inline 
//  		std::function<Real(Stencil*, const int i, const int j, GridFunction& gf)>
//  	getStencilOperatorFunctionSingle(const StencilOperator so)
//  		{
//  			switch(so)
//  			{
//  				case StencilOperator::Fx_forward:
//  					return &Stencil::Fx_forward;
//  					break;
//  				case StencilOperator::Fy_forward:
//  					return &Stencil::Fy_forward;
//  					break;
//  				case StencilOperator::Fx_backward:
//  					return &Stencil::Fx_backward;
//  					break;
//  				case StencilOperator::Fy_backward:
//  					return &Stencil::Fy_backward;
//  					break;
//  				case StencilOperator::Fxx:
//  					return &Stencil::Fxx;
//  					break;
//  				case StencilOperator::Fyy:
//  					return &Stencil::Fyy;
//  					break;
//  				case StencilOperator::FFx:
//  					return &Stencil::FFx;
//  					break;
//  				case StencilOperator::FFy:
//  					return &Stencil::FFy;
//  					break;
//  				default: /* uuhh... */
//  					break;
//  			};
//  			return 0;
//  		}
//  	inline std::function<Real(Stencil*, const int i, const int j, GridFunction& gf, GridFunction& gf2)>
//  	getStencilOperatorFunctionDouble(const StencilOperator so)
//  		{
//  			switch(so)
//  			{
//  				case StencilOperator::FGx:
//  					return &Stencil::FGx;
//  					break;
//  				case StencilOperator::FGy:
//  					return &Stencil::FGy;
//  					break;
//  				default: /* uuhh... */
//  					break;
//  			};
//  			return 0;
//  		}
//  
//  public:
//  	Stencil(const Point h) : h(h) {};
//  	virtual ~Stencil();
//  
//  #define forall(F,S,B,E) for(int F=B[0];F<=E[0];F++)for(int S=B[1];S<=E[1];S++)
//  
//  	inline void applyStencilOperator(
//  		const StencilOperator stencil,
//  		const Index& gridreadbegin, 
//  		const Index& gridreadend,
//  		const Index& gridwritebegin,
//  		const Index& gridwriteend, 
//  		GridFunction& sourcegrid,
//  		GridFunction& targetgrid)
//  	{
//  		auto stenciloperator = 
//  			getStencilOperatorFunctionSingle(stencil);
//  		forall(i,j,gridreadbegin,gridreadend)
//  		{
//  			targetgrid(
//  					gridwritebegin[0]+i-gridreadbegin[0],
//  					gridwritebegin[1]+j-gridreadbegin[1])
//  				= stenciloperator(this,i,j,sourcegrid);
//  		}
//  	}
//  
//  	inline void applyStencilOperator(
//  		const StencilOperator stencil,
//  		const Index& gridreadbegin, 
//  		const Index& gridreadend,
//  		const Index& gridwritebegin,
//  		const Index& gridwriteend, 
//  		GridFunction& sourcegrid,
//  		GridFunction& sourcegrid2,
//  		GridFunction& targetgrid)
//  	{
//  		auto stenciloperator = 
//  			getStencilOperatorFunctionDouble(stencil);
//  		forall(i,j,gridreadbegin,gridreadend)
//  		{
//  			targetgrid(
//  					gridwritebegin[0]+i-gridreadbegin[0],
//  					gridwritebegin[1]+j-gridreadbegin[1])
//  				= stenciloperator(this,i,j,sourcegrid,sourcegrid2);
//  		}
//  	}
//  
//  	enum class StencilOperator {
//  		Fx_forward,
//  		Fy_forward,
//  		Fx_backward,
//  		Fy_backward,
//  		Fxx,
//  		Fyy,
//  		FFx,
//  		FFy,
//  		FGx,
//  		FGy
//  	};
//  };

#endif

