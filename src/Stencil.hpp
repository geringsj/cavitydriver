//! This class implements the stencil
/*!
 * @author ...
 * @date 2014
 */


#ifndef Stencil_HPP_
#define Stencil_HPP_

#include "typedef.hpp"
#include "GridFunction.hpp"


class Stencil
{
public:
	/** ctor */
	Stencil(int stencilwidth, const PointType h);
	/** dtor */
	~Stencil();
	/**
	 *
	 */
	void ApplyStencilOperator(const MultiIndexType& gridbegin,
		const MultiIndexType& gridreadend,
		const MultiIndexType& gridwritebegin,
		const MultiIndexType& gridwriteend, 
		const GridFunction sourcegridfunction, GridFunction imagegridfunction);
	/**
	*
	*/
	void setFxStencil();
	/**
	*
	*/
	void setFyStencil();
	/**
	 *
	 */
	void setFFxStencil();
	void setFFyStencil();
	void setFxxStencil();
	/**
	*
	*/
	void setFyyStencil();
private:
	RealType rawmemory[9];
	StencilType stencil;
	const int stencilwidth;
	const PointType h;
};

#endif
