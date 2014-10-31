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
	void setFxxStencil();
	/**
	*
	*/
	void setFxyStencil();
	/**
	*
	*/
	void setFyxStencil();
	/**
	*
	*/
	void setFyyStencil();
private:
	StencilType stencil;
	int stencilwidth;
	const PointType h;
};

#endif
