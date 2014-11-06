//! Typedefs for the apllication 
/*!
 * @author diehlpk
 * @date 2012
 */

#ifndef TYPEDEF
#define TYPEDEF

#include "template.hpp"

typedef unsigned int uint;

/*! Creates a type name for RealType */
typedef double RealType;

/*! Creates a type name for IndexType */
typedef int IndexType;

/*! Creates a type name for MultiIndexType */
typedef Array < IndexType, 2 > MultiIndexType;

/*! Creates a type name for GridFunctionType */
typedef RealType **GridFunctionType;

/*! Creates a type name for StencilType */
typedef GridFunctionType StencilType;

/*! Creates a type name for PointType */
typedef Array < RealType, 2 > PointType;

/**
 * Additional typedefs
 */

#include "GridFunction.hpp"

typedef uint DimensionType;

typedef Array < IndexType, 3 > Index3D;
typedef Array < RealType, 3 > Point3D;

#endif
