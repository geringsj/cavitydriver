
/** 
 * Header containing globally used definitions, types and structs. 
 * 
 * @file Structs.hpp
 * @author becherml, friesfn, geringsj
 * @date 11/2014
 *
 * We threw away (some of) the initially given typedefs and structs, 
 * and instead use our own set of helpful definitions.
 *
 * Most important: we define the makro DIMENSIONS as a compile-time global variable to
 * set the dimensionality of the problem. We intend to later expand to 3D.
 */


#ifndef Structs_hpp
#define Structs_hpp


/** Define in how many dimensions the simulation should be computed.
 * @def DIMENSIONS
 *
 * At the moment we only support 2D. We intend to later expand to 3D.
 */
#define DIMENSIONS 2


/** Define a c++11 compatible NULL pointer.
 * @def NULL
 */
#ifndef NULL
#define NULL nullptr
#endif


/** For convenience, define uint.
 * @var typedef unsigned int uint
 */
typedef unsigned int uint;


/** The Real type is (should be!) used for all floating point operations.
 * @var typedef double Real
 *
 * Hard-coding float or double seems to not be very popular.
 * Also, appending 'Type' to the name of a type looks very ugly and does not
 * improve readability of code, so we don't do this. 
 */
typedef double Real;


/** The Index struct allows easy handling and access of 2D (and 3D) signed integral indices.
 * Indices can be accessed via .i , .j ( .k ) or via [0], [1] ([2]), respectively.
 * Note: no array templates are used here.
 */
struct Index
{
	Index() : i(0), j(0), k(0) {}
	Index(int i, int j, int k = 0)
		: i(i), j(j), k(k) {} 

	int i; /**< Index for x dimension. */
	int j; /**< Index for y dimension. */
	int k; /**< Index for z dimension. Is set to 0 when not specified otherwise. */

	int& operator[](const uint index)
	{
		switch (index)
		{
		case 0:
			return i;
			break;
		case 1:
			return j;
			break;
		case 2:
			return k;
			break;
		default:
			return i;
			break;
		}
	}
	int operator[](const uint index) const
	{
		switch (index)
		{
		case 0:
			return i;
			break;
		case 1:
			return j;
			break;
		case 2:
			return k;
			break;
		default:
			return i;
			break;
		}
	}
};
/** 
 * For convenience.
 * @var typedef Index Dimension
 * @see Index
 *
 * When passing around dimensions of grids we want to pass a type that gives us
 * the look an feel of a dimension. 
 * Of course the underlying structure is simply what we implemented in the Index struct.
 */
typedef Index Dimension;

struct IndexRange {
	Index begin;
	Index end;

	IndexRange() : begin(), end() {};
	IndexRange(Index b, Index e) : begin(b), end(e) {};

	Index operator[](const uint index) const 
	{
		switch(index)
		{
			case 0:
				return begin;
				break;
			case 1:
				return end;
				break;
			default:
				return end;
				break;
		}
	}
};
typedef IndexRange Range;


/** The Point struct allows easy handling and access of 2D (and 3D) points.
 * Point dimensions can be accessed via .x , .y ( .z ) or via [0], [1] ([2]), respectively.o
 * Internally the Real
 * type is used as floating point data type of the x/y/z members.
 * Note: no array templates are used here.
 */
struct Point 
{
	Point() : x(0), y(0), z(0) {}
	Point(Real x, Real y, Real z = 0)
		: x(x), y(y), z(z) {} 

	Real x;
	Real y;
	Real z;

	Real& operator[](const uint index)
	{
		switch (index)
		{
		case 0:
			return x;
			break;
		case 1:
			return y;
			break;
		case 2:
			return z;
			break;
		default:
			return x;
			break;
		}
	}
	Real operator[](const uint index) const
	{
		switch (index)
		{
		case 0:
			return x;
			break;
		case 1:
			return y;
			break;
		case 2:
			return z;
			break;
		default:
			return x;
			break;
		}
	}
};
/** For convenience.
 * @var typedef Point Delta;
 * @see Point 
 *
 * When passing around floating point distances (deltas) we want to pass a type 
 * that gives us the look an feel of a delta. 
 * Of course the underlying structure is simply what we implemented in the Point struct.
 */
typedef Point Delta;


/**
 * We use this for the Checkerboard-Solving Pattern in SOR, 
 * and also in Communication for exchanging pressures.
 */
enum class Color {
	Red,
	Black
};

struct Gridvertex {
	Gridvertex() : x(0), y(0), z(0), w(0) {}
	Gridvertex(float x_coord, float y_coord, float z_coord, float w_coord) :
		x(x_coord), y(y_coord), z(z_coord), w(w_coord) {}
	float x, y, z, w;
};

struct VertexUV
{
	VertexUV() : m_x(0), m_y(0), m_z(0), m_u(0), m_v(0) {}
	VertexUV(float x, float y, float z, float u, float v)
		: m_x(x), m_y(y), m_z(z), m_u(u), m_v(v) {}

	float m_x,m_y,m_z;
	float m_u,m_v;
};

#endif
