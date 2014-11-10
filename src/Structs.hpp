//! Used sane typedefs 
/*!
 * @author becherml, friesfn, geringsj
 * @date 2014
 */

#ifndef Structs_h
#define Structs_h

// struct Dimension3D
// {
// 	Dimension3D() : m_x(0), m_y(0), m_z(0) {}
// 	Dimension3D(DimensionType x, DimensionType y, DimensionType z = 1)
// 		: m_x(x), m_y(y), m_z(z) {} 
// 
// 	DimensionType m_x,m_y,m_z;
// 
// 	DimensionType operator[](uint index)
// 	{
// 		switch (index)
// 		{
// 		case 0:
// 			return m_x;
// 		case 1:
// 			return m_y;
// 		case 2:
// 			return m_z;
// 		default:
// 			return 0;
// 		}
// 	}
// };


#define NULL nullptr

#define DIMENSIONS 2

typedef unsigned int uint;

typedef double Real;

struct Index
{
	Index() : i(0), j(0), k(1) {}
	Index(int i, int j, int k = 1)
		: i(i), j(j), k(k) {} 

	int i,j,k;

	int& operator[](const uint index)
	{
		switch (index)
		{
		case 0:
			return i;
		case 1:
			return j;
		case 2:
			return k;
		default:
			return i;
		}
	}
	int operator[](const uint index) const
	{
		switch (index)
		{
		case 0:
			return i;
		case 1:
			return j;
		case 2:
			return k;
		default:
			return i;
		}
	}
};
typedef Index Dimension;

struct Point 
{
	Point() : x(0), y(0), z(0) {}
	Point(Real x, Real y, Real z = 0)
		: x(x), y(y), z(z) {} 

	Real x,y,z;

	Real& operator[](const uint index)
	{
		switch (index)
		{
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		default:
			return x;
		}
	}
	Real operator[](const uint index) const
	{
		switch (index)
		{
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		default:
			return x;
		}
	}
};
typedef Point Delta;


#endif
