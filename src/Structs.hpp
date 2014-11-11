//! Here we define our sane, well-designed and useful types and structs!
/*!
 * @author becherml, friesfn, geringsj
 * @date 2014
 */

#ifndef Structs_hpp
#define Structs_hpp

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
