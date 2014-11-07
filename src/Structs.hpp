#ifndef Structs_h
#define Structs_h

#include "typedef.hpp"

struct Index3D
{
	Index3D() : m_i(0), m_j(0), m_k(0) {}
	Index3D(IndexType i, IndexType j, IndexType k) : m_i(i), m_j(j), m_k(k) {}

	IndexType m_i, m_j, m_k;

	IndexType operator[](uint index)
	{
		switch (index)
		{
		case 0:
			return m_i;
		case 1:
			return m_j;
		case 2:
			return m_k;
		default:
			exit(-1);
		}
	}
};

struct Point3D
{
	RealType m_x, m_y, m_z;

	Point3D() : m_x(0), m_y(0), m_z(0) {}
	Point3D(RealType x, RealType y, RealType z = 1)
		: m_x(x), m_y(y), m_z(z) {} 

	RealType m_x,m_y,m_z;

	RealType operator[](uint index)
	{
		switch (index)
		{
		case 0:
			return m_x;
		case 1:
			return m_y;
		case 2:
			return m_z;
		default:
			exit(-1);
		}
	}
};

struct Dimension3D
{
	Dimension3D() : m_x(0), m_y(0), m_z(0) {}
	Dimension3D(DimensionType x, DimensionType y, DimensionType z = 1)
		: m_x(x), m_y(y), m_z(z) {} 

	DimensionType m_x,m_y,m_z;

	DimensionType operator[](uint index)
	{
		switch (index)
		{
		case 0:
			return m_x;
		case 1:
			return m_y;
		case 2:
			return m_z;
		default:
			exit(-1);
		}
	}
};

struct Grid3D
{
	Grid3D(Dimension3D dimension)
		: m_u(dimension.m_x,dimension.m_y),
		m_v(dimension.m_x,dimension.m_y),
		m_w(dimension.m_x,dimension.m_y) {}

	GridFunction m_u;
	GridFunction m_v;
	GridFunction m_w;

	GridFunction& operator[](uint index)
	{
		switch (index)
		{
		case 1:
			return m_u;
		case 2:
			return m_v;
		case 3:
			return m_w;
		default:
			exit(-1); //not to sure about this
		}
	}
};

#endif