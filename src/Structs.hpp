#ifndef Structs_h
#define Structs_h

#include "typedef.hpp"

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
			return 0;
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