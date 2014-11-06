#include "Domain.hpp"


Domain::Domain(Dimension3D dimension, Point3D delta)
	: m_p(dimension[0],dimension[1]), m_velocity(dimension),
	m_F(dimension[0],dimension[1]), m_G(dimension[0],dimension[1]), m_H(dimension[0],dimension[1]),
	m_gx(dimension[0],dimension[1]), m_gy(dimension[0],dimension[1]), m_gz(dimension[0],dimension[1]),
	m_rhs(dimension[0],dimension[1])
{
}

Domain::~Domain()
{
}