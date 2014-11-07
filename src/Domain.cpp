#include "Domain.hpp"


Domain::Domain(Dimension dimension, Point delta) : 
		m_dimension(dimension),
		m_delta(delta),
		m_p(dimension[0]+2,dimension[1]+2), 

		m_velocity(dimension),
		m_preliminary_velocities_FGH(dimension),

		m_rhs(dimension[0],dimension[1]),
		m_gx(0.0), m_gy(0.0), m_gz(0.0)
{
}

Domain::Domain(Dimension dimension, Point delta, Point extForces) :
		m_dimension(dimension),
		m_delta(delta),
		m_p(dimension[0],dimension[1]), 

		m_velocity(dimension),
		m_preliminary_velocities_FGH(dimension),

		m_rhs(dimension[0],dimension[1]),
		m_gx(extForces[0]), m_gy(extForces[2]), m_gz(extForces[3])
{
}

Domain::~Domain()
{
}
