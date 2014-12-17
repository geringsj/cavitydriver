#ifndef CavityPainter_hpp
#define CavityPainter_hpp

#include "glew.h"
#include "glfw3.h"

#include "glowl\glowl.h"

class CavityPainter
{
public:
	CavityPainter();
	~CavityPainter();

	void paint(unsigned int window_width, unsigned int window_height);

private:
	unsigned int m_window_width;
	unsigned int m_window_height;

	unsigned int m_grid_size_x;
	unsigned int m_grid_size_y;
	Mesh m_grid_mesh;

	/** Draw overlay grid */
	void drawGridOverlay();

	/** Draw active field data, i.e. velocity, pressure,.. */
	void drawField();

	/** Draw velocity field using line integral convolution */
	void drawLIC();

	/** Draw geometric shapes */
	void drawGeometry();
};

#endif