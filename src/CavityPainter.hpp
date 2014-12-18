#ifndef CavityPainter_hpp
#define CavityPainter_hpp

#include <iostream>
#include <vector>

#include "glew.h"
#include "glfw3.h"

#include "glowl\glowl.h"

#include "Structs.hpp"

class CavityPainter
{
public:
	CavityPainter();
	~CavityPainter();

	void paint(unsigned int window_width, unsigned int window_height);

	bool createGrid(Range p, Range u, Range v);

private:
	enum Grid { P,U,V};
	Grid m_show_grid;
	Mesh m_p_grid, m_u_grid, m_v_grid;

	unsigned int m_window_width;
	unsigned int m_window_height;

	Gridvertex* createSingleGrid(Range innerRange, unsigned int& data_size);

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