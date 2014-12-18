#ifndef CavityPainter_hpp
#define CavityPainter_hpp

#include <fstream>
#include <sstream>
#include <vector>

#include "glew.h"
#include "glfw3.h"

#include "glowl/glowl.h"
#include "AntTweakBar.h"

#include "Structs.hpp"
#include "CameraSystem.h"

class CavityPainter
{
public:
	CavityPainter();
	~CavityPainter();

	/**
	 * Function to create the OpenGL context and creata a 
	 * window. Also initializes the camera system.
	 */
	bool init(
		unsigned int window_width,
		/**< Width of the window in pixels [>=0]. */
		unsigned int window_height
		/**< Height of the window in pixels [>=0]. */
		);
	/**
	 * The main render loop. Render a frame, swap buffer
	 * and poll events until the window gets closed.
	 */
	void paint();
	/**
	 * Create the three grids for pressure, u and v by 
	 * calling createSingleGrid for each range. As the
	 * data is buffered here as well call this function
	 * after calling init.
	 */
	bool createGrids(
		Range p,
		/**< Range that defines the pressure grid. */
		Range u,
		/**< Range that defines the u grid.*/
		Range v
		/**< Range that defines the v grid. */
		);

private:
	TwBar* bar;
	GLFWwindow* m_window; /**< Pointer to the window that glfw will use. */
	CameraSystem m_cam_sys; /**< The camera system. Stores the camera data and can perform translation and rotation. */
	enum Grid {
		P,
		/**< Render the pressure grid. */
		U,
		/**< Render the u grid. */
		V
		/**< Render the v grid */
	};
	Grid m_show_grid; /**< Based on the Grid enumerator render the associated grid. */
	Mesh m_p_grid; /**< Store the mesh data of the pressure grid. */
	Mesh m_u_grid; /**< Store the mesh data of the u grid. */
	Mesh m_v_grid; /**< Store the mesh data of the v grid. */
	GLSLProgram m_grid_prgm; /**< Store the programm data for the grid rendering programm. */
	unsigned int m_window_width; /**< Store the window width in pixel. */
	unsigned int m_window_height; /**< Store the window height in pixel. */
	/** 
	 * Function to create the data and index array for a grid
	 * of the given Range.
	 */
	void createSingleGrid(
		Range innerRange,
		/**< Range that defines the grid. */
		unsigned int& data_size, 
		/**< Refrence of the size of the data [>= 0]. */
		std::vector<unsigned int>& index,
		/**< Reference to the index array. */
		std::vector<Gridvertex>& data
		/**< Reference to the data array. */
		);

	/** Draw overlay grid */
	void drawGridOverlay();

	/** Draw active field data, i.e. velocity, pressure,.. */
	void drawField();

	/** Draw velocity field using line integral convolution */
	void drawLIC();

	/** Draw geometric shapes */
	void drawGeometry();

	/* Read a shader source file */
	const std::string readShaderFile(const char* const path);
};

#endif