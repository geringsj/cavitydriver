#ifndef CavityRenderer_hpp
#define CavityRenderer_hpp

#include <fstream>
#include <sstream>
#include <vector>
#include <queue>

#include "GL/glew.h"
#include "GLFW/glfw3.h"

#include "glowl/glowl.h"
#include "AntTweakBar.h"

#include "Structs.hpp"
#include "CameraSystem.h"
#include "SimulationParameters.hpp"

class CavityRenderer
{
public:
	CavityRenderer();
	~CavityRenderer();

	/**
	 * Function to create the OpenGL context and creata a 
	 * window. Also initializes the camera system.
	 */
	bool initVis(
		unsigned int window_width, /**< Width of the window in pixels [>=0]. */
		unsigned int window_height, /**< Height of the window in pixels [>=0]. */
		SimulationParameters& sim_params
		);

	bool initBakeryVis(
		unsigned int window_width, /**< Width of the window in pixels [>=0]. */
		unsigned int window_height, /**< Height of the window in pixels [>=0]. */
		SimulationParameters& sim_params
		);

	/**
	 * The main render loop. Render a frame, swap buffer
	 * and poll events until the window gets closed.
	 */
	void paint();

	void loadFieldData();

	/**
	 * Create the three grids for pressure, u and v by 
	 * calling createSingleGrid for each range. As the
	 * data is buffered here as well call this function
	 * after calling init.
	 */
	bool createGrid(
		Range grid_size
		/**< Range that defines the grid. */
		);

	void setWindowSize(int width, int height)
		{ m_window_width=width; m_window_height = height; }

	void getijMax(int& p_iMax, int& p_jMax){ p_iMax = m_simparams.iMax; p_jMax = m_simparams.jMax; }

private:
	GLFWwindow* m_window; /**< Pointer to the window that glfw will use. */
	unsigned int m_window_width; /**< Store the window width in pixel. */
	unsigned int m_window_height; /**< Store the window height in pixel. */
	GLfloat m_window_background_colour[3];
	TwBar* bar;
	float m_zoom;

	SimulationParameters m_simparams;

	CameraSystem m_cam_sys; /**< The camera system. Stores the camera data and can perform translation and rotation. */

	Mesh m_grid;
	GLSLProgram m_grid_prgm; /**< Store the programm data for the grid rendering programm. */
	bool m_show_grid;

	Texture2D m_arrow;
	GLSLProgram m_arrow_prgm;
	Mesh m_arrow_quad;

	Mesh m_field_quad;
	GLSLProgram m_field_prgm;
	Texture2D m_pressure_tx;
	Texture2D m_velocity_tx;
	bool m_show_field;

	/** 
	 * Function to create the data and index array for a grid
	 * of the given Range.
	 */
	void createSingleGrid(
		Range innerRange,
		/**< Range that defines the grid. */
		std::vector<unsigned int>& index,
		/**< Reference to the index array. */
		std::vector<Gridvertex>& data
		/**< Reference to the data array. */
		);

	void addFloatParam(const char* name, const char* def, void* var,
		std::string mode = "RW",
		float min = std::numeric_limits<float>::min(),
		float max = std::numeric_limits<float>::max()
		);
	void addIntParam(const char* name, const char* def, void* var, 
		std::string mode = "RW",
		int min = std::numeric_limits<int>::min(),
		int max = std::numeric_limits<int>::max()
		);
	void addBoolParam(const char* name, const char* def, void* var, 
		std::string mode = "RW");
	void addVec3Param(const char* name, const char* def, void* var, 
		std::string mode = "RW");
	void addButtonParam(const char* name, const char* def, TwButtonCallback callback);
	void addStringParam(const char* name, const char* def, void* var, 
		std::string mode = "RW");
	void removeParam(const char* name);

	/** Upload field data of a specified timestep to texture ojects */
	void updateTextures(unsigned int timestep);

	/** Draw overlay grid */
	void drawGridOverlay();

	/** Draw active field data, i.e. velocity, pressure,.. */
	void drawField();

	/** Draw velocity field using line integral convolution */
	void drawLIC();

	/** Draw geometric shapes */
	void drawGeometry();

	void drawBoundaryCondition(Boundary::BoundaryPiece boundarypiece);
		/* BoundaryPiece holds all of the following... */
		//Boundary::Direction dir,
		//Boundary::Condition cond,
		//Boundary::Grid grid_type, 
		//Real condition_value, 
		//Range range

	void reloadSimParams(SimulationParameters& sim_params);

	/* Static GLFW callback functions - Primarily calls AntTweakBar functions */
	inline static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
	{ 
		/* get rid of stupid warnings: */
		if( 1 || window || mods)
		TwEventMouseButtonGLFW(button, action);
	}
	inline static void mousePositionCallback(GLFWwindow* window, double x, double y)
	{
		/* get rid of stupid warnings: */
		if( 1 || window)
		TwMouseMotion(int(x), int(y));
	}
	inline static void mouseWheelCallback(GLFWwindow* window, double x_offset, double y_offset)
	{
		/* get rid of stupid warnings: */
		if( 1 || window || x_offset)
		TwEventMouseWheelGLFW((int)y_offset);
	}
	inline static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		/* get rid of stupid warnings: */
		if(1 || window || scancode || mods)
		TwEventKeyGLFW(key, action);
	}
	inline static void charCallback(GLFWwindow* window, int codepoint)
	{
		/* get rid of stupid warnings: */
		if(1 || window)
		TwEventCharGLFW(codepoint, GLFW_PRESS);
	}
	inline static void windowResizeCallback(GLFWwindow* window, int width, int height)
	{
		CavityRenderer* painter = reinterpret_cast<CavityRenderer*>(glfwGetWindowUserPointer(window));
		painter->setWindowSize(width,height);
		TwWindowSize(width, height);
	}

	/** Read a shader source file */
	const std::string readShaderFile(const char* const path);

	bool createGLSLProgramms();

	/**
	* \brief Read a the header of a ppm image file. Courtesy to the computer vision lecture I attended.
	* \param filename Location of the image file
	* \param headerEndPos Out parameter, marks the point where the header of the ppm file ends
	* \param imgDimX Out parameter, containing the dimension of the image in X direction in pixels
	* \param imgDimY Out parameter, containing the dimension of the image in Y direction in pixels
	* \return Returns true if the ppm header was succesfully read, false otherwise
	*/
	bool readPpmHeader(const char* filename, unsigned long& headerEndPos, int& imgDimX, int& imgDimY);
	/**
	* \brief Read a the data of a ppm image file. Courtesy to the computer vision lecture I attended.
	* \param filename Location of the image file
	* \param imageData Pointer to the data buffer, that the image data will be written to
	* \param dataBegin Marks the location within the ppm file, where the data block begins
	* \param imgDimX Dimension of the image in X direction in pixels
	* \param imgDimY Dimension of the image in Y direction in pixels
	* \return Returns true if the ppm header was succesfully read, false otherwise
	*/
	bool readPpmData(const char* filename, char* imageData, unsigned long dataBegin, int imgDimX, int imgDimY);
};

#endif
