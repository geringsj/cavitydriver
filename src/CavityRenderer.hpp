#ifndef CavityRenderer_hpp
#define CavityRenderer_hpp

#include <fstream>
#include <sstream>
#include <vector>

#include "glew.h"
#include "glfw3.h"

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
	bool init(
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

	int* getijMax(){ int ij_max[2]; ij_max[0] = m_iMax; ij_max[1] = m_jMax; return ij_max; }

private:
	GLFWwindow* m_window; /**< Pointer to the window that glfw will use. */
	unsigned int m_window_width; /**< Store the window width in pixel. */
	unsigned int m_window_height; /**< Store the window height in pixel. */
	GLfloat m_window_background_colour[3];
	TwBar* bar;
	float m_zoom;

	/**
	 * I'm very sorry but we can't do this:
	 * SimulationParameters m_sim_params;
	 * because SimulationParameters doesn't have
	 * a default constructor (so far?). So we will
	 * create a copy of each simparam value in here.
	 * Just because it's late and i'm hungry.
	 */
	float m_alpha, m_deltaT, m_eps, m_gx, m_gy, m_KarmanAngle, 
		m_KarmanObjectWidth,m_pi, m_re, m_tau, m_tEnd, 
		m_ui, m_vi, m_xLength, m_yLength, m_omg;
	union{
		float m_tDeltaWriteVTK;
		float m_deltaVec;
	};
	int m_iterMax;
	union{
		int m_iMax;
		int m_xCells;
	};
	union{
		int m_jMax;
		int m_yCells;
	};
	std::string m_name;
	/**
	 * And this also doesn't work because AntTweakBar,
	 * apparently hates pointers and crashes for (no?)
	 * good reason.
	 */
	//SimulationParameters* m_sim_params;


	CameraSystem m_cam_sys; /**< The camera system. Stores the camera data and can perform translation and rotation. */

	Mesh m_grid;
	GLSLProgram m_grid_prgm; /**< Store the programm data for the grid rendering programm. */
	bool m_show_grid;

	Mesh m_field_quad;
	GLSLProgram m_field_prgm;
	Texture2D m_pressure_tx;
	Texture2D m_velocity_tx;
	bool m_show_field;
	bool m_grid_resize;

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
		float min = std::numeric_limits<float>::min(), 
		float max = std::numeric_limits<float>::max());
	void addIntParam(const char* name, const char* def, void* var, 
		int min = std::numeric_limits<int>::min(),
		int max = std::numeric_limits<int>::max());
	void addBoolParam(const char* name, const char* def, void* var);
	void addVec3Param(const char* name, const char* def, void* var);
	void addButtonParam(const char* name, const char* def, TwButtonCallback callback);
	void addStringParam(const char* name, const char* def, void* var);
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

	/* Static GLFW callback functions - Primarily calls AntTweakBar functions */
	inline static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
	{TwEventMouseButtonGLFW(button, action);}
	inline static void mousePositionCallback(GLFWwindow* window, double x, double y)
	{TwMouseMotion(int(x), int(y));}
	inline static void mouseWheelCallback(GLFWwindow* window, double x_offset, double y_offset)
	{TwEventMouseWheelGLFW(y_offset);}
	inline static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
	{TwEventKeyGLFW(key, action);}
	inline static void charCallback(GLFWwindow* window, int codepoint)
	{TwEventCharGLFW(codepoint, GLFW_PRESS);}
	inline static void windowResizeCallback(GLFWwindow* window, int width, int height)
	{
		CavityRenderer* painter = reinterpret_cast<CavityRenderer*>(glfwGetWindowUserPointer(window));
		painter->setWindowSize(width,height);
		TwWindowSize(width, height);
	}

	/** Read a shader source file */
	const std::string readShaderFile(const char* const path);
};

#endif