#ifndef CavityRenderer_hpp
#define CavityRenderer_hpp

#include <fstream>
#include <sstream>
#include <vector>
#include <typeinfo>

#include "GL/glew.h"
#include "GLFW/glfw3.h"

#include "glowl/glowl.h"
#include "AntTweakBar.h"

#include "Structs.hpp"
#include "SimulationParameters.hpp"
#include "MTQueue.hpp"

class CavityRenderer
{
private:
	class CameraSystem
	{
	public:
		//! constructor creates a camera system with the camera in position (0,0,0)
		//! and up-vector (0,1,0) forward-vector (0,0,1) and right-vector (1,0,0)
		CameraSystem()
			: cam_pos(), up_vector(0.0f,1.0f,0.0f), forward_vector(0.0f,0.0f,1.0f),
			right_vector(1.0f,0.0f,0.0f), center(), translation(1.0), rotation(1.0),
			m_aspect_ratio(16.0/9.0), m_field_of_view(60.0/m_aspect_ratio) {}
		//! constructor creates a camera system with the given parameters
		CameraSystem(glm::vec3 p_cam_pos,glm::vec3 p_up_vector,glm::vec3 p_forward_vector,glm::vec3 p_right_vector)
			: cam_pos(p_cam_pos), up_vector(p_up_vector), forward_vector(p_forward_vector),
			right_vector(p_right_vector), center(), translation(1.0), rotation(1.0),
			m_aspect_ratio(16.0/9.0), m_field_of_view(60.0/m_aspect_ratio) {}
		~CameraSystem() {};
		//! rotation around (up,forward,right)-vector v with angle alpha (in degree)
		void Rotation(glm::vec3 v,float alpha)
		{
			rotation = glm::rotate(glm::mat4(1.0f),alpha,v);
			up_vector = glm::vec3(rotation * glm::vec4(up_vector,1.0f));
			forward_vector = glm::vec3(rotation * glm::vec4(forward_vector,1.0));
			right_vector =  glm::vec3(rotation * glm::vec4(right_vector,1.0));
			center = glm::vec3(rotation * glm::vec4(center,1.0));
		}
		//! translation with a step in (up,forward,right)-vector direction
		void Translation(glm::vec3 t, float step_size)
		{
			translation = glm::translate(glm::mat4(1.0f),step_size*t);
			cam_pos = glm::vec3(translation * glm::vec4(cam_pos,1.0f));
			center = glm::vec3(translation * glm::vec4(center,1.0));
		}
		void zoom(float factor) { m_field_of_view = (60.0/m_aspect_ratio) * factor; }
		//! return the view matrix
		glm::mat4 GetViewMatrix() { return glm::lookAt(cam_pos,center,up_vector); }
		//! returns the camera position
		glm::vec3 GetCamPos() { return cam_pos; }
		//! returns the up-vector
		glm::vec3 GetUpVector() { return up_vector; }
		//! returns the forwards vector
		glm::vec3 GetForwardVector() { return forward_vector; }
		//! returns the right vector
		glm::vec3 GetRightVector() { return right_vector; }
		//! returns the center vector
		glm::vec3 GetCenterVector() { return center; }
		float getFieldOfView() { return m_field_of_view; }
		float& accessFieldOfView() { return m_field_of_view; }
		float getAspectRatio() { return m_aspect_ratio; }
		void setAspectRatio(float aspect_ratio) { m_aspect_ratio = aspect_ratio; }
	private:
		//! stores the camera position
		glm::vec3 cam_pos;
		//! stores the up-vector
		glm::vec3 up_vector;
		//! stores the forward vector
		glm::vec3 forward_vector;
		//! stores the right vector
		glm::vec3 right_vector;
		//! stores the center vector, i.e. the positin we look at in world space
		glm::vec3 center;
		//! stores the translation
		glm::mat4 translation;
		//! stores the rotation
		glm::mat4 rotation;

		float m_aspect_ratio;
		float m_field_of_view;
	};

public:
	CavityRenderer(MTQueue<SimulationParameters>& inbox, MTQueue<SimulationParameters>& outbox);
	~CavityRenderer();

	/**
	 * Function to create the OpenGL context and creata a 
	 * window. Also initializes the camera system.
	 */
	bool initPainterVis(
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

	void pushSimParams();


	void setWindowSize(int width, int height)
	{
		m_window_width=width;
		m_window_height = height;

		m_cam_sys.setAspectRatio((float)m_window_width/(float)m_window_height);
	}


	//void setMaxBoundaryPiece(int max_boundary_piece) { m_max_boundary_piece = max_boundary_piece; }
	//int getMaxBoundaryPiece() { return m_max_boundary_piece; }
	//int getBoundaryPieceIndex() { return m_nmbr_boundary_piece; }
	//void showBoundaryPiece(unsigned int index);
	//void addBoundaryPiece(Boundary::BoundaryPiece piece) { m_boundary_conditions.push_back(piece); }
	//void deleteBoundaryPiece(unsigned int index) { m_boundary_conditions.erase(m_boundary_conditions.begin() + index); }
	//void getBoundaryPieceParams(Boundary::Direction& dir, Boundary::Condition& cond,
	//	Boundary::Grid& grid, Real& value, int& i_begin, int& i_end, int& j_begin, int& j_end)
	//{
	//	dir = m_direction_enum; cond = m_condition_enum; grid = m_grid_enum; 
	//	value = m_condition_value; i_begin = m_i_begin; i_end = m_i_end;
	//	j_begin = m_j_begin; j_end = m_j_end;
	//}
	//void modifyBoundaryPieceParams(unsigned int index);

private:
	GLFWwindow* m_window; /**< Pointer to the window that glfw will use. */
	unsigned int m_window_width; /**< Store the window width in pixel. */
	unsigned int m_window_height; /**< Store the window height in pixel. */
	GLfloat m_window_background_colour[3];
	TwBar* bar;

	/* FBOs for visualization layers */
	FramebufferObject m_field_fbo;
	FramebufferObject m_grid_fbo;
	FramebufferObject m_boundary_gylphs_fbo;
	FramebufferObject m_boundary_cells_fbo;
	FramebufferObject m_geometry_fbo;

	/* local simparams copy */
	SimulationParameters m_simparams;

	MTQueue<SimulationParameters>& m_inbox;
	MTQueue<SimulationParameters>& m_outbox;

	/** Active camera */
	CameraSystem m_cam_sys;

	/* Grid related variables*/
	Mesh m_grid;
	GLSLProgram m_grid_prgm; /**< Store the programm data for the grid rendering programm. */
	GLfloat m_grid_colour[3];

	/* Boundary glyph variables*/
	Texture2D m_arrow;
	GLSLProgram m_arrow_prgm;
	Mesh m_arrow_quad;

	/* Boundary cells variables*/
	Texture2D m_boundary_cell_tx;
	Texture2D m_boundary_cell_positions_tx;
	GLSLProgram m_boundary_cell_prgm;
	Mesh m_boundary_cell;

	/* Field related variables */
	Mesh m_field_quad;
	GLSLProgram m_field_prgm;
	Texture2D m_pressure_tx;
	Texture2D m_velocity_tx;

	/*****************
	 *Visibility flags
	 *****************/
	bool m_show_field;
	bool m_show_grid;
	bool m_show_boundary_glyphs;
	bool m_show_boundary_cells;
	bool m_show_geometry;


	//int m_nmbr_boundary_piece, m_max_boundary_piece;
	//Boundary::Direction m_direction_enum;
	//Boundary::Condition m_condition_enum;
	//Boundary::Grid m_grid_enum;
	//Real m_condition_value;
	//Range m_range;
	//int m_i_begin, m_i_end, m_j_begin, m_j_end;
	//std::vector<Boundary::BoundaryPiece> m_boundary_conditions;
	//bool m_modify_cond;


	/** Upload field data of a specified timestep to texture ojects */
	void updateTextures(unsigned int timestep);

	/***********************************************************************
	 * Initialize graphics resources, e.g. shader programs, textures, meshes
	 **********************************************************************/

	bool createGLSLPrograms();
	/** Create overlay grid (really just a grid made of lines) */
	bool createOverlayGrid();

	bool createTextures();

	bool createFramebuffers();


	/****************************************************
	 * Drawing method for the single visualization layers
	 ***************************************************/

	/** Draw overlay grid */
	void drawOverlayGrid();

	/** Draw active field data or visualization, i.e. velocity, pressure, LIC,.. etc. */
	void drawField();

	/** Draw geometric shapes */
	void drawGeometry();

	void drawBoundaryGlyphs();

	void drawBoundaryCells();

	void postProcessing();



	//void drawBoundaryCondition(Boundary::BoundaryPiece boundarypiece);
		/* BoundaryPiece holds all of the following... */
		//Boundary::Direction dir,
		//Boundary::Condition cond,
		//Boundary::Grid grid_type, 
		//Real condition_value, 
		//Range range


	//void addBoundaryPieceToBar(std::string mode = "RW");

	//void reloadSimParams(SimulationParameters& sim_params);


	/****************************
	 * AntTweakBar helper methods
	 ***************************/
	void addFloatParam(const char* name, const char* def, void* var,
		std::string mode = "RW",
		float min = std::numeric_limits<float>::min(),
		float max = std::numeric_limits<float>::max()
		);
	void addDoubleParam(const char* name, const char* def, void* var,
		std::string mode = "RW",
		double min = std::numeric_limits<double>::min(),
		double max = std::numeric_limits<double>::max()
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
	void addEnumParam(const char* name, const char* def, void* var, TwEnumVal* _enum,
		int size, std::string mode = "RW");
	void modifyIntParam(const char* name, int min, int max);
	void removeParam(const char* name);


	/************************************************************************
	 * Static GLFW callback functions - Primarily calls AntTweakBar functions
	 ***********************************************************************/
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
		CavityRenderer* renderer = reinterpret_cast<CavityRenderer*>(glfwGetWindowUserPointer(window));
		renderer->setWindowSize(width,height);
		TwWindowSize(width, height);
	}


	/***************************************
	 * Utility methods e.g. for file loading
	 **************************************/

	/** Read a shader source file */
	const std::string readShaderFile(const char* const path);

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
