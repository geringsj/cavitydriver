#ifndef CavityRenderer_hpp
#define CavityRenderer_hpp

#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <memory>
#include <random>
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
			right_vector(1.0f,0.0f,0.0f), translation(1.0f), rotation(1.0),
			m_aspect_ratio(16.0f/9.0f), m_field_of_view(60.0f/m_aspect_ratio) {}
		//! constructor creates a camera system with the given parameters
		CameraSystem(glm::vec3 p_cam_pos,glm::vec3 p_up_vector,glm::vec3 p_forward_vector,glm::vec3 p_right_vector)
			: cam_pos(p_cam_pos), up_vector(p_up_vector), forward_vector(p_forward_vector),
			right_vector(p_right_vector), translation(1.0), rotation(1.0),
			m_aspect_ratio(16.0f/9.0f), m_field_of_view(60.0f/m_aspect_ratio) {}
		~CameraSystem() {};
		//! rotation around (up,forward,right)-vector v with angle alpha (in degree)
		void Rotation(glm::vec3 v,float alpha)
		{
			rotation = glm::rotate(glm::mat4(1.0f),alpha,v);
			up_vector = glm::vec3(rotation * glm::vec4(up_vector,1.0f));
			forward_vector = glm::vec3(rotation * glm::vec4(forward_vector,1.0));
			right_vector =  glm::vec3(rotation * glm::vec4(right_vector,1.0));
		}
		//! translation with a step in (up,forward,right)-vector direction
		void Translation(glm::vec3 t, float step_size)
		{
			translation = glm::translate(glm::mat4(1.0f),step_size*t);
			cam_pos = glm::vec3(translation * glm::vec4(cam_pos,1.0f));
		}
		void zoom(float factor) { m_field_of_view = (60.0f/m_aspect_ratio) * factor; }
		//! return the view matrix
		glm::mat4 GetViewMatrix() { return glm::lookAt(cam_pos,cam_pos+forward_vector,up_vector); }
		glm::mat4 GetProjectionMatrix() { return glm::perspective(m_field_of_view*3.14f/180.0f,m_aspect_ratio, 0.1f, 100.0f); }
		//! returns the camera position
		glm::vec3 GetCamPos() { return cam_pos; }
		glm::vec3& accessCamPos() { return cam_pos; }
		//! returns the up-vector
		glm::vec3 GetUpVector() { return up_vector; }
		//! returns the forwards vector
		glm::vec3 GetForwardVector() { return forward_vector; }
		//! returns the right vector
		glm::vec3 GetRightVector() { return right_vector; }
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
		//! stores the translation
		glm::mat4 translation;
		//! stores the rotation
		glm::mat4 rotation;

		float m_aspect_ratio;
		float m_field_of_view;
	};

	struct Layer
	{
		bool m_show;
		GLSLProgram m_prgm;

		virtual bool createResources(SimulationParameters& simparams) = 0;
		virtual void draw(CameraSystem& camera) = 0;
	};

	struct FieldLayer : public Layer
	{
		GLSLProgram m_ibfvAdvection_prgm;
		GLSLProgram m_ibfvMerge_prgm;
		GLSLProgram m_dyeInjection_prgm;
		GLSLProgram m_fieldPicking_prgm;
		GLSLProgram m_streamline_prgm;

		std::shared_ptr<FramebufferObject> m_ibfv_fbo0;
		std::shared_ptr<FramebufferObject> m_ibfv_fbo1;
		Mesh m_field_quad;
		Mesh m_ibfv_grid;
		Mesh m_dye_blob;
		Mesh m_fullscreen_quad;
		std::vector<std::shared_ptr<Mesh>> m_streamlines;
		std::shared_ptr<Texture2D> m_field_tx;
		std::shared_ptr<Texture2D> m_ibfvBackground_tx;
		std::shared_ptr<Texture2D> m_dyeBlob_tx;

		// coordinates given in uv space of field
		std::vector<Point> m_dye_seedpoints;
		std::vector<Point> m_streamline_seedpoints;

		Dimension m_field_dimension;
		std::vector<std::vector<float>> m_field_data;
		std::vector<glm::vec3> m_field_max_values;
		std::vector<glm::vec3> m_field_min_values;
		int m_num_fields;
		int m_current_field;
		int m_display_mode;
		Real m_dx, m_dy;
		GLfloat m_stream_colour[3];

		bool m_show_streamlines = false;
		bool m_play_animation = false;
		double m_requested_frametime = 0.033;
		double m_time_tracker = 0.0;
		double m_elapsed_time;

		bool createResources(SimulationParameters& simparams);
		void draw(CameraSystem& camera);
		void drawFieldPicking(CameraSystem& camera);
		void setFieldData(std::string path);
		bool setFieldTexture(unsigned int requested_frame);
		void updateFieldTexture(double current_time);
		void addDyeSeedpoint(float x, float y);
		void clearDye();
		void interpolateUV(PointVertex p, float& u, float& v);
		void addStreamlineSeedpoint(float x, float y);
		void clearStreamline();
	};

	struct OverlayGridLayer : public Layer
	{
		GLfloat m_colour[3];
		Mesh m_grid;

		bool createResources(SimulationParameters& simparams);
		void draw(CameraSystem& camera);

		bool updateGridMesh(SimulationParameters& simparams);
	};
	
	struct BoundaryCellsLayer : public Layer
	{
		std::vector<Point> m_cell_positions;
		Mesh m_cell;
		std::shared_ptr<Texture2D> m_cell_tx;
		//std::shared_ptr<Texture2D> m_cell_positions_tx;

		bool createResources(SimulationParameters& simparams);
		void draw(CameraSystem& camera);

		void setCellPositions(SimulationParameters& simparams);
		bool updateCellMesh(SimulationParameters& simparams);
	};
	
	struct BoundaryGlyphLayer : public Layer
	{
		int m_glyph_mode;
		std::vector<std::pair<Point,Point>> m_velocity_glyphs;
		std::vector<std::pair<Point,Boundary::Condition>> m_pressure_glyphs;
		Mesh m_glyph;

		std::shared_ptr<Texture2D> m_velocity_glyph_tx;
		std::shared_ptr<Texture2D> m_pdlt_glyph_tx;
		std::shared_ptr<Texture2D> m_pnm_glyph_tx;

		bool createResources(SimulationParameters& simparams);
		void draw(CameraSystem& camera);

		void setGlyphs(SimulationParameters& simparams);
		bool updateGlyphMesh(SimulationParameters& simparams);
	};
	
	struct GeometryLayer : public Layer
	{
		bool createResources(SimulationParameters& simparams);
		void draw(CameraSystem& camera);
	};

	struct InterfaceLayer : public Layer
	{
		Mesh m_domainIndicators;

		bool createResources(SimulationParameters& simparams);
		void draw(CameraSystem& camera);

		bool updateResources(SimulationParameters& simparams);
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
		SimulationParameters& sim_params,
		std::string fields_filename
		);

	bool initBakeryVis(
		unsigned int window_width, /**< Width of the window in pixels [>=0]. */
		unsigned int window_height, /**< Height of the window in pixels [>=0]. */
		SimulationParameters& sim_params
		);

	/**
	 * The main render loop. Render a frame, swap buffer and poll events
	 * until the window gets closed.
	 */
	void paint();
	void paintBakery();

	void pushSimParams();

	void addDyeSeedpoint(int x, int y)
	{
		if(m_fieldPicking_fbo != nullptr)
		{
			GLfloat* data = new GLfloat[4];
			m_fieldPicking_fbo->bind();
			glReadBuffer(GL_COLOR_ATTACHMENT0);
			glReadPixels(x, y, 1, 1, GL_RGBA, GL_FLOAT, data);

			// check if click was inside field
			if(data[3]>0.5)
				m_field_layer.addDyeSeedpoint(data[0],data[1]);

			delete[] data;
		}
	}
	void addStreamline(int x, int y)
	{
		if(m_fieldPicking_fbo != nullptr)
		{
			GLfloat* data = new GLfloat[4];
			m_fieldPicking_fbo->bind();
			glReadBuffer(GL_COLOR_ATTACHMENT0);
			glReadPixels(x, y, 1, 1, GL_RGBA, GL_FLOAT, data);

			// check if click was inside field
			if(data[3]>0.5)
				m_field_layer.addStreamlineSeedpoint(data[0],data[1]);

			delete[] data;
		}
	}

	void clearDye()
		{ m_field_layer.clearDye(); }
	void clearStreamline()
		{ m_field_layer.clearStreamline(); }

	void setWindowSize(int width, int height)
	{
		m_window_width=width;
		m_window_height = height;

		m_activeCamera.setAspectRatio((float)m_window_width/(float)m_window_height);

		if(m_fieldPicking_fbo != nullptr)
			m_fieldPicking_fbo->resize(width,height);
	}

	void zoomCamera(float factor)
	{
		m_activeCamera.accessFieldOfView() -= factor;
	}

	void moveCamera(float dx, float dy, float dz, float speed)
	{
		m_activeCamera.Translation(glm::vec3(dx,dy,dz),speed);
	}

	CameraSystem& getCamera()
	{
		return m_activeCamera;
	}
	unsigned int getWindowWidth()
		{ return m_window_width; }
	unsigned int getWindowHeight()
		{ return m_window_height; }

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
	bool m_dye;
	bool m_stream;
	TwBar* bar;

	/* Individual layers of the visualization */
	FieldLayer m_field_layer;
	OverlayGridLayer m_overlayGrid_layer;
	BoundaryCellsLayer m_boundaryCells_layer;
	BoundaryGlyphLayer m_boundaryGlyph_layer;
	GeometryLayer m_geometry_layer;
	InterfaceLayer m_interface_layer;

	GLSLProgram m_postProc_prgm;
	std::shared_ptr<FramebufferObject> m_fieldPicking_fbo;
	Mesh m_screen_quad;

	/** Local simparams copy */
	SimulationParameters m_simparams;

	/** FIFO for incomming simulation parameters */
	MTQueue<SimulationParameters>& m_inbox;
	/** FIFO for outgoing simulation parameters */
	MTQueue<SimulationParameters>& m_outbox;

	/** Active camera */
	CameraSystem m_activeCamera;
	CameraSystem m_perspectiveCamera_tgt;
	CameraSystem m_topCamera_tgt;

	//int m_nmbr_boundary_piece, m_max_boundary_piece;
	//Boundary::Direction m_direction_enum;
	//Boundary::Condition m_condition_enum;
	//Boundary::Grid m_grid_enum;
	//Real m_condition_value;
	//Range m_range;
	//int m_i_begin, m_i_end, m_j_begin, m_j_end;
	//std::vector<Boundary::BoundaryPiece> m_boundary_conditions;
	//bool m_modify_cond;
	
	/** Merge visualization layers */
	void postProcessing();

	/***********************************************************************
	 * Initialize graphics resources, e.g. shader programs, textures, meshes
	 **********************************************************************/

	bool createGLSLPrograms();

	bool createMeshes();

	//void addBoundaryPieceToBar(std::string mode = "RW");
	//void reloadSimParams(SimulationParameters& sim_params);

	/****************************
	 * AntTweakBar helper methods
	 ***************************/

	void initBakeryTweakBar();
	void initPainterTweakBar();

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
	static double last_mouse_x;
	static double last_mouse_y;

	inline static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
	{
		/* get rid of stupid warnings: */
		if( 1 || window || mods)
		if( !TwEventMouseButtonGLFW(button, action))
		{

			if(button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
			{
				CavityRenderer* renderer = reinterpret_cast<CavityRenderer*>(glfwGetWindowUserPointer(window));

				double xPos, yPos;
				glfwGetCursorPos(window, &xPos, &yPos);

				yPos = renderer->getWindowHeight() - yPos;

				renderer->addDyeSeedpoint(xPos,yPos);
				renderer->addStreamline(xPos,yPos);
			}
		}
	}
	inline static void mousePositionCallback(GLFWwindow* window, double x, double y)
	{
		/* get rid of stupid warnings: */
		if( 1 || window)
		if( !TwMouseMotion(int(x), int(y)))
		{
			CavityRenderer* renderer = reinterpret_cast<CavityRenderer*>(glfwGetWindowUserPointer(window));

			double dx = last_mouse_x - x;
			double dy = last_mouse_y - y;
			
			last_mouse_x = x;
			last_mouse_y = y;

			if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_2) == GLFW_PRESS)
			{
				//TODO decent zoom depented movement speed
				float speed = std::pow(renderer->getCamera().getFieldOfView()*0.001f,2.0f);
				renderer->moveCamera((float)dx,(float)-dy,0.0f,speed);
			}
		}
	}
	inline static void mouseWheelCallback(GLFWwindow* window, double x_offset, double y_offset)
	{
		/* get rid of stupid warnings: */
		if( 1 || window || x_offset)
		if( !TwEventMouseWheelGLFW((int)y_offset) )
		{
			CavityRenderer* renderer = reinterpret_cast<CavityRenderer*>(glfwGetWindowUserPointer(window));
			renderer->zoomCamera(y_offset);
		}
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

};

namespace Renderer
{
namespace IO
{
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

	bool readPfmHeader(const char* filename, unsigned long& headerEndPos, int& imgDimX, int& imgDimY);
	bool readPfmData(const char* filename, float* imageData, unsigned long dataBegin);
}
}

#endif
