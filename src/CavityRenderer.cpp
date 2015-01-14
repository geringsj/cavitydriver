#include "CavityRenderer.hpp"

#include <iostream>

void CavityRenderer::FieldLayer::draw(CameraSystem& camera, float* background_colour)
{

}

void CavityRenderer::OverlayGridLayer::draw(CameraSystem& camera, float* background_colour)
{
	m_fbo->bind();
	glClearColor(
		background_colour[0],
		background_colour[1],
		background_colour[2], 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, m_fbo->getWidth(), m_fbo->getHeight());

	if(m_show)
	{
		m_prgm.use();
		glm::mat4 proj_mat = glm::perspective(camera.getFieldOfView(), camera.getAspectRatio(), 0.1f, 100.0f);
		glm::mat4 model_mat = glm::mat4(1.0f);
		glm::mat4 view_mat = camera.GetViewMatrix();
		glm::mat4 mvp_mat = proj_mat * view_mat * model_mat;
		m_prgm.setUniform("mvp_matrix", mvp_mat);
		m_prgm.setUniform("colour", glm::vec3(m_colour[0],m_colour[1],m_colour[2]));
		
		m_grid.draw();
	}
}

bool CavityRenderer::OverlayGridLayer::updateGridMesh(SimulationParameters& simparams)
{
	std::vector<unsigned int> grid_index_array;
	std::vector<Gridvertex> grid_vertex_array;

	float x_length = (float)simparams.xLength / (float)simparams.iMax;
	float y_length = (float)simparams.yLength / (float)simparams.jMax;

	float bottom_left_i = 0.0f;
	float bottom_left_j = 0.0f;

	float bottom_right_i = bottom_left_i + (simparams.xCells + 2) * x_length;
	float bottom_right_j = bottom_left_j;

	float top_right_i = bottom_right_i;
	float top_right_j = bottom_right_j + (simparams.yCells + 2) * y_length;

	float top_left_i = bottom_left_i;
	float top_left_j = bottom_left_j + (simparams.yCells + 2) * y_length;

	float right_shift = top_right_i / 2.0f;
	float up_shift = top_right_j / 2.0f;

	grid_vertex_array.push_back(Gridvertex(bottom_left_i, bottom_left_j, -1.0f, 1.0f));
	grid_index_array.push_back(0);
	grid_vertex_array.push_back(Gridvertex(bottom_right_i, bottom_right_j, -1.0f, 1.0f));
	grid_index_array.push_back(1);
	grid_vertex_array.push_back(Gridvertex(top_right_i, top_right_j, -1.0f, 1.0f));
	grid_index_array.push_back(2);
	grid_vertex_array.push_back(Gridvertex(top_left_i, top_left_j, -1.0f, 1.0f));
	grid_index_array.push_back(3);

	grid_index_array.push_back(0);
	grid_index_array.push_back(3);
	grid_index_array.push_back(1);
	grid_index_array.push_back(2);

	//LEFT RIGHT
	float left = bottom_left_i;
	float right = bottom_right_i;
	unsigned int index_value = 4;
	for (float j = bottom_left_j; j <= top_left_j; j=j+y_length)
	{
		grid_vertex_array.push_back(Gridvertex(left, j, -1.0f, 1.0f));
		grid_index_array.push_back(index_value); index_value++;
		grid_vertex_array.push_back(Gridvertex(right, j, -1.0f, 1.0f));
		grid_index_array.push_back(index_value); index_value++;
	}
	//TOP BOTTOM
	float bottom = bottom_right_j;
	float top = top_right_j;
	for (float i = bottom_left_i; i <= bottom_right_i; i=i+x_length)
	{
		grid_vertex_array.push_back(Gridvertex(i, bottom, -1.0f, 1.0f));
		grid_index_array.push_back(index_value); index_value++;
		grid_vertex_array.push_back(Gridvertex(i, top, -1.0f, 1.0f));
		grid_index_array.push_back(index_value); index_value++;
	}
	
	if(!m_grid.bufferDataFromArray(grid_vertex_array.data(),grid_index_array.data(),
		(GLsizei)(grid_vertex_array.size()*sizeof(Gridvertex)),(GLsizei)(grid_index_array.size()*sizeof(unsigned int)),GL_LINES))
		return false;
	m_grid.setVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(Gridvertex), 0);
}

void CavityRenderer::BoundaryCellsLayer::draw(CameraSystem& camera, float* background_colour)
{
	m_fbo->bind();
	glClearColor(
		background_colour[0],
		background_colour[1],
		background_colour[2], 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, m_fbo->getWidth(), m_fbo->getHeight());

	if(m_show)
	{
		m_prgm.use();
		
		glEnable(GL_TEXTURE_2D);
		glActiveTexture(GL_TEXTURE0);
		m_prgm.setUniform("boundary_tx2D",0);
		m_cell_tx->bindTexture();

		
		for(auto& cell : m_cell_positions)
		{
			glm::mat4 proj_mat = glm::perspective(camera.getFieldOfView(), camera.getAspectRatio(), 0.1f, 100.0f);
			glm::mat4 model_mat = glm::translate(glm::mat4(1.0),glm::vec3(cell.x,cell.y,0.0));
			glm::mat4 view_mat = camera.GetViewMatrix();
			glm::mat4 mvp_mat = proj_mat * view_mat * model_mat;
			m_prgm.setUniform("mvp_matrix", mvp_mat);
			
			m_cell.draw();
		}

	}
}

void CavityRenderer::BoundaryCellsLayer::setCellPositions(SimulationParameters& simparams)
{
	m_cell_positions.clear();

	for(auto& boundary_piece : simparams.boundary_conditions)
	{
		Range range(boundary_piece.range);
		float x_length = (float)simparams.xLength / (float)simparams.iMax;
		float y_length = (float)simparams.yLength / (float)simparams.jMax;

		for_range(i, j, range)
		{
			float pos[] = { (float)(i) * x_length/1.0 + x_length / 2.0f, (float)(j) * y_length/1.0 + y_length/2.0f };
			switch(boundary_piece.direction)
			{
			case Boundary::Direction::Up:
				pos[1] += y_length;
				break;
			case Boundary::Direction::Down:
				pos[1] -= y_length;
				break;
			case Boundary::Direction::Left:
				pos[0] -= x_length;
				break;
			case Boundary::Direction::Right:
				pos[0] += x_length;
				break;
			}

			m_cell_positions.push_back(Point(pos[0],pos[1],0.0));
		}
	}
}

bool CavityRenderer::BoundaryCellsLayer::updateCellMesh(SimulationParameters& simparams)
{
	// Create mesh for boundary cells
	float dx = (float)simparams.xLength / (float)simparams.iMax;
	float dy = (float)simparams.yLength / (float)simparams.jMax;

	dx /= 2.0;
	dy /= 2.0;

	std::array< VertexUV, 4 > cell_vertex_array = {{ VertexUV(-dx,-dy,-1.0,0.0,0.0),
											VertexUV(-dx,dy,-1.0,0.0,1.0),
											VertexUV(dx,dy,-1.0,1.0,1.0),
											VertexUV(dx,-dy,-1.0,1.0,0.0) }};

	std::array< GLuint, 6 > cell_index_array = {{ 0,2,1,2,0,3 }};

	if(!(m_cell.bufferDataFromArray(cell_vertex_array.data(),cell_index_array.data(),sizeof(VertexUV)*4,sizeof(GLuint)*6,GL_TRIANGLES))) return false;
	m_cell.setVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VertexUV),0);
	m_cell.setVertexAttribPointer(1,2,GL_FLOAT,GL_FALSE,sizeof(VertexUV),(GLvoid*) (sizeof(float)*3));
}

void CavityRenderer::BoundaryGlyphLayer::draw(CameraSystem& camera, float* background_colour)
{
	m_fbo->bind();
	glClearColor(
		background_colour[0],
		background_colour[1],
		background_colour[2], 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, m_fbo->getWidth(), m_fbo->getHeight());

	if(m_show)
	{
		m_prgm.use();

		if(m_glyph_mode == 0)
		{
			for(auto& velocity_glyph : m_velocity_glyphs)
			{
				//if(boundary_piece.gridtype == Boundary::Grid::P &&  boundary_piece.condition==Boundary::Condition::INFLOW)
				//	m_boundary_pdlt_glyph_tx->bindTexture();
				//
				//if(boundary_piece.gridtype == Boundary::Grid::P &&  boundary_piece.condition==Boundary::Condition::OUTFLOW)
				//	m_boundary_pnm_glyph_tx->bindTexture();

				glm::mat4 proj_mat = glm::perspective(camera.getFieldOfView(), camera.getAspectRatio(), 0.1f, 100.0f);
				glm::mat4 model_mat = glm::translate(glm::mat4(1.0),glm::vec3(velocity_glyph.first.x,velocity_glyph.first.y,velocity_glyph.first.z));
				glm::mat4 view_mat = camera.GetViewMatrix();
				glm::mat4 mvp_mat = proj_mat * view_mat * model_mat;
				m_prgm.setUniform("mvp_matrix", mvp_mat);
				
				m_glyph.draw();
			}
		}
		else if(m_glyph_mode == 1)
		{
			for(auto& pressure_glyph : m_pressure_glyphs)
			{
				if( pressure_glyph.second==Boundary::Condition::INFLOW )
					m_pdlt_glyph_tx->bindTexture();

				if( pressure_glyph.second==Boundary::Condition::OUTFLOW )
					m_pnm_glyph_tx->bindTexture();

				glm::mat4 proj_mat = glm::perspective(camera.getFieldOfView(), camera.getAspectRatio(), 0.1f, 100.0f);
				glm::mat4 model_mat = glm::translate(glm::mat4(1.0),glm::vec3(pressure_glyph.first.x,pressure_glyph.first.y,pressure_glyph.first.z));
				glm::mat4 view_mat = camera.GetViewMatrix();
				glm::mat4 mvp_mat = proj_mat * view_mat * model_mat;
				m_prgm.setUniform("mvp_matrix", mvp_mat);
				
				m_glyph.draw();
			}
		}
	}
}

void CavityRenderer::BoundaryGlyphLayer::setGlyphs(SimulationParameters& simparams)
{
	m_velocity_glyphs.clear();
	m_pressure_glyphs.clear();

	for(auto& boundary_piece : simparams.boundary_conditions)
	{
		Range range(boundary_piece.range);
		float x_length = (float)simparams.xLength / (float)simparams.iMax;
		float y_length = (float)simparams.yLength / (float)simparams.jMax;

		for_range(i, j, range)
		{
			float pos[] = { (float)(i) * x_length/1.0 + x_length / 2.0f, (float)(j) * y_length/1.0 + y_length/2.0f };
			switch(boundary_piece.direction)
			{
			case Boundary::Direction::Up:
				pos[1] += y_length;
				break;
			case Boundary::Direction::Down:
				pos[1] -= y_length;
				break;
			case Boundary::Direction::Left:
				pos[0] -= x_length;
				break;
			case Boundary::Direction::Right:
				pos[0] += x_length;
				break;
			}

			if( boundary_piece.gridtype==Boundary::Grid::P )
				m_pressure_glyphs.push_back(std::pair<Point,Boundary::Condition>(
				Point(pos[0],pos[1],0.0),
				boundary_piece.condition));

			//if( boundary_piece.gridtype==Boundary::Grid::U || boundary_piece.gridtype==Boundary::Grid::V )
			//	m_velocity_glyphs.push_back(std::pair<Point,Point>(
			//	Point(pos[0],pos[1],0.0),
			//	Point(boundary_piece.)));
		}
	}
}

bool CavityRenderer::BoundaryGlyphLayer::updateGlyphMesh(SimulationParameters& simparams)
{
	float dx = (float)simparams.xLength / (float)simparams.iMax;
	float dy = (float)simparams.yLength / (float)simparams.jMax;

	dx /= 2.0;
	dy /= 2.0;

	std::array< VertexUV, 4 > glyph_vertex_array = {{ VertexUV(-dx,-dy,-1.0,0.0,0.0),
											VertexUV(-dx,dy,-1.0,0.0,1.0),
											VertexUV(dx,dy,-1.0,1.0,1.0),
											VertexUV(dx,-dy,-1.0,1.0,0.0) }};

	std::array< GLuint, 6 > glyph_index_array = {{ 0,2,1,2,0,3 }};

	if(!(m_glyph.bufferDataFromArray(glyph_vertex_array.data(),glyph_index_array.data(),sizeof(VertexUV)*4,sizeof(GLuint)*6,GL_TRIANGLES))) return false;
	m_glyph.setVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VertexUV),0);
	m_glyph.setVertexAttribPointer(1,2,GL_FLOAT,GL_FALSE,sizeof(VertexUV),(GLvoid*) (sizeof(float)*3));

}

void CavityRenderer::GeometryLayer::draw(CameraSystem& camera, float* background_colour)
{

}

/* Static tweak bar callback definitions. Implementation see further down. */
void TW_CALL Bake(void* clientData);

CavityRenderer::CavityRenderer(MTQueue<SimulationParameters>& inbox, MTQueue<SimulationParameters>& outbox)
	: m_inbox(inbox), m_outbox(outbox)
{
}

CavityRenderer::~CavityRenderer()
{
}

bool CavityRenderer::initPainterVis(unsigned int window_width, unsigned int window_height, SimulationParameters& sim_params)
{
	m_window_width = window_width;
	m_window_height = window_height;
	m_window_background_colour[0] = 0.2f;
	m_window_background_colour[1] = 0.2f;
	m_window_background_colour[2] = 0.2f;

	m_simparams = sim_params;

	m_overlayGrid_layer.m_show = true;

	/* Initialize the library */
	if (!glfwInit()) return false;

	/* Create a windowed mode window and its OpenGL context */
	glfwWindowHint(GLFW_SAMPLES, 8);
	m_window = glfwCreateWindow(window_width, window_height, "Cavity", NULL, NULL);
	if (!m_window)
	{
		//std::cout<<"Couldn't create glfw window."<<std::endl;
		glfwTerminate();
		return false;
	}
	glfwMakeContextCurrent(m_window);
	glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_NORMAL); // can be GLFW_CURSOR_HIDDEN

	// Initialize AntTweakBar
	TwInit(TW_OPENGL_CORE, NULL); // TwInit(TW_OPENGL, NULL);

	// Create a tweak bar
	bar = TwNewBar("Cavity-TweakBar");
	TwWindowSize(window_width, window_height);
	// TwDefine(" GLOBAL help='This example shows how to integrate AntTweakBar with GLFW and OpenGL.' "); // Message added to the help bar.
	// Add 'bgColor' to 'bar': it is a modifable variable of type TW_TYPE_COLOR3F (3 floats color)
	TwAddVarRW(bar, "m_window_background_colour", TW_TYPE_COLOR3F, &m_window_background_colour, " label='Background color' ");
	//addFloatParam("m_zoom", " label='Zoom' ", &m_zoom, "RW",0.0f, 999.0f);
	addBoolParam("m_show_grid", " label='Show grid' ", &m_overlayGrid_layer.m_show);
	TwAddSeparator(bar, "SimulationParameters", " label='SimulationParameters' ");

	Real test; float __float; double __double;
	const char* _double = typeid(__double).name();
	const char* _float = typeid(__float).name();
	if (strcmp(typeid(test).name(), _double) == 0)
	{
		addDoubleParam("m_alpha", " step=0.1 label='alpha' ", &m_simparams.alpha, "RO");
		addDoubleParam("m_deltaT", " step=0.1 label='deltaT' ", &m_simparams.deltaT, "RO");
		addDoubleParam("m_deltaVec", " step=0.1 label='deltaVec' ", &m_simparams.deltaVec, "RO");
		addDoubleParam("m_eps", " step=0.001 label='eps' ", &m_simparams.eps, "RO");
		addDoubleParam("m_gx", " step=0.1 label='gx' ", &m_simparams.gx, "RO");
		addDoubleParam("m_gy", " step=0.1 label='gy' ", &m_simparams.gy, "RO");
		addDoubleParam("m_KarmanAngle", " step=0.1 label='KarmanAngle' ", &m_simparams.KarmanAngle, "RO");
		addDoubleParam("m_KarmanObjectWidth", " step=0.1 label='KarmanObjectWidth' ", &m_simparams.KarmanObjectWidth, "RO");
		addDoubleParam("m_pi", " step=0.1 label='pi' ", &m_simparams.pi, "RO");
		addDoubleParam("m_re", " step=0.1 label='re' ", &m_simparams.re, "RO");
		addDoubleParam("m_tau", " step=0.1 label='tau' ", &m_simparams.tau, "RO");
		addDoubleParam("m_tEnd", " step=0.1 label='tEnd' ", &m_simparams.tEnd, "RO");
		addDoubleParam("m_ui", " step=0.1 label='ui' ", &m_simparams.ui, "RO");
		addDoubleParam("m_vi", " step=0.1 label='vi' ", &m_simparams.vi, "RO");
		addDoubleParam("m_xLength", " step=0.1 label='xLength' ", &m_simparams.xLength, "RO");
		addDoubleParam("m_yLength", " step=0.1 label='yLength' ", &m_simparams.yLength, "RO");
		addDoubleParam("m_omg", " step=0.1 label='omega' ", &m_simparams.omg, "RO");
	}
	if (strcmp(typeid(test).name(), _float) == 0)
	{
		addFloatParam("m_alpha", " step=0.1 label='alpha' ", &m_simparams.alpha, "RO");
		addFloatParam("m_deltaT", " step=0.1 label='deltaT' ", &m_simparams.deltaT, "RO");
		addFloatParam("m_deltaVec", " step=0.1 label='deltaVec' ", &m_simparams.deltaVec, "RO");
		addFloatParam("m_eps", " step=0.001 label='eps' ", &m_simparams.eps, "RO");
		addFloatParam("m_gx", " step=0.1 label='gx' ", &m_simparams.gx, "RO");
		addFloatParam("m_gy", " step=0.1 label='gy' ", &m_simparams.gy, "RO");
		addFloatParam("m_KarmanAngle", " step=0.1 label='KarmanAngle' ", &m_simparams.KarmanAngle, "RO");
		addFloatParam("m_KarmanObjectWidth", " step=0.1 label='KarmanObjectWidth' ", &m_simparams.KarmanObjectWidth, "RO");
		addFloatParam("m_pi", " step=0.1 label='pi' ", &m_simparams.pi, "RO");
		addFloatParam("m_re", " step=0.1 label='re' ", &m_simparams.re, "RO");
		addFloatParam("m_tau", " step=0.1 label='tau' ", &m_simparams.tau, "RO");
		addFloatParam("m_tEnd", " step=0.1 label='tEnd' ", &m_simparams.tEnd, "RO");
		addFloatParam("m_ui", " step=0.1 label='ui' ", &m_simparams.ui, "RO");
		addFloatParam("m_vi", " step=0.1 label='vi' ", &m_simparams.vi, "RO");
		addFloatParam("m_xLength", " step=0.1 label='xLength' ", &m_simparams.xLength, "RO");
		addFloatParam("m_yLength", " step=0.1 label='yLength' ", &m_simparams.yLength, "RO");
		addFloatParam("m_omg", " step=0.1 label='omega' ", &m_simparams.omg, "RO");
	}
	addIntParam("m_iterMax", " label='iterMax' ", &m_simparams.iterMax, "RO");
	addIntParam("m_iMax", " label='iMax' ", &m_simparams.iMax, "RO");
	addIntParam("m_jMax", " label='jMax' ", &m_simparams.jMax, "RO");
	addIntParam("m_xCells", " label='xCells' ", &m_simparams.xCells, "RO");
	addIntParam("m_yCells", " label='yCells' ", &m_simparams.yCells, "RO");
	addStringParam("m_name", " label='name' ", &m_simparams.name, "RO");

	TwAddSeparator(bar, "BoundaryConditions", " label='BoundaryConditions' ");
	//for (auto b : sim_params.boundary_conditions)
	//	m_boundary_conditions.push_back(b);
	//m_max_boundary_piece = (int)sim_params.boundary_conditions.size();
	//addIntParam("m_nmbr_boundary_piece", " label='Show boundary piece: ' ", &m_nmbr_boundary_piece, "RW", 0, m_max_boundary_piece);
	//addBoundaryPieceToBar("RO");

	// Set GLFW event callbacks
	glfwSetWindowUserPointer(m_window,this);
	glfwSetWindowSizeCallback(m_window, (GLFWwindowposfun)windowResizeCallback);
	glfwSetMouseButtonCallback(m_window, (GLFWmousebuttonfun)mouseButtonCallback);
	glfwSetCursorPosCallback(m_window, (GLFWcursorposfun)mousePositionCallback);
	glfwSetScrollCallback(m_window, (GLFWscrollfun)mouseWheelCallback);
	glfwSetKeyCallback(m_window, (GLFWkeyfun)keyCallback);
	glfwSetCharCallback(m_window, (GLFWcharfun)charCallback);

	/*	Initialize glew */
	//glewExperimental = GL_TRUE;
	GLenum error = glewInit();
	if (GLEW_OK != error)
	{
		std::cout << "-----\n"
			<< "The time is out of joint - O cursed spite,\n"
			<< "That ever I was born to set it right!\n"
			<< "-----\n"
			<< "Error: " << glewGetErrorString(error);
		return false;
	}
	/* Apparently glewInit() causes a GL ERROR 1280, so let's just catch that... */
	glGetError();


	m_cam_sys = CameraSystem(
		glm::vec3(0.0f, 0.0f, 1.0),
		glm::vec3(0.0f, 1.0f, 0.0f),
		glm::vec3(0.0f, 0.0f, -1.0f),
		glm::vec3(1.0f, 0.0f, 0.0f));

	createGLSLPrograms();

	return true;
}

bool CavityRenderer::initBakeryVis(unsigned int window_width, unsigned int window_height, SimulationParameters& sim_params)
{
	m_window_width = window_width;
	m_window_height = window_height;
	m_window_background_colour[0] = 0.2f; m_window_background_colour[1] = 0.2f; m_window_background_colour[2] = 0.2f;

	m_cam_sys = CameraSystem(
		glm::vec3(0.0f, 0.0f, 1.0),
		glm::vec3(0.0f, 1.0f, 0.0f),
		glm::vec3(0.0f, 0.0f, -1.0f),
		glm::vec3(1.0f, 0.0f, 0.0f));
	m_cam_sys.setAspectRatio((float)window_width/(float)window_height);
	m_cam_sys.accessFieldOfView() = (60.0 * ((float)window_height/(float)window_width));

	m_overlayGrid_layer.m_show = true;
	m_boundaryCells_layer.m_show = true;

	m_simparams = sim_params;

	/* Initialize the library */
	if (!glfwInit()) return false;

	/* Create a windowed mode window and its OpenGL context */
	glfwWindowHint(GLFW_SAMPLES, 8);
	m_window = glfwCreateWindow(window_width, window_height, "CavityBaker", NULL, NULL);
	if (!m_window)
	{
		//std::cout<<"Couldn't create glfw window."<<std::endl;
		glfwTerminate();
		return false;
	}
	glfwMakeContextCurrent(m_window);
	glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_NORMAL); // can be GLFW_CURSOR_HIDDEN

	// Initialize AntTweakBar
	TwInit(TW_OPENGL_CORE, NULL); // TwInit(TW_OPENGL, NULL);

	// Create a tweak bar
	bar = TwNewBar("CavityBaker-Settings");
	TwWindowSize(window_width, window_height);
	addIntParam("m_window_width", "label='Window width' group='Window' ", &m_window_width, "RO");
	addIntParam("m_window_height", "label='Window height' group='Window' ", &m_window_height, "RO");
	TwAddVarRW(bar, "m_window_background_colour", TW_TYPE_COLOR3F, &m_window_background_colour, " label='Background color' group='Window' ");
	addFloatParam("m_fieldOfView", " step=0.1 label='Field of View' group='Camera' ", &m_cam_sys.accessFieldOfView(), "RW", 1.0f, 180.0f);
	addFloatParam("x", " step=0.01 label='X postion' group='Camera' ", &m_cam_sys.accessCamPos().x, "RW", 0.0f, (float)m_simparams.xLength);
	addFloatParam("y", " step=0.01 label='Y position' group='Camera' ", &m_cam_sys.accessCamPos().y, "RW", 0.0f, (float)m_simparams.yLength);
	addBoolParam("m_show_grid", " label='Show grid' group='Grid' ", &m_overlayGrid_layer.m_show);
	TwAddVarRW(bar, "m_grid_colour", TW_TYPE_COLOR3F, &m_overlayGrid_layer.m_colour, " label='Grid color' group='Grid' ");

	addBoolParam("m_show_boundary_cells", " label='Show boundary cells' group='Boundary' ", &m_boundaryCells_layer.m_show);
	addBoolParam("m_show_boundary_glyphs", " label='Show boundary glyphs' group='Boundary' ", &m_boundaryGlyph_layer.m_show);
	addIntParam("m_boundary_glyph_mode", " label='Glyph display mode' group='Boundary' ", &m_boundaryGlyph_layer.m_glyph_mode, "RW", 0, 1);


	
	addIntParam("useComplexGeometry", " label='Scenario' group='Simulation Parameters' ", &m_simparams.useComplexGeometry, "RW", 0, 4);

	Real test; float __float; double __double;
	const char* _double = typeid(__double).name();
	const char* _float = typeid(__float).name();
	if (strcmp(typeid(test).name(), _double) == 0)
	{
		addDoubleParam("alpha", " step=0.1 label='alpha' group='Simulation Parameters' ", &m_simparams.alpha);
		addDoubleParam("deltaT", " step=0.1 label='deltaT' group='Simulation Parameters' ", &m_simparams.deltaT);
		addDoubleParam("deltaVec", " step=0.1 label='deltaVec' group='Simulation Parameters' ", &m_simparams.deltaVec);
		addDoubleParam("eps", " step=0.001 label='eps' group='Simulation Parameters' ", &m_simparams.eps);
		addDoubleParam("gx", " step=0.1 label='gx' group='Simulation Parameters' ", &m_simparams.gx);
		addDoubleParam("gy", " step=0.1 label='gy' group='Simulation Parameters' ", &m_simparams.gy);
		addDoubleParam("KarmanAngle", " step=0.1 label='KarmanAngle' group='Simulation Parameters' ", &m_simparams.KarmanAngle);
		addDoubleParam("KarmanObjectWidth", " step=0.1 label='KarmanObjectWidth' group='Simulation Parameters' ", &m_simparams.KarmanObjectWidth);
		addDoubleParam("pi", " step=0.1 label='pi' group='Simulation Parameters' ", &m_simparams.pi);
		addDoubleParam("re", " step=0.1 label='re' group='Simulation Parameters' ", &m_simparams.re);
		addDoubleParam("tau", " step=0.1 label='tau' group='Simulation Parameters' ", &m_simparams.tau);
		addDoubleParam("tEnd", " step=0.1 label='tEnd' group='Simulation Parameters' ", &m_simparams.tEnd);
		addDoubleParam("ui", " step=0.1 label='ui' group='Simulation Parameters' ", &m_simparams.ui);
		addDoubleParam("vi", " step=0.1 label='vi' group='Simulation Parameters' ", &m_simparams.vi);
		addDoubleParam("xLength", " step=0.1 label='xLength' group='Simulation Parameters' ", &m_simparams.xLength);
		addDoubleParam("yLength", " step=0.1 label='yLength' group='Simulation Parameters' ", &m_simparams.yLength);
		addDoubleParam("omg", " step=0.1 label='omega' group='Simulation Parameters' ", &m_simparams.omg);
	}
	if (strcmp(typeid(test).name(), _float) == 0)
	{
		addFloatParam("alpha", " step=0.1 label='alpha' ", &m_simparams.alpha);
		addFloatParam("deltaT", " step=0.1 label='deltaT' ", &m_simparams.deltaT);
		addFloatParam("deltaVec", " step=0.1 label='deltaVec' ", &m_simparams.deltaVec);
		addFloatParam("eps", " step=0.001 label='eps' ", &m_simparams.eps);
		addFloatParam("gx", " step=0.1 label='gx' ", &m_simparams.gx);
		addFloatParam("gy", " step=0.1 label='gy' ", &m_simparams.gy);
		addFloatParam("KarmanAngle", " step=0.1 label='KarmanAngle' ", &m_simparams.KarmanAngle);
		addFloatParam("KarmanObjectWidth", " step=0.1 label='KarmanObjectWidth' ", &m_simparams.KarmanObjectWidth);
		addFloatParam("pi", " step=0.1 label='pi' ", &m_simparams.pi);
		addFloatParam("re", " step=0.1 label='re' ", &m_simparams.re);
		addFloatParam("tau", " step=0.1 label='tau' ", &m_simparams.tau);
		addFloatParam("tEnd", " step=0.1 label='tEnd' ", &m_simparams.tEnd);
		addFloatParam("ui", " step=0.1 label='ui' ", &m_simparams.ui);
		addFloatParam("vi", " step=0.1 label='vi' ", &m_simparams.vi);
		addFloatParam("xLength", " step=0.1 label='xLength' ", &m_simparams.xLength);
		addFloatParam("yLength", " step=0.1 label='yLength' ", &m_simparams.yLength);
		addFloatParam("omg", " step=0.1 label='omega' ", &m_simparams.omg);
	}
	addIntParam("iterMax", " label='iterMax' group='Simulation Parameters' ", &m_simparams.iterMax);
	addIntParam("iMax", " label='iMax' group='Simulation Parameters' ", &m_simparams.iMax);
	addIntParam("jMax", " label='jMax' group='Simulation Parameters' ", &m_simparams.jMax);
	addIntParam("xCells", " label='xCells' group='Simulation Parameters' ", &m_simparams.xCells);
	addIntParam("yCells", " label='yCells' group='Simulation Parameters' ", &m_simparams.yCells);
	addStringParam("name", " label='name' group='Simulation Parameters' ", &m_simparams.name);

	//TwAddSeparator(bar, "BoundaryConditions", " label='BoundaryConditions' ");
	//for (auto b : sim_params.boundary_conditions)
	//	m_boundary_conditions.push_back(b);
	//m_max_boundary_piece = (int)sim_params.boundary_conditions.size();
	//addBoolParam("m_modify_cond", " label='Modify boundary piece' ", &m_modify_cond);
	//addIntParam("m_nmbr_boundary_piece", " label='Show boundary piece: ' ", &m_nmbr_boundary_piece, "RW", 0, m_max_boundary_piece);
	//addBoundaryPieceToBar("RW");

	//addButtonParam("m_boundarypiece", " label='Add boundary condition' ", BoundaryPiece);
	//addButtonParam("m_boundarypiece_mod", " label='Modify boundary condition' ", ModifyBoundaryPiece);
	//addButtonParam("m_boundarypiece_del", " label='Delete boundary condition' ", RemoveBoundaryPiece);

	addButtonParam("m_bake", " label='bake parameters' group='Simulation Parameters' ", Bake);

	// Set GLFW event callbacks
	glfwSetWindowUserPointer(m_window, this);
	glfwSetWindowSizeCallback(m_window, (GLFWwindowposfun)windowResizeCallback);
	glfwSetMouseButtonCallback(m_window, (GLFWmousebuttonfun)mouseButtonCallback);
	glfwSetCursorPosCallback(m_window, (GLFWcursorposfun)mousePositionCallback);
	glfwSetScrollCallback(m_window, (GLFWscrollfun)mouseWheelCallback);
	glfwSetKeyCallback(m_window, (GLFWkeyfun)keyCallback);
	glfwSetCharCallback(m_window, (GLFWcharfun)charCallback);

	/*	Initialize glew */
	//glewExperimental = GL_TRUE;
	GLenum error = glewInit();
	if (GLEW_OK != error)
	{
		std::cout << "-----\n"
			<< "The time is out of joint - O cursed spite,\n"
			<< "That ever I was born to set it right!\n"
			<< "-----\n"
			<< "Error: " << glewGetErrorString(error);
		return false;
	}
	/* Apparently glewInit() causes a GL ERROR 1280, so let's just catch that... */
	glGetError();

	// Create resources
	if(!createGLSLPrograms()) { return false; }
	if(!createMeshes()) { return false; }
	if(!createTextures()) { return false; }
	if(!createFramebuffers()) { return false; }

	// Init layers where necessary
	m_boundaryCells_layer.setCellPositions(m_simparams);
	m_boundaryGlyph_layer.setGlyphs(m_simparams);

	return true;
}

bool CavityRenderer::createGLSLPrograms()
{
	/* Arrow texture programm */
	m_boundaryGlyph_layer.m_prgm.init();
	
	std::string arrow_vertex = readShaderFile("./shader/arrowVertex.glsl");
	if (!m_boundaryGlyph_layer.m_prgm.compileShaderFromString(&arrow_vertex, GL_VERTEX_SHADER)) { std::cout << m_boundaryGlyph_layer.m_prgm.getLog(); return false; };
	
	std::string arrow_fragment = readShaderFile("./shader/arrowFragment.glsl");
	if (!m_boundaryGlyph_layer.m_prgm.compileShaderFromString(&arrow_fragment, GL_FRAGMENT_SHADER)) { std::cout << m_boundaryGlyph_layer.m_prgm.getLog(); return false; };
	
	m_boundaryGlyph_layer.m_prgm.bindAttribLocation(0, "v_position");
	m_boundaryGlyph_layer.m_prgm.bindAttribLocation(1, "v_uv");
	
	m_boundaryGlyph_layer.m_prgm.link();

	/* Grid programm */
	m_overlayGrid_layer.m_prgm.init();
	
	std::string grid_vertex = readShaderFile("./shader/gridVertex.glsl");
	if (!m_overlayGrid_layer.m_prgm.compileShaderFromString(&grid_vertex, GL_VERTEX_SHADER)) { std::cout << m_overlayGrid_layer.m_prgm.getLog(); return false; };
	
	std::string grid_fragment = readShaderFile("./shader/gridFragment.glsl");
	if (!m_overlayGrid_layer.m_prgm.compileShaderFromString(&grid_fragment, GL_FRAGMENT_SHADER)) { std::cout << m_overlayGrid_layer.m_prgm.getLog(); return false; };
	
	m_overlayGrid_layer.m_prgm.bindAttribLocation(0, "in_position");
	
	m_overlayGrid_layer.m_prgm.link();

	/* Boundary cells program */
	m_boundaryCells_layer.m_prgm.init();

	std::string boundary_vertex_shdr = readShaderFile("./shader/boundaryCellVertex.glsl");
	if (!m_boundaryCells_layer.m_prgm.compileShaderFromString(&boundary_vertex_shdr, GL_VERTEX_SHADER)) { std::cout << m_boundaryCells_layer.m_prgm.getLog(); return false; };

	std::string boundary_fragment_shdr = readShaderFile("./shader/boundaryCellFragment.glsl");
	if (!m_boundaryCells_layer.m_prgm.compileShaderFromString(&boundary_fragment_shdr, GL_FRAGMENT_SHADER)) { std::cout << m_boundaryCells_layer.m_prgm.getLog(); return false; };

	m_boundaryCells_layer.m_prgm.bindAttribLocation(0, "in_position");
	m_boundaryCells_layer.m_prgm.bindAttribLocation(1, "in_uv");

	m_boundaryCells_layer.m_prgm.link();

	/* Post processing program */
	m_postProc_prgm.init();

	std::string postProc_vertex_shdr = readShaderFile("./shader/postProcVertex.glsl");
	if (!m_postProc_prgm.compileShaderFromString(&postProc_vertex_shdr, GL_VERTEX_SHADER)) { std::cout << m_postProc_prgm.getLog(); return false; };

	std::string postProc_fragment_shdr = readShaderFile("./shader/postProcFragment.glsl");
	if (!m_postProc_prgm.compileShaderFromString(&postProc_fragment_shdr, GL_FRAGMENT_SHADER)) { std::cout << m_postProc_prgm.getLog(); return false; };

	m_postProc_prgm.bindAttribLocation(0, "v_position");
	m_postProc_prgm.bindAttribLocation(1, "v_uv");

	if (!m_postProc_prgm.link()) { std::cout << m_postProc_prgm.getLog(); return false; };

	return true; /* return with great success */
}

bool CavityRenderer::createMeshes()
{
	std::vector<unsigned int> grid_index_array;
	std::vector<Gridvertex> grid_vertex_array;

	float x_length = (float)m_simparams.xLength / (float)m_simparams.iMax;
	float y_length = (float)m_simparams.yLength / (float)m_simparams.jMax;

	float bottom_left_i = 0.0f;
	float bottom_left_j = 0.0f;

	float bottom_right_i = bottom_left_i + (m_simparams.xCells + 2) * x_length;
	float bottom_right_j = bottom_left_j;

	float top_right_i = bottom_right_i;
	float top_right_j = bottom_right_j + (m_simparams.yCells + 2) * y_length;

	float top_left_i = bottom_left_i;
	float top_left_j = bottom_left_j + (m_simparams.yCells + 2) * y_length;

	float right_shift = top_right_i / 2.0f;
	float up_shift = top_right_j / 2.0f;

	//why is this in here ?
	//m_cam_sys.Translation(m_cam_sys.GetRightVector(), 0.0f - m_cam_sys.GetCamPos().x);
	m_cam_sys.Translation(m_cam_sys.GetRightVector(), right_shift);
	//m_cam_sys.Translation(m_cam_sys.GetUpVector(), 0.0f - m_cam_sys.GetCamPos().y);
	m_cam_sys.Translation(m_cam_sys.GetUpVector(), up_shift);

	grid_vertex_array.push_back(Gridvertex(bottom_left_i, bottom_left_j, -1.0f, 1.0f));
	grid_index_array.push_back(0);
	grid_vertex_array.push_back(Gridvertex(bottom_right_i, bottom_right_j, -1.0f, 1.0f));
	grid_index_array.push_back(1);
	grid_vertex_array.push_back(Gridvertex(top_right_i, top_right_j, -1.0f, 1.0f));
	grid_index_array.push_back(2);
	grid_vertex_array.push_back(Gridvertex(top_left_i, top_left_j, -1.0f, 1.0f));
	grid_index_array.push_back(3);

	grid_index_array.push_back(0);
	grid_index_array.push_back(3);
	grid_index_array.push_back(1);
	grid_index_array.push_back(2);

	//LEFT RIGHT
	float left = bottom_left_i;
	float right = bottom_right_i;
	unsigned int index_value = 4;
	for (float j = bottom_left_j; j <= top_left_j; j=j+y_length)
	{
		grid_vertex_array.push_back(Gridvertex(left, j, -1.0f, 1.0f));
		grid_index_array.push_back(index_value); index_value++;
		grid_vertex_array.push_back(Gridvertex(right, j, -1.0f, 1.0f));
		grid_index_array.push_back(index_value); index_value++;
	}
	//TOP BOTTOM
	float bottom = bottom_right_j;
	float top = top_right_j;
	for (float i = bottom_left_i; i <= bottom_right_i; i=i+x_length)
	{
		grid_vertex_array.push_back(Gridvertex(i, bottom, -1.0f, 1.0f));
		grid_index_array.push_back(index_value); index_value++;
		grid_vertex_array.push_back(Gridvertex(i, top, -1.0f, 1.0f));
		grid_index_array.push_back(index_value); index_value++;
	}
	
	if(!m_overlayGrid_layer.m_grid.bufferDataFromArray(grid_vertex_array.data(),grid_index_array.data(),
		(GLsizei)(grid_vertex_array.size()*sizeof(Gridvertex)),(GLsizei)(grid_index_array.size()*sizeof(unsigned int)),GL_LINES))
		return false;
	m_overlayGrid_layer.m_grid.setVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(Gridvertex), 0);


	// Create mesh for boundary cells
	float dx = (float)m_simparams.xLength / (float)m_simparams.iMax;
	float dy = (float)m_simparams.yLength / (float)m_simparams.jMax;

	dx /= 2.0;
	dy /= 2.0;

	std::array< VertexUV, 4 > cell_vertex_array = {{ VertexUV(-dx,-dy,-1.0,0.0,0.0),
											VertexUV(-dx,dy,-1.0,0.0,1.0),
											VertexUV(dx,dy,-1.0,1.0,1.0),
											VertexUV(dx,-dy,-1.0,1.0,0.0) }};

	std::array< GLuint, 6 > cell_index_array = {{ 0,2,1,2,0,3 }};

	if(!(m_boundaryCells_layer.m_cell.bufferDataFromArray(cell_vertex_array.data(),cell_index_array.data(),sizeof(VertexUV)*4,sizeof(GLuint)*6,GL_TRIANGLES))) return false;
	m_boundaryCells_layer.m_cell.setVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VertexUV),0);
	m_boundaryCells_layer.m_cell.setVertexAttribPointer(1,2,GL_FLOAT,GL_FALSE,sizeof(VertexUV),(GLvoid*) (sizeof(float)*3));

	
	// Create mesh for boundary glyphs
	std::array< VertexUV, 4 > glyph_vertex_array = {{ VertexUV(-dx,-dy,-1.0,0.0,0.0),
											VertexUV(-dx,dy,-1.0,0.0,1.0),
											VertexUV(dx,dy,-1.0,1.0,1.0),
											VertexUV(dx,-dy,-1.0,1.0,0.0) }};

	std::array< GLuint, 6 > glyph_index_array = {{ 0,2,1,2,0,3 }};

	if(!(m_boundaryGlyph_layer.m_glyph.bufferDataFromArray(glyph_vertex_array.data(),glyph_index_array.data(),sizeof(VertexUV)*4,sizeof(GLuint)*6,GL_TRIANGLES))) return false;
	m_boundaryGlyph_layer.m_glyph.setVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VertexUV),0);
	m_boundaryGlyph_layer.m_glyph.setVertexAttribPointer(1,2,GL_FLOAT,GL_FALSE,sizeof(VertexUV),(GLvoid*) (sizeof(float)*3));


	// create mesh for screen filling quad
	std::array< VertexUV, 4 > vertex_array = {{ VertexUV(-1.0,-1.0,0.0,0.0,0.0),
											VertexUV(-1.0,1.0,0.0,0.0,1.0),
											VertexUV(1.0,1.0,0.0,1.0,1.0),
											VertexUV(1.0,-1.0,0.0,1.0,0.0) }};

	std::array< GLuint, 6 > index_array = {{ 0,2,1,2,0,3 }};

	if(!(m_screen_quad.bufferDataFromArray(vertex_array.data(),index_array.data(),sizeof(VertexUV)*4,sizeof(GLuint)*6,GL_TRIANGLES))) return false;
	m_screen_quad.setVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VertexUV),0);
	m_screen_quad.setVertexAttribPointer(1,2,GL_FLOAT,GL_FALSE,sizeof(VertexUV),(GLvoid*) (sizeof(float)*3));

	return true;
}

bool CavityRenderer::createTextures()
{
	unsigned long begin_pos;
	int x_dim, y_dim;
	char* m_img_data;
	readPpmHeader("arrow.ppm", begin_pos, x_dim, y_dim);
	m_img_data = new char[x_dim * y_dim * 3];
	readPpmData("arrow.ppm", m_img_data, begin_pos, x_dim, y_dim);
	m_boundaryGlyph_layer.m_velocity_glyph_tx = std::make_shared<Texture2D>("m_velocity_glyph_tx",GL_RGB, x_dim, y_dim, GL_RGB, GL_UNSIGNED_BYTE, m_img_data);
	delete[] m_img_data;

	readPpmHeader("boundary_cell.ppm", begin_pos, x_dim, y_dim);
	m_img_data = new char[x_dim * y_dim * 3];
	readPpmData("boundary_cell.ppm", m_img_data, begin_pos, x_dim, y_dim);
	m_boundaryCells_layer.m_cell_tx = std::make_shared<Texture2D>("m_cell_tx",GL_RGB, x_dim, y_dim, GL_RGB, GL_UNSIGNED_BYTE, m_img_data);
	delete[] m_img_data;

	readPpmHeader("boundary_p_dlt_glyph.ppm", begin_pos, x_dim, y_dim);
	m_img_data = new char[x_dim * y_dim * 3];
	readPpmData("boundary_p_dlt_glyph.ppm", m_img_data, begin_pos, x_dim, y_dim);
	m_boundaryGlyph_layer.m_pdlt_glyph_tx = std::make_shared<Texture2D>("m_pdlt_glyph_tx",GL_RGB, x_dim, y_dim, GL_RGB, GL_UNSIGNED_BYTE, m_img_data);
	delete[] m_img_data;

	readPpmHeader("boundary_p_nm_glyph.ppm", begin_pos, x_dim, y_dim);
	m_img_data = new char[x_dim * y_dim * 3];
	readPpmData("boundary_p_nm_glyph.ppm", m_img_data, begin_pos, x_dim, y_dim);
	m_boundaryGlyph_layer.m_pnm_glyph_tx = std::make_shared<Texture2D>("m_pnm_glyph_tx",GL_RGB, x_dim, y_dim, GL_RGB, GL_UNSIGNED_BYTE, m_img_data);
	delete[] m_img_data;

	return true;
}

bool CavityRenderer::createFramebuffers()
{
	//TODO create all framebuffers based on window width and height
	m_field_layer.m_fbo = std::make_shared<FramebufferObject>(m_window_width,m_window_height);
	m_field_layer.m_fbo->createColorAttachment(GL_RGBA32F,GL_RGBA,GL_FLOAT);

	m_overlayGrid_layer.m_fbo = std::make_shared<FramebufferObject>(m_window_width,m_window_height);
	m_overlayGrid_layer.m_fbo->createColorAttachment(GL_RGBA32F,GL_RGBA,GL_FLOAT);

	m_boundaryGlyph_layer.m_fbo = std::make_shared<FramebufferObject>(m_window_width,m_window_height);
	m_boundaryGlyph_layer.m_fbo->createColorAttachment(GL_RGBA32F,GL_RGBA,GL_FLOAT);

	m_boundaryCells_layer.m_fbo = std::make_shared<FramebufferObject>(m_window_width,m_window_height);
	m_boundaryCells_layer.m_fbo->createColorAttachment(GL_RGBA32F,GL_RGBA,GL_FLOAT);

	m_geometry_layer.m_fbo = std::make_shared<FramebufferObject>(m_window_width,m_window_height);
	m_geometry_layer.m_fbo->createColorAttachment(GL_RGBA32F,GL_RGBA,GL_FLOAT);

	return true;
}

//void CavityRenderer::reloadSimParams(SimulationParameters& sim_params)
//{
//	m_simparams = sim_params;
//	m_max_boundary_piece = (int)sim_params.boundary_conditions.size();
//}

void CavityRenderer::paint()
{
	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(m_window))
	{
		m_field_layer.draw(m_cam_sys,m_window_background_colour);
		m_boundaryGlyph_layer.draw(m_cam_sys,m_window_background_colour);
		m_boundaryCells_layer.draw(m_cam_sys,m_window_background_colour);
		m_geometry_layer.draw(m_cam_sys,m_window_background_colour);
		m_overlayGrid_layer.draw(m_cam_sys,m_window_background_colour);

		// Merge layers, do additional post processing. Output to primary framebuffer
		postProcessing();

		// Draw TB
		TwDraw();

		/* Swap front and back buffers */
		glfwSwapBuffers(m_window);

        /* Poll for and process events */
        glfwPollEvents();

		/* Check for new simparams */
		SimulationParameters received_params;
		if(m_inbox.tryPop(received_params))
		{
			// overwrite simparams
			m_simparams = received_params;
			// update grid etc.
			m_boundaryCells_layer.setCellPositions(m_simparams);
			m_boundaryCells_layer.updateCellMesh(m_simparams);
			m_boundaryGlyph_layer.setGlyphs(m_simparams);
			m_boundaryGlyph_layer.updateGlyphMesh(m_simparams);
			m_overlayGrid_layer.updateGridMesh(m_simparams);
		}
    }

	//TODO cleanup
	TwTerminate();
	glfwTerminate();
}

void CavityRenderer::postProcessing()
{
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glClearColor(
		m_window_background_colour[0],
		m_window_background_colour[1],
		m_window_background_colour[2], 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	int width, height;
	glfwGetFramebufferSize(m_window, &width, &height);
	glViewport(0, 0, width, height);

	m_postProc_prgm.use();

	glEnable(GL_TEXTURE_2D);
	glActiveTexture(GL_TEXTURE0);
	m_postProc_prgm.setUniform("grid_tx2D",0);
	m_overlayGrid_layer.m_fbo->bindColorbuffer(0);

	glActiveTexture(GL_TEXTURE1);
	m_postProc_prgm.setUniform("boundary_cells_tx2D",1);
	m_boundaryCells_layer.m_fbo->bindColorbuffer(0);

	glActiveTexture(GL_TEXTURE2);
	m_postProc_prgm.setUniform("boundary_glyphs_tx2D",2);
	m_boundaryGlyph_layer.m_fbo->bindColorbuffer(0);

	m_postProc_prgm.setUniform("background_colour",glm::vec3(m_window_background_colour[0],
															m_window_background_colour[1],
															m_window_background_colour[2]));

	m_screen_quad.draw();
}

void CavityRenderer::pushSimParams()
{
	m_outbox.push(m_simparams);
}

/**
 * Example for a callback function that is used in addButtonParam
 */
void TW_CALL Callback(void *clientData)
{
	if(clientData)
		clientData=NULL;
	// do something
}

//void TW_CALL RemoveBoundaryPiece(void* clientData)
//{
//	CavityRenderer* cr = (CavityRenderer*)clientData;
//	cr->setMaxBoundaryPiece(cr->getMaxBoundaryPiece() - 1);
//	cr->deleteBoundaryPiece(cr->getBoundaryPieceIndex());
//}

//void TW_CALL ModifyBoundaryPiece(void* clientData)
//{
//	CavityRenderer* cr = (CavityRenderer*)clientData;
//	cr->modifyBoundaryPieceParams(cr->getBoundaryPieceIndex());
//}

//void TW_CALL BoundaryPiece(void* clientData)
//{
//	CavityRenderer* cr = (CavityRenderer*)clientData;
//	cr->setMaxBoundaryPiece(cr->getMaxBoundaryPiece() + 1);
//	Boundary::Direction dir;
//	Boundary::Condition cond;
//	Boundary::Grid grid;
//	Real value;
//	int i_begin;
//	int i_end;
//	int j_begin;
//	int j_end;
//	cr->getBoundaryPieceParams(dir, cond, grid, value, i_begin, i_end, j_begin, j_end);
//	Index begin = Index(i_begin, j_begin);
//	Index end = Index(i_end, j_end);
//	Range range = Range(begin, end);
//	cr->addBoundaryPiece(Boundary::BoundaryPiece(dir, cond, grid, value, range));
//	cr->showBoundaryPiece(cr->getMaxBoundaryPiece() - 1);
//}

void TW_CALL Bake(void* clientData)
{
	CavityRenderer* cr = (CavityRenderer*)clientData;

	// this ones job will be to push updates simulation parameters to communication queue
	cr->pushSimParams();
}

//void CavityRenderer::addBoundaryPieceToBar(std::string mode)
//{
//	TwEnumVal direction_enum[] = {
//		{ (int)Boundary::Direction::Down, "Down" },
//		{ (int)Boundary::Direction::Left, "Left" },
//		{ (int)Boundary::Direction::Right, "Right" },
//		{ (int)Boundary::Direction::Up, "Up" }
//	};
//	m_direction_enum = Boundary::Direction::Down;
//	addEnumParam("Direction", "DirectionType", &m_direction_enum, direction_enum, 4, mode);
//
//	TwEnumVal condition_enum[] = {
//		{ (int)Boundary::Condition::NOSLIP, "NOSLIP" },
//		{ (int)Boundary::Condition::INFLOW, "INFLOW" },
//		{ (int)Boundary::Condition::OUTFLOW, "OUTFLOW" },
//		{ (int)Boundary::Condition::SLIP, "SLIP" }
//	};
//	m_condition_enum = Boundary::Condition::NOSLIP;
//	addEnumParam("Condition", "ConditionType", &m_condition_enum, condition_enum, 4, mode);
//
//	TwEnumVal grid_enum[] = {
//		{ (int)Boundary::Grid::U, "U" },
//		{ (int)Boundary::Grid::V, "V" },
//		{ (int)Boundary::Grid::P, "P" },
//		{ (int)Boundary::Grid::F, "F" },
//		{ (int)Boundary::Grid::G, "G" }
//	};
//	m_grid_enum = Boundary::Grid::U;
//	addEnumParam("Grid", "GridType", &m_grid_enum, grid_enum, 5, mode);
//	Real test; float __float; double __double;
//	const char* _double = typeid(__double).name();
//	const char* _float = typeid(__float).name();
//	if (strcmp(typeid(test).name(), _double) == 0)
//		addDoubleParam("m_condition_value", " step=0.1 label='Condition Value' ", &m_condition_value, mode);
//	if (strcmp(typeid(test).name(), _float) == 0)
//		addFloatParam("m_condition_value", " step=0.1 label='Condition Value' ", &m_condition_value, mode);
//	addIntParam("m_i_begin", " label='i begin' ", &m_i_begin, mode);
//	addIntParam("m_i_end", " label='i end' ", &m_i_end, mode);
//	addIntParam("m_j_begin", " label='j begin' ", &m_j_begin, mode);
//	addIntParam("m_j_end", " label='j end' ", &m_j_end, mode);
//}

//void CavityRenderer::showBoundaryPiece(unsigned int index)
//{
//	if (!m_boundary_conditions.empty() && !m_modify_cond)
//	{
//		m_direction_enum = m_boundary_conditions.at(index).direction;
//		m_condition_enum = m_boundary_conditions.at(index).condition;
//		m_grid_enum = m_boundary_conditions.at(index).gridtype;
//		m_condition_value = m_boundary_conditions.at(index).condition_value;
//		m_i_begin = m_boundary_conditions.at(index).range.begin[0];
//		m_j_begin = m_boundary_conditions.at(index).range.begin[1];
//		m_i_end = m_boundary_conditions.at(index).range.end[0];
//		m_j_end = m_boundary_conditions.at(index).range.end[1];
//
//		modifyIntParam("m_nmbr_boundary_piece", 0, m_max_boundary_piece - 1);
//		drawBoundaryCondition(m_boundary_conditions.at(index));
//	}
//}

//void CavityRenderer::modifyBoundaryPieceParams(unsigned int index)
//{
//	if (!m_boundary_conditions.empty())
//	{
//		m_boundary_conditions.at(index).direction = m_direction_enum;
//		m_boundary_conditions.at(index).condition = m_condition_enum;
//		m_boundary_conditions.at(index).gridtype = m_grid_enum;
//		m_boundary_conditions.at(index).condition_value = m_condition_value;
//		m_boundary_conditions.at(index).range.begin[0] = m_i_begin;
//		m_boundary_conditions.at(index).range.begin[1] = m_j_begin;
//		m_boundary_conditions.at(index).range.end[0] = m_i_end;
//		m_boundary_conditions.at(index).range.end[1] = m_j_end;
//	}
//}

void CavityRenderer::addFloatParam(const char* name, const char* def, void* var, std::string mode, float min, float max)
{
	if (mode == "RW")
	{
		TwAddVarRW(bar, name, TW_TYPE_FLOAT, var, def);
		TwSetParam(bar, name, "min", TW_PARAM_FLOAT, 1, &min);
		TwSetParam(bar, name, "max", TW_PARAM_FLOAT, 1, &max);
	}
	if (mode == "RO")
	{
		TwAddVarRO(bar, name, TW_TYPE_FLOAT, var, def);
	}
}

void CavityRenderer::addDoubleParam(const char* name, const char* def, void* var, std::string mode, double min, double max)
{
	if (mode == "RW")
	{
		TwAddVarRW(bar, name, TW_TYPE_DOUBLE, var, def);
		TwSetParam(bar, name, "min", TW_PARAM_DOUBLE, 1, &min);
		TwSetParam(bar, name, "max", TW_PARAM_DOUBLE, 1, &max);
	}
	if (mode == "RO")
	{
		TwAddVarRO(bar, name, TW_TYPE_DOUBLE, var, def);
	}
}

void CavityRenderer::addIntParam(const char* name, const char* def, void* var, std::string mode, int min, int max)
{
	if (mode == "RW")
	{
		TwAddVarRW(bar, name, TW_TYPE_INT32, var, def);
		TwSetParam(bar, name, "min", TW_PARAM_INT32, 1, &min);
		TwSetParam(bar, name, "max", TW_PARAM_INT32, 1, &max);
	}
	if (mode == "RO")
	{
		TwAddVarRO(bar, name, TW_TYPE_INT32, var, def);
	}
}

void CavityRenderer::addBoolParam(const char* name, const char* def, void* var, std::string mode)
{
	if (mode == "RW")
	{
		TwAddVarRW(bar, name, TW_TYPE_BOOLCPP, var, def);
	}
	if (mode == "RO")
	{
		TwAddVarRO(bar, name, TW_TYPE_BOOLCPP, var, def);
	}
}

void CavityRenderer::addVec3Param(const char* name, const char* def, void* var, std::string mode)
{
	if (mode == "RW")
	{
		TwAddVarRW(bar, name, TW_TYPE_DIR3F, var, def);
	}
	if (mode == "RO")
	{
		TwAddVarRO(bar, name, TW_TYPE_DIR3F, var, def);
	}
}

void CavityRenderer::addButtonParam(const char* name, const char* def, TwButtonCallback callback)
{
	TwAddButton(bar, name, callback, this, def);
}

void CavityRenderer::addStringParam(const char* name, const char* def, void* var, std::string mode)
{
	if (mode == "RW")
	{
		TwAddVarRW(bar, name, TW_TYPE_CDSTRING, var, def);
	}
	if (mode == "RO")
	{
		TwAddVarRO(bar, name, TW_TYPE_CDSTRING, var, def);
	}
}

void CavityRenderer::addEnumParam(const char* name, const char* def, void* var, TwEnumVal* _enum, int size, std::string mode)
{
	TwType enumType = TwDefineEnum(def, _enum, size);
	if (mode == "RW") TwAddVarRW(bar, name, enumType, var, NULL);
	if (mode == "RO") TwAddVarRO(bar, name, enumType, var, NULL);
}

void CavityRenderer::modifyIntParam(const char* name, int min, int max)
{
	TwSetParam(bar, name, "min", TW_PARAM_INT32, 1, &min);
	TwSetParam(bar, name, "max", TW_PARAM_INT32, 1, &max);
}

void CavityRenderer::removeParam(const char* name)
{
	TwRemoveVar(bar, name);
}


const std::string CavityRenderer::readShaderFile(const char* const path)
{
	std::ifstream inFile(path, std::ios::in);
	std::ostringstream source;
	while (inFile.good()) {
		int c = inFile.get();
		if (!inFile.eof()) source << (char)c;
	}
	inFile.close();
	return source.str();
}

bool CavityRenderer::readPpmHeader(const char* filename, unsigned long& headerEndPos, int& imgDimX, int& imgDimY)
{
	int currentComponent = 0;
	bool firstline = false;
	std::string::iterator itr1;
	std::string::iterator itr2;
	std::string buffer;
	std::string compBuffer;
	std::ifstream file(filename, std::ios::in | std::ios::binary);
	/*
	/ Check if the file could be opened.
	*/
	if (!(file.is_open()))return false;
	/*
	/ Go to the beginning of the file and read the first line.
	*/
	file.seekg(0, file.beg);
	std::getline(file, buffer, '\n');
	itr1 = buffer.begin();
	for (itr2 = buffer.begin(); itr2 != buffer.end(); itr2++)
	{
		/*
		/ Check if the first line contains more than just ppm's magic number.
		/ If it does, it should look like this:
		/ "magic_number image_dimension_x image_dimension_y maximum_value"
		/ Therefore we scan the string for a space character and start parsing it.
		*/
		if (*itr2 == ' ')
		{
			if (currentComponent == 0)
			{
				/* The first component is the magic number. We don't need it. */
				currentComponent++;
				firstline = true;
				itr1 = (itr2 + 1);
			}
			else if (currentComponent == 1)
			{
				/* Get the image dimension in x. */
				compBuffer.assign(itr1, itr2);
				imgDimX = atoi(compBuffer.c_str());
				currentComponent++;
				itr1 = (itr2 + 1);
			}
			else if (currentComponent == 2)
			{
				/* Get the image dimension in y. */
				compBuffer.assign(itr1, itr2);
				imgDimY = atoi(compBuffer.c_str());
				currentComponent++;
				itr1 = (itr2 + 1);
			}
		}
	}
	/*
	/ If the information we were looking for was inside the first line, we are done here.
	/ Note the position where we left off and exit with return true after closing the file.
	*/
	if (firstline)
	{
		headerEndPos = static_cast<long>(file.tellg());
		file.close();
		return true;
	}
	/*
	/ If the information wasn't inside the first line we have to keep reading lines.
	/ Skip all comment lines (first character = '#').
	*/
	std::getline(file, buffer, '\n');
	while (buffer[0] == '#' || (buffer.size() < 1))
	{
		std::getline(file, buffer, '\n');
	}
	/*
	/ Now we should have a string containing the image dimensions and can extract them.
	*/
	itr1 = buffer.begin();
	for (itr2 = buffer.begin(); itr2 != buffer.end(); itr2++)
	{
		/* Get the image dimension in x. */
		if (*itr2 == ' ')
		{
			compBuffer.assign(itr1, itr2);
			imgDimX = atoi(compBuffer.c_str());
			currentComponent++;
			itr1 = (itr2 + 1);
		}
	}
	/*
	/ The last component of a line can't be parsed within the loop since it isn't followed by
	/ a space character, but an end-of-line.
	/
	/ Get the image dimension in x.
	*/
	compBuffer.assign(itr1, itr2);
	imgDimY = atoi(compBuffer.c_str());
	/*
	/ Read one more line. This should contain the maximum value of the image, but we don't need
	/ that.
	/ Note down the position after this line and exit with return true after closing the file.
	*/
	std::getline(file, buffer, '\n');
	headerEndPos = static_cast<unsigned long>(file.tellg());
	file.close();
	return true;
}

bool CavityRenderer::readPpmData(const char* filename, char* imageData, unsigned long dataBegin, int imgDimX, int imgDimY)
{
	std::ifstream file(filename, std::ios::in | std::ios::binary);
	/*
	/ Check if the file could be opened.
	*/
	if (!(file.is_open()))return false;
	/*
	/ Determine the length from the beginning of the image data to the end of the file.
	*/
	file.seekg(0, file.end);
	unsigned long length = static_cast<unsigned long>(file.tellg());
	length = length - dataBegin;
	char* buffer = new char[length];
	file.seekg(dataBegin, std::ios::beg);
	file.read(buffer, length);
	/*
	/ Rearrange the image information so that the data begins with the lower left corner.
	*/
	int k = 0;
	for (int i = 0; i < imgDimY; i++)
	{
		int dataLoc = (imgDimY - 1 - i)*imgDimX * 3;
		for (int j = 0; j < imgDimX; j++)
		{
			imageData[k] = buffer[dataLoc + (j * 3)];
			k++;
			imageData[k] = buffer[dataLoc + (j * 3) + 1];
			k++;
			imageData[k] = buffer[dataLoc + (j * 3) + 2];
			k++;
		}
	}
	file.close();
	delete[] buffer;
	return true;
}
