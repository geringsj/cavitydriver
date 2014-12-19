#include "CavityRenderer.hpp"

/**
 * Example for a callback function that is used in addButtonParam
 */
void TW_CALL Callback(void *clientData)
{
	// do something
}

void TW_CALL Bake(void* clientData)
{
	// start baking
}

CavityRenderer::CavityRenderer()
{
	/**
	 * And this also doesn't work because AntTweakBar,
	 * apparently hates pointers and crashes for (no?)
	 * good reason.
	 */
	//m_sim_params = NULL;
}

CavityRenderer::~CavityRenderer()
{
}

bool CavityRenderer::init(unsigned int window_width, unsigned int window_height, SimulationParameters& sim_params)
{
	m_window_width = window_width;
	m_window_height = window_height;
	m_window_background_colour[0] = 0.2;m_window_background_colour[1] = 0.2;m_window_background_colour[2] = 0.2;
	m_zoom = 80.0f;

	m_show_grid = true;

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
	bar = TwNewBar("TweakBar");
	TwWindowSize(window_width, window_height);
	TwDefine(" GLOBAL help='This example shows how to integrate AntTweakBar with GLFW and OpenGL.' "); // Message added to the help bar.
	// Add 'bgColor' to 'bar': it is a modifable variable of type TW_TYPE_COLOR3F (3 floats color)
	TwAddVarRW(bar, "m_window_background_colour", TW_TYPE_COLOR3F, &m_window_background_colour, " label='Background color' ");
	addFloatParam("m_zoom", " label='Zoom' ", &m_zoom, 1.0f, 999.0f);
	addBoolParam("m_show_grid", " label='Show grid' ", &m_show_grid);

	// Set GLFW event callbacks
	glfwSetWindowUserPointer(m_window,this);
	glfwSetWindowSizeCallback(m_window, (GLFWwindowposfun)windowResizeCallback);
    glfwSetMouseButtonCallback(m_window, (GLFWmousebuttonfun)mouseButtonCallback);
    glfwSetCursorPosCallback(m_window, (GLFWcursorposfun)mousePositionCallback);
    glfwSetScrollCallback(m_window, (GLFWscrollfun)mouseWheelCallback);
    glfwSetKeyCallback(m_window, (GLFWkeyfun)keyCallback);
    glfwSetCharCallback(m_window, (GLFWcharfun)glfwSetCharCallback);

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


	m_cam_sys = CameraSystem(glm::vec3(0.0f, 0.0f, m_zoom), glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(0.0f, 0.0f, -1.0f), glm::vec3(1.0f, 0.0f, 0.0f));
	return true;
}

bool CavityRenderer::initBakeryVis(unsigned int window_width, unsigned int window_height, SimulationParameters& sim_params)
{
	/**
	 * And this also doesn't work because AntTweakBar,
	 * apparently hates pointers and crashes for (no?)
	 * good reason.
	 */
	//m_sim_params = new SimulationParameters(sim_params);
	m_window_width = window_width;
	m_window_height = window_height;
	m_window_background_colour[0] = 0.2; m_window_background_colour[1] = 0.2; m_window_background_colour[2] = 0.2;
	m_zoom = 80.0f;

	m_show_grid = true;

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
	bar = TwNewBar("TweakBar");
	TwWindowSize(window_width, window_height);
	TwDefine(" GLOBAL help='This example shows how to integrate AntTweakBar with GLFW and OpenGL.' "); // Message added to the help bar.
	// Add 'bgColor' to 'bar': it is a modifable variable of type TW_TYPE_COLOR3F (3 floats color)
	TwAddVarRW(bar, "m_window_background_colour", TW_TYPE_COLOR3F, &m_window_background_colour, " label='Background color' ");
	addFloatParam("m_zoom", " label='Zoom' ", &m_zoom, 1.0f, 999.0f);
	addBoolParam("m_show_grid", " label='Show grid' ", &m_show_grid);
	TwAddSeparator(bar, "SimulationParameters", " label='SimulationParameters' ");
	m_alpha = sim_params.alpha;
	m_deltaT = sim_params.deltaT; 
	m_deltaVec = sim_params.deltaVec; 
	m_eps = sim_params.eps; 
	m_gx = sim_params.gx; 
	m_gy = sim_params.gy; 
	m_iMax = sim_params.iMax;
	m_iterMax = sim_params.iterMax; 
	m_jMax = sim_params.jMax; 
	m_KarmanAngle = sim_params.KarmanAngle; 
	m_KarmanObjectWidth = sim_params.KarmanObjectWidth; 
	m_name = sim_params.name; 
	m_omg = sim_params.omg;
	m_pi = sim_params.pi; 
	m_re = sim_params.re; 
	m_tau = sim_params.tau; 
	m_tDeltaWriteVTK = sim_params.tDeltaWriteVTK; 
	m_tEnd = sim_params.tEnd; 
	m_ui = sim_params.ui; 
	m_vi = sim_params.vi;
	m_xCells = sim_params.xCells; 
	m_xLength = sim_params.xLength; 
	m_yCells = sim_params.yCells; 
	m_yLength = sim_params.yLength;
	addFloatParam("m_alpha", " label='alpha' ", &m_alpha);
	addFloatParam("m_deltaT", " label='deltaT' ", &m_deltaT);
	addFloatParam("m_deltaVec", " label='deltaVec' ", &m_deltaVec);
	addFloatParam("m_eps", " label='eps' ", &m_eps);
	addFloatParam("m_gx", " label='gx' ", &m_gx);
	addFloatParam("m_gy", " label='gy' ", &m_gy);
	addFloatParam("m_KarmanAngle", " label='KarmanAngle' ", &m_KarmanAngle);
	addFloatParam("m_KarmanObjectWidth", " label='KarmanObjectWidth' ", &m_KarmanObjectWidth);
	addFloatParam("m_pi", " label='pi' ", &m_pi);
	addFloatParam("m_re", " label='re' ", &m_re);
	addFloatParam("m_tau", " label='tau' ", &m_tau);
	addFloatParam("m_tDeltaWriteVTK", " label='tDeltaWriteVTK' ", &m_tDeltaWriteVTK);
	addFloatParam("m_tEnd", " label='tEnd' ", &m_tEnd);
	addFloatParam("m_ui", " label='ui' ", &m_ui);
	addFloatParam("m_vi", " label='vi' ", &m_vi);
	addFloatParam("m_xLength", " label='xLength' ", &m_xLength);
	addFloatParam("m_yLength", " label='yLength' ", &m_yLength);
	addFloatParam("m_omg", " label='omega' ", &m_omg);
	addIntParam("m_iterMax", " label='iterMax' ", &m_iterMax);
	addIntParam("m_iMax", " label='iMax' ", &m_iMax);
	addIntParam("m_jMax", " label='jMax' ", &m_jMax);
	addIntParam("m_xCells", " label='xCells' ", &m_xCells);
	addIntParam("m_yCells", " label='yCells' ", &m_yCells);
	addStringParam("m_name", " label='name' ", &m_name);

	addButtonParam("m_bake", " label='bake the parameter' ", Bake);

	/**
	 * And this also doesn't work because AntTweakBar,
	 * apparently hates pointers and crashes for (no?)
	 * good reason.
	 */
	//addFloatParam("m_alpha", " label='alpha' ", &m_sim_params->alpha);
	//addFloatParam("m_deltaT", " label='deltaT' ", &m_sim_params->deltaT);
	//addFloatParam("m_deltaVec", " label='deltaVec' ", &m_sim_params->deltaVec);
	//addFloatParam("m_eps", " label='eps' ", &m_sim_params->eps);
	//addFloatParam("m_gx", " label='gx' ", &m_sim_params->gx);
	//addFloatParam("m_gy", " label='gy' ", &m_sim_params->gy);
	//addFloatParam("m_KarmanAngle", " label='KarmanAngle' ", &m_sim_params->KarmanAngle);
	//addFloatParam("m_KarmanObjectWidth", " label='KarmanObjectWidth' ", &m_sim_params->KarmanObjectWidth);
	//addFloatParam("m_pi", " label='pi' ", &m_sim_params->pi);
	//addFloatParam("m_re", " label='re' ", &m_sim_params->re);
	//addFloatParam("m_tau", " label='tau' ", &m_sim_params->tau);
	//addFloatParam("m_tDeltaWriteVTK", " label='tDeltaWriteVTK' ", &m_sim_params->tDeltaWriteVTK);
	//addFloatParam("m_tEnd", " label='tEnd' ", &m_sim_params->tEnd);
	//addFloatParam("m_ui", " label='ui' ", &m_sim_params->ui);
	//addFloatParam("m_vi", " label='vi' ", &m_sim_params->vi);
	//addFloatParam("m_xLength", " label='xLength' ", &m_sim_params->xLength);
	//addFloatParam("m_yLength", " label='yLength' ", &m_sim_params->yLength);
	//addFloatParam("m_omg", " label='omega' ", &m_sim_params->omg);
	//addIntParam("m_iMax", " label='iMax' ", &m_sim_params->iMax);
	//addIntParam("m_iterMax", " label='iterMax' ", &m_sim_params->iterMax);
	//addIntParam("m_jMax", " label='jMax' ", &m_sim_params->jMax);
	//addIntParam("m_xCells", " label='xCells' ", &m_sim_params->xCells);
	//addIntParam("m_yCells", " label='yCells' ", &m_sim_params->yCells);
	//addStringParam("m_name", " label='name' ", &m_sim_params->name);

	// Set GLFW event callbacks
	glfwSetWindowUserPointer(m_window, this);
	glfwSetWindowSizeCallback(m_window, (GLFWwindowposfun)windowResizeCallback);
	glfwSetMouseButtonCallback(m_window, (GLFWmousebuttonfun)mouseButtonCallback);
	glfwSetCursorPosCallback(m_window, (GLFWcursorposfun)mousePositionCallback);
	glfwSetScrollCallback(m_window, (GLFWscrollfun)mouseWheelCallback);
	glfwSetKeyCallback(m_window, (GLFWkeyfun)keyCallback);
	glfwSetCharCallback(m_window, (GLFWcharfun)glfwSetCharCallback);

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


	m_cam_sys = CameraSystem(glm::vec3(0.0f, 0.0f, m_zoom), glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(0.0f, 0.0f, -1.0f), glm::vec3(1.0f, 0.0f, 0.0f));
	return true;
}

void CavityRenderer::paint()
{

	/* Loop until the user closes the window */
    while (!glfwWindowShouldClose(m_window))
    {
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glClearColor(m_window_background_colour[0], m_window_background_colour[1], m_window_background_colour[2], 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		int width, height;
		glfwGetFramebufferSize(m_window, &width, &height);
		glViewport(0, 0, width, height);
		
		drawGridOverlay();

		// Draw TB
		TwDraw();
		// Draw TB

        /* Swap front and back buffers */
        glfwSwapBuffers(m_window);

        /* Poll for and process events */
        glfwPollEvents();
    }

	//TODO cleanup
	TwTerminate();
	glfwTerminate();
}

bool CavityRenderer::createGrid(Range grid_size)
{
	std::vector<unsigned int> index_array;
	std::vector<Gridvertex> vertex_array;
	createSingleGrid(grid_size, index_array, vertex_array);

	float right_shift = (float)((grid_size.end[0] + 1) - (grid_size.begin[0] - 1)) / 2.0f;
	float up_shift = (float)((grid_size.end[1] + 1) - (grid_size.begin[1] - 1)) / 2.0f;

	m_cam_sys.Translation(glm::vec3(1.0f, 0.0f, 0.0f), right_shift*1.0f);
	m_cam_sys.Translation(glm::vec3(0.0f, 1.0f, 0.0f), up_shift*1.0f);

	m_grid.bufferDataFromArray(vertex_array.data(),index_array.data(),
		(GLsizei)(vertex_array.size()*sizeof(Gridvertex)),(GLsizei)(index_array.size()*sizeof(unsigned int)),GL_LINES);
	m_grid.setVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(Gridvertex), 0);

	m_grid_prgm.init();

	std::string grid_vertex = readShaderFile("./shader/gridVertex.glsl");
	if (!m_grid_prgm.compileShaderFromString(&grid_vertex, GL_VERTEX_SHADER)) { std::cout << m_grid_prgm.getLog(); return false; };

	std::string grid_fragment = readShaderFile("./shader/gridFragment.glsl");
	if (!m_grid_prgm.compileShaderFromString(&grid_fragment, GL_FRAGMENT_SHADER)) { std::cout << m_grid_prgm.getLog(); return false; };

	m_grid_prgm.bindAttribLocation(0, "in_position");

	m_grid_prgm.link();

	return true;
}

void CavityRenderer::createSingleGrid(Range innerRange, std::vector<unsigned int>& index, std::vector<Gridvertex>& data)
{
	unsigned int bottom_left_i = innerRange.begin[0] - 1;
	unsigned int bottom_left_j = innerRange.begin[1] - 1;

	unsigned int bottom_right_i = innerRange.end[0] + 1;
	unsigned int bottom_right_j = innerRange.begin[1] - 1;

	unsigned int top_right_i = innerRange.end[0] + 1;
	unsigned int top_right_j = innerRange.end[1] + 1;

	unsigned int top_left_i = innerRange.begin[0] - 1;
	unsigned int top_left_j = innerRange.end[1] + 1;

	data.push_back(Gridvertex((float)bottom_left_i, (float)bottom_left_j, -1.0f, 1.0f));
	index.push_back(0);
	data.push_back(Gridvertex((float)bottom_right_i, (float)bottom_right_j, -1.0f, 1.0f));
	index.push_back(1);
	data.push_back(Gridvertex((float)top_right_i, (float)top_right_j, -1.0f, 1.0f));
	index.push_back(2);
	data.push_back(Gridvertex((float)top_left_i, (float)top_left_j, -1.0f, 1.0f));
	index.push_back(3);

	index.push_back(0);
	index.push_back(3);
	index.push_back(1);
	index.push_back(2);

	//LEFT RIGHT
	unsigned int left = bottom_left_i;
	unsigned int right = bottom_right_i;
	unsigned int index_value = 4;
	for (unsigned int j = bottom_left_j; j <= top_left_j; j++)
	{
		data.push_back(Gridvertex((float)left, (float)j, -1.0f, 1.0f));
		index.push_back(index_value); index_value++;
		data.push_back(Gridvertex((float)right, (float)j, -1.0f, 1.0f));
		index.push_back(index_value); index_value++;
	}
	//TOP BOTTOM
	unsigned int bottom = bottom_right_j;
	unsigned int top = top_left_j;
	for (unsigned int i = bottom_left_i; i <= bottom_right_i; i++)
	{
		data.push_back(Gridvertex((float)i, (float)bottom, -1.0f, 1.0f));
		index.push_back(index_value); index_value++;
		data.push_back(Gridvertex((float)i, (float)top, -1.0f, 1.0f));
		index.push_back(index_value); index_value++;
	}
}

void CavityRenderer::drawGridOverlay()
{
	m_cam_sys.Translation(glm::vec3(0.0f, 0.0f, -1.0f), m_cam_sys.GetCamPos().z - m_zoom);
	m_grid_prgm.use();
	glm::mat4 proj_mat = glm::perspective(45.0f, (float)m_window_width / (float)m_window_height, 0.1f, 1000.0f);
	glm::mat4 model_mat = glm::mat4(1.0f);
	glm::mat4 view_mat = m_cam_sys.GetViewMatrix();
	glm::mat4 mvp_mat = proj_mat * view_mat * model_mat;
	m_grid_prgm.setUniform("mvp_matrix", mvp_mat);

	if(m_show_grid)
		m_grid.draw();
}

void CavityRenderer::drawField()
{
}

void CavityRenderer::drawLIC()
{
}

void CavityRenderer::drawGeometry()
{
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

void CavityRenderer::addFloatParam(const char* name, const char* def, void* var, float min, float max)
{
	TwAddVarRW(bar, name, TW_TYPE_FLOAT, var, def);
	TwSetParam(bar, name, "min", TW_PARAM_FLOAT, 1, &min);
	TwSetParam(bar, name, "max", TW_PARAM_FLOAT, 1, &max);
}

void CavityRenderer::addIntParam(const char* name, const char* def, void* var, int min, int max)
{
	TwAddVarRW(bar, name, TW_TYPE_INT32, var, def);
	TwSetParam(bar, name, "min", TW_PARAM_INT32, 1, &min);
	TwSetParam(bar, name, "max", TW_PARAM_INT32, 1, &max);
}

void CavityRenderer::addBoolParam(const char* name, const char* def, void* var)
{
	TwAddVarRW(bar, name, TW_TYPE_BOOL32, var, def);
}

void CavityRenderer::addVec3Param(const char* name, const char* def, void* var)
{
	TwAddVarRW(bar, name, TW_TYPE_DIR3F, var, def);
}

void CavityRenderer::addButtonParam(const char* name, const char* def, TwButtonCallback callback)
{
	TwAddButton(bar, name, callback, this, def);
}

void CavityRenderer::addStringParam(const char* name, const char* def, void* var)
{
	TwAddVarRW(bar, name, TW_TYPE_CDSTRING, var, def);
}

void CavityRenderer::removeParam(const char* name)
{
	TwRemoveVar(bar, name);
}