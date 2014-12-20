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
	CavityRenderer* cr = (CavityRenderer*)clientData;
	Index begin = Index(1, 1);
	int iMax, jMax;
	cr->getijMax(iMax, jMax);
	Index end = Index(iMax,jMax);
	Range new_range = Range(begin, end);
	cr->createGrid(new_range);
}

CavityRenderer::CavityRenderer()
{
	/**
	 * And this also doesn't work because AntTweakBar,
	 * apparently hates pointers and crashes for (no?)
	 * good reason.
	 */
	//m_sim_params = NULL;
	m_grid_resize = false;
}

CavityRenderer::~CavityRenderer()
{
}

bool CavityRenderer::init(unsigned int window_width, unsigned int window_height, SimulationParameters& sim_params)
{
	m_window_width = window_width;
	m_window_height = window_height;
	m_window_background_colour[0] = 0.2f;m_window_background_colour[1] = 0.2f;m_window_background_colour[2] = 0.2f;
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
	m_window_background_colour[0] = 0.2f; m_window_background_colour[1] = 0.2f; m_window_background_colour[2] = 0.2f;
	m_zoom = 1.1f;

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
	addFloatParam("m_zoom", " step=0.1 label='Zoom' ", &m_zoom, 1.0f, 999.0f);
	addBoolParam("m_show_grid", " label='Show grid' ", &m_show_grid);
	TwAddSeparator(bar, "SimulationParameters", " label='SimulationParameters' ");
	m_alpha = (float)sim_params.alpha;
	m_deltaT = (float)sim_params.deltaT;
	m_deltaVec = (float)sim_params.deltaVec;
	m_eps = (float)sim_params.eps;
	m_gx = (float)sim_params.gx;
	m_gy = (float)sim_params.gy;
	m_iMax = (int)sim_params.iMax;
	m_iterMax = (int)sim_params.iterMax;
	m_jMax = (int)sim_params.jMax;
	m_KarmanAngle = (float)sim_params.KarmanAngle;
	m_KarmanObjectWidth = (float)sim_params.KarmanObjectWidth;
	m_name = sim_params.name;
	m_omg = (float)sim_params.omg;
	m_pi = (float)sim_params.pi;
	m_re = (float)sim_params.re;
	m_tau = (float)sim_params.tau;
	m_tDeltaWriteVTK = (float)sim_params.tDeltaWriteVTK;
	m_tEnd = (float)sim_params.tEnd;
	m_ui = (float)sim_params.ui;
	m_vi = (float)sim_params.vi;
	m_xCells = (int)sim_params.xCells;
	m_xLength = (float)sim_params.xLength;
	m_yCells = (int)sim_params.yCells;
	m_yLength = (float)sim_params.yLength;
	printf("m_xLength: %f m_yLength: %f \n", m_xLength, m_yLength);
	addFloatParam("m_alpha", " step=0.1 label='alpha' ", &m_alpha);
	addFloatParam("m_deltaT", " step=0.1 label='deltaT' ", &m_deltaT);
	addFloatParam("m_deltaVec", " step=0.1 label='deltaVec' ", &m_deltaVec);
	addFloatParam("m_eps", " step=0.001 label='eps' ", &m_eps);
	addFloatParam("m_gx", " step=0.1 label='gx' ", &m_gx);
	addFloatParam("m_gy", " step=0.1 label='gy' ", &m_gy);
	addFloatParam("m_KarmanAngle", " step=0.1 label='KarmanAngle' ", &m_KarmanAngle);
	addFloatParam("m_KarmanObjectWidth", " step=0.1 label='KarmanObjectWidth' ", &m_KarmanObjectWidth);
	addFloatParam("m_pi", " step=0.1 label='pi' ", &m_pi);
	addFloatParam("m_re", " step=0.1 label='re' ", &m_re);
	addFloatParam("m_tau", " step=0.1 label='tau' ", &m_tau);
	addFloatParam("m_tDeltaWriteVTK", " step=0.1 label='tDeltaWriteVTK' ", &m_tDeltaWriteVTK);
	addFloatParam("m_tEnd", " step=0.1 label='tEnd' ", &m_tEnd);
	addFloatParam("m_ui", " step=0.1 label='ui' ", &m_ui);
	addFloatParam("m_vi", " step=0.1 label='vi' ", &m_vi);
	addFloatParam("m_xLength", " step=0.1 label='xLength' ", &m_xLength);
	addFloatParam("m_yLength", " step=0.1 label='yLength' ", &m_yLength);
	addFloatParam("m_omg", " step=0.1 label='omega' ", &m_omg);
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
	
	m_grid.bufferDataFromArray(vertex_array.data(),index_array.data(),
		(GLsizei)(vertex_array.size()*sizeof(Gridvertex)),(GLsizei)(index_array.size()*sizeof(unsigned int)),GL_LINES);
	m_grid.setVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(Gridvertex), 0);

	if (!m_grid_resize)
	{
		m_grid_resize = true;
		m_grid_prgm.init();

		std::string grid_vertex = readShaderFile("./shader/gridVertex.glsl");
		if (!m_grid_prgm.compileShaderFromString(&grid_vertex, GL_VERTEX_SHADER)) { std::cout << m_grid_prgm.getLog(); return false; };

		std::string grid_fragment = readShaderFile("./shader/gridFragment.glsl");
		if (!m_grid_prgm.compileShaderFromString(&grid_fragment, GL_FRAGMENT_SHADER)) { std::cout << m_grid_prgm.getLog(); return false; };

		m_grid_prgm.bindAttribLocation(0, "in_position");

		m_grid_prgm.link();
	}
	return true;
}

void CavityRenderer::createSingleGrid(Range innerRange, std::vector<unsigned int>& index, std::vector<Gridvertex>& data)
{
	float x_length = m_xLength / (float)m_iMax;
	float y_length = m_yLength / (float)m_jMax;

	printf("m_xLength: %f m_yLength: %f \n", m_xLength, m_yLength);

	float bottom_left_i = (float)innerRange.begin[0] - 1.0f;
	float bottom_left_j = (float)innerRange.begin[1] - 1.0f;

	float bottom_right_i = bottom_left_i + m_xCells * x_length;
	float bottom_right_j = bottom_left_j;

	float top_right_i = bottom_right_i;
	float top_right_j = bottom_right_j + m_yCells * y_length;

	float top_left_i = bottom_left_i;
	float top_left_j = bottom_left_j + m_yCells * y_length;

	float right_shift = (m_xCells * x_length) / 2.0f;
	float up_shift = (m_yCells * y_length) / 2.0f;

	m_cam_sys.Translation(m_cam_sys.GetRightVector(), 0.0f - m_cam_sys.GetCamPos().x);
	m_cam_sys.Translation(m_cam_sys.GetRightVector(), right_shift);
	m_cam_sys.Translation(m_cam_sys.GetUpVector(), 0.0f - m_cam_sys.GetCamPos().y);
	m_cam_sys.Translation(m_cam_sys.GetUpVector(), up_shift);

	data.push_back(Gridvertex(bottom_left_i, bottom_left_j, -1.0f, 1.0f));
	index.push_back(0);
	data.push_back(Gridvertex(bottom_right_i, bottom_right_j, -1.0f, 1.0f));
	index.push_back(1);
	data.push_back(Gridvertex(top_right_i, top_right_j, -1.0f, 1.0f));
	index.push_back(2);
	data.push_back(Gridvertex(top_left_i, top_left_j, -1.0f, 1.0f));
	index.push_back(3);

	index.push_back(0);
	index.push_back(3);
	index.push_back(1);
	index.push_back(2);

	//LEFT RIGHT
	float left = bottom_left_i;
	float right = bottom_right_i;
	unsigned int index_value = 4;
	for (unsigned int j = innerRange.begin[1] - 1; j <= (unsigned int)innerRange.end[1]; j++)
	{
		data.push_back(Gridvertex(left, (float)j * y_length, -1.0f, 1.0f));
		index.push_back(index_value); index_value++;
		data.push_back(Gridvertex(right, (float)j * y_length, -1.0f, 1.0f));
		index.push_back(index_value); index_value++;
	}
	//TOP BOTTOM
	float bottom = bottom_right_j;
	float top = top_right_j;
	for (unsigned int i = innerRange.begin[0] - 1; i <= (unsigned int)innerRange.end[0]; i++)
	{
		data.push_back(Gridvertex((float)i * x_length, bottom, -1.0f, 1.0f));
		index.push_back(index_value); index_value++;
		data.push_back(Gridvertex((float)i * x_length, top, -1.0f, 1.0f));
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