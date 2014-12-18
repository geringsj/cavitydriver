#include "CavityPainter.hpp"

CavityPainter::CavityPainter()
{
	m_show_grid = P;
}
CavityPainter::~CavityPainter()
{
}

bool CavityPainter::init(unsigned int window_width, unsigned int window_height)
{
	m_window_width = window_width;
	m_window_height = window_height;

	/* Initialize the library */
	if (!glfwInit()) return false;

	/* Create a windowed mode window and its OpenGL context */
	glfwWindowHint(GLFW_SAMPLES, 8);
	m_window = glfwCreateWindow(window_width, window_height, "Lasagne", NULL, NULL);
	if (!m_window)
	{
		//std::cout<<"Couldn't create glfw window."<<std::endl;
		glfwTerminate();
		return false;
	}
	glfwMakeContextCurrent(m_window);

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

	m_cam_sys = CameraSystem(glm::vec3(0.0f, 0.0f, 80.0f), glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(0.0f, 0.0f, -1.0f), glm::vec3(1.0f, 0.0f, 0.0f));
	return true;
}

void CavityPainter::paint()
{

	/* Loop until the user closes the window */
    while (!glfwWindowShouldClose(m_window))
    {
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		int width, height;
		glfwGetFramebufferSize(m_window, &width, &height);
		glViewport(0, 0, width, height);
		
		drawGridOverlay();

        /* Swap front and back buffers */
        glfwSwapBuffers(m_window);

        /* Poll for and process events */
        glfwPollEvents();
    }

	//TODO cleanup
	glfwTerminate();
}

bool CavityPainter::createGrids(Range p, Range u, Range v)
{
	unsigned int p_size, u_size, v_size;
	std::vector<unsigned int> p_index, u_index, v_index;
	std::vector<Gridvertex> p_data, u_data, v_data;
	createSingleGrid(p, p_size, p_index, p_data);
	createSingleGrid(u, u_size, u_index, u_data);
	createSingleGrid(v, v_size, v_index, v_data);

	float right_shift = (float)((p.end[0] + 1) - (p.begin[0] - 1)) / 2.0f;
	float up_shift = (float)((p.end[1] + 1) - (p.begin[1] - 1)) / 2.0f;

	m_cam_sys.Translation(glm::vec3(1.0f, 0.0f, 0.0f), right_shift);
	m_cam_sys.Translation(glm::vec3(0.0f, 1.0f, 0.0f), up_shift);

	m_p_grid.bufferDataFromArray(p_data.data(), p_index.data(), 
		(GLsizei)(p_size*sizeof(Gridvertex)), (GLsizei)(p_index.size()*sizeof(unsigned int)), 
		GL_LINES);
	m_p_grid.setVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(Gridvertex), 0);

	m_u_grid.bufferDataFromArray(u_data.data(), u_index.data(), 
		(GLsizei)(u_size*sizeof(Gridvertex)), (GLsizei)(u_index.size()*sizeof(unsigned int)), 
		GL_LINES);
	m_u_grid.setVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(Gridvertex), 0);

	m_v_grid.bufferDataFromArray(v_data.data(), v_index.data(), 
		(GLsizei)(v_size*sizeof(Gridvertex)), (GLsizei)(v_index.size()*sizeof(unsigned int)), 
		GL_LINES);
	m_v_grid.setVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(Gridvertex), 0);

	m_grid_prgm.init();

	std::string grid_vertex = readShaderFile("./shader/gridVertex.glsl");
	if (!m_grid_prgm.compileShaderFromString(&grid_vertex, GL_VERTEX_SHADER)) { std::cout << m_grid_prgm.getLog(); return false; };

	std::string grid_fragment = readShaderFile("./shader/gridFragment.glsl");
	if (!m_grid_prgm.compileShaderFromString(&grid_fragment, GL_FRAGMENT_SHADER)) { std::cout << m_grid_prgm.getLog(); return false; };

	m_grid_prgm.bindAttribLocation(0, "in_position");

	m_grid_prgm.link();

	return true;
}

void CavityPainter::createSingleGrid(Range innerRange, unsigned int& data_size, std::vector<unsigned int>& index, std::vector<Gridvertex>& data)
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

	data_size = (unsigned int)data.size();
}

void CavityPainter::drawGridOverlay()
{
	m_grid_prgm.use();
	glm::mat4 proj_mat = glm::perspective(45.0f, (float)m_window_width / (float)m_window_height, 0.1f, 100.0f);
	glm::mat4 model_mat = glm::mat4(1.0f);
	glm::mat4 view_mat = m_cam_sys.GetViewMatrix();
	glm::mat4 mvp_mat = proj_mat * view_mat * model_mat;
	m_grid_prgm.setUniform("mvp_matrix", mvp_mat);
	switch (m_show_grid)
	{
	case CavityPainter::P:
		m_p_grid.draw();
		break;
	case CavityPainter::U:
		m_u_grid.draw();
		break;
	case CavityPainter::V:
		m_v_grid.draw();
		break;
	default:
		break;
	}
}

void CavityPainter::drawField()
{
}

void CavityPainter::drawLIC()
{
}

void CavityPainter::drawGeometry()
{
}

const std::string CavityPainter::readShaderFile(const char* const path)
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