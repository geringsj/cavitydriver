#include "CavityPainter.hpp"

CavityPainter::CavityPainter()
{
	m_show_grid = P;
}
CavityPainter::~CavityPainter()
{
}

void CavityPainter::paint(unsigned int window_width, unsigned int window_height)
{
	m_window_width = window_width;
	m_window_height = window_height;

	GLFWwindow* window;

    /* Initialize the library */
    if (!glfwInit())

    /* Create a windowed mode window and its OpenGL context */
	glfwWindowHint(GLFW_SAMPLES, 8);
    window = glfwCreateWindow(window_width, window_height, "Lasagne", NULL, NULL);
    if (!window)
    {
		//std::cout<<"Couldn't create glfw window."<<std::endl;
        glfwTerminate();
    }

	/*	Initialize glew */
	//glewExperimental = GL_TRUE;
	GLenum error = glewInit();
	if( GLEW_OK != error)
	{
		std::cout<<"-----\n"
				<<"The time is out of joint - O cursed spite,\n"
				<<"That ever I was born to set it right!\n"
				<<"-----\n"
				<<"Error: "<<glewGetErrorString(error);
	}
	/* Apparently glewInit() causes a GL ERROR 1280, so let's just catch that... */
	glGetError();

	/* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window))
    {
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		int width, height;
		glfwGetFramebufferSize(window, &width, &height);
		glViewport(0, 0, width, height);
		
		drawGridOverlay();

        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events */
        glfwPollEvents();
    }

	//TODO cleanup
}

bool CavityPainter::createGrid(Range p, Range u, Range v)
{
	unsigned int p_size, u_size, v_size;
	Gridvertex* p_data = createSingleGrid(p, p_size);
	Gridvertex* u_data = createSingleGrid(u, u_size);
	Gridvertex* v_data = createSingleGrid(v, v_size);
}

Gridvertex* CavityPainter::createSingleGrid(Range innerRange, unsigned int& data_size)
{
	std::vector<Gridvertex> data;
	unsigned int bottom_left_i = innerRange.begin[0] - 1;
	unsigned int bottom_left_j = innerRange.begin[1] - 1;

	unsigned int bottom_right_i = innerRange.end[0] + 1;
	unsigned int bottom_right_j = innerRange.begin[1] - 1;

	unsigned int top_right_i = innerRange.end[0] + 1;
	unsigned int top_right_j = innerRange.end[1] + 1;

	unsigned int top_left_i = innerRange.begin[0] - 1;
	unsigned int top_left_j = innerRange.end[1] + 1;

	data.push_back(Gridvertex((float)bottom_left_i, (float)bottom_left_j, -1.0f, 1.0f));
	data.push_back(Gridvertex((float)bottom_right_i, (float)bottom_right_j, -1.0f, 1.0f));
	data.push_back(Gridvertex((float)top_right_i, (float)top_right_j, -1.0f, 1.0f));
	data.push_back(Gridvertex((float)top_left_i, (float)top_left_j, -1.0f, 1.0f));

	data.push_back(Gridvertex((float)bottom_left_i, (float)bottom_left_j, -1.0f, 1.0f));
	data.push_back(Gridvertex((float)top_left_i, (float)top_left_j, -1.0f, 1.0f));
	data.push_back(Gridvertex((float)bottom_right_i, (float)bottom_right_j, -1.0f, 1.0f));
	data.push_back(Gridvertex((float)top_right_i, (float)top_right_j, -1.0f, 1.0f));

	//LEFT RIGHT
	unsigned int left = bottom_left_i;
	unsigned int right = bottom_right_i;
	for (unsigned int j = bottom_left_j; j <= top_left_j; j++)
	{
		data.push_back(Gridvertex((float)left, (float)j, -1.0f, 1.0f));
		data.push_back(Gridvertex((float)right, (float)j, -1.0f, 1.0f));
	}
	//TOP BOTTOM
	unsigned int bottom = bottom_right_j;
	unsigned int top = top_left_j;
	for (unsigned int i = bottom_left_i; i <= bottom_right_i; i++)
	{
		data.push_back(Gridvertex((float)i, (float)bottom, -1.0f, 1.0f));
		data.push_back(Gridvertex((float)i, (float)top, -1.0f, 1.0f));
	}

	data_size = (unsigned int)data.size();
	return data.data();
}

void CavityPainter::drawGridOverlay()
{
	switch (m_show_grid)
	{
	case CavityPainter::P:
		break;
	case CavityPainter::U:
		break;
	case CavityPainter::V:
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