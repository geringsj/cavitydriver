#include "CavityPainter.hpp"

CavityPainter::CavityPainter()
{
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
		
        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events */
        glfwPollEvents();
    }

	//TODO cleanup
}

void CavityPainter::drawGridOverlay()
{
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