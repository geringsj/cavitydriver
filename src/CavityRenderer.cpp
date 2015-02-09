#include "CavityRenderer.hpp"

#include <iostream>

double CavityRenderer::last_mouse_x = 0.0;
double CavityRenderer::last_mouse_y = 0.0;

bool CavityRenderer::FieldLayer::createResources(SimulationParameters& simparams)
{
	// Create field quad
	float x_length = (float)simparams.xLength / (float)simparams.iMax; m_dx = x_length;
	float y_length = (float)simparams.yLength / (float)simparams.jMax; m_dy = y_length;


	float left_i = 0.0f;
	float right_i = left_i +(simparams.xCells + 2) * x_length;

	float bottom_j =  0.0f;
	float top_j = bottom_j + (simparams.yCells + 2) * y_length;
	
	std::array< VertexUV, 4 > vertex_array = {{ VertexUV(left_i,bottom_j,-1.0,0.0,0.0),
											VertexUV(left_i,top_j,-1.0,0.0,1.0),
											VertexUV(right_i,top_j,-1.0,1.0,1.0),
											VertexUV(right_i,bottom_j,-1.0,1.0,0.0) }};

	std::array< GLuint, 6 > index_array = {{ 0,2,1,2,0,3 }};

	if(!(m_field_quad.bufferDataFromArray(vertex_array.data(),index_array.data(),sizeof(VertexUV)*4,sizeof(GLuint)*6,GL_TRIANGLES))) return false;
	m_field_quad.setVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VertexUV), 0);
	m_field_quad.setVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(VertexUV), (GLvoid*) (sizeof(float)*3));

	// Create grid mesh for ibfv
	std::vector<VertexUV> ibfv_vertex_array;
	std::vector<GLuint> ibfv_index_array;

	int base_subdivisions = std::min(simparams.iMax,simparams.jMax)/2;
	int iSubdivisions, jSubdivisions;

	// Build vertex array
	float iStep, jStep;
	if( right_i > top_j)
	{
		jSubdivisions =  base_subdivisions;
		iSubdivisions = (int)std::floor(base_subdivisions * (right_i/top_j));
	}
	else
	{
		iSubdivisions = base_subdivisions;
		jSubdivisions = (int)std::floor(base_subdivisions * (top_j/right_i));
	}

	iStep = 1.0f/(float)iSubdivisions;
	jStep = 1.0f/(float)jSubdivisions;

	for(float j = 0.0; j<1.0+(jStep*0.1); j=j+jStep) // use (jStep*0.1) as temporary epsilon
	{
		for(float i = 0.0; i<1.0+(iStep*0.1); i=i+iStep)
		{
			ibfv_vertex_array.push_back(VertexUV(i,j,-1.0,i/1.0,j/1.0));
		}
	}

	// Build index array
	for(GLuint j=0; j<static_cast<GLuint>(jSubdivisions); j++)
	{
		for(GLuint i=0; i<static_cast<GLuint>(iSubdivisions); i++)
		{
			ibfv_index_array.push_back(i + j*(iSubdivisions+1));
			ibfv_index_array.push_back(i+iSubdivisions+1 + j*(iSubdivisions+1));
			ibfv_index_array.push_back(i+1 + j*(iSubdivisions+1));

			ibfv_index_array.push_back(i+1 + j*(iSubdivisions+1));
			ibfv_index_array.push_back(i+iSubdivisions+1 + j*(iSubdivisions+1));
			ibfv_index_array.push_back(i+1+iSubdivisions+1 + j*(iSubdivisions+1));
			
		}
	}

	if(!(m_ibfv_grid.bufferDataFromArray(ibfv_vertex_array.data(),ibfv_index_array.data(),sizeof(VertexUV)*ibfv_vertex_array.size(),sizeof(GLuint)*ibfv_index_array.size(),GL_TRIANGLES))) return false;
	m_ibfv_grid.setVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VertexUV), 0);
	m_ibfv_grid.setVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(VertexUV), (GLvoid*) (sizeof(float)*3));

	// Fullscreen quad
	std::array< VertexUV, 4 > fullscreenQuad_va = {{ VertexUV(-1.0,-1.0,-1.0,0.0,0.0),
											VertexUV(-1.0,1.0,-1.0,0.0,1.0),
											VertexUV(1.0,1.0,-1.0,1.0,1.0),
											VertexUV(1.0,-1.0,-1.0,1.0,0.0) }};
	std::array< GLuint, 6 > fullscreenQuad_ia = {{ 0,2,1,2,0,3 }};

	if(!(m_fullscreen_quad.bufferDataFromArray(fullscreenQuad_va.data(),fullscreenQuad_ia.data(),sizeof(VertexUV)*4,sizeof(GLuint)*6,GL_TRIANGLES))) return false;
	m_fullscreen_quad.setVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VertexUV), 0);
	m_fullscreen_quad.setVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(VertexUV), (GLvoid*) (sizeof(float)*3));

	// Dye injection sprite
	float blob_size = (1.0f/(float)simparams.xLength);
	blob_size = 0.05f;
	std::array< VertexUV, 4 > dyeBlobSprite_va = {{ VertexUV(-blob_size,-blob_size,-1.0,0.0,0.0),
											VertexUV(-blob_size,blob_size,-1.0,0.0,1.0),
											VertexUV(blob_size,blob_size,-1.0,1.0,1.0),
											VertexUV(blob_size,-blob_size,-1.0,1.0,0.0) }};
	std::array< GLuint, 6 > dyeBlobSprite_ia = {{ 0,2,1,2,0,3 }};

	if(!(m_dye_blob.bufferDataFromArray(dyeBlobSprite_va.data(),dyeBlobSprite_ia.data(),sizeof(VertexUV)*4,sizeof(GLuint)*6,GL_TRIANGLES))) return false;
	m_dye_blob.setVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VertexUV), 0);
	m_dye_blob.setVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(VertexUV), (GLvoid*) (sizeof(float)*3));

	// Shader program for rendering field quad
	m_prgm.init();
	
	std::string field_vertex = Renderer::IO::readShaderFile("./shader/fieldVertex.glsl");
	if (!m_prgm.compileShaderFromString(&field_vertex, GL_VERTEX_SHADER))
		{ std::cout << m_prgm.getLog(); return false; };
	
	std::string field_fragment = Renderer::IO::readShaderFile("./shader/fieldFragment.glsl");
	if (!m_prgm.compileShaderFromString(&field_fragment, GL_FRAGMENT_SHADER))
		{ std::cout << m_prgm.getLog(); return false; };
	
	m_prgm.bindAttribLocation(0, "v_position");
	m_prgm.bindAttribLocation(1, "v_uv");
	
	if (!m_prgm.link()) { std::cout << m_prgm.getLog(); return false; };

	// Shader program for image based flow visualition
	m_ibfvAdvection_prgm.init();
	
	std::string ibfv_vertex = Renderer::IO::readShaderFile("./shader/ibfvAdvectionVertex.glsl");
	if (!m_ibfvAdvection_prgm.compileShaderFromString(&ibfv_vertex, GL_VERTEX_SHADER))
		{ std::cout << m_ibfvAdvection_prgm.getLog(); return false; };
	
	std::string ibfv_fragment = Renderer::IO::readShaderFile("./shader/ibfvAdvectionFragment.glsl");
	if (!m_ibfvAdvection_prgm.compileShaderFromString(&ibfv_fragment, GL_FRAGMENT_SHADER))
		{ std::cout << m_ibfvAdvection_prgm.getLog(); return false; };
	
	m_ibfvAdvection_prgm.bindAttribLocation(0, "v_position");
	m_ibfvAdvection_prgm.bindAttribLocation(1, "v_uv");
	
	if (!m_ibfvAdvection_prgm.link())
		{ std::cout << m_ibfvAdvection_prgm.getLog(); return false; };

	m_ibfvMerge_prgm.init();
	
	std::string ibfvMerge_vertex = Renderer::IO::readShaderFile("./shader/ibfvMergeVertex.glsl");
	if (!m_ibfvMerge_prgm.compileShaderFromString(&ibfvMerge_vertex, GL_VERTEX_SHADER))
		{ std::cout << m_ibfvMerge_prgm.getLog(); return false; };
	
	std::string ibfvMerge_fragment = Renderer::IO::readShaderFile("./shader/ibfvMergeFragment.glsl");
	if (!m_ibfvMerge_prgm.compileShaderFromString(&ibfvMerge_fragment, GL_FRAGMENT_SHADER))
		{ std::cout << m_ibfvMerge_prgm.getLog(); return false; };
	
	m_ibfvMerge_prgm.bindAttribLocation(0, "v_position");
	m_ibfvMerge_prgm.bindAttribLocation(1, "v_uv");
	
	if (!m_ibfvMerge_prgm.link())
		{ std::cout << m_ibfvMerge_prgm.getLog(); return false; };

	// Shader program for dye injection
	m_dyeInjection_prgm.init();
	std::string dyeInjection_vertex = Renderer::IO::readShaderFile("./shader/dyeInjectionVertex.glsl");
	if (!m_dyeInjection_prgm.compileShaderFromString(&dyeInjection_vertex, GL_VERTEX_SHADER))
		{ std::cout << m_dyeInjection_prgm.getLog(); return false; };
	
	std::string dyeInjection_fragment = Renderer::IO::readShaderFile("./shader/dyeInjectionFragment.glsl");
	if (!m_dyeInjection_prgm.compileShaderFromString(&dyeInjection_fragment, GL_FRAGMENT_SHADER))
		{ std::cout << m_dyeInjection_prgm.getLog(); return false; };
	
	m_dyeInjection_prgm.bindAttribLocation(0, "v_position");
	m_dyeInjection_prgm.bindAttribLocation(1, "v_uv");
	
	if (!m_dyeInjection_prgm.link())
		{ std::cout << m_dyeInjection_prgm.getLog(); return false; };

	// Shader program for picking within the field
	m_fieldPicking_prgm.init();
	std::string fieldPicking_vertex = Renderer::IO::readShaderFile("./shader/fieldPickingVertex.glsl");
	if(!m_fieldPicking_prgm.compileShaderFromString(&fieldPicking_vertex, GL_VERTEX_SHADER))
		{ std::cout << m_fieldPicking_prgm.getLog(); return false; };
	
	std::string fieldPicking_fragment = Renderer::IO::readShaderFile("./shader/fieldPickingFragment.glsl");
	if (!m_fieldPicking_prgm.compileShaderFromString(&fieldPicking_fragment, GL_FRAGMENT_SHADER))
		{ std::cout << m_fieldPicking_prgm.getLog(); return false; };
	
	m_fieldPicking_prgm.bindAttribLocation(0, "v_position");
	m_fieldPicking_prgm.bindAttribLocation(1, "v_uv");
	
	if (!m_fieldPicking_prgm.link())
		{ std::cout << m_fieldPicking_prgm.getLog(); return false; };

	// Shader program for streamlines
	m_streamline_prgm.init();
	std::string streamline_vertex = Renderer::IO::readShaderFile("./shader/streamlinesVertex.glsl");
	if (!m_streamline_prgm.compileShaderFromString(&streamline_vertex, GL_VERTEX_SHADER))
		{ std::cout << m_fieldPicking_prgm.getLog(); return false; };

	std::string streamline_fragment = Renderer::IO::readShaderFile("./shader/streamlinesFragment.glsl");
	if (!m_streamline_prgm.compileShaderFromString(&streamline_fragment, GL_FRAGMENT_SHADER))
		{ std::cout << m_fieldPicking_prgm.getLog(); return false; };

	m_streamline_prgm.bindAttribLocation(0, "v_position");
	

	// Create framebuffers "ping-pong" ibfv rendering
	m_ibfv_fbo0 = std::make_shared<FramebufferObject>(10*(float)simparams.iMax,10*(float)simparams.jMax,false,false);
	m_ibfv_fbo0->createColorAttachment(GL_RGBA32F,GL_RGBA,GL_FLOAT);
	m_ibfv_fbo0->bindColorbuffer(0);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	m_ibfv_fbo1 = std::make_shared<FramebufferObject>(10*(float)simparams.iMax,10*(float)simparams.jMax,false,false);
	m_ibfv_fbo1->createColorAttachment(GL_RGBA32F,GL_RGBA,GL_FLOAT);
	m_ibfv_fbo1->bindColorbuffer(0);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);

	// Create empty texture for holding ibfv background/source image
	m_ibfvBackground_tx = std::make_shared<Texture2D>("m_ibfvBackground_tx",GL_RGB, 512, 512, GL_RGB, GL_FLOAT, nullptr);

	// Load mask texture for dye injection blob
	unsigned long begin_pos;
	int x_dim, y_dim;
	char* m_img_data;
	Renderer::IO::readPpmHeader("dye_blob.ppm", begin_pos, x_dim, y_dim);
	m_img_data = new char[x_dim * y_dim * 3];
	Renderer::IO::readPpmData("dye_blob.ppm", m_img_data, begin_pos, x_dim, y_dim);
	m_dyeBlob_tx = std::make_shared<Texture2D>("m_dyeBlob_tx",GL_RGB, x_dim, y_dim, GL_RGB, GL_UNSIGNED_BYTE, m_img_data);
	delete[] m_img_data;

	return true;
}

void CavityRenderer::FieldLayer::draw(CameraSystem& camera)
{
	if(m_show)
	{
		m_prgm.use();

		glEnable(GL_TEXTURE_2D);
		glActiveTexture(GL_TEXTURE0);
		m_prgm.setUniform("field_tx2D",0);
		m_field_tx->bindTexture();

		glActiveTexture(GL_TEXTURE1);
		m_prgm.setUniform("ibfv_tx2D",1);
		m_ibfv_fbo1->bindColorbuffer(0);

		m_prgm.setUniform("min_values",m_field_min_values[m_current_field]);
		m_prgm.setUniform("max_values",m_field_max_values[m_current_field]);

		m_prgm.setUniform("mode",m_display_mode);

		glm::mat4 proj_mat = camera.GetProjectionMatrix();
		glm::mat4 model_mat = glm::mat4(1.0f);
		glm::mat4 view_mat = camera.GetViewMatrix();
		glm::mat4 mvp_mat = proj_mat * view_mat * model_mat;
		m_prgm.setUniform("mvp_matrix", mvp_mat);

		m_field_quad.draw();

		if(m_show_streamlines)
		{
			//TODO draw streamlines
			m_streamline_prgm.use();
			glm::vec4 col = glm::vec4(m_stream_colour[0], m_stream_colour[1], m_stream_colour[2], 1.0);
			m_streamline_prgm.setUniform("mvp_matrix", mvp_mat);
			m_streamline_prgm.setUniform("colour", col);
			for (auto& streamline : m_streamlines)
			{
				streamline->draw();
			}
		}
	}
}

void CavityRenderer::FieldLayer::drawFieldPicking(CameraSystem& camera)
{
	if(m_show)
	{
		m_fieldPicking_prgm.use();

		glm::mat4 proj_mat = camera.GetProjectionMatrix();
		glm::mat4 model_mat = glm::mat4(1.0f);
		glm::mat4 view_mat = camera.GetViewMatrix();
		glm::mat4 mvp_mat = proj_mat * view_mat * model_mat;
		m_fieldPicking_prgm.setUniform("mvp_matrix", mvp_mat);

		m_field_quad.draw();
	}
}

void CavityRenderer::FieldLayer::setFieldData(std::string path)
{
	// variables for field/image attributes
	int dim_x;
	int dim_y;

	// construct first filename
	unsigned int i = 0;
	std::string filename = path;
	filename.append(std::to_string(i));
	filename.append(".pfm");

	std::ifstream file(filename);

	while(file.is_open())
	{
		file.close();

		unsigned long header_end_pos;
		if(!Renderer::IO::readPfmHeader(filename.c_str(), header_end_pos, dim_x, dim_y)) {std::cout<<"Failed to open pfm header of file \""<<filename<<"\""<<std::endl;}

		//TODO do something if image size stays not the same
		m_field_data.push_back(std::vector<float>( dim_x * dim_y * 3));

		if(!Renderer::IO::readPfmData(filename.c_str(), m_field_data.back().data(), header_end_pos)) {std::cout<<"Failed to read pfm data of file \""<<filename<<"\""<<std::endl;}

		//TODO compute min/max
		m_field_max_values.push_back(glm::vec3(std::numeric_limits<float>::min()));
		m_field_min_values.push_back(glm::vec3(std::numeric_limits<float>::max()));

		//m_field_max_values.push_back(glm::vec3(0.0));
		//m_field_min_values.push_back(glm::vec3(1.0));

		unsigned int index = 0;
		for(int x=0; x<dim_x; x++)
		{
			for(int y=0; y<dim_y; y++)
			{
				for(int c=0; c<3; c++)
				{
					m_field_max_values.back()[c] = std::max(m_field_max_values.back()[c],m_field_data.back()[index]);
					m_field_min_values.back()[c] = std::min(m_field_max_values.back()[c],m_field_data.back()[index]);
					index++;
				}
			}
		}

		i++;
		filename = path;
		filename.append(std::to_string(i));
		filename.append(".pfm");
		file.open(filename);
	}

	m_field_dimension.i = dim_x;
	m_field_dimension.j = dim_y;
	m_num_fields = m_field_data.size();
}

bool CavityRenderer::FieldLayer::setFieldTexture(unsigned int requested_frame)
{
	if(requested_frame < m_field_data.size())
	{
		m_current_field = requested_frame;

		m_field_tx = std::make_shared<Texture2D>( "m_field_tx",
												GL_RGB32F,
												m_field_dimension.i,
												m_field_dimension.j,
												GL_RGB, GL_FLOAT,
												m_field_data[requested_frame].data());

		m_field_tx->texParameteri(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		m_field_tx->texParameteri(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

		return true;
	}

	return false;
}

void CavityRenderer::FieldLayer::updateFieldTexture(double current_time)
{
	if(m_play_animation)
	{
		m_elapsed_time += current_time - m_time_tracker;

		if(m_elapsed_time > m_requested_frametime)
		{
			m_current_field = (m_current_field + (int)std::floor(m_elapsed_time / m_requested_frametime)) % m_field_data.size();
			m_elapsed_time = 0.0;
		}

		//Streamlines begin
		int dim_x = m_field_dimension.i;
		int dim_y = m_field_dimension.j;
		Real cell_size = std::sqrt(m_dx*m_dx + m_dy*m_dy);
		Real domain_size_x = (Real)dim_x * m_dx;
		Real domain_size_y = (Real)dim_y * m_dy;

		// build streamlines
		Real start_x = m_dx;
		Real start_y = 0.0;
		unsigned int index = 0;
		for (auto& seedpoint : m_streamline_seedpoints)
		{
			std::vector<PointVertex> points;
			std::vector<unsigned int> points_index;
			unsigned int vi_index = 0;
			while (0.0 < seedpoint.x && seedpoint.x < domain_size_x &&
				0.0 < seedpoint.y && seedpoint.y < domain_size_y)
			{
				Point p = seedpoint;
				Real u, v;
				interpolateUV(p, u, v);

				if (!((u + v)>0.0))
					break;

				Real scaling = cell_size / std::sqrt(u*u + v*v);
				PointVertex p_new(p.x + scaling*u, p.y + scaling * v, -1.0f);
				points.push_back(p_new);
				points_index.push_back(vi_index);
				vi_index++;

				if ((int)points.size() > std::max(dim_x, dim_y))
					break;
			}
			m_streamlines[index]->bufferDataFromArray(points.data(), points_index.data(),sizeof(PointVertex)*points.size(),sizeof(unsigned int) * points_index.size(),GL_LINE_STRIP);
			index++;
		}
		// Streamlines end
	}
	else
	{
		m_elapsed_time = 0.0;
	}

	m_time_tracker = current_time;

	//TODO use reload function
	m_field_tx = std::make_shared<Texture2D>( "m_field_tx",
											GL_RGB32F,
											m_field_dimension.i,
											m_field_dimension.j,
											GL_RGB, GL_FLOAT,
											m_field_data[m_current_field].data());

	m_field_tx->texParameteri(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	m_field_tx->texParameteri(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);


	// Generate noise texture
	std::mt19937 rg(std::chrono::system_clock::now().time_since_epoch().count());
	float scale = 1.0f/(float)rg.max();
	int size = m_ibfvBackground_tx->getWidth() * m_ibfvBackground_tx->getHeight();
	GLfloat* noise_data = new GLfloat[size * 3];
	for(int i=0;i<size;i++)
	{
		float value = rg()*scale;
		noise_data[i*3] = value;
		noise_data[i*3 +1] = value;
		noise_data[i*3 +2] = value;
	}
	m_ibfvBackground_tx->reload(m_ibfvBackground_tx->getWidth(),m_ibfvBackground_tx->getHeight(),noise_data);
	delete[] noise_data;

	m_ibfv_fbo0->bind();
	//glClearColor(0.0f,0.0f,0.0f,0.0f);
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, m_ibfv_fbo0->getWidth(), m_ibfv_fbo0->getHeight());

	m_ibfvAdvection_prgm.use();

	glEnable(GL_TEXTURE_2D);
	glActiveTexture(GL_TEXTURE0);
	m_ibfvAdvection_prgm.setUniform("field_tx2D",0);
	m_field_tx->bindTexture();
	glActiveTexture(GL_TEXTURE1);
	m_ibfvAdvection_prgm.setUniform("previous_frame_tx2d",1);
	m_ibfv_fbo1->bindColorbuffer(0);

	m_ibfv_grid.draw();

	m_ibfv_fbo1->bind();
	//glClearColor(0.0f,0.0f,0.0f,0.0f);
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, m_ibfv_fbo1->getWidth(), m_ibfv_fbo1->getHeight());

	m_ibfvMerge_prgm.use();

	glActiveTexture(GL_TEXTURE0);
	m_ibfvMerge_prgm.setUniform("background_tx2D",0);
	m_ibfvBackground_tx->bindTexture();
	glActiveTexture(GL_TEXTURE1);
	m_ibfvMerge_prgm.setUniform("advected_tx2D",1);
	m_ibfv_fbo0->bindColorbuffer(0);

	m_fullscreen_quad.draw();

	m_dyeInjection_prgm.use();
	glActiveTexture(GL_TEXTURE0);
	m_dyeInjection_prgm.setUniform("dyeBlob_tx2D",0);
	m_dyeBlob_tx->bindTexture();

	//TODO upload position matrices
	std::string uniform_name;
	for(int instance=0; instance<m_dye_seedpoints.size(); instance++)
	{
		uniform_name = "dye_location[" + std::to_string(instance) + "]";
		m_dyeInjection_prgm.setUniform(uniform_name.c_str(),glm::vec2((float)m_dye_seedpoints[instance].x,
																(float)m_dye_seedpoints[instance].y));
	}

	//TODO instanced draw call
	m_dye_blob.draw(m_dye_seedpoints.size());
}

void CavityRenderer::FieldLayer::addDyeSeedpoint(float x, float y)
{
	m_dye_seedpoints.push_back(Point(x,y,0.0));
}

void CavityRenderer::FieldLayer::clearDye()
{
	m_dye_seedpoints.clear();
}

void CavityRenderer::FieldLayer::interpolateUV(Point p, Real& u, Real& v)
{
	// step1: find grid cell
	int i = (int)(p.x / m_dx) + 1;
	int j = (int)((p.y + (m_dy / 2.0)) / m_dy) + 1;

	// step2: find neighbour cells
	Real x_1 = (Real)(i - 1) * m_dx;
	Real x_2 = (Real)i * m_dx;
	Real y_1 = (Real)(j - 1)*m_dy - (m_dy / 2.0);
	Real y_2 = (Real)(j)*m_dy - (m_dy / 2.0);

	// step3: bilinear interpolation
	u = 1.0 / (m_dx*m_dy) * ((x_2 - p.x)*(y_2 - p.y) * m_field_data[m_current_field][(i - 1) * 3 + (j - 1)*m_field_dimension.i] +
		(p.x - x_1)*(y_2 - p.y) * m_field_data[m_current_field][(i)* 3 + (j - 1)*m_field_dimension.i] +
		(x_2 - p.x)*(p.y - y_1) * m_field_data[m_current_field][(i - 1) * 3 + (j)*m_field_dimension.i] +
		(p.x - x_1)*(p.y - y_1) * m_field_data[m_current_field][(i)* 3 + (j)*m_field_dimension.i]);

	v = 1.0 / (m_dx*m_dy) * ((x_2 - p.x)*(y_2 - p.y) * m_field_data[m_current_field][((i - 1) * 3) + 1 + (j - 1)*m_field_dimension.i] +
		(p.x - x_1)*(y_2 - p.y) * m_field_data[m_current_field][((i)* 3) + 1 + (j - 1)*m_field_dimension.i] +
		(x_2 - p.x)*(p.y - y_1) * m_field_data[m_current_field][((i - 1) * 3) + 1 + (j)*m_field_dimension.i] +
		(p.x - x_1)*(p.y - y_1) * m_field_data[m_current_field][((i)* 3) + 1 + (j)*m_field_dimension.i]);
}

void CavityRenderer::FieldLayer::addStreamlineSeedpoint(float x, float y)
{
	m_streamline_seedpoints.push_back(Point(x, y, 0.0));
	m_streamlines.push_back(std::make_shared<Mesh> ());
}

void CavityRenderer::FieldLayer::clearStreamline()
{
	m_streamlines.clear();
}

bool CavityRenderer::OverlayGridLayer::createResources(SimulationParameters& simparams)
{
	updateGridMesh(simparams);

	m_prgm.init();
	
	std::string grid_vertex = Renderer::IO::readShaderFile("./shader/gridVertex.glsl");
	if (!m_prgm.compileShaderFromString(&grid_vertex, GL_VERTEX_SHADER)) { std::cout << m_prgm.getLog(); return false; };
	
	std::string grid_fragment = Renderer::IO::readShaderFile("./shader/gridFragment.glsl");
	if (!m_prgm.compileShaderFromString(&grid_fragment, GL_FRAGMENT_SHADER)) { std::cout << m_prgm.getLog(); return false; };
	
	m_prgm.bindAttribLocation(0, "in_position");
	
	m_prgm.link();

	return true;
}

void CavityRenderer::OverlayGridLayer::draw(CameraSystem& camera)
{
	if(m_show)
	{
		m_prgm.use();
		glm::mat4 proj_mat = camera.GetProjectionMatrix();
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
	up_shift = up_shift*right_shift; /* get rid of compiler warning unused variables */

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

	return true;
}

bool CavityRenderer::BoundaryCellsLayer::createResources(SimulationParameters& simparams)
{
	updateCellMesh(simparams);

	m_prgm.init();

	std::string boundary_vertex_shdr = Renderer::IO::readShaderFile("./shader/boundaryCellVertex.glsl");
	if (!m_prgm.compileShaderFromString(&boundary_vertex_shdr, GL_VERTEX_SHADER)) { std::cout << m_prgm.getLog(); return false; };

	std::string boundary_fragment_shdr = Renderer::IO::readShaderFile("./shader/boundaryCellFragment.glsl");
	if (!m_prgm.compileShaderFromString(&boundary_fragment_shdr, GL_FRAGMENT_SHADER)) { std::cout << m_prgm.getLog(); return false; };

	m_prgm.bindAttribLocation(0, "in_position");
	m_prgm.bindAttribLocation(1, "in_uv");

	m_prgm.link();

	unsigned long begin_pos;
	int x_dim, y_dim;
	char* m_img_data;
	Renderer::IO::readPpmHeader("boundary_cell.ppm", begin_pos, x_dim, y_dim);
	m_img_data = new char[x_dim * y_dim * 3];
	Renderer::IO::readPpmData("boundary_cell.ppm", m_img_data, begin_pos, x_dim, y_dim);
	m_cell_tx = std::make_shared<Texture2D>("m_cell_tx",GL_RGB, x_dim, y_dim, GL_RGB, GL_UNSIGNED_BYTE, m_img_data);
	delete[] m_img_data;

	return true;
}

void CavityRenderer::BoundaryCellsLayer::draw(CameraSystem& camera)
{
	if(m_show)
	{
		m_prgm.use();
		
		glEnable(GL_TEXTURE_2D);
		glActiveTexture(GL_TEXTURE0);
		m_prgm.setUniform("boundary_tx2D",0);
		m_cell_tx->bindTexture();

		
		for(auto& cell : m_cell_positions)
		{
			glm::mat4 proj_mat = camera.GetProjectionMatrix();
			glm::mat4 model_mat = glm::translate(glm::mat4(1.0),glm::vec3(cell.x,cell.y,0.0));
			glm::mat4 view_mat = camera.GetViewMatrix();
			glm::mat4 mvp_mat = proj_mat * view_mat * model_mat;
			m_prgm.setUniform("mvp_matrix", mvp_mat);
			
			m_cell.draw();
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

	return true;
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
			float pos[2] = 
				{ float((float)(i) * x_length + x_length / 2.0), float((float)(j) * y_length + y_length/2.0) };
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

bool CavityRenderer::BoundaryGlyphLayer::createResources(SimulationParameters& simparams)
{
	updateGlyphMesh(simparams);

	m_prgm.init();
	
	std::string arrow_vertex = Renderer::IO::readShaderFile("./shader/arrowVertex.glsl");
	if (!m_prgm.compileShaderFromString(&arrow_vertex, GL_VERTEX_SHADER)) { std::cout << m_prgm.getLog(); return false; };
	
	std::string arrow_fragment = Renderer::IO::readShaderFile("./shader/arrowFragment.glsl");
	if (!m_prgm.compileShaderFromString(&arrow_fragment, GL_FRAGMENT_SHADER)) { std::cout << m_prgm.getLog(); return false; };
	
	m_prgm.bindAttribLocation(0, "v_position");
	m_prgm.bindAttribLocation(1, "v_uv");
	
	m_prgm.link();

	unsigned long begin_pos;
	int x_dim, y_dim;
	char* m_img_data;
	Renderer::IO::readPpmHeader("arrow.ppm", begin_pos, x_dim, y_dim);
	m_img_data = new char[x_dim * y_dim * 3];
	Renderer::IO::readPpmData("arrow.ppm", m_img_data, begin_pos, x_dim, y_dim);
	m_velocity_glyph_tx = std::make_shared<Texture2D>("m_velocity_glyph_tx",GL_RGB, x_dim, y_dim, GL_RGB, GL_UNSIGNED_BYTE, m_img_data);
	delete[] m_img_data;

	Renderer::IO::readPpmHeader("boundary_p_dlt_glyph.ppm", begin_pos, x_dim, y_dim);
	m_img_data = new char[x_dim * y_dim * 3];
	Renderer::IO::readPpmData("boundary_p_dlt_glyph.ppm", m_img_data, begin_pos, x_dim, y_dim);
	m_pdlt_glyph_tx = std::make_shared<Texture2D>("m_pdlt_glyph_tx",GL_RGB, x_dim, y_dim, GL_RGB, GL_UNSIGNED_BYTE, m_img_data);
	delete[] m_img_data;

	Renderer::IO::readPpmHeader("boundary_p_nm_glyph.ppm", begin_pos, x_dim, y_dim);
	m_img_data = new char[x_dim * y_dim * 3];
	Renderer::IO::readPpmData("boundary_p_nm_glyph.ppm", m_img_data, begin_pos, x_dim, y_dim);
	m_pnm_glyph_tx = std::make_shared<Texture2D>("m_pnm_glyph_tx",GL_RGB, x_dim, y_dim, GL_RGB, GL_UNSIGNED_BYTE, m_img_data);
	delete[] m_img_data;

	return true;
}

void CavityRenderer::BoundaryGlyphLayer::draw(CameraSystem& camera)
{
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
				m_velocity_glyph_tx->bindTexture();

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

				glm::mat4 proj_mat = camera.GetProjectionMatrix();
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
			float pos[] = { float((float)(i) * x_length/1.0 + x_length / 2.0f), float((float)(j) * y_length/1.0 + y_length/2.0f) };
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

	return true;
}

bool CavityRenderer::GeometryLayer::createResources(SimulationParameters& simparams)
{
	if( 1 || (&simparams + 1)) /* get rid of compiler warning */
	return true;
}

void CavityRenderer::GeometryLayer::draw(CameraSystem& camera)
{
}

bool CavityRenderer::InterfaceLayer::createResources(SimulationParameters& simparams)
{
	std::vector<unsigned int> domainIndicator_ia;
	std::vector<VertexUV> domainIndicator_va;

	// compute size of domain indicators
	float domain_sizeX = (float)simparams.xLength * (1.0f + 2.0f/(float)simparams.iMax);
	float domain_sizeY = (float)simparams.yLength * (1.0f + 2.0f/(float)simparams.jMax); 
	float domainIndicator_size = std::min( std::max(domain_sizeX,domain_sizeY)/100.0f ,
											std::min(domain_sizeX,domain_sizeY) );

	// southwest corner
	domainIndicator_va.push_back(VertexUV(-domainIndicator_size/2.0f,-domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_va.push_back(VertexUV(-domainIndicator_size/2.0f,domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_va.push_back(VertexUV(domainIndicator_size/2.0f,-domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_ia.push_back(0); domainIndicator_ia.push_back(1); // first line
	domainIndicator_ia.push_back(0); domainIndicator_ia.push_back(2); // second line

	// northwest corner
	domainIndicator_va.push_back(VertexUV(-domainIndicator_size/2.0f,domain_sizeY+domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_va.push_back(VertexUV(-domainIndicator_size/2.0f,domain_sizeY-domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_va.push_back(VertexUV(domainIndicator_size/2.0f,domain_sizeY+domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_ia.push_back(3); domainIndicator_ia.push_back(4); // first line
	domainIndicator_ia.push_back(3); domainIndicator_ia.push_back(5); // second line

	// northeast corner
	domainIndicator_va.push_back(VertexUV(domain_sizeX+domainIndicator_size/2.0f,domain_sizeY+domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_va.push_back(VertexUV(domain_sizeX+domainIndicator_size/2.0f,domain_sizeY-domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_va.push_back(VertexUV(domain_sizeX-domainIndicator_size/2.0f,domain_sizeY+domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_ia.push_back(6); domainIndicator_ia.push_back(7); // first line
	domainIndicator_ia.push_back(6); domainIndicator_ia.push_back(8); // second line

	// southeast corner
	domainIndicator_va.push_back(VertexUV(domain_sizeX+domainIndicator_size/2.0f,-domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_va.push_back(VertexUV(domain_sizeX+domainIndicator_size/2.0f,domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_va.push_back(VertexUV(domain_sizeX-domainIndicator_size/2.0f,-domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_ia.push_back(9); domainIndicator_ia.push_back(10); // first line
	domainIndicator_ia.push_back(9); domainIndicator_ia.push_back(11); // second line
	

	if(!m_domainIndicators.bufferDataFromArray(domainIndicator_va.data(),domainIndicator_ia.data(),
		(GLsizei)(domainIndicator_va.size()*sizeof(VertexUV)),(GLsizei)(domainIndicator_ia.size()*sizeof(unsigned int)),GL_LINES))
		return false;
	m_domainIndicators.setVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VertexUV), 0);
	m_domainIndicators.setVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(VertexUV), (GLvoid*) (sizeof(GLfloat)*3));

	// Create shader program
	m_prgm.init();
	
	std::string arrow_vertex = Renderer::IO::readShaderFile("./shader/interfaceVertex.glsl");
	if (!m_prgm.compileShaderFromString(&arrow_vertex, GL_VERTEX_SHADER)) { std::cout << m_prgm.getLog(); return false; };
	
	std::string arrow_fragment = Renderer::IO::readShaderFile("./shader/interfaceFragment.glsl");
	if (!m_prgm.compileShaderFromString(&arrow_fragment, GL_FRAGMENT_SHADER)) { std::cout << m_prgm.getLog(); return false; };
	
	m_prgm.bindAttribLocation(0, "v_position");
	m_prgm.bindAttribLocation(1, "v_uv");
	
	m_prgm.link();

	return true;
}

void CavityRenderer::InterfaceLayer::draw(CameraSystem& camera)
{
	m_prgm.use();
	glm::mat4 proj_mat = camera.GetProjectionMatrix();
	glm::mat4 model_mat = glm::mat4(1.0f);
	glm::mat4 view_mat = camera.GetViewMatrix();
	glm::mat4 mvp_mat = proj_mat * view_mat * model_mat;
	m_prgm.setUniform("mvp_matrix", mvp_mat);
	
	m_domainIndicators.draw();
}

bool CavityRenderer::InterfaceLayer::updateResources(SimulationParameters& simparams)
{
	std::vector<unsigned int> domainIndicator_ia;
	std::vector<VertexUV> domainIndicator_va;

	// compute size of domain indicators
	float domain_sizeX = (float)simparams.xLength * (1.0 + 2.0/(float)simparams.iMax);
	float domain_sizeY = (float)simparams.yLength * (1.0 + 2.0/(float)simparams.jMax); 
	float domainIndicator_size = std::min( std::max(domain_sizeX,domain_sizeY)/100.0f ,
											std::min(domain_sizeX,domain_sizeY) );

	// southwest corner
	domainIndicator_va.push_back(VertexUV(-domainIndicator_size/2.0f,-domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_va.push_back(VertexUV(-domainIndicator_size/2.0f,domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_va.push_back(VertexUV(domainIndicator_size/2.0f,-domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_ia.push_back(0); domainIndicator_ia.push_back(1); // first line
	domainIndicator_ia.push_back(0); domainIndicator_ia.push_back(2); // second line

	// northwest corner
	domainIndicator_va.push_back(VertexUV(-domainIndicator_size/2.0f,domain_sizeY+domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_va.push_back(VertexUV(-domainIndicator_size/2.0f,domain_sizeY-domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_va.push_back(VertexUV(domainIndicator_size/2.0f,domain_sizeY+domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_ia.push_back(3); domainIndicator_ia.push_back(4); // first line
	domainIndicator_ia.push_back(3); domainIndicator_ia.push_back(5); // second line

	// northeast corner
	domainIndicator_va.push_back(VertexUV(domain_sizeX+domainIndicator_size/2.0f,domain_sizeY+domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_va.push_back(VertexUV(domain_sizeX+domainIndicator_size/2.0f,domain_sizeY-domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_va.push_back(VertexUV(domain_sizeX-domainIndicator_size/2.0f,domain_sizeY+domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_ia.push_back(6); domainIndicator_ia.push_back(7); // first line
	domainIndicator_ia.push_back(6); domainIndicator_ia.push_back(8); // second line

	// southeast corner
	domainIndicator_va.push_back(VertexUV(domain_sizeX+domainIndicator_size/2.0f,-domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_va.push_back(VertexUV(domain_sizeX+domainIndicator_size/2.0f,domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_va.push_back(VertexUV(domain_sizeX-domainIndicator_size/2.0f,-domainIndicator_size/2.0f,-1.0,0.0,0.0));
	domainIndicator_ia.push_back(9); domainIndicator_ia.push_back(10); // first line
	domainIndicator_ia.push_back(9); domainIndicator_ia.push_back(11); // second line
	

	if(!m_domainIndicators.bufferDataFromArray(domainIndicator_va.data(),domainIndicator_ia.data(),
		(GLsizei)(domainIndicator_va.size()*sizeof(VertexUV)),(GLsizei)(domainIndicator_ia.size()*sizeof(unsigned int)),GL_LINES))
		return false;
	m_domainIndicators.setVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VertexUV), 0);
	m_domainIndicators.setVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(VertexUV), (GLvoid*) (sizeof(GLfloat)*3));

	return true;
}


/* Static tweak bar callback definitions. Implementation see further down. */
void TW_CALL Bake(void* clientData);

void TW_CALL ClearDye(void* clientData);

void TW_CALL ClearStream(void* clientData);

void TW_CALL CopyStdStringToClient(std::string& destinationClientString, const std::string& sourceLibraryString);

CavityRenderer::CavityRenderer(MTQueue<SimulationParameters>& inbox, MTQueue<SimulationParameters>& outbox)
	: m_inbox(inbox), m_outbox(outbox)
{
}

CavityRenderer::~CavityRenderer()
{
}

bool CavityRenderer::initPainterVis(unsigned int window_width, unsigned int window_height, SimulationParameters& sim_params, std::string fields_filename)
{
	m_window_width = window_width;
	m_window_height = window_height;
	m_window_background_colour[0] = 0.2f; m_window_background_colour[1] = 0.2f; m_window_background_colour[2] = 0.2f;

	m_activeCamera = CameraSystem(
		glm::vec3(0.0f, 0.0f, 1.0),
		glm::vec3(0.0f, 1.0f, 0.0f),
		glm::vec3(0.0f, 0.0f, -1.0f),
		glm::vec3(1.0f, 0.0f, 0.0f));
	m_activeCamera.setAspectRatio((float)window_width/(float)window_height);
	m_activeCamera.accessFieldOfView() = (60.0f * ((float)window_height/(float)window_width));

	m_field_layer.m_show = true;
	m_overlayGrid_layer.m_show = true;
	m_boundaryCells_layer.m_show = true;

	m_simparams = sim_params;

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

	// Init layers
	m_field_layer.setFieldData(fields_filename);
	if(!m_field_layer.createResources(sim_params)) { return false; }
	if(!m_field_layer.setFieldTexture(0)) { return false; }

	if(!m_overlayGrid_layer.createResources(sim_params)) { return false; }

	if(!m_boundaryCells_layer.createResources(sim_params)) { return false; }
	m_boundaryCells_layer.setCellPositions(m_simparams);

	if(!m_boundaryGlyph_layer.createResources(sim_params)) {return false; }
	m_boundaryGlyph_layer.setGlyphs(m_simparams);

	// TODO geometry layer

	if(!m_interface_layer.createResources(sim_params)) { return false; }

	// some additional resource intialization
	m_fieldPicking_fbo = std::make_shared<FramebufferObject>(m_window_width,m_window_height,false,false);
	m_fieldPicking_fbo->createColorAttachment(GL_RGBA32F,GL_RGBA,GL_FLOAT);

	// Init tweak bar entries
	initPainterTweakBar();

	return true;
}

bool CavityRenderer::initBakeryVis(unsigned int window_width, unsigned int window_height, SimulationParameters& sim_params)
{
	m_window_width = window_width;
	m_window_height = window_height;
	m_window_background_colour[0] = 0.2f; m_window_background_colour[1] = 0.2f; m_window_background_colour[2] = 0.2f;

	m_activeCamera = CameraSystem(
		glm::vec3(0.0f, 0.0f, 1.0),
		glm::vec3(0.0f, 1.0f, 0.0f),
		glm::vec3(0.0f, 0.0f, -1.0f),
		glm::vec3(1.0f, 0.0f, 0.0f));
	m_activeCamera.setAspectRatio((float)window_width/(float)window_height);
	m_activeCamera.accessFieldOfView() = (60.0 * ((float)window_height/(float)window_width));

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
	TwCopyStdStringToClientFunc(CopyStdStringToClient);

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

	// Init layers
	if(!m_overlayGrid_layer.createResources(sim_params)) { return false; }

	if(!m_boundaryCells_layer.createResources(sim_params)) { return false; }
	m_boundaryCells_layer.setCellPositions(m_simparams);

	if(!m_boundaryGlyph_layer.createResources(sim_params)) {return false; }
	m_boundaryGlyph_layer.setGlyphs(m_simparams);

	// TODO geometry layer

	if(!m_interface_layer.createResources(sim_params)) { return false; }

	initBakeryTweakBar();

	return true;
}

bool CavityRenderer::createGLSLPrograms()
{
	/* Post processing program */
	m_postProc_prgm.init();

	std::string postProc_vertex_shdr = Renderer::IO::readShaderFile("./shader/postProcVertex.glsl");
	if (!m_postProc_prgm.compileShaderFromString(&postProc_vertex_shdr, GL_VERTEX_SHADER)) { std::cout << m_postProc_prgm.getLog(); return false; };

	std::string postProc_fragment_shdr = Renderer::IO::readShaderFile("./shader/postProcFragment.glsl");
	if (!m_postProc_prgm.compileShaderFromString(&postProc_fragment_shdr, GL_FRAGMENT_SHADER)) { std::cout << m_postProc_prgm.getLog(); return false; };

	m_postProc_prgm.bindAttribLocation(0, "v_position");
	m_postProc_prgm.bindAttribLocation(1, "v_uv");

	if (!m_postProc_prgm.link()) { std::cout << m_postProc_prgm.getLog(); return false; };

	return true; /* return with great success */
}

bool CavityRenderer::createMeshes()
{
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

void CavityRenderer::initBakeryTweakBar()
{
	// Create a tweak bar
	bar = TwNewBar("CavityBaker-Settings");
	TwWindowSize(m_window_width, m_window_height);
	addIntParam("m_window_width", "label='Window width' group='Window' ", &m_window_width, "RO");
	addIntParam("m_window_height", "label='Window height' group='Window' ", &m_window_height, "RO");
	TwAddVarRW(bar, "m_window_background_colour", TW_TYPE_COLOR3F, &m_window_background_colour, " label='Background color' group='Window' ");
	addFloatParam("m_fieldOfView", " step=0.1 label='Field of View' group='Camera' ", &m_activeCamera.accessFieldOfView(), "RW", 1.0f, 180.0f);
	addFloatParam("x", " step=0.01 label='X postion' group='Camera' ", &m_activeCamera.accessCamPos().x, "RW", 0.0f, (float)m_simparams.xLength);
	addFloatParam("y", " step=0.01 label='Y position' group='Camera' ", &m_activeCamera.accessCamPos().y, "RW", 0.0f, (float)m_simparams.yLength);
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
	//TODO
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

}

void CavityRenderer::initPainterTweakBar()
{
	bar = TwNewBar("CavityBaker-Settings");
	TwWindowSize(m_window_width, m_window_height);
	addIntParam("m_window_width", "label='Window width' group='Window' ", &m_window_width, "RO");
	addIntParam("m_window_height", "label='Window height' group='Window' ", &m_window_height, "RO");
	TwAddVarRW(bar, "m_window_background_colour", TW_TYPE_COLOR3F, &m_window_background_colour, " label='Background color' group='Window' ");

	addFloatParam("m_fieldOfView", " step=0.1 label='Field of View' group='Camera' ", &m_activeCamera.accessFieldOfView(), "RW", 1.0f, 180.0f);
	addFloatParam("x", " step=0.01 label='X postion' group='Camera' ", &m_activeCamera.accessCamPos().x, "RW", 0.0f, (float)m_simparams.xLength);
	addFloatParam("y", " step=0.01 label='Y position' group='Camera' ", &m_activeCamera.accessCamPos().y, "RW", 0.0f, (float)m_simparams.yLength);

	addBoolParam("m_show_grid", " label='Show grid' group='Grid' ", &m_overlayGrid_layer.m_show);
	TwAddVarRW(bar, "m_grid_colour", TW_TYPE_COLOR3F, &m_overlayGrid_layer.m_colour, " label='Grid color' group='Grid' ");

	addBoolParam("m_show_field", " label='Show field' group='Field' ", &m_field_layer.m_show);
	addIntParam("m_current_field", " label='Frame' group='Field' ", &m_field_layer.m_current_field, "RW", 0, m_field_layer.m_num_fields);
	addIntParam("m_display_mode", " label='Mode' group='Field' ", &m_field_layer.m_display_mode, "RW", 0, 5);
	addBoolParam("m_play_animation", " label='Play animation' group='Field' ", &m_field_layer.m_play_animation);
	addDoubleParam("m_requested_frametime", " step=0.001 label='Frametime' group='Field' ", &m_field_layer.m_requested_frametime, "RW");
	addButtonParam("m_clearDye", " label='Clear dye' group='Field' ", ClearDye);
	addButtonParam("m_clearStream", "label = 'Clear streamlines' group = 'Field' ", ClearStream);
	addBoolParam("m_dye", "label = 'Place dye' group = 'Field' ", &m_dye, "RW"); 
	addBoolParam("m_stream", "label = 'Place streamline' group = 'Field' ", &m_stream, "RW");
	TwAddVarRW(bar, "m_stream_colour", TW_TYPE_COLOR3F, &m_field_layer.m_stream_colour, " label='Streamline color' group='Field' ");

	addBoolParam("m_show_boundary_cells", " label='Show boundary cells' group='Boundary' ", &m_boundaryCells_layer.m_show);
	addBoolParam("m_show_boundary_glyphs", " label='Show boundary glyphs' group='Boundary' ", &m_boundaryGlyph_layer.m_show);
	addIntParam("m_boundary_glyph_mode", " label='Glyph display mode' group='Boundary' ", &m_boundaryGlyph_layer.m_glyph_mode, "RW", 0, 1);
	
	addIntParam("useComplexGeometry", " label='Scenario' group='Simulation Parameters' ", &m_simparams.useComplexGeometry, "RO", 0, 4);

	Real test; float __float; double __double;
	const char* _double = typeid(__double).name();
	const char* _float = typeid(__float).name();
	if (strcmp(typeid(test).name(), _double) == 0)
	{
		addDoubleParam("alpha", " step=0.1 label='alpha' group='Simulation Parameters' ", &m_simparams.alpha, "RO");
		addDoubleParam("deltaT", " step=0.1 label='deltaT' group='Simulation Parameters' ", &m_simparams.deltaT, "RO");
		addDoubleParam("deltaVec", " step=0.1 label='deltaVec' group='Simulation Parameters' ", &m_simparams.deltaVec, "RO");
		addDoubleParam("eps", " step=0.001 label='eps' group='Simulation Parameters' ", &m_simparams.eps, "RO");
		addDoubleParam("gx", " step=0.1 label='gx' group='Simulation Parameters' ", &m_simparams.gx, "RO");
		addDoubleParam("gy", " step=0.1 label='gy' group='Simulation Parameters' ", &m_simparams.gy, "RO");
		addDoubleParam("KarmanAngle", " step=0.1 label='KarmanAngle' group='Simulation Parameters' ", &m_simparams.KarmanAngle, "RO");
		addDoubleParam("KarmanObjectWidth", " step=0.1 label='KarmanObjectWidth' group='Simulation Parameters' ", &m_simparams.KarmanObjectWidth, "RO");
		addDoubleParam("pi", " step=0.1 label='pi' group='Simulation Parameters' ", &m_simparams.pi, "RO");
		addDoubleParam("re", " step=0.1 label='re' group='Simulation Parameters' ", &m_simparams.re, "RO");
		addDoubleParam("tau", " step=0.1 label='tau' group='Simulation Parameters' ", &m_simparams.tau, "RO");
		addDoubleParam("tEnd", " step=0.1 label='tEnd' group='Simulation Parameters' ", &m_simparams.tEnd, "RO");
		addDoubleParam("ui", " step=0.1 label='ui' group='Simulation Parameters' ", &m_simparams.ui, "RO");
		addDoubleParam("vi", " step=0.1 label='vi' group='Simulation Parameters' ", &m_simparams.vi, "RO");
		addDoubleParam("xLength", " step=0.1 label='xLength' group='Simulation Parameters' ", &m_simparams.xLength, "RO");
		addDoubleParam("yLength", " step=0.1 label='yLength' group='Simulation Parameters' ", &m_simparams.yLength, "RO");
		addDoubleParam("omg", " step=0.1 label='omega' group='Simulation Parameters' ", &m_simparams.omg, "RO");
	}
	if (strcmp(typeid(test).name(), _float) == 0)
	{
		addFloatParam("alpha", " step=0.1 label='alpha' ", &m_simparams.alpha, "RO");
		addFloatParam("deltaT", " step=0.1 label='deltaT' ", &m_simparams.deltaT, "RO");
		addFloatParam("deltaVec", " step=0.1 label='deltaVec' ", &m_simparams.deltaVec, "RO");
		addFloatParam("eps", " step=0.001 label='eps' ", &m_simparams.eps, "RO");
		addFloatParam("gx", " step=0.1 label='gx' ", &m_simparams.gx, "RO");
		addFloatParam("gy", " step=0.1 label='gy' ", &m_simparams.gy, "RO");
		addFloatParam("KarmanAngle", " step=0.1 label='KarmanAngle' ", &m_simparams.KarmanAngle, "RO");
		addFloatParam("KarmanObjectWidth", " step=0.1 label='KarmanObjectWidth' ", &m_simparams.KarmanObjectWidth, "RO");
		addFloatParam("pi", " step=0.1 label='pi' ", &m_simparams.pi, "RO");
		addFloatParam("re", " step=0.1 label='re' ", &m_simparams.re, "RO");
		addFloatParam("tau", " step=0.1 label='tau' ", &m_simparams.tau, "RO");
		addFloatParam("tEnd", " step=0.1 label='tEnd' ", &m_simparams.tEnd, "RO");
		addFloatParam("ui", " step=0.1 label='ui' ", &m_simparams.ui, "RO");
		addFloatParam("vi", " step=0.1 label='vi' ", &m_simparams.vi, "RO");
		addFloatParam("xLength", " step=0.1 label='xLength' ", &m_simparams.xLength, "RO");
		addFloatParam("yLength", " step=0.1 label='yLength' ", &m_simparams.yLength, "RO");
		addFloatParam("omg", " step=0.1 label='omega' ", &m_simparams.omg, "RO");
	}
	addIntParam("iterMax", " label='iterMax' group='Simulation Parameters' ", &m_simparams.iterMax, "RO");
	addIntParam("iMax", " label='iMax' group='Simulation Parameters' ", &m_simparams.iMax, "RO");
	addIntParam("jMax", " label='jMax' group='Simulation Parameters' ", &m_simparams.jMax, "RO");
	addIntParam("xCells", " label='xCells' group='Simulation Parameters' ", &m_simparams.xCells, "RO");
	addIntParam("yCells", " label='yCells' group='Simulation Parameters' ", &m_simparams.yCells, "RO");
	addStringParam("name", " label='name' group='Simulation Parameters' ", &m_simparams.name, "RO");

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

	//addButtonParam("m_bake", " label='bake parameters' group='Simulation Parameters' ", Bake);
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
		// Field picking render pass
		m_fieldPicking_fbo->bind();
		glClearColor(0.0f,0.0f,0.0f,0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glViewport(0, 0, m_fieldPicking_fbo->getWidth(), m_fieldPicking_fbo->getHeight());
		m_field_layer.drawFieldPicking(m_activeCamera);

		// Update field
		m_field_layer.updateFieldTexture(glfwGetTime());

		// Render all layers to screen
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glClearColor(
			m_window_background_colour[0],
			m_window_background_colour[1],
			m_window_background_colour[2], 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		int width, height;
		glfwGetFramebufferSize(m_window, &width, &height);
		glViewport(0, 0, width, height);

		m_field_layer.draw(m_activeCamera);
		m_boundaryGlyph_layer.draw(m_activeCamera);
		m_boundaryCells_layer.draw(m_activeCamera);
		m_geometry_layer.draw(m_activeCamera);
		m_overlayGrid_layer.draw(m_activeCamera);
		m_interface_layer.draw(m_activeCamera);

		// Draw TB
		TwRefreshBar(bar);
		TwDraw();

		/* Swap front and back buffers */
		glfwSwapBuffers(m_window);

        /* Poll for and process events */
        glfwPollEvents();
    }

	//TODO cleanup
	TwTerminate();
	glfwTerminate();
}

void CavityRenderer::paintBakery()
{
	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(m_window))
	{
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glClearColor(
			m_window_background_colour[0],
			m_window_background_colour[1],
			m_window_background_colour[2], 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		int width, height;
		glfwGetFramebufferSize(m_window, &width, &height);
		glViewport(0, 0, width, height);

		m_boundaryGlyph_layer.draw(m_activeCamera);
		m_boundaryCells_layer.draw(m_activeCamera);
		m_geometry_layer.draw(m_activeCamera);
		m_overlayGrid_layer.draw(m_activeCamera);
		m_interface_layer.draw(m_activeCamera);

		// Draw TB
		TwRefreshBar(bar);
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
			m_interface_layer.updateResources(received_params);
		}
    }

	//TODO cleanup
	TwTerminate();
	glfwTerminate();
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

void TW_CALL ClearDye(void* clientData)
{
	CavityRenderer* cr = (CavityRenderer*)clientData;

	cr->clearDye();
}

void TW_CALL ClearStream(void* clientData)
{
	CavityRenderer* cr = (CavityRenderer*)clientData;

	cr->clearStreamline();
}

void TW_CALL CopyStdStringToClient(std::string& destinationClientString, const std::string& sourceLibraryString)
{
	destinationClientString = sourceLibraryString;
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
	if (mode.compare("RW") == 0)
	{
		TwAddVarRW(bar, name, TW_TYPE_FLOAT, var, def);
		TwSetParam(bar, name, "min", TW_PARAM_FLOAT, 1, &min);
		TwSetParam(bar, name, "max", TW_PARAM_FLOAT, 1, &max);
	}
	if (mode.compare("RO") == 0)
	{
		TwAddVarRO(bar, name, TW_TYPE_FLOAT, var, def);
	}
}

void CavityRenderer::addDoubleParam(const char* name, const char* def, void* var, std::string mode, double min, double max)
{
	if (mode.compare("RW") == 0)
	{
		TwAddVarRW(bar, name, TW_TYPE_DOUBLE, var, def);
		TwSetParam(bar, name, "min", TW_PARAM_DOUBLE, 1, &min);
		TwSetParam(bar, name, "max", TW_PARAM_DOUBLE, 1, &max);
	}
	if (mode.compare("RO") == 0)
	{
		TwAddVarRO(bar, name, TW_TYPE_DOUBLE, var, def);
	}
}

void CavityRenderer::addIntParam(const char* name, const char* def, void* var, std::string mode, int min, int max)
{
	if (mode.compare("RW") == 0)
	{
		TwAddVarRW(bar, name, TW_TYPE_INT32, var, def);
		TwSetParam(bar, name, "min", TW_PARAM_INT32, 1, &min);
		TwSetParam(bar, name, "max", TW_PARAM_INT32, 1, &max);
	}
	if (mode.compare("RO") == 0)
	{
		TwAddVarRO(bar, name, TW_TYPE_INT32, var, def);
	}
}

void CavityRenderer::addBoolParam(const char* name, const char* def, void* var, std::string mode)
{
	if (mode.compare("RW") == 0)
	{
		TwAddVarRW(bar, name, TW_TYPE_BOOLCPP, var, def);
	}
	if (mode.compare("RO") == 0)
	{
		TwAddVarRO(bar, name, TW_TYPE_BOOLCPP, var, def);
	}
}

void CavityRenderer::addVec3Param(const char* name, const char* def, void* var, std::string mode)
{
	if (mode.compare("RW") == 0)
	{
		TwAddVarRW(bar, name, TW_TYPE_DIR3F, var, def);
	}
	if (mode.compare("RO") == 0)
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
	if (mode.compare("RW") == 0)
	{
		TwAddVarRW(bar, name, TW_TYPE_STDSTRING, var, def);
	}
	if (mode.compare("RO") == 0)
	{
		TwAddVarRO(bar, name, TW_TYPE_STDSTRING, var, def);
	}
}

void CavityRenderer::addEnumParam(const char* name, const char* def, void* var, TwEnumVal* _enum, int size, std::string mode)
{
	TwType enumType = TwDefineEnum(def, _enum, size);
	if (mode.compare("RW") == 0) TwAddVarRW(bar, name, enumType, var, NULL);
	if (mode.compare("RO") == 0) TwAddVarRO(bar, name, enumType, var, NULL);
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


const std::string Renderer::IO::readShaderFile(const char* const path)
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

bool Renderer::IO::readPpmHeader(const char* filename, unsigned long& headerEndPos, int& imgDimX, int& imgDimY)
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

bool Renderer::IO::readPpmData(const char* filename, char* imageData, unsigned long dataBegin, int imgDimX, int imgDimY)
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

//#include "Debug.hpp"
bool Renderer::IO::readPfmHeader(const char* filename, unsigned long& headerEndPos, int& imgDimX, int& imgDimY)
{
	std::string buffer;
	std::ifstream file(filename, std::ios::in | std::ios::binary);

	if( !file.is_open() )
		return false;

	file.seekg(0, std::ios::beg);
	std::getline(file, buffer, '\n');
	if(buffer.compare("PF\n") == 0)
	{
		//debug("failing: buffer=\"%s\"", buffer.c_str());
		return false;
	}

	std::getline(file, buffer, ' ');
	imgDimX = std::stoi(buffer);
	std::getline(file, buffer, '\n');
	imgDimY = std::stoi(buffer);

	std::getline(file, buffer);

	headerEndPos = static_cast<unsigned long>(file.tellg());

	file.close();

	return true;
}

bool Renderer::IO::readPfmData(const char* filename, float* imageData, unsigned long dataBegin)
{
	std::ifstream file(filename, std::ios::in | std::ios::binary);

	if (!file.is_open()) return false;
	/*
	/ Determine the length from the beginning of the image data to the end of the file.
	*/
	file.seekg(0, file.end);
	unsigned long length = static_cast<unsigned long>(file.tellg());
	length = length - dataBegin;
	file.seekg(dataBegin, std::ios::beg);
	file.read( reinterpret_cast<char*>(imageData), length);

	file.close();

	return true;
}
