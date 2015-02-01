#version 330

uniform mat4 mvp_matrix;

in vec3 v_position;
in vec2 v_uv;

out vec2 uv;

void main()
{
	uv = v_uv;
	
	gl_Position = mvp_matrix * vec4(v_position,1.0);
}