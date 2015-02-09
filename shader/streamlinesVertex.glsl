#version 330

uniform mat4 mvp_matrix;

in vec3 v_position;

void main()
{
	gl_Position = mvp_matrix * vec4(v_position,1.0);
}