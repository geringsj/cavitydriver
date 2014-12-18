#version 430

in vec4 in_position;

uniform mat4 mvp_matrix;

void main(void)
{
	gl_Position = mvp_matrix * in_position;
};