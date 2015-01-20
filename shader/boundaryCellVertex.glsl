#version 330

uniform mat4 mvp_matrix;

in vec3 in_position;
in vec2 in_uv;

out vec2 uv_coord;

void main() 
{
	uv_coord = in_uv;
	
	vec4 pos = vec4(in_position,1.0);
	gl_Position = mvp_matrix * pos;;
};
