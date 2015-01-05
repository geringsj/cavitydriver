#version 430

uniform mat4 mvp_matrix;

in vec4 in_position;

out vec2 tex_coord;

void main() 
{
	// The texture coords are saved in the z and w component of the in_position
	tex_coord = in_position.zw;
	// Reset the z and w component
	vec4 pos = vec4(in_position.x,in_position.y,-1.0,1.0);
	gl_Position = mvp_matrix * pos;
};
