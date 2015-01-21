#version 330

uniform mat4 mvp_matrix;

in vec3 v_position;
in vec2 v_uv;

out vec2 tex_coord;

void main() 
{
	// The texture coords are saved in the z and w component of the in_position
	tex_coord = v_uv;
	// Reset the z and w component
	vec4 pos = vec4(v_position.x,v_position.y,-1.0,1.0);
	gl_Position = mvp_matrix * pos;
};
