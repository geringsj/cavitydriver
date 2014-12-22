#version 430

uniform mat4 mvp_matrix;

in vec4 in_position;

out vec2 tex_coord;

void main() {
	vec2 madd = vec2(0.5f, 0.5f);
	gl_Position = mvp_matrix * in_position;
	tex_coord = in_position.xy * madd + madd;
};