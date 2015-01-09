#version 430

uniform sampler2D grid_tx2D;

in vec2 uv;

out vec4 frag_colour;

void main()
{
	frag_colour = texture(grid_tx2D,uv);
}