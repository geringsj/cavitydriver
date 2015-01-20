#version 330

uniform sampler2D boundary_tx2D;

in vec2 uv_coord;

out vec4 frag_colour;

void main(void) 
{
	vec4 colour_rgba = texture(boundary_tx2D,uv_coord);
	
	if(colour_rgba.r > 0.7 && colour_rgba.g < 0.3 && colour_rgba.a > 0.7)
		colour_rgba.a = 0.0;

	frag_colour = colour_rgba;
}