#version 430

uniform sampler2D arrowTexture;

in vec2 tex_coord;

out vec4 frag_colour;

void main(void) 
{		
	vec4 colour_rgba = texture(arrowTexture,tex_coord);
	
	if(colour_rgba.r > 0.7 && colour_rgba.g < 0.3 && colour_rgba.a > 0.7)
		colour_rgba.a = 0.0;
	
	frag_colour = colour_rgba;
};