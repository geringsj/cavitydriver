#version 430

uniform sampler2D arrowTexture;

in vec2 tex_coord;

void main(void) 
{		
	vec4 colour = texture2D(arrowTexture,tex_coord);
	if(colour.r > 0.5 && colour.g > 0.5 && colour.b > 0.5)
		colour.w = 0.0;
	gl_FragColor = colour;
};