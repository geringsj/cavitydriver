#version 430

uniform sampler2D arrowTexture;

in vec2 tex_coord;

void main(void) 
{		
	gl_FragColor = texture2D(arrowTexture,tex_coord);
};