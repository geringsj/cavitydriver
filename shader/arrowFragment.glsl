#version 430

uniform sampler2D arrowTexture;

in vec2 tex_coord;

void main(void) 
{		
	vec4 colour = texture2D(arrowTexture,tex_coord);
	if(colour == vec4(1.0,1.0,1.0,1.0)) colour = vec4(1.0,1.0,1.0,0.0);
	gl_FragColor = colour;
};