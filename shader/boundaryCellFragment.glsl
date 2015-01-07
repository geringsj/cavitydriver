#version 430

uniform sampler2D boundary_tx2D;

in vec2 uv_coord;

void main(void) 
{
	gl_FragColor =  texture(boundary_tx2D,uv_coord);
}