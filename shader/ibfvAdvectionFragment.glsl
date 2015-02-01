#version 330

uniform sampler2D previous_frame_tx2d;

in vec2 uv;

out vec4 frag_colour;

void main()
{
	frag_colour = texture(previous_frame_tx2d,uv);
}