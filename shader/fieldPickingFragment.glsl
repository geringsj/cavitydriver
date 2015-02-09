#version 330

in vec2 uv;

out vec4 fragColour;

void main()
{
	// 1.0 in third entry is used to identify field pixel from background pixels
	fragColour = vec4(uv,1.0,1.0);
}