#version 330

uniform vec4 colour;

out vec4 fragColour;

void main()
{
	fragColour = colour;
	fragColour = vec4(1.0);
}