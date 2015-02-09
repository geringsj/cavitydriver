#version 330

uniform sampler2D background_tx2D;
uniform sampler2D advected_tx2D;

in vec2 uv;

out vec4 fragColour;

void main()
{
	vec4 background_rgba = texture(background_tx2D,uv);
	vec4 advected_rgba = texture(advected_tx2D,uv);
	float alpha = 0.97;

	fragColour = advected_rgba*alpha + (1.0-alpha)*background_rgba;
}