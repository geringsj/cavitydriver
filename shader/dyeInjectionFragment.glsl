#version 330

uniform sampler2D dyeBlob_tx2D;
uniform vec4 dye_colour;

in vec2 uv;

out vec4 fragColour;

void main()
{
	float mask_value = texture(dyeBlob_tx2D,uv).x;
	
	if(mask_value < 0.5)
		discard;
	else
		fragColour = dye_colour;
}