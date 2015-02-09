#version 330

uniform sampler2D dyeBlob_tx2D;

in vec2 uv;

out vec4 fragColour;

void main()
{
	float mask_value = texture(dyeBlob_tx2D,uv).x;
	
	if(mask_value < 0.5)
		discard;
	else
		fragColour = vec4(0.8,0.0,0.6,1.0);
		
}