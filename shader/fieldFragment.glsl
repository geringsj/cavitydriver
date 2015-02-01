#version 330

uniform int mode;

uniform sampler2D field_tx2D;
uniform sampler2D ibfv_tx2D;
uniform vec3 min_values;
uniform vec3 max_values;

in vec2 uv;

out vec4 frag_colour;

//TODO mapping functions R -> R^3
float scaleValue(float value, float min, float max)
{
	return (value - min) / (max-min);
}

// Paraview style blue to red color mapping
vec3 transferFuncton(float value)
{
	vec3 blue = vec3(0.2,0.2,0.7);
	vec3 mintgreen = vec3(0.6,0.6,0.55);
	vec3 red = vec3(0.7,0.0,0.0);
	
	vec3 color_value = (value < 0.5) ? mix(blue,mintgreen,value*2.0) : mix(mintgreen,red,(value*2.0)-1.0);
	
	return color_value;
}

void main()
{
	frag_colour = vec4(0.0,0.0,0.0,1.0);
	
	vec3 field_uvp = texture(field_tx2D,uv).rgb;
	
	if(mode == 0) // u velocity
	{
		frag_colour.rgb = transferFuncton(scaleValue(field_uvp.r,min_values.x,max_values.x));
	}
	else if(mode == 1) // v velocity
	{
		frag_colour.rgb = transferFuncton(scaleValue(field_uvp.g,min_values.y,max_values.y));
	}
	else if(mode == 2) // uv magnitude
	{
		float magnitude = sqrt(field_uvp.r * field_uvp.r + field_uvp.g * field_uvp.g);
		frag_colour.rgb = transferFuncton(scaleValue(magnitude,min_values.z,max_values.z));
	}
	else if(mode == 3) // pressure
	{
		frag_colour.rgb = transferFuncton(scaleValue(field_uvp.b,min_values.z,max_values.z));
	}
	else if(mode == 4)
	{
		frag_colour.rgb = texture(ibfv_tx2D,uv).xyz;
	}
	else if(mode == 5) //combine 2 and 4
	{
		float magnitude = sqrt(field_uvp.r * field_uvp.r + field_uvp.g * field_uvp.g);
		frag_colour.rgb = transferFuncton(scaleValue(magnitude,min_values.z,max_values.z))/2.0;
		frag_colour.rgb += texture(ibfv_tx2D,uv).xyz/2.0;
	}
}