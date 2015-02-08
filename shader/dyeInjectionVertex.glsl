#version 330

uniform vec2 dye_location[128];

in vec3 v_position;
in vec2 v_uv;

out vec2 uv;

void main()
{
	uv = v_uv;

	vec3 shifted_position = v_position;
	shifted_position.xy += + dye_location[gl_InstanceID];
	
	gl_Position = vec4(shifted_position,1.0);
}