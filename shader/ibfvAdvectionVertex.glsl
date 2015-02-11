#version 330

uniform sampler2D field_tx2D;
uniform float field_global_maxU;

in vec3 v_position;
in vec2 v_uv;

out vec2 uv;

vec2 scaleVelocity(vec2 velocity)
{
	// Adaptive scaling for the mesh advection
	return 0.1 * (velocity)/field_global_maxU;
}

void main()
{
	uv = v_uv;
	
	// Adevected vertex position
	vec3 position = vec3(v_position.xy*2.0 - vec2(1.0),v_position.z);
	// Adaptive scaling for the mesh advection
	vec2 field_uv = scaleVelocity(texture(field_tx2D,uv).xy);
	
	position.xy = position.xy + field_uv;
	
	//gl_Position = mvp_matrix * vec4(position,1.0);
	gl_Position = vec4(position,1.0);
}