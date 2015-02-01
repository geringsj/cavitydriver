#version 330

uniform mat4 mvp_matrix;
uniform sampler2D field_tx2D;

in vec3 v_position;
in vec2 v_uv;

out vec2 uv;

void main()
{
	uv = v_uv;
	
	// Adevected vertex position
	vec3 position = vec3(v_position.xy*2.0 - vec2(1.0),v_position.z);
	position.xy = position.xy + 0.05 * texture(field_tx2D,uv).xy;
	
	//gl_Position = mvp_matrix * vec4(position,1.0);
	gl_Position = vec4(position,1.0);
}