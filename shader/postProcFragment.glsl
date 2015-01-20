#version 330

uniform sampler2D grid_tx2D;
uniform sampler2D boundary_cells_tx2D;
uniform sampler2D boundary_glyphs_tx2D;
uniform sampler2D field_tx2D;

uniform vec3 background_colour;

in vec2 uv;

out vec4 frag_colour;

void main()
{	
	frag_colour = vec4(background_colour,1.0);
	
	vec4 grid_rgba = texture(grid_tx2D,uv);
	vec4 boundary_cells_rgba = texture(boundary_cells_tx2D,uv);
	vec4 boundary_glyphs_rgba = texture(boundary_glyphs_tx2D,uv);
	vec4 field_rgba = texture(field_tx2D,uv);
	
	frag_colour.rgb = mix(frag_colour.rgb,field_rgba.rgb,field_rgba.a);
	frag_colour.rgb = mix(frag_colour.rgb,grid_rgba.rgb,grid_rgba.a);
	frag_colour.rgb = mix(frag_colour.rgb,boundary_cells_rgba.rgb,boundary_cells_rgba.a);
	frag_colour.rgb = mix(frag_colour.rgb,boundary_glyphs_rgba.rgb,boundary_glyphs_rgba.a);
}