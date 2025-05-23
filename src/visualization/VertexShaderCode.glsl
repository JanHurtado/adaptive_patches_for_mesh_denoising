#version 430

in layout(location=0) vec3 position;
in layout(location=1) vec3 vertexColor;

out vec3 theColor;
uniform mat4 mvp;

void main()
{
	gl_Position = mvp*vec4(position,1.0);
	theColor = vertexColor;
}