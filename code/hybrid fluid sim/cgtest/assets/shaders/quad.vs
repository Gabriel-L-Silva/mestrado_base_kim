#version 430 core

void main()
{
  const vec4 vertices[] = vec4[](
    vec4(-0.5, -0.5, 0.0, 1.0),
    vec4( 0.5, -0.5, 0.0, 1.0),
    vec4(-0.5,  0.5, 0.0, 1.0),
    vec4( 0.5,  0.5, 0.0, 1.0));

  gl_Position = vertices[gl_VertexID];
}
