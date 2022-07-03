#include "ExTriangle.h"

#define STRINGIFY(A) "#version 400\n"#A

static const char* vs = STRINGIFY(
uniform vec4 vertices[3];
uniform vec4 vertexColors[3];
out vec4 color;

void main()
{
  gl_Position = vertices[gl_VertexID];
  color = vertexColors[gl_VertexID];
}
);

static const char* fs = STRINGIFY(
in vec4 color;
out vec4 fragmentColor;

void main()
{
  fragmentColor = color;
}
);

ExTriangle::ExTriangle(int width, int height):
  cg::GLWindow("Exercise: Triangle", width, height)
{
  _vertices[0] = {-0.75f, -0.75f, 0, 1};
  _vertices[1] = {+0.75f, -0.75f, 0, 1};
  _vertices[2] = {0, +0.75f, 0, 1};
  _vertexColors[0] = cg::Color::red;
  _vertexColors[1] = cg::Color::green;
  _vertexColors[2] = cg::Color::blue;
}

void
setShader(GLuint program, const char* source, GLenum type)
{
  auto shader = glCreateShader(type);

  glShaderSource(shader, 1, &source, 0);
  glCompileShader(shader);

  GLint success;

  glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
  assert(success == GL_TRUE);
  glAttachShader(program, shader);
  glDeleteShader(shader);
}

void
ExTriangle::initialize()
{
  _program = glCreateProgram();
  setShader(_program, vs, GL_VERTEX_SHADER);
  setShader(_program, fs, GL_FRAGMENT_SHADER);
  glLinkProgram(_program);

  GLint success;

  glGetProgramiv(_program, GL_LINK_STATUS, &success);
  assert(success == GL_TRUE);
  glUseProgram(_program);
  _vLoc[0] = glGetUniformLocation(_program, "vertices[0]");
  _vLoc[1] = glGetUniformLocation(_program, "vertices[1]");
  _vLoc[2] = glGetUniformLocation(_program, "vertices[2]");
  _cLoc[0] = glGetUniformLocation(_program, "vertexColors[0]");
  _cLoc[1] = glGetUniformLocation(_program, "vertexColors[1]");
  _cLoc[2] = glGetUniformLocation(_program, "vertexColors[2]");
  glGenVertexArrays(1, &_vao);
  glBindVertexArray(_vao);
}

void
ExTriangle::render()
{
  clear(cg::Color::darkGray);
  glUniform4fv(_vLoc[0], 1, (float*)(_vertices + 0));
  glUniform4fv(_vLoc[1], 1, (float*)(_vertices + 1));
  glUniform4fv(_vLoc[2], 1, (float*)(_vertices + 2));
  glUniform4fv(_cLoc[0], 1, (float*)(_vertexColors + 0));
  glUniform4fv(_cLoc[1], 1, (float*)(_vertexColors + 1));
  glUniform4fv(_cLoc[2], 1, (float*)(_vertexColors + 2));
  glDrawArrays(GL_TRIANGLES, 0, 3);
}

void
ExTriangle::terminate()
{
  glDeleteVertexArrays(1, &_vao);
  glDeleteProgram(_program);
}
