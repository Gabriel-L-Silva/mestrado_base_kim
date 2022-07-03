#include "ExRandomTriangles.h"
#include "Utils.h"

#define STRINGIFY(A) "#version 400\n"#A

static const char* vs = STRINGIFY(
layout(location = 0) in vec4 vertex;
layout(location = 1) in vec4 vertexColor;
out vec4 color;

void main()
{
  gl_Position = vertex;
  color = vertexColor;
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

ExRandomTriangles::ExRandomTriangles(int width, int height, int nt):
  cg::GLWindow("Exercise: Random Triangles", width, height),
  _program{"ExRandomTriangles"},
  _nv{3 * nt}
{
  _vertices = new cg::vec4f[_nv];
  _vertexColors = new cg::Color[_nv];
}

ExRandomTriangles::~ExRandomTriangles()
{
  delete[] _vertexColors;
  delete[] _vertices;
}

template <typename T>
inline int
size(int count)
{
  return count * sizeof(T);
}

void
ExRandomTriangles::initialize()
{
  _program.setShaders(vs, fs).use();
  glGenVertexArrays(1, &_vao);
  glBindVertexArray(_vao);
  glGenBuffers(2, _buffers);
  glBindBuffer(GL_ARRAY_BUFFER, _buffers[0]);
  glBufferData(GL_ARRAY_BUFFER, size<cg::vec4f>(_nv), NULL, GL_DYNAMIC_DRAW);
  glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, _buffers[1]);
  glBufferData(GL_ARRAY_BUFFER, size<cg::Color>(_nv), NULL, GL_DYNAMIC_DRAW);
  glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(1);
}

inline void
ExRandomTriangles::randTriangles()
{
  for (int i = 0; i < _nv; ++i)
  {
    _vertices[i] = {frand(-1, 1), frand(-1, 1), 0, 1};
    _vertexColors[i] = cg::Color{frand(0, 1), frand(0, 1), frand(0, 1)};
  }
  glBindBuffer(GL_ARRAY_BUFFER, _buffers[0]);
  glBufferSubData(GL_ARRAY_BUFFER, 0, size<cg::vec4f>(_nv), _vertices);
  glBindBuffer(GL_ARRAY_BUFFER, _buffers[1]);
  glBufferSubData(GL_ARRAY_BUFFER, 0, size<cg::Color>(_nv), _vertexColors);
}

void
ExRandomTriangles::render()
{
  constexpr int drawCount = 60;
  static int count = drawCount;

  clear(cg::Color::darkGray);
  if (count++ == drawCount)
  {
    randTriangles();
    count = 0;
  }
  glDrawArrays(GL_TRIANGLES, 0, _nv);
}

void
ExRandomTriangles::terminate()
{
  glDeleteBuffers(2, _buffers);
  glDeleteVertexArrays(1, &_vao);
}
