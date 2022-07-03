#include "ExRGBCube.h"
#include "Utils.h"
#include <algorithm>

#define STRINGIFY(A) "#version 400\n"#A

static const char* vs = STRINGIFY(
layout(location = 0) in vec4 vertex;
layout(location = 1) in vec4 vertexColor;
uniform mat4 m;
out vec4 color;

void main()
{
  gl_Position = m * vertex;
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

struct uint3
{
  unsigned int x;
  unsigned int y;
  unsigned int z;

}; // int3

namespace cube
{ // begin namespace cube

static const cg::vec4f v[]
{
  {-0.5f, -0.5f, -0.5f, 1}, // 0
  {+0.5f, -0.5f, -0.5f, 1}, // 1
  {+0.5f, +0.5f, -0.5f, 1}, // 2
  {-0.5f, +0.5f, -0.5f, 1}, // 3
  {-0.5f, -0.5f, +0.5f, 1}, // 4
  {+0.5f, -0.5f, +0.5f, 1}, // 5
  {+0.5f, +0.5f, +0.5f, 1}, // 6
  {-0.5f, +0.5f, +0.5f, 1}  // 7
};

static const cg::Color c[]
{
  cg::Color{0.0f, 0.0f, 0.0f}, // 0
  cg::Color{1.0f, 0.0f, 0.0f}, // 1
  cg::Color{1.0f, 1.0f, 0.0f}, // 2
  cg::Color{0.0f, 1.0f, 0.0f}, // 3
  cg::Color{0.0f, 0.0f, 1.0f}, // 4
  cg::Color{1.0f, 0.0f, 1.0f}, // 5
  cg::Color{1.0f, 1.0f, 1.0f}, // 6
  cg::Color{0.0f, 1.0f, 1.0f}  // 7
};

static const uint3 t[]
{
  {0, 3, 1}, {1, 3, 2}, // back
  {4, 5, 7}, {5, 6, 7}, // front
  {0, 4, 3}, {4, 7, 3}, // left
  {1, 2, 5}, {5, 2, 6}, // right
  {0, 1, 4}, {1, 5, 4}, // bottom
  {7, 6, 3}, {6, 2, 3}  // top
};

} // end namespace cube

ExRGBCube::ExRGBCube(int width, int height):
  cg::GLWindow("Exercise: RGB Cube", width, height),
  _program{"ExRGBCube"}
{
  // do nothing
}

template <typename T>
inline int
size(int count)
{
  return count * sizeof(T);
}

inline void
ExRGBCube::setViewport(int w, int h)
{
  _viewport.size = std::min(w, h);
  _viewport.x = (w - _viewport.size) >> 1;
  _viewport.y = (h - _viewport.size) >> 1;
  glViewport(_viewport.x, _viewport.y, _viewport.size, _viewport.size);
}

void
ExRGBCube::initialize()
{
  _program.setShaders(vs, fs).use();
  glGenVertexArrays(1, &_vao);
  glBindVertexArray(_vao);
  glGenBuffers(3, _buffers);
  glBindBuffer(GL_ARRAY_BUFFER, _buffers[0]);
  glBufferData(GL_ARRAY_BUFFER, size<cg::vec4f>(8), cube::v, GL_STATIC_DRAW);
  glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, _buffers[1]);
  glBufferData(GL_ARRAY_BUFFER, size<cg::Color>(8), cube::c, GL_STATIC_DRAW);
  glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(1);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _buffers[2]);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
    size<uint3>(12),
    cube::t,
    GL_STATIC_DRAW);
  glEnable(GL_DEPTH_TEST);
  setViewport(width(), height());
}

void
ExRGBCube::gui()
{
  using namespace cg;

  ImGui::SetNextWindowSize(ImVec2(240, 240), ImGuiCond_FirstUseEver);
  ImGui::Begin("Inspector");
  ImGui::DragFloat3("Rotation", (float*)&_rotation, 1.0f, 0, 0, "%.2f");
  
  auto m = mat4f::rotation(quatf::eulerAngles(_rotation), vec3f::null());

  _program.setUniformMat4("m", m);
  ImGui::End();
}

void
ExRGBCube::render()
{
  clear(cg::Color::darkGray);
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
}

void
ExRGBCube::terminate()
{
  glDeleteBuffers(3, _buffers);
  glDeleteVertexArrays(1, &_vao);
}

bool
ExRGBCube::windowResizeEvent(int width, int height)
{
  setViewport(width, height);
  return true;
}