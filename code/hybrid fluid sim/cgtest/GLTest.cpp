#include "GLTest.h"

GLTest::GLTest(const char* title,
  const char* vs,
  const char* fs,
  const char* runCode):
  _program{title},
  _vs{vs},
  _fs{fs},
  _runCode{runCode}
{
  // do nothing
}

void
GLTest::initialize()
{
  _program.setShaders(_vs, _fs);
  glGenVertexArrays(1, &_vao);
}

void
GLTest::setCurrent()
{
  _program.use();
  glBindVertexArray(_vao);
}

void
GLTest::programText()
{
  ImGui::TextColored({1, 1, 0, 1}, "Vertex Shader");
  ImGui::Text(_vs);
  ImGui::Separator();
  ImGui::TextColored({0, 1, 0, 1}, "Fragment Shader");
  ImGui::Text(_fs);
}

void
GLTest::runCodeText()
{
  ImGui::Text(_runCode);
}

void
GLTest::gui()
{
  // do nothing
}

void
GLTest::terminate()
{
  glDeleteVertexArrays(1, &_vao);
}

bool
GLTest::onResizeWindow(int, int)
{
  return false;
}

bool
GLTest::onScroll(double, double)
{
  return false;
}

bool
GLTest::onMouseDown(int, int)
{
  return false;
}

bool
GLTest::onMouseDrag(int, int)
{
  return false;
}
