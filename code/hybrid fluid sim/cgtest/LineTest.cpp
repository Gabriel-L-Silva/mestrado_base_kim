#include "LineTest.h"

namespace
{

static const char* vs =
  "#version 430 core\n\n"
  "uniform vec4 points[2];\n\n"
  "void main()\n{\n"
  "  gl_Position = points[gl_VertexID];\n"
  "}";

static const char* fs =
  "#version 430 core\n\n"
  "uniform vec4 color = vec4(1.0);\n"
  "out vec4 fragmentColor;\n\n"
  "void main()\n{\n"
  "  fragmentColor = color;\n"
  "}";

static const char* runCode =
  "clearWindow(cg::Color::darkGray);\n"
  "program().setUniformVec4(\"points[0]\", cg::vec4f{_p1.x, _p1.y, 0, 1});\n"
  "program().setUniformVec4(\"points[1]\", cg::vec4f{_p2.x, _p2.y, 0, 1});\n"
  "program().setUniformVec4(\"color\", _color);\n"
  "glDrawArrays(GL_LINES, 0, 2);\n";

}

LineTest::LineTest():
  GLTest{"Line", vs, fs, runCode}
{
  // do nothing
}

void
LineTest::run()
{
  clearWindow(cg::Color::darkGray);
  program().setUniformVec4("points[0]", cg::vec4f{_p1.x, _p1.y, 0, 1});
  program().setUniformVec4("points[1]", cg::vec4f{_p2.x, _p2.y, 0, 1});
  program().setUniformVec4("color", _color);
  glDrawArrays(GL_LINES, 0, 2);
}

void
LineTest::gui()
{
  ImGui::dragVec2("P1", _p1, {-1.0, 1.0});
  ImGui::dragVec2("P2", _p2, {-1.0, 1.0});
  ImGui::colorEdit3("Color", _color);
}
