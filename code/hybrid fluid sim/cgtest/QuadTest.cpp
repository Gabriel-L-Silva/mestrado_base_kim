#include "QuadTest.h"

namespace
{

static const char* vs =
  "#version 430 core\n\n"
  "uniform vec2 scale = vec2(1.0);\n\n"
  "void main()\n{\n"
  "  const vec4 vertices[] = vec4[](\n"
  "    vec4(-1.0, -1.0, 0.0, 1.0),\n"
  "    vec4( 1.0, -1.0, 0.0, 1.0),\n"
  "    vec4( 1.0,  1.0, 0.0, 1.0),\n"
  "    vec4(-1.0,  1.0, 0.0, 1.0));\n"
  "  vec4 p = vertices[gl_VertexID];\n\n"
  "  p.x *= scale.x;\n"
  "  p.y *= scale.y;\n"
  "  gl_Position = p;\n"
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
  "program().setUniformVec2(\"scale\", _scale);\n"
  "program().setUniformVec4(\"color\", _color);\n"
  "glDrawArrays(GL_TRIANGLE_FAN, 0, 4);\n";

}

QuadTest::QuadTest():
  GLTest{"Quad", vs, fs, runCode}
{
  // do nothing
}

void
QuadTest::run()
{
  clearWindow(cg::Color::darkGray);
  program().setUniformVec2("scale", _scale);
  program().setUniformVec4("color", _color);
  glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
}

void
QuadTest::gui()
{
  ImGui::dragVec2("Scale", _scale, {0.0, 1.0});
  ImGui::colorEdit3("Color", _color);
}
