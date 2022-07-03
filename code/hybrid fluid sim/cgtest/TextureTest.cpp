#include "graphics/GLImage.h"
#include "TextureTest.h"

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
  "uniform sampler2D tex;\n"
  "out vec4 fragmentColor;\n\n"
  "void main()\n{\n"
  "  fragmentColor = texture(tex, gl_FragCoord.xy / textureSize(tex, 0));\n"
  "}";

static const char* runCode =
  "clearWindow(cg::Color::darkGray);\n"
  "program().setUniformVec2(\"scale\", _scale);\n"
  "glDrawArrays(GL_TRIANGLE_FAN, 0, 4);\n";

}

TextureTest::TextureTest(const char* title):
  Base{title, vs, fs, runCode}
{
  // do nothing
}

void
TextureTest::initialize()
{
  Base::initialize();
  _texSelector.initialize();
}

void
TextureTest::setCurrent()
{
  Base::setCurrent();
  program().setUniform("tex", _texSelector.textureUnitIndex());
  _texSelector.bindTexture();
}

void
TextureTest::run()
{
  clearWindow(cg::Color::darkGray);
  program().setUniformVec2("scale", _scale);
  glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
}

void
TextureTest::gui()
{
  ImGui::dragVec2("Scale", _scale, {0.0, 1.0});
  _texSelector.inspect();
}
