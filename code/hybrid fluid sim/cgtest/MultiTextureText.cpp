#include "MultiTextureTest.h"

namespace
{

static const char* vs =
  "#version 430 core\n\n"
  "uniform float time;\n"
  "uniform vec2 scale = vec2(1.0);\n"
  "out vec2 uv1;\n"
  "out vec2 uv2;\n\n"
  "void main()\n{\n"
  "  const vec4 vertices[] = vec4[](\n"
  "    vec4(-1.0, -1.0, 0.0, 1.0),\n"
  "    vec4( 1.0, -1.0, 0.0, 1.0),\n"
  "    vec4( 1.0,  1.0, 0.0, 1.0),\n"
  "    vec4(-1.0,  1.0, 0.0, 1.0));\n"
  "  const vec2 uv[] = vec2[](\n"
  "    vec2(0.0, 0.0),\n"
  "    vec2(1.0, 0.0),\n"
  "    vec2(1.0, 1.0),\n"
  "    vec2(0.0, 1.0));\n"
  "  mat2 m = mat2(\n"
  "    vec2(+cos(time), +sin(time)),\n"
  "    vec2(-sin(time), +cos(time)));\n"
  "  vec4 p = vertices[gl_VertexID];\n\n"
  "  p.x *= scale.x;\n"
  "  p.y *= scale.y;\n"
  "  gl_Position = p;\n"
  "  uv1 = uv[gl_VertexID] * m;\n"
  "  uv2 = uv[gl_VertexID] * transpose(m);\n"
  "}";

static const char* fs =
  "#version 430 core\n\n"
  "in vec2 uv1;\n"
  "in vec2 uv2;\n"
  "uniform sampler2D tex1;\n"
  "uniform sampler2D tex2;\n"
  "out vec4 fragmentColor;\n\n"
  "void main()\n{\n"
  "  fragmentColor = texture(tex1, uv1) + texture(tex2, uv2);\n"
  "}";

static const char* runCode =
  "clearWindow(cg::Color::darkGray);\n"
  "program().setUniform(\"time\", _time);\n"
  "program().setUniformVec2(\"scale\", _scale);\n"
  "glDrawArrays(GL_TRIANGLE_FAN, 0, 4);\n"
  "if (_rotationFlag)\n"
  "  _time += 0.1f / ImGui::GetIO().Framerate;\n";

}

MultiTextureTest::MultiTextureTest(const char* title):
  Base{title, vs, fs, runCode}
{
  // do nothing
}

void
MultiTextureTest::initialize()
{
  Base::initialize();
  _texSelector1.initialize();
  _texSelector2.initialize();
}

void
MultiTextureTest::setCurrent()
{
  Base::setCurrent();
  program().setUniform("tex1", _texSelector1.textureUnitIndex());
  program().setUniform("tex2", _texSelector2.textureUnitIndex());
  _texSelector1.bindTexture();
  _texSelector2.bindTexture();
}

void
MultiTextureTest::run()
{
  clearWindow(cg::Color::darkGray);
  program().setUniform("time", _time);
  program().setUniformVec2("scale", _scale);
  glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
  if (_rotationFlag)
    _time += 0.1f / ImGui::GetIO().Framerate;
}

void
MultiTextureTest::gui()
{
  ImGui::dragVec2("Scale", _scale, {0.0, 1.0});
  _texSelector1.inspect();
  _texSelector2.inspect();
  ImGui::Checkbox("Rotate Textures", &_rotationFlag);
}
