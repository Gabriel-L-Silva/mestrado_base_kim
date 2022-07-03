#include "GLTestSuite.h"
#include "FramebufferTest.h"

namespace
{

static const char* runCode =
  "clearWindow(cg::Color::darkGray);\n"
  "program().setUniform(\"time\", _time);\n"
  "program().setUniformVec2(\"scale\", _scale);\n"
  "glDrawArrays(GL_TRIANGLE_FAN, 0, 4);\n"
  "if (_rotationFlag)\n"
  "  _time += 0.1f / ImGui::GetIO().Framerate;\n";

}

FramebufferTest::FramebufferTest():
  Base{"Multiple Textures + Mesh Renderer"}
{
  _material.spot = cg::Color::gray;
  _material.shine = 20;
}

void
FramebufferTest::initialize()
{
  constexpr float d = 3;

  Base::initialize();
  _fb = new cg::GLFramebuffer{1024, 1024};
  _camera = new cg::Camera;
  _camera->setPosition({0, 0, d});
  _camera->setDistance(d);
  _meshRenderer = new cg::GLMeshRenderer{_camera};
  _mesh = cg::GLGraphics3::sphere();

  cg::Light lights[2];

  lights[0].position = {+d, +d, 0};
  lights[0].flags.set(cg::Light::LightBits::Camera);
  lights[1].position = {-d, -d, 0};
  lights[1].flags.set(cg::Light::LightBits::Camera);
  _meshRenderer->begin();
  _meshRenderer->setLights(lights, lights + 2);
  _meshRenderer->end();
  glEnable(GL_LINE_SMOOTH);
}

void
FramebufferTest::setCurrent()
{
  _camera->setAspectRatio((float)parent()->width() / parent()->height());
  Base::setCurrent();
}

void
FramebufferTest::run()
{
  _fb->use();
  Base::run();
  _fb->disuse();
  _material.texture = (void*)(intptr_t)(_fb->texture());
  clearWindow(cg::Color::darkGray);
  _meshRenderer->begin();
  _meshRenderer->setMaterial(_material);
  _meshRenderer->render(*_mesh);
  _meshRenderer->end();
}

bool
FramebufferTest::onResizeWindow(int width, int height)
{
  _camera->setAspectRatio((float)width / height);
  glViewport(0, 0, width, height);
  return true;
}

constexpr auto CAMERA_RES = 0.01f;
constexpr auto ZOOM_SCALE = 1.01f;

bool
FramebufferTest::onScroll(double, double dy)
{
  _camera->zoom(dy < 0 ? 1.0f / ZOOM_SCALE : ZOOM_SCALE);
  return true;
}

bool
FramebufferTest::onMouseDrag(int dx, int dy)
{
  const auto da = -_camera->viewAngle() * CAMERA_RES;
  _camera->rotateYX(dx * da, dy * da, true);
  return true;
}
