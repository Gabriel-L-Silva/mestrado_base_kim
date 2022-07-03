#ifndef __ExRGBCube_h
#define __ExRGBCube_h

#include "graphics/GLWindow.h"
#include "math/Matrix4x4.h"

class ExRGBCube final: public cg::GLWindow
{
public:
  ExRGBCube(int width, int height);

private:
  cg::GLSL::Program _program;
  GLuint _vao{};
  GLuint _buffers[3];
  struct
  {
    int x;
    int y;
    int size;
  } _viewport;
  cg::vec3f _rotation{0.0f};

  void initialize() override;
  void gui() override;
  void render() override;
  void terminate() override;
  bool windowResizeEvent(int, int) override;

  void setViewport(int, int);

}; // ExRGBCube

#endif // __ExRGBCube_h
