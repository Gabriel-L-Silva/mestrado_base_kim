#ifndef __ExTexture_h
#define __ExTexture_h

#include "graphics/GLWindow.h"

class ExTexture final: public cg::GLWindow
{
public:
  ExTexture(int width, int height);

private:
  cg::GLSL::Program _program;
  GLuint _vao{};
  GLuint _texture{};

  void initialize() override;
  void render() override;
  void terminate() override;

}; // ExTexture

#endif // __ExTexture_h
