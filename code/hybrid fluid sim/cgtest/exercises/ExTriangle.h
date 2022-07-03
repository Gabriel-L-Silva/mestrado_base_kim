#ifndef __ExTriangle_h
#define __ExTriangle_h

#include "graphics/GLWindow.h"

class ExTriangle final: public cg::GLWindow
{
public:
  ExTriangle(int width, int height);

private:
  GLuint _program{};
  cg::vec4f _vertices[3];
  cg::Color _vertexColors[3];
  GLint _vLoc[3];
  GLint _cLoc[3];
  GLuint _vao{};

  void initialize() override;
  void render() override;
  void terminate() override;

}; // ExTriangle

#endif // __ExTriangle_h
