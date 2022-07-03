#ifndef __ExRandomTriangles_h
#define __ExRandomTriangles_h

#include "graphics/GLWindow.h"

class ExRandomTriangles final: public cg::GLWindow
{
public:
  ExRandomTriangles(int width, int height, int nt);

  ~ExRandomTriangles();

private:
  int _nv;
  cg::GLSL::Program _program;
  cg::vec4f* _vertices;
  cg::Color* _vertexColors;
  GLuint _vao{};
  GLuint _buffers[2];

  void initialize() override;
  void render() override;
  void terminate() override;

  void randTriangles();

}; // ExRandomTriangles

#endif // __ExRandomTriangles_h
