#ifndef __FramebufferTest_h
#define __FramebufferTest_h

#include "graphics/GLMeshRenderer.h"
#include "GLFramebuffer.h"
#include "MultiTextureTest.h"

class FramebufferTest final: public MultiTextureTest
{
public:
  FramebufferTest();

private:
  using Base = MultiTextureTest;

  cg::Reference<cg::GLFramebuffer> _fb;
  cg::Reference<cg::TriangleMesh> _mesh;
  cg::Reference<cg::GLMeshRenderer> _meshRenderer;
  cg::Reference<cg::Camera> _camera;
  cg::Material _material;

  void setCurrent() override;
  void initialize() override;
  void run() override;

  bool onResizeWindow(int, int) override;
  bool onScroll(double, double) override;
  bool onMouseDrag(int, int) override;

}; // FramebufferTest

#endif // __FramebufferTest_h
