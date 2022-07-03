#ifndef __TriangleTest_h
#define __TriangleTest_h

#include "GLTest.h"

class TriangleTest final: public GLTest
{
public:
  TriangleTest();

private:
  cg::vec2f _p1{-1.0f, -1.0f};
  cg::vec2f _p2{+1.0f, -1.0f};
  cg::vec2f _p3{+1.0f, +1.0f};
  cg::Color _color{0.1f, 0.4f, 0.8f};

  void run() override;
  void gui() override;

}; // TriangleTest

#endif // __TriangleTest_h
