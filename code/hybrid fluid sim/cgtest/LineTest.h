#ifndef __LineTest_h
#define __LineTest_h

#include "GLTest.h"

class LineTest final: public GLTest
{
public:
  LineTest();

private:
  cg::vec2f _p1{-1.0f};
  cg::vec2f _p2{+1.0f};
  cg::Color _color{cg::Color::black};

  void run() override;
  void gui() override;

}; // LineTest

#endif // __LineTest_h
