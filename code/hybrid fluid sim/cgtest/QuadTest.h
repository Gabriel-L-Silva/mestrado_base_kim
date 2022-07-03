#ifndef __QuadTest_h
#define __QuadTest_h

#include "GLTest.h"

class QuadTest final: public GLTest
{
public:
  QuadTest();

private:
  cg::vec2f _scale{1.0f};
  cg::Color _color{1.0f, 0.6f, 0.2f};

  void run() override;
  void gui() override;

}; // QuadTest

#endif // __QuadTest_h
