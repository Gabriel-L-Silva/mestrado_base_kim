#ifndef __GLTestSuite_h
#define __GLTestSuite_h

#include "GLTest.h"
#include <vector>

class GLTestSuite final: public cg::GLWindow
{
public:
  GLTestSuite(std::initializer_list<GLTest*>);

private:
  std::vector<cg::Reference<GLTest>> _tests;
  GLTest* _currentTest{};
  struct
  {
    int px, py;
    int cx, cy;
    bool dragging{};
  } _mouse;

  void initialize() override;
  void render() override;
  void gui() override;
  void terminate() override;

  bool windowResizeEvent(int, int) override;
  bool scrollEvent(double, double) override;
  bool mouseButtonInputEvent(int, int, int) override;
  bool mouseMoveEvent(double, double) override;

  void setCurrentTest(size_t);

  void inspectTests();
  void inspectTestCode();

}; // GLTestSuite

#endif // __GLTestSuite
