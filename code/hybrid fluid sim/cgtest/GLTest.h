#ifndef __GLTest_h
#define __GLTest_h

#include "core/SharedObject.h"
#include "graphics/GLWindow.h"
#include <string>

//
// Forward definition
//
class GLTestSuite;

class GLPipeline; // TODO

class GLTest abstract: public cg::SharedObject
{
public:
  auto title() const
  {
    return _program.name();
  }

protected:
  GLTest(const char*, const char*, const char*, const char*);

  virtual void initialize();
  virtual void setCurrent();
  virtual void programText();
  virtual void runCodeText();
  virtual void run() = 0;
  virtual void gui();
  virtual void terminate();

  virtual bool onResizeWindow(int, int);
  virtual bool onScroll(double, double);
  virtual bool onMouseDown(int, int);
  virtual bool onMouseDrag(int, int);

  auto& program()
  {
    return _program;
  }

  auto parent() const
  {
    return _parent;
  }

private:
  const char* _vs;
  const char* _fs;
  const char* _runCode;
  GLTestSuite* _parent{};
  cg::GLSL::Program _program;
  GLuint _vao{};

  friend GLTestSuite;

}; // GLTest

inline void
clearWindow(const cg::Color& color)
{
  glClearColor(color.r, color.g, color.b, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

namespace ImGui
{

inline bool
dragVec2(const char* label, cg::vec2f& v)
{
  return ImGui::DragFloat2(label, &v.x, 0.01f, 0, 0, "%.2g");
}

inline bool
dragVec2(const char* label, cg::vec2f& v, cg::vec2f r)
{
  return ImGui::SliderFloat2(label, &v.x, r.x, r.y, "%.2g");
}

inline bool
dragVec3(const char* label, cg::vec3f& v)
{
  return ImGui::DragFloat3(label, &v.x, 0.01f, 0, 0, "%.2g");
}

inline bool
dragVec3(const char* label, cg::vec3f& v, cg::vec2f r)
{
  return ImGui::SliderFloat3(label, &v.x, r.x, r.y, "%.2g");
}

inline bool
colorEdit3(const char* label, cg::Color& color)
{
  return ImGui::ColorEdit3(label, &color.r);
}

}

#endif // __GLTest_h
