#include "GLTestSuite.h"

GLTestSuite::GLTestSuite(std::initializer_list<GLTest*> tests):
  cg::GLWindow{"OpenGL Test Suite", 1280, 720}
{
  for (auto t : tests)
    _tests.push_back(t);
}

void
GLTestSuite::initialize()
{
  if (_tests.empty())
    return;
  for (GLTest* t : _tests)
  {
    t->_parent = this;
    t->initialize();
  }
  setCurrentTest(0);
}

void
GLTestSuite::render()
{
  _currentTest ? _currentTest->run() : void();
}

inline void
GLTestSuite::setCurrentTest(size_t index)
{
  if (_currentTest == _tests[index])
    return;
  (_currentTest = _tests[index])->setCurrent();
}

inline void
GLTestSuite::inspectTestCode()
{
  if (ImGui::CollapsingHeader("GLSL Program Code"))
  {
    ImGui::BeginChild("PC",
      ImVec2(0, 200),
      true,
      ImGuiWindowFlags_HorizontalScrollbar);
    _currentTest->programText();
    ImGui::EndChild();
  }
  if (ImGui::CollapsingHeader("Run Code"))
  {
    ImGui::BeginChild("RC",
      ImVec2(0, 100),
      true,
      ImGuiWindowFlags_HorizontalScrollbar);
    _currentTest->runCodeText();
    ImGui::EndChild();
  }
}

inline void
GLTestSuite::inspectTests()
{
  static size_t testIndex;

  if (ImGui::BeginCombo("Current Test", _currentTest->title()))
  {
    for (size_t i = 0; i < _tests.size(); ++i)
    {
      bool selected = testIndex == i;

      if (ImGui::Selectable(_tests[i]->title(), selected))
        testIndex = i;
      if (selected)
        ImGui::SetItemDefaultFocus();
    }
    ImGui::EndCombo();
    setCurrentTest(testIndex);
  }
}

void
GLTestSuite::gui()
{
  ImGui::SetNextWindowSize(ImVec2(240, 680), ImGuiCond_FirstUseEver);
  ImGui::Begin("Inspector");
  if (!_tests.empty())
  {
    inspectTests();
    inspectTestCode();
    _currentTest->gui();
  }
  ImGui::End();
}

void
GLTestSuite::terminate()
{
  for (GLTest* t : _tests)
  {
    t->terminate();
    t->_parent = nullptr;
  }
}

bool
GLTestSuite::windowResizeEvent(int width, int height)
{
  return _currentTest && _currentTest->onResizeWindow(width, height) ? true :
    GLWindow::windowResizeEvent(width, height);
}

bool
GLTestSuite::scrollEvent(double dx, double dy)
{
  if (ImGui::GetIO().WantCaptureMouse)
    return false;
  return _currentTest ? _currentTest->onScroll(dx, dy) : true;
}

bool
GLTestSuite::mouseButtonInputEvent(int button, int actions, int)
{
  if (ImGui::GetIO().WantCaptureMouse)
    return false;
  if (button == GLFW_MOUSE_BUTTON_LEFT)
    _mouse.dragging = actions == GLFW_PRESS;
  if (_mouse.dragging)
  {
    cursorPosition(_mouse.px, _mouse.py);
    if (_currentTest)
      return _currentTest->onMouseDown(_mouse.px, _mouse.py);
  }
  return true;
}

bool
GLTestSuite::mouseMoveEvent(double xPos, double yPos)
{
  if (!_mouse.dragging)
    return false;
  _mouse.cx = (int)xPos;
  _mouse.cy = (int)yPos;

  const auto dx = (_mouse.cx - _mouse.px);
  const auto dy = (_mouse.cy - _mouse.py);

  _mouse.px = _mouse.cx;
  _mouse.py = _mouse.cy;
  if (dx != 0 || dy != 0 && _currentTest)
    return _currentTest->onMouseDrag(dx, dy);
  return true;
}
