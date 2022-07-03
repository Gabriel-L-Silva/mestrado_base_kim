//[]---------------------------------------------------------------[]
//|                                                                 |
//| Copyright (C) 2018, 2019 Orthrus Group.                         |
//|                                                                 |
//| This software is provided 'as-is', without any express or       |
//| implied warranty. In no event will the authors be held liable   |
//| for any damages arising from the use of this software.          |
//|                                                                 |
//| Permission is granted to anyone to use this software for any    |
//| purpose, including commercial applications, and to alter it and |
//| redistribute it freely, subject to the following restrictions:  |
//|                                                                 |
//| 1. The origin of this software must not be misrepresented; you  |
//| must not claim that you wrote the original software. If you use |
//| this software in a product, an acknowledgment in the product    |
//| documentation would be appreciated but is not required.         |
//|                                                                 |
//| 2. Altered source versions must be plainly marked as such, and  |
//| must not be misrepresented as being the original software.      |
//|                                                                 |
//| 3. This notice may not be removed or altered from any source    |
//| distribution.                                                   |
//|                                                                 |
//[]---------------------------------------------------------------[]
//
// OVERVIEW: GLRenderWindow2.cpp
// ========
// Source file for OpenGL 2D render window.
//
// Author: Paulo Pagliosa
// Last revision: 16/02/2019

#include "graphics/GLRenderWindow2.h"

namespace cg
{ // begin namespace cg


//////////////////////////////////////////////////////////
//
// GLRenderWindow2 implementation
// ===============
void
GLRenderWindow2::initialize()
{
  _g2 = new GLGraphics2();
  _g2->setAspectRatio(float(width()) / height());
}

constexpr auto MOVE_SCALE = 0.01f;
constexpr auto ZOOM_SCALE = 1.01f;

void
GLRenderWindow2::render()
{
  if (_moveFlags)
  {
    const auto delta = _g2->bounds().size() * MOVE_SCALE;
    auto d = vec2f::null();

    if (_moveFlags.isSet(MoveBits::Left))
      d.x -= delta.x;
    if (_moveFlags.isSet(MoveBits::Right))
      d.x += delta.x;
    if (_moveFlags.isSet(MoveBits::Up))
      d.y += delta.y;
    if (_moveFlags.isSet(MoveBits::Down))
      d.y -= delta.y;
    _g2->pan(d.x, d.y);
  }
  GLWindow::render();
  _g2->updateView();
  renderScene();
}

void
GLRenderWindow2::renderScene()
{
  // do nothing
}

bool
GLRenderWindow2::windowResizeEvent(int width, int height)
{
  _g2->setAspectRatio((float)width / height);
  return true;
}

bool
GLRenderWindow2::keyInputEvent(int key, int action, int mods)
{
  auto active = action != GLFW_RELEASE && mods == GLFW_MOD_ALT;

  switch (key)
  {
    case GLFW_KEY_A:
      _moveFlags.enable(MoveBits::Left, active);
      break;
    case GLFW_KEY_D:
      _moveFlags.enable(MoveBits::Right, active);
      break;
    case GLFW_KEY_Q:
      _moveFlags.enable(MoveBits::Up, active);
      break;
    case GLFW_KEY_Z:
      _moveFlags.enable(MoveBits::Down, active);
      break;
  }
  return false;
}

bool
GLRenderWindow2::scrollEvent(double, double yOffset)
{
  if (ImGui::GetIO().WantCaptureMouse)
    return false;
  _g2->zoom(yOffset < 0 ? 1.0f / ZOOM_SCALE : ZOOM_SCALE);
  return true;
}

bool
GLRenderWindow2::mouseButtonInputEvent(int button, int actions, int mods)
{
  if (ImGui::GetIO().WantCaptureMouse)
    return false;
  (void)mods;

  auto active = actions == GLFW_PRESS;

  if (button == GLFW_MOUSE_BUTTON_RIGHT)
    _dragFlags.enable(DragBits::Force, active);
  else if (button == GLFW_MOUSE_BUTTON_LEFT)
    _dragFlags.enable(DragBits::Source, active);
  if (_dragFlags)
    cursorPosition(_pivotX, _pivotY);
  return true;
}

bool
GLRenderWindow2::mouseMoveEvent(double xPos, double yPos)
{
  if (!_dragFlags)
    return false;
  _mouseX = (int)xPos;
  _mouseY = (int)yPos;

  const auto dx = (_pivotX - _mouseX);
  const auto dy = (_pivotY - _mouseY);

  _pivotX = _mouseX;
  _pivotY = _mouseY;
  if (dx != 0 || dy != 0)
  {
    if (_dragFlags.isSet(DragBits::Source))
    {
      // TODO: rotate
      return true;
    }
    if (_dragFlags.isSet(DragBits::Force))
    {
      // TODO: pan
    }
  }
  return true;
}

} // end namespace cg
