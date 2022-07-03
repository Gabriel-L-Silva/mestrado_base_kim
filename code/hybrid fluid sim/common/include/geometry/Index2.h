//[]---------------------------------------------------------------[]
//|                                                                 |
//| Copyright (C) 2014, 2020 Orthrus Group.                         |
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
// OVERVIEW: Index2.h
// ========
// Class definition for 2D index.
//
// Author: Paulo Pagliosa
// Last revision: 28/01/2020

#ifndef __Index2_h
#define __Index2_h

#include "math/Real.h"
#include <iostream>

namespace cg
{ // begin namespace cg

template <int D> struct Index;


/////////////////////////////////////////////////////////////////////
//
// Index2: 2D index class
// ======
template <>
struct Index<2>
{
  using base_type = int64_t;

  base_type x;
  base_type y;

  HOST DEVICE
  Index()
  {
    // do nothing
  }

  HOST DEVICE
  explicit Index(base_type i, base_type j)
  {
    set(i, j);
  }

  HOST DEVICE
  explicit Index(base_type i)
  {
    set(i, i);
  }

  template <typename V>
  HOST DEVICE
  explicit Index(const V& v)
  {
    set(base_type(v.x), base_type(v.y));
  }

  HOST DEVICE
  void set(base_type i, base_type j)
  {
    x = i;
    y = j;
  }

  HOST DEVICE
  auto& operator =(base_type i)
  {
    set(i, i);
    return *this;
  }

  HOST DEVICE
  auto operator +(const Index<2>& other) const
  {
    return Index<2>{x + other.x, y + other.y};
  }

  HOST DEVICE
  auto operator +(base_type i) const
  {
    return operator +(Index<2>{i});
  }

  HOST DEVICE
  auto operator -(const Index<2>& other) const
  {
    return Index<2>{x - other.x, y - other.y};
  }

  HOST DEVICE
  auto operator -(base_type i) const
  {
    return operator -(Index<2>{i});
  }

  HOST DEVICE
  const auto& operator [](int i) const
  {
    return (&x)[i];
  }

  HOST DEVICE
  auto& operator [](int i)
  {
    return (&x)[i];
  }

  HOST DEVICE
  bool operator ==(const Index<2>& other) const
  {
    return x == other.x && y == other.y;
  }

  HOST DEVICE
  bool operator !=(const Index<2>& other) const
  {
    return !operator ==(other);
  }

  HOST DEVICE
  auto min() const
  {
    return math::min(x, y);
  }

  HOST DEVICE
  auto max() const
  {
    return math::max(x, y);
  }

  HOST DEVICE
  auto& clamp(const Index<2>& s)
  {
    x = x < 0 ? 0 : math::min(x, s.x - 1);
    y = y < 0 ? 0 : math::min(y, s.y - 1);
    return *this;
  }

  HOST DEVICE
  auto prod() const
  {
    return x * y;
  }

  void print(const char* s, std::ostream& f = std::cout) const
  {
    f << s << '(' << x << ',' << y << ")\n";
  }

}; // Index2

using Index2 = Index<2>;

} // end namespace cg

#endif // __Index2_h
