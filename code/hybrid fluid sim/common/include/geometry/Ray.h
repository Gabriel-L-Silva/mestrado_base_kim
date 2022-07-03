//[]---------------------------------------------------------------[]
//|                                                                 |
//| Copyright (C) 2019 Orthrus Group.                               |
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
//  OVERVIEW: Ray.h
//  ========
//  Class definition for 2D/3D ray.
//
// Author: Paulo Pagliosa
// Last revision: 15/05/2019

#ifndef __Ray_h
#define __Ray_h

#include "math/Matrix4x4.h"

namespace cg
{ // begin namespace cg


/////////////////////////////////////////////////////////////////////
//
// Ray: 2D/3D ray class
// ===
template <typename real, int D>
class Ray
{
public:
  ASSERT_REAL(real, "Ray: floating point type expected");

  using vec_type = Vector<real, D>;
  using mat_type = Matrix<real, D + 1, D + 1>;

  vec_type origin;
  vec_type direction;
  real tMin;
  real tMax;

  /// Constructs an empty Ray object.
  HOST DEVICE
  Ray()
  {
    // do nothing
  }

  HOST DEVICE
  Ray(const vec_type& origin, const vec_type& direction):
    tMin{real(0)},
    tMax{math::Limits<real>::inf()}
  {
    set(origin, direction);
  }

  HOST DEVICE
  Ray(const Ray<real, D>& ray, const mat_type& m):
    tMin{ray.tMin},
    tMax{ray.tMax}
  {
    set(m.transform(ray.origin), m.transformVector(ray.direction));
  }

  HOST DEVICE
  void set(const vec_type& origin, const vec_type& direction)
  {
    this->origin = origin;
    this->direction = direction.versor();
  }

  HOST DEVICE
  void transform(const mat_type& m)
  {
    origin = m.transform(origin);
    direction = m.transformVector(direction).versor();
  }

  HOST DEVICE
  vec_type operator ()(real t) const
  {
    return origin + direction * t;
  }

}; // Ray

template <typename real> using Ray2 = Ray<real, 2>;
template <typename real> using Ray3 = Ray<real, 3>;

} // end namespace cg

#endif // __Ray_h
