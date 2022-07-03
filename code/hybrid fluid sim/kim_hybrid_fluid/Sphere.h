#ifndef __Sphere_h
#define __Sphere_h

#include "math/Surface.h"

namespace cg
{

template <size_t D, typename real>
class Sphere : public math::Surface<D, real>
{
public:
  using Base = math::Surface<D, real>;
  using vec = Vector<real, D>;
  using bounds_type = Bounds<real, D>;
  using math::Surface<D, real>::Surface;

  // Sphere center
  vec center;

  // Sphere radius
  real radius = 1.0f;
  
  Sphere(const vec& center, real radius, const Transform<D, real>& t = Transform<D, real>())
    : center(center), radius(radius), Base(t)
  {
    // do nothing
  }

  Sphere(const Sphere<D, real>& other)
    : center(other.center), radius(other.radius), Base(other.transform)
  {
    // do nothing
  }

private:
  // overriding surface virtual methods
  bounds_type localBounds() const override;

  vec localClosestNormal(const vec& p) const override;

  vec localClosestPoint(const vec& p) const override;

  real localClosestDistance(const vec& p) const override;

}; // Sphere<D, real>

template<size_t D, typename real>
inline Bounds<real, D>
Sphere<D, real>::localBounds() const
{
  vec r{ radius };
  return bounds_type{ center - r, center + r };
}

template<size_t D, typename real>
inline Vector<real, D>
Sphere<D, real>::localClosestNormal(const vec& p) const
{
  auto r = p - center;
  if (r.isNull())
  {
    r = vec::null();
    r.x = real(1.0f);
    return r;
  }
  else
  {
    return r.versor();
  }
}

template<size_t D, typename real>
inline Vector<real, D>
Sphere<D, real>::localClosestPoint(const vec& p) const
{
  return localClosestNormal(p) * radius + center;
}

template<size_t D, typename real>
inline real
Sphere<D, real>::localClosestDistance(const vec& p) const
{
  return math::abs<real>((p - center).length() - radius);
}

} // end namespace cg

#endif // __Sphere_h
