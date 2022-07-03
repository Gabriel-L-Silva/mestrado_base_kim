#ifndef __Surface_h
#define __Surface_h

#include "core/SharedObject.h"
#include "math/Matrix3x3.h"
#include "geometry/Ray.h"
#include "geometry/Bounds2.h"
#include "geometry/Bounds3.h"
#include "C:\Users\gabriel\source\repos\fluid-sanches\code\hybrid fluid sim\kim_hybrid_fluid\Transform.h"

namespace cg::math
{

template <size_t D, typename real>
class Surface : public SharedObject
{
public:
  using vec_type = Vector<real, D>;
  using bounds_type = Bounds<real, D>;

  Transform<D, real> transform;

  Surface(const Transform<D, real>& t) :
    transform(t)
  {
    // do nothing
  }

  virtual ~Surface()
  {
    // do nothing
  }

  bounds_type bounds() const
  {
    auto b = localBounds();
    return bounds_type{ transform.transform(b.min()), transform.transform(b.max()) };
  }

  /// Returns the closest point from given point p to the surface.
  vec_type closestPoint(const vec_type& p) const
  {
    return transform.transform(localClosestPoint(transform.inverseTransform(p)));
  }

  /// Returns the closest distance from the given point to the surface
  real closestDistance(const vec_type& p) const
  {
    return localClosestDistance(transform.inverseTransform(p));
  }

  /// Returns the normal to the closest point on the surface from given point p
  vec_type closestNormal(const vec_type& p) const
  {
    return transform.transformDirection(localClosestNormal(transform.inverseTransform(p)));
  }

  bool intesects(const Ray<real, D>& ray) const
  {
    // TODO
  }

  /// Returns true if bounding box can be defined.
  virtual bool isBounded() const
  {
    return true;
  }
  
  /// Returns true if p is inside the surface.
  bool isInside(const vec_type& p) const
  {
    return localIsInside(transform.inverseTransform(p));
  }

protected:
  /// Returns the bounding box of this surface object in local coordinates
  virtual bounds_type localBounds() const = 0;

  /// Returns the closest point at the surface from the point p in local coordinates
  virtual vec_type localClosestPoint(const vec_type& p) const = 0;

  /// Returns the normal to the closest point on the surface from p in local coordinates
  virtual vec_type localClosestNormal(const vec_type& p) const = 0;

  /// Returns the closest distance from given point to the point on the surface in local coordinates
  virtual real localClosestDistance(const vec_type& p) const
  {
    return (p - localClosestPoint(p)).length();
  }

  /// Returns true if p is inside the surface
  virtual bool localIsInside(const vec_type& p) const
  {
    auto _localCP = localClosestPoint(p);
    auto _localN = localClosestNormal(p);
    return (p - _localCP).dot(_localN) < 0.0f;
  }

}; // Surface

} // end namespace cg::math

#endif // __Surface_h
