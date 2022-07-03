#ifndef __Box_h
#define __Box_h

#include "math/Surface.h"

namespace cg
{

/**
* D dimensional Box geometry.
* 
* This class represents a D dimensional Box geometry which extends Surface
* by overriding its queries. This implementation is an axis-aligned box
* that wraps the Bounds type.
* 
* \tparam D Defines the number of dimensions.
* \tparam real A floating point type.
*/
template <size_t D, typename real>
class Box: public math::Surface<D, real>
{
public:
  using Base = math::Surface<D, real>;  ///< Base class alias.
  using vec_type = Vector<real, D>;     ///< Vector type alias.
  using bounds_type = Bounds<real, D>;  ///< Bounding box type alias.

  /** Bounding box of this box. */
  bounds_type boxBounds = bounds_type{ (real)0.0f, (real)1.0f };

  /** Constructs default box */
  Box(const Transform<D, real>& t = Transform<D, real>())
    : Base(t)
  {
    // do nothing
  }

  /** Constructs box with bounds from \p min to \p max */
  Box(const vec_type& min, const vec_type& max, const Transform<D, real>& t = Transform<D, real>())
    : boxBounds(min, max), Base(t)
  {
    // do nothing
  }

  /** Constructs box with given bounds instance \p boundingBox */
  explicit Box(const bounds_type& boundingBox, const Transform<D, real>& t = Transform<D, real>())
    : boxBounds(boundingBox), Base(t)
  {
    // do nothing
  }

  /** Copy constructor */
  Box(const Box<D, real>& other)
    : boxBounds(other.boxBounds), Base(other.transform)
  {
    // do nothing
  }

  virtual ~Box()
  {
    // do nothing
  }

protected:
  // Surface methods

  bounds_type localBounds() const override { return boxBounds; }

  vec_type localClosestPoint(const vec_type& p) const override;

  vec_type localClosestNormal(const vec_type& p) const override;

}; // Box

template<size_t D, typename real>
inline Vector<real, D>
Box<D, real>::localClosestPoint(const vec_type& p) const
{
  auto min = boxBounds.min();
  auto max = boxBounds.max();
  vec_type ret{ real(0.0f) };
  if (boxBounds.contains(p))
  {
    auto r = p - max;
    auto distanceSquared = math::Limits<real>::inf();
    for (size_t i = 0; i < D; ++i)
    {
      vec_type normal{ real(0.0f) };
      normal[i] = 1.0f;
      auto result = p - max;
      result = result - normal.dot(result) * normal + max;
      auto t = (result - p).squaredNorm();

      if (t < distanceSquared)
      {
        ret = result;
        distanceSquared = t;
      }

      normal[i] = -1.0f;
      result = p - min;
      result = result - normal.dot(result) * normal + min;
      t = (result - p).squaredNorm();

      if (t < distanceSquared)
      {
        ret = result;
        distanceSquared = t;
      }
    }
  }
  else
  {
    for (size_t i = 0; i < D; ++i)
      ret[i] = math::clamp<real>(p[i], min[i], max[i]);
  }
  return ret;
}

template<size_t D, typename real>
inline Vector<real, D>
Box<D, real>::localClosestNormal(const vec_type& p) const
{
  vec_type closestNormal{ static_cast<real>(0.0f) };
  auto min = boxBounds.min();
  auto max = boxBounds.max();
  if (boxBounds.contains(p))
  {
    real minDistanceSquared = math::Limits<real>::inf();
    for (size_t i = 0; i < D; ++i)
    {
      vec_type normal{ real(0.0f) };
      normal[i] = 1.0f;
      auto closestPoint = p - max;
      closestPoint = closestPoint - normal.dot(closestPoint) * normal + max;
      auto t = (closestPoint - p).squaredNorm();

      if (t < minDistanceSquared)
      {
        closestNormal = normal;
        minDistanceSquared = t;
      }

      normal[i] = -1.0f;
      closestPoint = p - min;
      closestPoint = closestPoint - normal.dot(closestPoint) * normal + min;
      t = (closestPoint - p).squaredNorm();

      if (t < minDistanceSquared)
      {
        closestNormal = normal;
        minDistanceSquared = t;
      }
    }
  }
  else
  {
    vec_type q{ real(0.0f) };
    for (size_t i = 0; i < D; ++i)
      q[i] = math::clamp<real>(p[i], min[i], max[i]);
    auto r = p - q;
    real maxCosine = -1.0f;

    // go through box planes
    for (size_t i = 0; i < D; ++i)
    {
      // planes on bounds max
      vec_type normal{ real(0.0f) };
      normal[i] = 1.0f;
      real cosineAngle = normal.dot(r);

      if (cosineAngle > maxCosine)
      {
        closestNormal = normal;
        maxCosine = cosineAngle;
      }

      // planes on bounds min
      normal[i] = -1.0f;
      cosineAngle = normal.dot(r);

      if (cosineAngle > maxCosine)
      {
        closestNormal = normal;
        maxCosine = cosineAngle;
      }
    }
  }
  return closestNormal;
}

} // end namespace cg

#endif // __Box_h
