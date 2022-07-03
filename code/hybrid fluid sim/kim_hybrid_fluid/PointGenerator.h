#ifndef __PointGenerator_h
#define __PointGenerator_h

#include <vector>
#include <functional>
#include "math/Vector4.h"
#include "geometry/Bounds3.h"

namespace cg
{

// Abstract template for point generator
template <size_t D, typename real>
class PointGenerator
{
public:
  using vec_type = Vector<real, D>;
  using bounds_type = Bounds<real, D>;
  using callback_type = std::function<bool(const vec_type&)>;

  virtual ~PointGenerator()
  {
    // do nothing
  }

  // Generates points to output array points inside given bounds
  // with target point spacing.
  void generate(const bounds_type& bounds, real spacing, std::vector<vec_type>& points) const
  {
    forEachPoint(bounds, spacing, [&points](const vec_type& v) {
      points.push_back(v);
      return true;
    });
  }


  // Iterates every point within the bounding box with specified
  // point pattern and invokes the callback function.
  //
  // This function iterates every point within the bounding box and invokes
  // the callback function. The position of the point is specified by the
  // actual implementation. The suggested spacing between the points are
  // given by \p spacing. The input parameter of the callback function is
  // the position of the point and the return value tells whether the
  // iteration should stop or not.
  virtual void forEachPoint(const bounds_type& bounds, real spacing, const callback_type& callback) const = 0;
}; // PointGenerator

} // end namespace cg

#endif // __PointGenerator_h
