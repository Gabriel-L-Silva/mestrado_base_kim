#ifndef __TrianglePointGenerator_h
#define __TrianglePointGenerator_h

#include "PointGenerator.h"

namespace cg
{

template <typename real>
class TrianglePointGenerator final : public PointGenerator<2, real>
{
  ASSERT_REAL(real, "*TrianglePointGenerator: real should be float or double");
public:
  using vec_type = Vector<real, 2>;
  using bounds_type = Bounds<real, 2>;
  using callback_type = typename PointGenerator<2, real>::callback_type;

  // Invokes callback function for each right triangle points
  // inside boundingBox.
  //
  // This function iterates every right triangle points inside boundingBox
  // where spacing is the size of the right triangle structure.
  void forEachPoint(const bounds_type& bounds, real spacing, const callback_type& callback) const override
  {
    const auto halfSpacing = spacing * 0.5f;
    const auto ySpacing = spacing * real(std::sqrt(3.0f) / 2.0f);
    const auto boxSize = bounds.size();
    const auto& min = bounds.min();
    const auto& max = bounds.max();

    vec_type position;
    bool hasOffset = false;
    bool shouldQuit = false;
    for (int j = 0; j * ySpacing <= boxSize.y && !shouldQuit; ++j) {
      position.y = j * ySpacing + min.y;

      auto offset = (hasOffset) ? halfSpacing : real(0.0f);

      for (int i = 0; i * spacing + offset <= boxSize.x && !shouldQuit; ++i) {
        position.x = i * spacing + offset + min.x;
        if (!callback(position)) {
          shouldQuit = true;
          break;
        }
      }

      hasOffset = !hasOffset;
    }
  }

}; // TrianglePointGenerator

} // end namespace cg

#endif // __TrianglePointGenerator_h
