#ifndef __CellCenteredScalarGrid_h
#define __CellCenteredScalarGrid_h

#include "ScalarGrid.h"

namespace cg
{

/**
* D-dimensional Cell-centered scalar grid structure.
* 
* This class represents a D-dimensional cell-centered scalar grid, which
* extends ScalarGrid. The class defines the data point at the center of a grid cell.
* 
* \tparam D Defines the number of dimensions.
* \tparam real A floating point type.
*/
template <size_t D, typename real>
class CellCenteredScalarGrid final: public ScalarGrid<D, real>
{
public:
  using Base = ScalarGrid<D, real>;     ///< Base class alias.
  using vec_type = Vector<real, D>;     ///< Vector type alias.
  using bounds_type = Bounds<real, D>;  ///< Bounding box type alias.

  /** Constructs a grid with given size, cell spacing, origin and initial value */
  CellCenteredScalarGrid(const Index<D>& size, const vec_type& gridSpacing, const vec_type& origin, real initialValue = 0.0f) : Base(size, gridSpacing, origin, initialValue)
  {
    // do nothing
  }

  /** Returns the data size of this structure. */
  Index<D> dataSize() const override
  {
    return this->size();
  }

  const auto linearSize2() const
  {
    auto size = this->size();
    return size.x*size.y;
  }
  /**
  * Returns the data position for the first grid cell.
  * 
  * Returns the data position for the first grid cell, in this case is
  * the center of first cell.
  */
  vec_type dataOrigin() const override
  {
    // origin + 0.5 * gridSpacing
    return this->_untouched_bounds[0]
      + static_cast<real>(0.5f) * this->_untouched_cellSize;
  }

  auto const cellSize() const
  {
    return this->_untouched_cellSize;
  }

  auto const bounds() const
  {
    return this->_untouched_bounds;
  }
};

} // end namespace cg

#endif // __CellCenteredScalarGrid_h