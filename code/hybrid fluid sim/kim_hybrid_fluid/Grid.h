#ifndef __Grid_h
#define __Grid_h

#include "geometry/GridBase.h"
#include "math/Vector3.h"
#include "geometry/Bounds3.h"
#include "geometry/Index3.h"

namespace cg
{

/**
* Template class for D-dimensional axis-aligned grid structure.
* 
* This class represents a D-dimensional cartesian grid structure. This class
* stores nothing but the shape of the grid, which can have different cell
* sizes, or grid spacing, per axis. No public constructor is avaliable, thus,
* this class must be used by inheritance.
* 
* \tparam D Defines the number of dimensions.
* \tparam real A floating point type.
*/
template<size_t D, typename real>
class SpacialGrid : public SharedObject {
public:
  ASSERT_REAL(real, "*SpacialGrid: real must be float or double");
  static_assert(D == 2 || D == 3, "Not implemented, use with D = {2, 3}");

  using vec_type = Vector<real, D>; ///< Vector type alias.
  using index_type = Index<D>; ///< Index type alias.
  using bounds_type = Bounds<real, D>; ///< Bounds type alias.

  /** Default destructor. */
  virtual ~SpacialGrid()
  { /*do nothing*/ }

  /** Returns the grid size. */
  const index_type& size() const
  {
    return _size;
  }

  /** Returns the grid origin. */
  const vec_type& origin() const
  {
    return _origin;
  }

  /** Returns the grid spacing. */
  const vec_type& gridSpacing() const
  {
    return _gridSpacing;
  }

  /** Returns the bounding box of the grid. */
  const bounds_type& bounds() const
  {
    return _bounds;
  }

  /** Returns the function that maps grid index to the cell-center position. */
  auto cellCenterPosition() const
  {
    return [=](const index_type& index) -> vec_type {
      return _origin + _gridSpacing * (vec_type(index) + vec_type(static_cast<real>(0.5)));
    };
  }

protected:
  /**
  * Constructs a SpacialGrid with given parameters.
  */
  SpacialGrid(const index_type& size, const vec_type& spacing,
    const vec_type& origin):
    _size(size),
    _gridSpacing(spacing),
    _origin(origin),
    _bounds(bounds_type{ origin, origin + vec_type{size} * spacing  })
  {
    // do nothing 
  }

  /**
  * Sets the size parameters of this grid instance.
  * 
  * \param size The number of cells in each dimension.
  * \param spacing The cell size, or spacing.
  * \param origin The grid origin.
  */
  void set(const index_type& size, const vec_type& spacing,
    const vec_type& origin)
  {
    _size = size;
    _gridSpacing = spacing;
    _origin = origin;
    _bounds = bounds_type{ bounds_type{ origin, origin + vec_type{size} *spacing  } };
  }

  /**
  * Sets the size parameters of this grid instance using another grid.
  * 
  * \param other Another grid instance.
  */
  void set(const SpacialGrid& other)
  {
    _size = other._size;
    _gridSpacing = other._gridSpacing;
    _origin = other._origin;
    _bounds = other._bounds;
  }

private:
  /** Grid size. */
  index_type _size;
  /** Cell size. */
  vec_type _gridSpacing;
  /** Grid origin. */
  vec_type _origin;
  /** Grid bounds. */
  bounds_type _bounds;
}; // SpacialGrid

}  // end namespace cg

#endif // !__Grid_h