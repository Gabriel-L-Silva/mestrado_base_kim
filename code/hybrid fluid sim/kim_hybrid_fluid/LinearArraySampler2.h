#ifndef __LinearArraySampler2_h
#define __LinearArraySampler2_h

#include "geometry/Grid3.h"
#include "math/Vector3.h"
#include "geometry/Index3.h"
#include "MathUtils.h"

namespace cg
{

template <typename real, size_t D> class LinearArraySampler;

template <typename real>
class LinearArraySampler<real, 2> final
{
public:
  using vec_type = Vector2<real>;
  using int_type = Index2::base_type;
  ASSERT_REAL(real, "LinearArraySampler2: floating point type expected");

  explicit LinearArraySampler(
    const Grid2<real>* grid,
    const vec_type& gridSpacing,
    const vec_type& gridOrigin
  ) :
    _gridSpacing(gridSpacing),
    _origin(gridOrigin),
    _grid(grid)
  {
    _invGridSpacing = gridSpacing.inverse();
  }

  LinearArraySampler():
    LinearArraySampler(nullptr, vec_type(1, 1), vec_type::null())
  { }

  // copy constructor
  LinearArraySampler(const LinearArraySampler& other) {
    _gridSpacing = other._gridSpacing;
    _invGridSpacing = other._invGridSpacing;
    _origin = other._origin;
    _grid = other._grid;
  }

  real operator()(const Vector2<real>& x) const {
    int_type i, j;
    real fx, fy;

#ifdef _DEBUG
    assert(math::isPositive(_gridSpacing.x) && math::isPositive(_gridSpacing.y));
#endif

    Vector2<real> normalizedX = (x - _origin) * _invGridSpacing;

    auto size = _grid ? _grid->size() :_data->size();

    int_type iSize = static_cast<int_type>(size.x);
    int_type jSize = static_cast<int_type>(size.y);

    getBarycentric(normalizedX.x, iSize - 1, &i, &fx);
    getBarycentric(normalizedX.y, jSize - 1, &j, &fy);

    int_type iPlus1 = std::min(i + 1, iSize - 1);
    int_type jPlus1 = std::min(j + 1, jSize - 1);
    
    if (_grid)
      return bilerp<real, real>(
        (*_grid)[Index2{ i, j }],
        (*_grid)[Index2{ iPlus1, j }],
        (*_grid)[Index2{ i, jPlus1 }],
        (*_grid)[Index2{ iPlus1, jPlus1 }],
        fx,
        fy);
    else
      return bilerp<real, real>(
        (*_data)[_data->id(Index2{ i, j })],
        (*_data)[_data->id(Index2{ iPlus1, j })],
        (*_data)[_data->id(Index2{ i, jPlus1 })],
        (*_data)[_data->id(Index2{ iPlus1, jPlus1 })],
        fx, fx
        );
  };

  void getCoordinatesAndWeights(
    const Vector2<real>& p,
    std::array<Index2, 4>& indices,
    std::array<real, 4>& weights
  ) const
  {
    int64_t i, j;
    real fx, fy;
    vec_type normalizeP{ (p - _origin) * _invGridSpacing };
    Index2 size = _grid ? _grid->size() : _data->size();

    getBarycentric<real>(normalizeP.x, size.x - 1, &i, &fx);
    getBarycentric<real>(normalizeP.y, size.y - 1, &j, &fy);

    auto iPlus1 = std::min(i + 1, size.x - 1);
    auto jPlus1 = std::min(j + 1, size.y - 1);

    indices[0] = Index2{i, j};
    indices[1] = Index2{iPlus1, j};
    indices[2] = Index2{i, jPlus1};
    indices[3] = Index2{iPlus1, jPlus1};

    weights[0] = (1 - fx) * (1 - fy);
    weights[1] = fx * (1 - fy);
    weights[2] = (1 - fx) * fy;
    weights[3] = fx * fy;
  }

  void setGridData(const GridData<2, real>* data)
  {
    _data = data;
  }

private:
  vec_type _gridSpacing;
  vec_type _invGridSpacing;
  vec_type _origin;
  // Removed reference because if we have a Reference as class member
  // the compiler does not generate the operator=() and the class becomes
  // non-assignable which becomes a problem in FaceCenteredGrid.ResetSampler
  // method
  // https://stackoverflow.com/questions/12387239/reference-member-variables-as-class-members#:~:text=There%20are%20a%20few%20important,the%20constructor%20member%20initializer%20list.
  const GridData<2, real>* _data{ nullptr };
  const Grid2<real>* _grid;
};

template <typename real>
using LinearArraySampler2 = LinearArraySampler<real, 2>;

}; // end namespace cg

#endif // __LinearArraySampler2_h