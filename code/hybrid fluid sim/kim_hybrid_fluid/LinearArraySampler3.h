#ifndef __LinearArraySampler3_h
#define __LinearArraySampler3_h

#include "LinearArraySampler2.h"

namespace cg
{

template <typename real>
class LinearArraySampler<real, 3> final
{
public:
  using vec_type = Vector3<real>;
  using int_type = Index3::base_type;
  ASSERT_REAL(real, "LinearArraySampler3: floating point type expected");

  explicit LinearArraySampler(
    const Grid3<real>* grid,
    const vec_type& gridSpacing,
    const vec_type& gridOrigin
  ) :
    _gridSpacing(gridSpacing),
    _origin(gridOrigin),
    _grid(grid)
  {
    _invGridSpacing = gridSpacing.inverse();
  }

  LinearArraySampler() :
    LinearArraySampler(nullptr, vec_type(1, 1, 1), vec_type::null())
  { }

  // copy constructor
  LinearArraySampler(const LinearArraySampler<real, 3>& other) {
    _gridSpacing = other._gridSpacing;
    _invGridSpacing = other._invGridSpacing;
    _origin = other._origin;
    _grid = other._grid;
  }

  real operator()(const vec_type& x) const
  {
    int_type i, j, k;
    real fx, fy, fz;

#ifdef _DEBUG
    assert(math::isPositive(_gridSpacing.x) && math::isPositive(_gridSpacing.y) && math::isPositive(_gridSpacing.z));
#endif

    vec_type normalizedX = (x - _origin) * _invGridSpacing;

    auto size = Index3{ _grid->size() };

    int_type iSize = int_type(size.x);
    int_type jSize = int_type(size.y);
    int_type kSize = int_type(size.z);

    getBarycentric(normalizedX.x, iSize - 1, &i, &fx);
    getBarycentric(normalizedX.y, jSize - 1, &j, &fy);
    getBarycentric(normalizedX.z, kSize - 1, &k, &fz);

    int_type iPlus1 = std::min(i + 1, iSize - 1);
    int_type jPlus1 = std::min(j + 1, jSize - 1);
    int_type kPlus1 = std::min(k + 1, kSize - 1);


    return trilerp(
      (*_grid)[Index3{ i, j, k }],
      (*_grid)[Index3{ iPlus1, j, k }],
      (*_grid)[Index3{ i, jPlus1, k }],
      (*_grid)[Index3{ iPlus1, jPlus1, k }],
      (*_grid)[Index3{ i, j, kPlus1 }],
      (*_grid)[Index3{ iPlus1, j, kPlus1 }],
      (*_grid)[Index3{ i, jPlus1, kPlus1 }],
      (*_grid)[Index3{ iPlus1, jPlus1, kPlus1 }],
      fx,
      fy,
      fz);
  };

  void getCoordinatesAndWeights(
    const vec_type& p,
    std::array<Index3, 8>& indices,
    std::array<real, 8>& weights
  ) const
  {
    int64_t i, j, k;
    real fx, fy, fz;
    vec_type normalizeP{ (p - _origin) * _invGridSpacing };
    Index3 size{ _grid->size() };

    getBarycentric<real>(normalizeP.x, size.x - 1, &i, &fx);
    getBarycentric<real>(normalizeP.y, size.y - 1, &j, &fy);
    getBarycentric<real>(normalizeP.z, size.z - 1, &k, &fz);

    auto iPlus1 = std::min(i + 1, size.x - 1);
    auto jPlus1 = std::min(j + 1, size.y - 1);
    auto kPlus1 = std::min(k + 1, size.z - 1);

    indices[0] = Index3{ i, j, k };
    indices[1] = Index3{ iPlus1, j, k };
    indices[2] = Index3{ i, jPlus1, k };
    indices[3] = Index3{ iPlus1, jPlus1, k };
    indices[4] = Index3{ i, j, kPlus1 };
    indices[5] = Index3{ iPlus1, j, kPlus1 };
    indices[6] = Index3{ i, jPlus1, kPlus1 };
    indices[7] = Index3{ iPlus1, jPlus1, kPlus1 };

    weights[0] = (1 - fx) * (1 - fy) * (1 - fz);
    weights[1] = fx * (1 - fy) * (1 - fz);
    weights[2] = (1 - fx) * fy * (1 - fz);
    weights[3] = fx * fy * (1 - fz);
    weights[4] = (1 - fx) * (1 - fy) * fz;
    weights[5] = fx * (1 - fy) * fz;
    weights[6] = (1 - fx) * fy * fz;
    weights[7] = fx * fy * fz;
  }

private:
  vec_type _gridSpacing;
  vec_type _invGridSpacing;
  vec_type _origin;
  const Grid3<real>* _grid;
}; // LinearArraySampler3

} // end namespace cg

#endif // __LinearArraySampler3_h