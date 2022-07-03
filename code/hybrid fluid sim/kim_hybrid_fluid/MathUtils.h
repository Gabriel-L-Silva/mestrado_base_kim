#ifndef __MathUtils_h
#define __MathUtils_h

#include "math/Vector4.h"
#include "geometry/Grid2.h"
#include "geometry/Grid3.h"
#include <iostream>
#include <array>
#include <functional>

namespace cg
{

template<typename T>
inline void getBarycentric(
  T x,
  int64_t iHigh,
  int64_t* i,
  T* f) {
  T s = std::floor(x);
  *i = static_cast<int64_t>(s);

  // the first 3 conditions treat boundary cases
  if (0 == iHigh) {
    *i = 0;
    *f = 0;
  }
  else if (*i < 0) {
    *i = 0;
    *f = 0;
  }
  else if (*i >= iHigh) {
    *i = iHigh - 1;
    *f = 1;
  }
  else {
    *f = static_cast<T>(x - s);
  }
}

template <typename S, typename T>
inline S lerp(const S& f0, const S& f1, T t)
{
  return (1 - t) * f0 + t * f1;
}

template <typename S, typename T>
inline S bilerp(
  const S& f00,
  const S& f10,
  const S& f01,
  const S& f11,
  T tx, T ty)
{
  return lerp(
    lerp(f00, f10, tx), // bottom x
    lerp(f01, f11, tx), // top x
    ty                  // interpolate on y
  );
}

template<typename S, typename T>
inline S trilerp(
  const S& f000,
  const S& f100,
  const S& f010,
  const S& f110,
  const S& f001,
  const S& f101,
  const S& f011,
  const S& f111,
  T tx,
  T ty,
  T tz) {
  return lerp(
    bilerp(f000, f100, f010, f110, tx, ty),
    bilerp(f001, f101, f011, f111, tx, ty),
    tz);
}

template <typename real>
inline Vector2<real>
gradient2(const Grid2<real>& grid, const Vector2<real>& spacing, Index2::base_type i, Index2::base_type j)
{
  auto size = grid.size();

  auto leftIndex  = Index2{ (i > 0) ? i - 1 : i, j };
  auto rightIndex = Index2{ (i + 1 < size.x) ? i + 1 : i, j };
  auto downIndex  = Index2{ i, (j > 0) ? j - 1 : j };
  auto upIndex    = Index2{ i, (j < size.y) ? j - 1 : j };

  real left = grid[leftIndex];
  real right = grid[rightIndex];
  real down = grid[downIndex];
  real up = grid[upIndex];

  return 0.5 * Vector2<real>{right - left, up - down} * spacing.inverse();
}

template <typename real>
inline Vector3<real>
gradient3(const Grid3<real>& grid, const Vector3<real>& spacing, Index3::base_type i, Index3::base_type j, Index3::base_type k)
{
  auto size = grid.size();

  auto leftIndex = Index3{ (i > 0) ? i - 1 : i, j, k };
  auto rightIndex = Index3{ (i + 1 < size.x) ? i + 1 : i, j, k };
  auto downIndex = Index3{ i, (j > 0) ? j - 1 : j, k };
  auto upIndex = Index3{ i, (j < size.y) ? j - 1 : j, k };
  auto backIndex = Index3{ i, j, (k > 0) ? k - 1 : k };
  auto frontIndex = Index3{ i, j, (k < size.z) ? k + 1 : k };

  real left = grid[leftIndex];
  real right = grid[rightIndex];
  real down = grid[downIndex];
  real up = grid[upIndex];
  real back = grid[backIndex];
  real front = grid[frontIndex];

  return 0.5 * Vector3<real>{right - left, up - down, front - back} * spacing.inverse();
}

template <typename real>
inline real
laplacian2(const Grid2<real>& grid, const Vector2<real>& spacing, Index2::base_type i, Index2::base_type j)
{
  auto size = grid.size();

  auto leftIndex = Index2{ (i > 0) ? i - 1 : i, j };
  auto rightIndex = Index2{ (i + 1 < size.x) ? i + 1 : i, j };
  auto downIndex = Index2{ i, (j > 0) ? j - 1 : j };
  auto upIndex = Index2{ i, (j < size.y) ? j - 1 : j };

  real left = grid[leftIndex];
  real right = grid[rightIndex];
  real down = grid[downIndex];
  real up = grid[upIndex];

  return (right - left)*spacing.inverse()*spacing.inverse() + (up - down)*spacing.inverse()*spacing.inverse();
}

template <typename real>
inline Vector3<real>
laplacian3(const Grid3<real>& grid, const Vector3<real>& spacing, Index3::base_type i, Index3::base_type j, Index3::base_type k)
{
  auto size = grid.size();

  auto leftIndex = Index3{ (i > 0) ? i - 1 : i, j, k };
  auto rightIndex = Index3{ (i + 1 < size.x) ? i + 1 : i, j, k };
  auto downIndex = Index3{ i, (j > 0) ? j - 1 : j, k };
  auto upIndex = Index3{ i, (j < size.y) ? j - 1 : j, k };
  auto backIndex = Index3{ i, j, (k > 0) ? k - 1 : k };
  auto frontIndex = Index3{ i, j, (k < size.z) ? k + 1 : k };

  real left = grid[leftIndex];
  real right = grid[rightIndex];
  real down = grid[downIndex];
  real up = grid[upIndex];
  real back = grid[backIndex];
  real front = grid[frontIndex];

  return 0.5 * Vector3<real>{right - left, up - down, front - back} *spacing.inverse();
}


template <typename Callback>
inline void forEachIndex2(const Index2& size, Callback func)
{
  Index2 index;
  for (index.y = 0; index.y < size.y; ++index.y)
    for (index.x = 0; index.x < size.x; ++index.x)
      func(index);
}

template <typename Callback>
inline void forEachIndex3(const Index3& size, Callback func)
{
  Index3 index;
  for(index.z = 0; index.z < size.z; ++index.z)
    for (index.y = 0; index.y < size.y; ++index.y)
      for (index.x = 0; index.x < size.x; ++index.x)
        func(index);
}

template <size_t D>
inline void
neighborCellIndexes(const Index<D>& index, std::array<Index<D>, 2*D>& indexes)
{
  // do nothing
}

template <>
inline void
neighborCellIndexes<2>(const Index2& index, std::array<Index2, 4>& indexes)
{
  indexes[0] = Index2{ index.x + 1, index.y };
  indexes[1] = Index2{ index.x - 1, index.y };
  indexes[2] = Index2{ index.x, index.y + 1 };
  indexes[3] = Index2{ index.x, index.y - 1 };
}

template <>
inline void
neighborCellIndexes<3>(const Index3& index, std::array<Index3, 6>& indexes)
{
  indexes[0] = Index3{ index.x + 1, index.y, index.z };
  indexes[1] = Index3{ index.x - 1, index.y, index.z };
  indexes[2] = Index3{ index.x, index.y + 1, index.z };
  indexes[3] = Index3{ index.x, index.y - 1, index.z };
  indexes[4] = Index3{ index.x, index.y, index.z + 1 };
  indexes[5] = Index3{ index.x, index.y, index.z - 1 };
}

template <size_t D, typename T>
inline void
neighborCellIds(const GridData<D, T>& grid, const Index<D>& index, std::array<typename Index<D>::base_type, 2 * D>& ids)
{
  std::array<Index<D>, 2 * D> _n;
  neighborCellIndexes<D>(index, _n);
  auto n = grid.size().prod();

  for (decltype(n) i = 0; i < 2 * D; ++i)
    ids[i] = grid.id(_n[i]);
}

template <size_t D, typename T>
inline void
neighborCellIds(const Grid<D, T>& grid, const Index<D>& index, std::array<typename Index<D>::base_type, 2 * D>& ids)
{
  std::array<Index<D>, 2 * D> _n;
  neighborCellIndexes<D>(index, _n);
  auto n = grid.size().prod();

  for (decltype(n) i = 0; i < 2 * D; ++i)
    ids[i] = grid.id(_n[i]);
}

template <bool checkValidity, size_t D, typename T>
inline void
neighborCellIds(T id, const Index<D>& size, std::array<T, 2 * D>& ids)
{
  static_assert(D > 0 && D < 4, "neighborCellIds can only be used with 0 < D < 4");
  static_assert(std::is_integral<T>::value, "Integral required");

  T inc = 1;
  T val = 1;

  for (size_t i = 0; i < D; ++i)
  {
    if constexpr (checkValidity) val *= size[int(i)];

    ids[i * 2] = id + inc;
    ids[i * 2 + 1] = id - inc;
    inc *= size[int(i)];

    if constexpr (checkValidity)
    {
      // in case checkValidity template argument is set to true
      // at each computation of a neighbor its validity is checked
      // if this neighbor is not valid its value in ids array
      // is set to the max value of type T
      if (std::floor((float(ids[i * 2]) / val)) != std::floor((float(id) / val)))
        ids[i * 2] = std::numeric_limits<T>::max();
      if (std::floor((float(ids[i * 2 + 1]) / val)) != std::floor((float(id) / val)))
        ids[i * 2 + 1] = std::numeric_limits<T>::max();
    }
  }
}

template <size_t D, typename Callback>
inline void forEachIndex(
  const Index<D>& size,
  Callback func
)
{
  Index<D> index;
  if constexpr (D == 2)
  {
    for (index.y = 0; index.y < size.y; index.y++)
      for (index.x = 0; index.x < size.x; index.x++)
        func(index);
  }
  else if constexpr (D == 3)
  {
    for (index.z = 0; index.z < size.z; index.z++)
      for (index.y = 0; index.y < size.y; index.y++)
        for (index.x = 0; index.x < size.x; index.x++)
          func(index);
  }
}

template <size_t D, typename real>
inline Vector<real, D>
projectAndApplyFriction(const Vector<real, D>& vel, const Vector<real, D>& normal, real frictionCoefficient) {
  // projecting velocity into surface using normal, this way velt is perpendicular to normal
  auto velt = vel - (vel.dot(normal)) * normal;
  if (velt.squaredNorm() > 0.0f) {
    double veln = math::max<real>(-vel.dot(normal), 0.0f);
    velt *= math::max<real>(1.0f - frictionCoefficient * veln / velt.length(), 0.0f);
  }

  return velt;
}

} // end namespace cg

#endif // __MathUtils_h