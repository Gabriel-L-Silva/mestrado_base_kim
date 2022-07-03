#ifndef __GridUtils_h
#define __GridUtils_h

#include "MathUtils.h"

namespace cg
{

template <typename real>
auto isInsideSdf(real x)
{
  return math::isNegative(x);
}

template <typename T>
T fractionInsideSdf(T phi0, T phi1) {
  if (isInsideSdf(phi0) && isInsideSdf(phi1)) {
    return 1;
  }
  else if (isInsideSdf(phi0) && !isInsideSdf(phi1)) {
    return phi0 / (phi0 - phi1);
  }
  else if (!isInsideSdf(phi0) && isInsideSdf(phi1)) {
    return phi1 / (phi1 - phi0);
  }
  else {
    return 0;
  }
}

template <typename T>
void cycleArray(T* arr, int size) {
  T t = arr[0];
  for (int i = 0; i < size - 1; ++i) arr[i] = arr[i + 1];
  arr[size - 1] = t;
}

template <typename T>
T fractionInside(T phiBottomLeft, T phiBottomRight, T phiTopLeft, T phiTopRight)
{
  int inside_count = (phiBottomLeft < 0 ? 1 : 0) + (phiTopLeft < 0 ? 1 : 0) +
    (phiBottomRight < 0 ? 1 : 0) + (phiTopRight < 0 ? 1 : 0);
  T list[] = { phiBottomLeft, phiBottomRight, phiTopRight, phiTopLeft };
  if (inside_count == 4) {
    return 1;
  }
  else if (inside_count == 3) {
    // rotate until the positive value is in the first position
    while (list[0] < 0) {
      cycleArray(list, 4);
    }

    // Work out the area of the exterior triangle
    T side0 = 1 - fractionInsideSdf(list[0], list[3]);
    T side1 = 1 - fractionInsideSdf(list[0], list[1]);
    return 1 - 0.5f * side0 * side1;
  }
  else if (inside_count == 2) {
    // rotate until a negative value is in the first position, and the next
    // negative is in either slot 1 or 2.
    while (list[0] >= 0 || !(list[1] < 0 || list[2] < 0)) {
      cycleArray(list, 4);
    }

    if (list[1] < 0) {  // the matching signs are adjacent
      T side_left = fractionInsideSdf(list[0], list[3]);
      T side_right = fractionInsideSdf(list[1], list[2]);
      return 0.5f * (side_left + side_right);
    }
    else {  // matching signs are diagonally opposite
     // determine the centre point's sign to disambiguate this case
      T middle_point = 0.25f * (list[0] + list[1] + list[2] + list[3]);
      if (middle_point < 0) {
        T area = 0;

        // first triangle (top left)
        T side1 = 1 - fractionInsideSdf(list[0], list[3]);
        T side3 = 1 - fractionInsideSdf(list[2], list[3]);

        area += 0.5f * side1 * side3;

        // second triangle (top right)
        T side2 = 1 - fractionInsideSdf(list[2], list[1]);
        T side0 = 1 - fractionInsideSdf(list[0], list[1]);
        area += 0.5f * side0 * side2;

        return 1 - area;
      }
      else {
        T area = 0;

        // first triangle (bottom left)
        T side0 = fractionInsideSdf(list[0], list[1]);
        T side1 = fractionInsideSdf(list[0], list[3]);
        area += 0.5f * side0 * side1;

        // second triangle (top right)
        T side2 = fractionInsideSdf(list[2], list[1]);
        T side3 = fractionInsideSdf(list[2], list[3]);
        area += 0.5f * side2 * side3;
        return area;
      }
    }
  }
  else if (inside_count == 1) {
    // rotate until the negative value is in the first position
    while (list[0] >= 0) {
      cycleArray(list, 4);
    }

    // Work out the area of the interior triangle, and subtract from 1.
    T side0 = fractionInsideSdf(list[0], list[3]);
    T side1 = fractionInsideSdf(list[0], list[1]);
    return 0.5f * side0 * side1;
  }
  else {
    return 0;
  }
}

/// <summary>
/// </summary>
/// <typeparam name="T"></typeparam>
/// <param name="grid"></param>

/// <summary>
/// Fills GridData with specified value
/// </summary>
/// <typeparam name="T">Data type to be stored</typeparam>
/// <param name="grid">GridData to fill</param>
/// <param name="value">Value to put in grid</param>
template <size_t D, typename T>
inline void
fill(GridData<D, T>& grid, const T& value)
{
  auto size = grid.length();
  for (decltype(size) i = 0; i < size; ++i)
    grid[i] = value;
}

template <size_t D, typename T>
inline void
extrapolateToRegion(const Grid<D, T>& input, const GridData<D, char>& valid, unsigned iterations, Grid<D, T>& output)
{
  const auto size = input.size();
  auto n = size.prod();

#ifdef _DEBUG
  assert(size == valid.size());
  assert(size == output.size());
#endif // _DEBUG

  std::vector<char> valid0, valid1;
  valid0.resize(n);
  valid1.resize(n);

  for (decltype(n) i = 0; i < n; ++i)
  {
    valid0[i] = valid[i];
    output[i] = input[i];
  }

  std::array<int64_t, 2*D> neighbors;

  for (unsigned iter = 0; iter < iterations; ++iter)
  {
    for (decltype(n) i = 0; i < n; ++i)
    {
      T sum = T(0);
      unsigned count = 0;

      if (!valid0[i])
      {
        auto index = input.index(i);
        neighborCellIds<D>(input, index, neighbors);
        // iterate through dimensions
        for (decltype(n) j = 0; j < D; ++j)
        {
          // check forward cell in dimension
          if (index[j] + 1 < size[j] && valid0[neighbors[2 * j]])
          {
            sum += output[neighbors[2 * j]];
            ++count;
          }
          // check backward cell in dimension
          if (index[j] > 0 && valid0[neighbors[2 * j + 1]])
          {
            sum += output[neighbors[2 * j + 1]];
            ++count;
          }
        }

        if (count > 0)
        {
          output[i] = sum / ((T)count);
          valid1[i] = 1;
        }
      }
      else
      {
        valid1[i] = 1;
      }
    }
    valid1.swap(valid0);
  }
}

} // end namespace cg

#endif // __GridUtils_h
