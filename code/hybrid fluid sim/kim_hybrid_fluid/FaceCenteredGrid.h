#ifndef __FaceCenteredGrid_h
#define __FaceCenteredGrid_h

#include <vector>
#include "geometry/Grid2.h"
#include "core/SoA.h"
#include "LinearArraySampler3.h"
#include "VectorField.h"
#include "Grid.h"

namespace cg
{

/**
* Face-centered D-dimensional (a.k.a MAC or staggered) grid.
* 
* This class implements face-centered grid, also known as marker-and-cell (MAC)
* or staggered grid. This vector grid stores each vector component at face
* center. Thus, u, v and possibly w components are not collocated.
* This class uses a SOA structure to store important data along with
* each scalar grid (representing a velocity component), such as the grid
* origin and a related ArraySampler<real, D> to operate with the grid data.
* 
* \tparam D Defines the number of dimensions.
* \tparam real A floating point type.
*/
template <size_t D, typename real>
class FaceCenteredGrid : public SpacialGrid<D, real>
{
public:
  ASSERT_REAL(real, "*FaceCenteredGrid: real must be float or double");
  /** Base class alias. */
  using Base = SpacialGrid<D, real>;
  /** Index type from Index<2> class. */
  using base_type = Index<2>::base_type;
  /** Vector type alias. */
  using vec_type = Vector<real, D>;
  /** Bounding box type alias. */
  using bounds_type = Bounds<real, D>;
  /** Struct-of-array type alias */
  using Data = SoA<ArrayAllocator, Reference<RegionGrid<D, real, real>>, vec_type, LinearArraySampler<real, D>>;
  using id_type = typename Index<D>::base_type;

  /** Construct a grid using given parameters. */
  FaceCenteredGrid(
    const Index<D>& size,
    const vec_type& gridSpacing,
    const vec_type& origin) :
    Base(size, gridSpacing, origin),
    _data(D)
  {
    if constexpr (D == 2)
    {
      auto uGrid = new RegionGrid<D, real, real>(bounds_type{ origin, gridSpacing * vec_type{size} }, size + Index<D>{1, 0});
      auto vGrid = new RegionGrid<D, real, real>(bounds_type{ origin, gridSpacing * vec_type{size} }, size + Index<D>{0, 1});
      auto uOrigin = origin + vec_type{ 0.0, gridSpacing.y * 0.5f };
      auto vOrigin = origin + vec_type{ gridSpacing.x * 0.5f, 0.0 };
      _data.set(
        0,
        uGrid,
        uOrigin,
        LinearArraySampler<real, D>(uGrid, gridSpacing, uOrigin)
      );
      _data.set(
        1,
        vGrid,
        vOrigin,
        LinearArraySampler<real, D>(vGrid, gridSpacing, vOrigin)
      );
    }
    else if constexpr (D == 3)
    {
      auto uGrid = new RegionGrid<D, real, real>(bounds_type{ origin, gridSpacing * vec_type{size} }, size + Index<D>{1, 0, 0});
      auto vGrid = new RegionGrid<D, real, real>(bounds_type{ origin, gridSpacing * vec_type{size} }, size + Index<D>{0, 1, 0});
      auto wGrid = new RegionGrid<D, real, real>(bounds_type{ origin, gridSpacing * vec_type{size} }, size + Index<D>{0, 0, 1});
      _data.set(
        0,
        uGrid,
        origin + (vec_type{ 0.0, gridSpacing.y, gridSpacing.z } * 0.5f),
        LinearArraySampler<real, D>(uGrid, gridSpacing, origin)
      );
      _data.set(
        1,
        vGrid,
        origin + (vec_type{ gridSpacing.x, 0.0, gridSpacing.z } * 0.5f),
        LinearArraySampler<real, D>(vGrid, gridSpacing, origin)
      );
      _data.set(
        2,
        wGrid,
        origin + (vec_type{ gridSpacing.x, gridSpacing.y, 0.0 } * 0.5f),
        LinearArraySampler<real, D>(wGrid, gridSpacing, origin)
      );
    }
  }

  /** Default destructor. */
  ~FaceCenteredGrid()
  {
    // delete
  }

  /**
  * Returns sampled vector for given point \p x.
  * 
  * This is done by using the ArraySampler<D, real> instances to sample
  * the velocity components.
  */
  vec_type sample(const vec_type& x) const
  {
    vec_type sampled = vec_type::null();
    for (size_t i = 0; i < D; ++i)
      sampled[(int)i] = _data.get<2>(i)(x);
    return sampled;
  }


  /**
  * Returns the I-th velocity component.
  * 
  * \param[in] i Grid cell index.
  * \tparam I One of the D dimensions.
  */
  template <size_t I>
  constexpr real& velocityAt(Index<D> i)
  {
    static_assert(I >= 0 && I < D);
    return (*_data.get<0>(I)).operator[](i);
  }

  /**
  * Returns the I-th velocity component.
  *
  * \param[in] i Grid cell index.
  * \tparam I One of the D dimensions.
  */
  template <size_t I>
  constexpr const real& velocityAt(Index<D> i) const
  {
    static_assert(I >= 0 && I < D);
    return (*_data.get<0>(I)).operator[](i);
  }

  /**
  * Returns the I-th velocity component.
  *
  * \param[in] id Grid cell id.
  * \tparam I One of the D dimensions.
  */
  template <size_t I>
  constexpr real& velocityAt(id_type id)
  {
    static_assert(I >= 0 && I < D);
    return (*_data.get<0>(I))[id];
  }

  /**
  * Returns the I-th velocity component.
  *
  * \param[in] id Grid cell id.
  * \tparam I One of the D dimensions.
  */
  template <size_t I>
  constexpr const real& velocityAt(id_type id) const
  {
    return (*_data.get<0>(I))[id];
  }

  /**
  * Returns the d-th velocity component.
  *
  * Equivalent to other FaceCenteredGrid<D, real>::velocityAt overloads,
  * the only difference is that this method doesn't take any template
  * parameters.
  * 
  * \param[in] d Selected dimension.
  * \param[in] i Grid cell index.
  */
  real& velocityAt(size_t d, const Index<D>& i)
  {
    if (d == 0)
      return velocityAt<0>(i);
    else if (d == 1)
      return velocityAt<1>(i);
    if constexpr (D == 3)
    {
      if (d == 2)
        return velocityAt<2>(i);
    }
    else
      return velocityAt<0>(i);
  }

  /**
  * Returns the d-th velocity component.
  *
  * Equivalent to other FaceCenteredGrid<D, real>::velocityAt overloads,
  * the only difference is that this method doesn't take any template
  * parameters.
  *
  * \param[in] d Selected dimension.
  * \param[in] i Grid cell index.
  */
  const real& velocityAt(size_t d, const Index<D>& i) const
  {
    if (d == 0)
      return velocityAt<0>(i);
    else if (d == 1)
      return velocityAt<1>(i);
    if constexpr (D == 3)
    {
      if (d == 2)
        return velocityAt<2>(i);
    }
  }

  /**
  * Returns the I-th grid size.
  * 
  * \tparam I One of the D dimensions.
  */
  template <size_t I>
  auto iSize() const
  {
    static_assert(I >= 0 && I < D);
    return (*_data.get<0>(I)).size();
  }


  /**
  * Returns the I-th grid data origin.
  * 
  * \tparam I One of the D dimensions.
  */
  template <size_t I>
  vec_type iOrigin() const
  {
    static_assert(I >= 0 && I < D);
    return _data.get<1>(I);
  }


  /**
  * Returns the I-th grid.
  *
  * \tparam I One of the D dimensions.
  */
  template <size_t I>
  auto data() const
  {
    static_assert(I >= 0 && I < D);
    return _data.get<0>(I);
  }

  /** Fills all D grids with \p value. */
  void fill(const vec_type& value)
  {
    for (size_t i = 0; i < D; ++i)
    {
      auto& grid = _data.get<0>(i);
      auto length = grid->length();
      for (int64_t j = 0; j < length; ++j)
        grid->operator[](j) = value[(int)i];
    }
  }

  /**
  * Returns a function that maps a grid cell to its actual position.
  * 
  * \tparam I One of the D dimensions.
  */
  template <size_t I>
  auto positionInSpace() const
  {
    static_assert(I >= 0 && I < D);

    auto o = _data.get<1>(I);
    auto s = this->gridSpacing();
    return [o, s](const Index<D>& i) -> vec_type {
      return o + s * vec_type { i };
    };
  }

  /**
  * Returns a function that maps a grid cell to its actual position.
  * 
  * This method is the same as the template
  * FaceCenteredGrid<D, real>::positionInSpace method, but taking \p I as
  * runtime parameter.
  * \param[in] I One of the D dimensions.
  */
  auto positionInSpace(size_t I) const
  {
#ifdef _DEBUG
    assert(I >= 0 && I < D);
#endif // _DEBUG

    auto o = _data.get<1>(I);
    auto s = this->gridSpacing();
    return [o, s](const Index<D>& i) -> vec_type {
      return o + s * vec_type { i };
    };
  }

  /**
  * Returns interpolated value at cell center.
  * \param[in] index Cell index.
  */
  auto valueAtCellCenter(const Index<D>& index) const
  {
#ifdef _DEBUG
    assert(index.x < this->size().x);
    assert(index.y < this->size().y);
    if constexpr(D == 3)
      assert(index.z < this->size().z);
#endif // _DEBUG

    vec_type v{};
    std::array<Index<D>, D> forward;
    for (size_t i = 0; i < D; ++i)
    {
      forward[i] = index;
      forward[i][int(i)]++;
    }

    v.x = velocityAt<0>(index) + velocityAt<0>(forward[0]);
    v.y = velocityAt<1>(index) + velocityAt<1>(forward[1]);
    if constexpr (D == 3)
      v.z = velocityAt<2>(index) + velocityAt<2>(forward[2]);
    return v * 0.5f;
  }

  /**
  * Invokes the given function \p func for each I-data point.
  * 
  * This method invokes the given function \p func for each I-data point.
  * The input parameter for the callback function is the cell index of a data
  * point.
  * 
  * \tparam I Specifies the I-th grid to consider.
  * \tparam Callback Generic callback type.
  */
  template <size_t I, typename Callback>
  void forEachIndex(Callback func) const
  {
    static_assert(I >= 0 && I < D);

    const auto& grid = _data.get<0>(I);
    auto& size = grid->size();
    auto index = Index<D>((int64_t)(0));
    if constexpr (D == 2)
    {
      for (index.y = 0; index.y < size.y; ++index.y)
        for (index.x = 0; index.x < size.x; ++index.x)
          func(index);
    }
    else if (D == 3)
    {
      for (index.z = 0; index.z < size.z; ++index.z)
        for (index.y = 0; index.y < size.y; ++index.y)
          for (index.x = 0; index.x < size.x; ++index.x)
            func(index);
    }
  }

  /**
  * Invokes the given function \p func for each I-data point.
  *
  * This method invokes the given function \p func for each I-data point.
  * The input parameter for the callback function is the value of a data point.
  *
  * \tparam I Specifies the I-th grid to consider.
  * \tparam Callback Generic callback type.
  */
  template <size_t I, typename Callback>
  void forEach(Callback func) const
  {
    static_assert(I >= 0 && I < D);

    const auto& grid = _data.get<0>(I);
    auto size = grid->size();
    auto index = Index<D>((int64_t)(0));
    if constexpr (D == 2)
    {
      for (index.y = 0; index.y < size.y; ++index.y)
        for (index.x = 0; index.x < size.x; ++index.x)
          func(grid->operator[](index));
    }
    else if (D == 3)
    {
      for (index.z = 0; index.z < size.z; ++index.z)
        for (index.y = 0; index.y < size.y; ++index.y)
          for (index.x = 0; index.x < size.x; ++index.x)
            func(grid->operator[](index));
    }
  }

  // IMPLEMENTAR OPERADORES DIFERENCIAIS

private:
  /** Grid spacing or cell size. */
  vec_type _gridSpacing;
  /** SOA instance that holds the data. */
  Data _data;
};

template <typename real>
using FaceCenteredGrid2 = FaceCenteredGrid<2, real>;

template <typename real>
using FaceCenteredGrid3 = FaceCenteredGrid<3, real>;

} // end namespace cg

#endif // !__FaceCenteredGrid_h

