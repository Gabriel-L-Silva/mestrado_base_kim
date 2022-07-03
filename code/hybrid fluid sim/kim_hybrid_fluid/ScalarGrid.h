#ifndef __CellCenteredGrid_h
#define __CellCenteredGrid_h

#include "geometry/Grid3.h"
#include "ScalarField.h"
#include "LinearArraySampler2.h"

namespace cg
{

/// <summary>
/// Modelo de classe derivada de RegionGrid para uma grade que armazena valores escalares.
/// </summary>
/// <typeparam name="D">Dimensão. Toma valor 2 ou 3</typeparam>
/// <typeparam name="real">Tipo de ponto flutuante</typeparam>
template <size_t D, typename real>
class ScalarGrid: public ScalarField<D, real>, public RegionGrid<D, real, real>
{
public:
  using Base = RegionGrid<D, real, real>;
  using bounds_type = typename RegionGrid<D, real, real>::bounds_type;
  using vec_type = typename ScalarField<D, real>::vec_type;
  // using id_type = Index<D>::base_type;

  ScalarGrid(
    const Index<D>& size,
    const vec_type& gridSpacing,
    const vec_type& origin,
    real initialValue = 0.0f):
    _untouched_bounds(bounds_type{ origin, origin + gridSpacing * vec_type{size} }),
    _untouched_cellSize(gridSpacing),
    Base(bounds_type{ origin, origin + gridSpacing * vec_type{size} }, size),
    _linearSampler(this, gridSpacing, origin)
  {
    for (auto& v : *this)
      v = initialValue;
  }

  virtual ~ScalarGrid()
  {
    // do nothing
  }

  virtual Index<D> dataSize() const = 0;

  virtual vec_type dataOrigin() const = 0;

  vec_type gradientAt(Index<D> index) const
  {
    if constexpr (D == 2)
      return gradient2(*this, _untouched_cellSize, index.x, index.y);
    else
      return gradient3(*this, _untouched_cellSize, index.x, index.y, index.z);

  }

  real laplacianAt(Index<D> index) const
  {
    // TODO
  }

  real sample(const vec_type& x) const override
  {
    return _linearSampler(x);
  }

  vec_type gradient(const vec_type& x) const override
  {
    constexpr int array_size = 1 << D;
    std::array<Index<D>, array_size> indices;
    std::array<real, array_size> weights;
    _linearSampler.getCoordinatesAndWeights(x, indices, weights);

    vec_type result{};
    for (int i = 0; i < array_size; ++i)
    {
      result += weights[i] * gradientAt(indices[i]);
    }
    return result;
  }

  real laplacian(const vec_type& x) const override
  {
    return 0.0f;
    // TODO
  }

  void clear()
  {
    for (auto& v : *this)
      v = static_cast<real>(0.0f);
  }

  vec_type dataPosition(const Index<D>& index) const
  {
    return dataOrigin() + vec_type{ index } *_untouched_cellSize;
  }

  void setData(std::vector<real> d)
  {
    auto it = this->begin();
    auto end = this->end();
    for (auto i = 0; it != end; it++, i++)
      (*it) = d[i];
  }
protected:
  // this instance of bounds_type, different from RegionGrid's bounds, will
  // not use any fatFactor
  bounds_type _untouched_bounds;
  // the same applies to cellSize below
  vec_type _untouched_cellSize;
  // NOTE: another option is to implement this class using
  // simple grid class and inherit from it

private:
  LinearArraySampler<real, D> _linearSampler;
};

} // end namespace cg

#endif // __CellCenteredGrid_h
