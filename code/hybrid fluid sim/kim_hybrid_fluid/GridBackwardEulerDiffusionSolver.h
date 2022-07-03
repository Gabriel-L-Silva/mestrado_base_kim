#ifndef __GridBackwardEulerDiffusionSolver_h
#define __GridBackwardEulerDiffusionSolver_h

#include <Eigen/Sparse>
#include "GridDiffusionSolver.h"
#include "GridUtils.h"
#include "MathUtils.h"

using namespace Eigen;
namespace cg
{

template <size_t D, typename real, bool isDirichlet = false>
class GridBackwardEulerDiffusionSolver : public GridDiffusionSolver<D, real>
{
public:
  using FCG = FaceCenteredGrid<D, real>;
  using FCGref = Reference<FCG>;
  using vec_type = Vector<real, D>;
  using id_type = Index<3>::base_type;

  explicit GridBackwardEulerDiffusionSolver()
  {
    // do nothing
  }

  virtual ~GridBackwardEulerDiffusionSolver()
  {
    // do nothing
  }

  void solve(
    const FCGref source,
    real diffusionCoefficient,
    double timeInterval,
    FCGref dest,
    const ScalarField<D, real>& boundarySdf
    = ConstantScalarField<D, real>(0.0),
    const ScalarField<D, real>& fluidSdf
    = ConstantScalarField<D, real>(0.0)
  );

  //void setLinearSystemSolver(const FdmLinearSystemSolver2Ptr& solver);

private:
  // system matrix
  SparseMatrix<real> A;
  // solution vector
  Eigen::Matrix<real, -1, 1> x;
  // RHS vector
  Eigen::Matrix<real, -1, 1> b;
  // system solver
  ConjugateGradient<SparseMatrix<real>, Lower | Upper> solver;
  /*FdmLinearSystemSolver2Ptr _systemSolver;*/
  GridData<D, char> _markers;

  void buildMarkers(
    const Index<D>& size,
    const std::function<vec_type(const Index<D>&)>& pos,
    const ScalarField<D, real>& boundarySdf,
    const ScalarField<D, real>& fluidSdf);

  void buildMatrix(
    const Index<D>& size,
    const vec_type& c);

  template <size_t I>
  void buildVectors(const FCGref source, const vec_type& c);

  void createCoefficient(const id_type row, const Index<D>& index, real& center, const real inc, std::vector<Eigen::Triplet<real, id_type>>& coefficients);

  enum kMarkers
  {
    Fluid,
    Air,
    Boundary
  };
};

template<size_t D, typename real, bool isDirichlet>
inline void
GridBackwardEulerDiffusionSolver<D, real, isDirichlet>::solve(const FCGref source, real diffusionCoefficient, double timeInterval, FCGref dest, const ScalarField<D, real>& boundarySdf, const ScalarField<D, real>& fluidSdf)
{
  auto h = (source->gridSpacing() * source->gridSpacing()).inverse();
  auto c = (timeInterval * diffusionCoefficient) * h;

  // u
  auto uPos = source->positionInSpace<0>();
  auto uSize = source->iSize<0>();
  buildMarkers(uSize, uPos, boundarySdf, fluidSdf);
  buildMatrix(uSize, c);
  buildVectors<0>(source, c);

  solver.compute(A);
  x = solver.solve(b);
  auto info = solver.info();
  assert(info == 0); // asserts success
  for (int64_t i = 0; i < uSize.prod(); ++i)
    dest->velocityAt<0>(i) = x[i];

  // v
  auto vPos = source->positionInSpace<1>();
  auto vSize = source->iSize<1>();
  buildMarkers(vSize, vPos, boundarySdf, fluidSdf);
  buildMatrix(vSize, c);
  buildVectors<1>(source, c);

  solver.compute(A);
  x = solver.solve(b);
  info = solver.info();
  assert(info == 0); // asserts success
  for (int64_t i = 0; i < vSize.prod(); ++i)
    dest->velocityAt<1>(i) = x[i];

  if constexpr (D == 3)
  {
    // w
    auto wPos = source->positionInSpace<2>();
    auto wSize = source->iSize<2>();
    buildMarkers(wSize, wPos, boundarySdf, fluidSdf);
    buildMatrix(wSize, c);
    buildVectors<2>(source, c);

    solver.compute(A);
    x = solver.solve(b);
    info = solver.info();
    assert(info == 0); // asserts success
    for (int64_t i = 0; i < wSize.prod(); ++i)
      dest->velocityAt<2>(i) = x[i];
  }
}

template<size_t D, typename real, bool isDirichlet>
inline void cg::GridBackwardEulerDiffusionSolver<D, real, isDirichlet>::buildMarkers(
  const Index<D>& size,
  const std::function<vec_type(const Index<D>&)>& pos,
  const ScalarField<D, real>& boundarySdf,
  const ScalarField<D, real>& fluidSdf)
{
  _markers.resize(size);
  // BUG: em GridData<3, T> temos alguns erros, um deles acontece ao efetuarmos resize()
  // o valor de _size_xy nao é alterado condizentemente, acarretando em erros em id() e index()
  // alem disso, temos erros em index(), em initialize e no construtor por copia
  if constexpr (D == 3)
    _markers.initialize(size);
  forEachIndex<D>(size, [&](const Index<D>& index) {
      if (isInsideSdf(boundarySdf.sample(pos(index))))
      {
        _markers[_markers.id(index)] = kMarkers::Boundary;
      }
      else if (isInsideSdf(fluidSdf.sample(pos(index))))
      {
        _markers[_markers.id(index)] = kMarkers::Fluid;
      }
      else
      {
        _markers[_markers.id(index)] = kMarkers::Air;
      }
    }
  );
}

template<size_t D, typename real, bool isDirichlet>
inline void cg::GridBackwardEulerDiffusionSolver<D, real, isDirichlet>::buildMatrix(
  const Index<D>& size, const vec_type& c)
{
  using Triplet = Eigen::Triplet<real, id_type>;
  auto numberOfCells = (Eigen::Index) size.prod();
  A.resize(numberOfCells, numberOfCells);
  Index<D> index;

  std::vector<Triplet> coefficients;

  // array to store cell neighbors in x+-, y+- and z+- respectively
  std::array<Index<D>, 2 * D> nbrs;
  //coefficients.reserve(); // estimate non zero values and reserve memory to speed things up
  if constexpr (D == 2)
  {
    forEachIndex<2>(size,
      [&](const Index2& index) {
        auto cellId = _markers.id(index);
        real centerValue = 1.0f;
        neighborCellIndexes(index, nbrs);

        // _markers[index] == Fluid
        if (_markers[cellId] == kMarkers::Fluid)
        {
          // has right neighbor
          if (index.x + 1 < size.x)
            createCoefficient(cellId, nbrs[0], centerValue, c.x, coefficients);

          // has left neighbor
          if (index.x > 0)
            createCoefficient(cellId, nbrs[1], centerValue, c.x, coefficients);

          // has up neighbor
          if (index.y + 1 < size.y)
            createCoefficient(cellId, nbrs[2], centerValue, c.y, coefficients);

          // has down neighbor
          if (index.y > 0)
            createCoefficient(cellId, nbrs[3], centerValue, c.y, coefficients);
        }

        // adding diagonal coefficient
        coefficients.push_back(Triplet(cellId, cellId, centerValue));
      }
    );
  }
  else if constexpr (D == 3)
  {
    forEachIndex<3>(size,
      [&](const Index3& index) {
        auto cellId = _markers.id(index);
        real centerValue = static_cast<real>(1.0);
        // _markers[index] == Fluid
        if (_markers[cellId] == kMarkers::Fluid)
        {
          // try to code a compile time for so you can use the index as template parameter
          // has right neighbor

          //TODO!!

          if (index.x + 1 < size.x)
            createCoefficient<0>(cellId, index + Index3(1, 0, 0), centerValue, coefficients, c);

          // has left neighbor
          if (index.x > 0)
            createCoefficient<0>(cellId, index - Index3(1, 0, 0), centerValue, coefficients, c);

          // has up neighbor
          if (index.y + 1 < size.y)
            createCoefficient<1>(cellId, index + Index3(0, 1, 0), centerValue, coefficients, c);

          // has down neighbor
          if (index.y > 0)
            createCoefficient<1>(cellId, index - Index3(0, 1, 0), centerValue, coefficients, c);

          // has front neighbor
          if (index.z + 1 < size.z)
            createCoefficient<2>(cellId, index + Index3(0, 0, 1), centerValue, coefficients, c);

          // has back neighbor
          if (index.z > 0)
            createCoefficient<2>(cellId, index - Index3(0, 0, 1), centerValue, coefficients, c);
        }
        // adding diagonal coefficient
        coefficients.push_back(Triplet(cellId, cellId, centerValue));
      }
    );
  }
  /*std::cout << "Coefficients:\n";
  for (auto& t : coefficients)
  {
    std::cout << "Row: " << t.row() << " Col: " << t.col() << "v: " << t.value() << '\n';
  }
  std::cout << "------------------------------------------------------------\n";*/
  // if this is taking too long, try assigning values directly into matrix
  A.setFromTriplets(coefficients.begin(), coefficients.end());
}

template<size_t D, typename real, bool isDirichlet>
template<size_t I>
inline void cg::GridBackwardEulerDiffusionSolver<D, real, isDirichlet>::buildVectors(const FCGref source, const vec_type& c)
{
  const auto f = source->data<I>();
  const auto& size = source->iSize<I>();
  auto numberOfCells = (Eigen::Index)size.prod();
  b.resize(numberOfCells);
  x.resize(numberOfCells); // is this necessary?

  std::array<id_type, 2 * D> neighborIds;

  // for each grid cell
  for (Eigen::Index i = 0; i < numberOfCells; ++i)
  {
    // why do we assign f(id) to x(id)?
    b(i) = x(i) = source->velocityAt<I>(i);

    // if boundary condition is Dirichlet
    if constexpr (isDirichlet)
    {
      if (_markers[i] == kMarkers::Fluid)
      {
        neighborCellIds<true, D, id_type>
          (i, size, neighborIds);
        for (size_t j = 0; j < 2 * D; ++j)
        {
          if (neighborIds[j] < std::numeric_limits<id_type>::max()
            && _markers[neighborIds[j]] == kMarkers::Boundary)
          {
            b(i) += c[j / 2] * source->velocityAt<I>(neighborIds[j]);
          }
        }
      }
    }
  }
}

/*
* Auxiliar method that checks grid cell and updates system matrix accordingly
*/
template<size_t D, typename real, bool isDirichlet>
inline void
GridBackwardEulerDiffusionSolver<D, real, isDirichlet>::createCoefficient(const id_type row, const Index<D>& index, real& center, const real inc, std::vector<Eigen::Triplet<real, id_type>>& coefficients)
{
  auto id = _markers.id(index);
  auto marker = _markers[id];
  if ((isDirichlet && marker != kMarkers::Air) || marker == kMarkers::Fluid)
    center += inc;
  if (_markers[id] == kMarkers::Fluid)
    coefficients.push_back(Eigen::Triplet<real, id_type>(row, id, -inc));
}

} // end namespace cg

#endif // __GridBackwardEulerDiffusionSolver_h
