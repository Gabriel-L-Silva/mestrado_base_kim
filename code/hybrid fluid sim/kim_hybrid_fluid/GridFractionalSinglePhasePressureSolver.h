#ifndef __GridFractionalSinglePhasePressureSolver_h
#define __GridFractionalSinglePhasePressureSolver_h

#include <Eigen/Sparse>
#include "GridPressureSolver.h"
#include "GridUtils.h"

using namespace Eigen;

namespace
{

using namespace cg;

template <size_t D>
auto __id(const cg::Index<D>& size, const cg::Index<D>& index)
{
  static_assert(D == 2 || D == 3, "__id can't be instantiated with D > 3 or D < 2");
  
  auto ret = index.x + index.y * size.x;
  if constexpr (D == 3)
  {
    ret += index.z * size.x * size.y;
  }
  return ret;
}

template <size_t D, typename T = typename cg::Index<D>::base_type>
auto __index(const cg::Index<D>& size, T id)
{
  static_assert(std::is_integral<T>::value, "Integral required.");

  cg::Index<D> i;

  if constexpr (D == 3)
  {
    auto size_xy = size.x * size.y;
    i.z = id / size_xy;
    id -= size_xy * i.z;
  }
  i.y = id / size.x;
  i.x = id - size.x * i.y;

  return i;
}

template <size_t D, typename real>
void buildSingleSystem(
  SparseMatrix<real>& A,
  Eigen::Matrix<real, -1, 1>& b,
  GridData<D, real>& fluidSdf,
  std::array<GridData<D, real>, D>& weights,
  const VectorField<D, real>& boundaryVel,
  const Reference<FaceCenteredGrid<D, real>> input)
{
  using id_type = typename cg::Index<D>::base_type;
  using Triplet = Eigen::Triplet<real, id_type>;

  const auto& size = input->size();
  const auto invH = input->gridSpacing().inverse();
  const auto invHSqr = invH * invH;
  
  // position in space functions
  
  std::array<std::function<cg::Vector<real, D>(const cg::Index<D>&)>, D> pos;
  pos[0] = input->positionInSpace<0>();
  pos[1] = input->positionInSpace<1>();
  if constexpr (D == 3)
    pos[2] = input->positionInSpace<2>();

  std::vector<Triplet> coefficients;

  forEachIndex<D>(
    size,
    [&](const cg::Index<D>& index) {
      auto id = __id<D>(size, index);
      auto centerPhi = fluidSdf[id];
      b(id) = 0.0f;

     real centerValue = static_cast<real>(0.0f);
      if (isInsideSdf(centerPhi))
      {
        real boundaryCondition = 0.0f;
        for (int k = 0; k < D; ++k) // iterate through D dimensions
        {
          real term;
          auto inc = cg::Index<D>((int64_t)0);
          inc[k] = 1;
          auto indexP1 = index + inc;
          auto indexM1 = index - inc;
          auto iP1 = weights[k].id(indexP1);
          auto iM1 = weights[k].id(indexM1);
          
          // dont even ask me....
          boundaryCondition += (1.0f - weights[k][iP1]) *
            boundaryVel.sample(pos[k](indexP1))[k] * invH[k] - 
            (1.0f - weights[k][id]) * boundaryVel.sample(pos[k](index))[k] * invH[k];

          if (index[k] + 1 < size[k])
          {
            term = weights[k][iP1] * invHSqr[k];
            auto phi = fluidSdf[fluidSdf.id(indexP1)];
            if (isInsideSdf(phi))
            {
              centerValue += term;
              // right neighbor in K direction
              coefficients.push_back(Triplet(id, __id(size, indexP1), -term));
            }
            else
            {
              real theta = math::max(
                fractionInsideSdf(centerPhi, phi), static_cast<real>(0.01f));
              centerValue += term / theta;
            }
            b(id) += weights[k][iP1] * input->velocityAt(k, indexP1) * invH[k];
          }
          else
          {
            b(id) += input->velocityAt(k, indexP1) * invH[k];
          }

          if (index[k] > 0)
          {
            auto wId = weights[k].id(index);
            term = weights[k][wId] * invHSqr[k];
            auto phi = fluidSdf[fluidSdf.id(indexM1)];
            if (isInsideSdf(phi))
            {
              centerValue += term;
              // left neighbor in K direction
              coefficients.push_back(Triplet(id, __id(size, indexM1), -term));
            }
            else
            {
              real theta = math::max(
                fractionInsideSdf(centerPhi, phi), static_cast<real>(0.01f));
              centerValue += term / theta;
            }
            b(id) -= weights[k][wId] * input->velocityAt(k, index) * invH[k];
          }
          else
          {
            b(id) -= input->velocityAt(k, index) * invH[k];
          }
        }
        
        // accumulate contributions from the moving boundary
        b(id) += boundaryCondition;

        // if centerValue is near-zero, the cell is likely inside a solid
        // boundary
        if (math::isZero(centerValue))
        {
          centerValue = 0.0f;
          b(id) = 0.0f;
        }

        coefficients.push_back(Triplet(id, id, centerValue));
      }
      else
        centerValue = 1.0f;
    }
  );

  A.setFromTriplets(coefficients.begin(), coefficients.end());
}

}; // end anonymous namespace

namespace cg
{

template <size_t D, typename real>
class GridFractionalSinglePhasePressureSolverBase : public GridPressureSolver<D, real>
{
public:
  using FCG = FaceCenteredGrid<D, real>;
  using FCGref = Reference<FCG>;
  using vec_type = Vector<real, D>;
  using id_type = typename Index<D>::base_type;
  using ScalarFieldType = ScalarField<D, real>;
  using VectorFieldType = VectorField<D, real>;

  GridFractionalSinglePhasePressureSolverBase()
  {
    // do nothing
  }

  virtual ~GridFractionalSinglePhasePressureSolverBase()
  {
    // do nothing
  }

  void solve(
    const FCGref input,
    double timeInterval,
    FCGref dest,
    const ScalarFieldType& boundarySdf = ConstantScalarField<D, real>(math::Limits<real>::inf()),
    const ScalarFieldType& fluidSdf = ConstantScalarField<D, real>(-math::Limits<real>::inf()),
    const VectorFieldType& boundaryVelocity = ConstantVectorField<D, real>(vec_type{ real(0.0f) })) override;

protected:
  // system matrix
  SparseMatrix<real> A;
  // solution vector
  Eigen::Matrix<real, -1, 1> x;
  // RHS vector
  Eigen::Matrix<real, -1, 1> b;
  // system solver
  ConjugateGradient<SparseMatrix<real>, Lower | Upper> solver;
  // array to hold the weights
  std::array<GridData<D, real>, D> _weights;

  GridData<D, real> _fluidSdf;

  virtual void buildWeights(
    const FCGref input,
    const ScalarFieldType& boundarySdf,
    const ScalarFieldType& fluidSdf,
    const VectorFieldType& boundaryVelocity
  ) = 0;

  void buildSystem(const FCGref input, const VectorFieldType& boundaryVelocity);

  void applyPressureGradient(FCGref input, FCGref dest);

  enum kMarkers
  {
    Fluid,
    Air,
    Boundary
  };

}; // GridFractionalSinglePhaseSolverBase

template<size_t D, typename real>
inline void
GridFractionalSinglePhasePressureSolverBase<D, real>::solve(
  const FCGref input,
  double timeInterval,
  FCGref dest,
  const ScalarFieldType& boundarySdf,
  const ScalarFieldType& fluidSdf,
  const VectorFieldType& boundaryVelocity)
{
  buildWeights(input, boundarySdf, fluidSdf, boundaryVelocity);
  buildSystem(input, boundaryVelocity);

  solver.compute(A);
  x = solver.solve(b);
  auto info = solver.info();
  if (info != Eigen::Success)
  {
    std::cout << "Pressure solver problem: " << info << "\n";
    std::cout << "Solver error: " << solver.error() << '\n';
    std::cout << "Solver max error: " << solver.tolerance() << '\n';
    std::cout << "Solver iterations: " << solver.iterations() << '\n';
    x.setZero();
    //throw std::runtime_error();
  }

  applyPressureGradient(input, dest);
}

template<size_t D, typename real>
inline void
GridFractionalSinglePhasePressureSolverBase<D, real>::buildSystem(const FCGref input, const VectorFieldType& boundaryVelocity)
{
  // Not considering the use of multi-grid solvers
  const auto& size = input->size();
  auto numberOfCells = (Eigen::Index) size.prod();
  A.resize(numberOfCells, numberOfCells);
  A.data().squeeze(); // release as much memory as possible
  b.resize(numberOfCells);

  buildSingleSystem(A, b, _fluidSdf, _weights, boundaryVelocity, input);
}

template<size_t D, typename real>
inline void
GridFractionalSinglePhasePressureSolverBase<D, real>::applyPressureGradient(FCGref input, FCGref dest)
{
  const auto& size = input->size();

  auto valueAt = [](const GridData<D, real>& grid, const Index<D>& index) {
    return grid[grid.id(index)];
  };

  std::vector<Index<D>> nbrs(D, Index<D>(0LL));
  auto invH = input->gridSpacing().inverse();
  auto _s = x.size();
  for (Eigen::Index i = 0; i < x.size(); ++i)
  {
    auto index = __index<D, Eigen::Index>(size, i);
    auto centerPhi = _fluidSdf[i];
    const auto centerPhiInside = isInsideSdf(centerPhi);

    // compute right, up and front neighbors
    for (int j = 0; j < D; ++j)
    {
      nbrs[j] = Index<D>(0LL);
      nbrs[j][j] = 1;
      nbrs[j] = nbrs[j] + index;
    }

    if (nbrs[0].x < size.x && valueAt(_weights[0], nbrs[0]) > 0.0f &&
      (centerPhiInside ||
        isInsideSdf(valueAt(_fluidSdf, nbrs[0]))))
    {
      auto rightPhi = valueAt(_fluidSdf, nbrs[0]);
      auto theta = fractionInsideSdf(centerPhi, rightPhi);
      theta = math::max(theta, (real)0.01f);

      dest->velocityAt<0>(nbrs[0]) = input->velocityAt<0>(nbrs[0]) + invH.x / theta * (x(__id(size, nbrs[0])) - x(i));
    }

    if (nbrs[1].y < size.y && valueAt(_weights[1], nbrs[1]) > 0.0f && (centerPhiInside || isInsideSdf(valueAt(_fluidSdf, nbrs[1]))))
    {
      auto upPhi = valueAt(_fluidSdf, nbrs[1]);
      auto theta = fractionInsideSdf(centerPhi, upPhi);
      theta = math::max(theta, (real)0.01f);

      dest->velocityAt<1>(nbrs[1]) = input->velocityAt<1>(nbrs[1]) + invH.y / theta * (x(__id(size, nbrs[1])) - x(i));
    }

    if constexpr (D == 3)
    {
      if (nbrs[2].z < size.z && valueAt(_weights[2], nbrs[2]) > 0.0f && (centerPhiInside || isInsideSdf(valueAt(_fluidSdf, nbrs[2]))))
      {
        auto frontPhi = valueAt(_fluidSdf, nbrs[2]);
        auto theta = fractionInsideSdf(centerPhi, frontPhi);
        theta = math::max(theta, (real)0.01f);

        dest->velocityAt<2>(nbrs[2]) = input->velocityAt<2>(nbrs[2]) + invH.y / theta * (x(__id(size, nbrs[2])) - x(i));
      }
    }
  }
}

template <size_t D, typename real> class GridFractionalSinglePhasePressureSolver;

template <typename real>
class GridFractionalSinglePhasePressureSolver<2, real>: public GridFractionalSinglePhasePressureSolverBase<2, real>
{
public:
  using Base = GridFractionalSinglePhasePressureSolverBase<2, real>;
  using FCG = typename Base::FCG;
  using FCGref = Reference<FCG>;
  using vec_type = typename Base::vec_type;
  using id_type = Index2::base_type;
  using ScalarFieldType = typename Base::ScalarFieldType;
  using VectorFieldType = typename Base::VectorFieldType;

  GridFractionalSinglePhasePressureSolver()
    : Base::GridFractionalSinglePhasePressureSolverBase()
  {
    // do nothing
  }

protected:
  void buildWeights(
    const FCGref input,
    const ScalarFieldType& boundarySdf,
    const ScalarFieldType& fluidSdf,
    const VectorFieldType& boundaryVelocity
  ) override
  {
    auto size = input->size();
    // @note we are EXCLUDING the multigrid functionality for the sake of simplicity
    // Build levels
    this->_fluidSdf.resize(size);
    this->_weights[0].resize(size + Index2{ 1, 0 });
    this->_weights[1].resize(size + Index2{ 0, 1 });

    auto cellPos = input->cellCenterPosition();
    auto h = input->gridSpacing();
    auto uPos = input->positionInSpace<0>();
    auto vPos = input->positionInSpace<1>();

    // build Top level grid
    forEachIndex<2>(size, [&](const Index2& index) {
      this->_fluidSdf[this->_fluidSdf.id(index)] = static_cast<real>(
        fluidSdf.sample(cellPos(index)));
    });

    // u
    forEachIndex<2>(
      this->_weights[0].size(),
      [&](const Index2& index) {
        auto id = this->_weights[0].id(index);
        auto pt = uPos(index);
        real phi0 = boundarySdf.sample(pt - vec_type(0.5f * h.x, 0.0f));
        real phi1 = boundarySdf.sample(pt + vec_type(0.5f * h.x, 0.0f));
        real frac = fractionInsideSdf(phi0, phi1);
        real weight = math::clamp<real>(1.0f - frac, 0.0f, 1.0f);

        // Clamp non-zero weight to kMinWeight. Having nearly-zero element
        // in the matrix can be an issue.
        if (weight < 0.01f && weight > 0.0f) {
          weight = 0.01f;
        }

        this->_weights[0][id] = weight;
      }
    );

    // v
    forEachIndex<2>(
      this->_weights[1].size(),
      [&](const Index2& index) {
        auto id = this->_weights[1].id(index);
        auto pt = uPos(index);
        real phi0 = boundarySdf.sample(pt - vec_type(0.0f, 0.5f * h.y));
        real phi1 = boundarySdf.sample(pt + vec_type(0.0f, 0.5f * h.y));
        real frac = fractionInsideSdf(phi0, phi1);
        real weight = math::clamp<real>(1.0f - frac, 0.0f, 1.0f);

        // Clamp non-zero weight to kMinWeight. Having nearly-zero element
        // in the matrix can be an issue.
        if (weight < 0.01f && weight > 0.0f) {
          weight = 0.01f;
        }

        this->_weights[1][id] = weight;
      }
    );
  }
}; // GridFractionalSinglePhasePressureSolver<2, real>

template <typename real>
class GridFractionalSinglePhasePressureSolver<3, real> : public GridFractionalSinglePhasePressureSolverBase<3, real>
{
public:
  using Base = GridFractionalSinglePhasePressureSolverBase<3, real>;
  using FCG = typename Base::FCG;
  using FCGref = Reference<FCG>;
  using vec_type = typename Base::vec_type;
  using id_type = Index2::base_type;
  using ScalarFieldType = typename Base::ScalarFieldType;
  using VectorFieldType = typename Base::VectorFieldType;

  GridFractionalSinglePhasePressureSolver()
    : Base::GridFractionalSinglePhasePressureSolverBase()
  {
    // do nothing
  }
protected:
  void buildWeights(
    const FCGref input,
    const ScalarFieldType& boundarySdf,
    const ScalarFieldType& fluidSdf,
    const VectorFieldType& boundaryVelocity
  ) override
  {
    auto size = input->size();
      // @note we are EXCLUDING the multigrid functionality for the sake of simplicity
      // Build levels
    this->_fluidSdf.resize(size);
    this->_fluidSdf.initialize(size);
    this->_weights[0].resize(size + Index3{ 1, 0, 0 });
    this->_weights[1].resize(size + Index3{ 0, 1, 0 });
    this->_weights[2].resize(size + Index3{ 0, 0, 1 });
    
    auto cellPos = input->cellCenterPosition();
    auto h = input->gridSpacing();
    auto uPos = input->positionInSpace<0>();
    auto vPos = input->positionInSpace<1>();
    auto wPos = input->positionInSpace<2>();

    // build Top level grid
    forEachIndex<3>(size, [&](const Index3& index) {
      this->_fluidSdf[this->_fluidSdf.id(index)] = static_cast<real>(
        fluidSdf.sample(cellPos(index)));
    });

    // u
    forEachIndex<3>(
      this->_weights[0].size(),
      [&](const Index3& index) {
        auto id = this->_weights[0].id(index);
        auto pt = uPos(index);
        real phi0 = boundarySdf.sample(pt + vec_type(0.0f, -0.5f * h.y, -0.5f * h.z));
        real phi1 = boundarySdf.sample(pt + vec_type(0.0f,  0.5f * h.y, -0.5f * h.z));
        real phi2 = boundarySdf.sample(pt + vec_type(0.0f, -0.5f * h.y,  0.5f * h.z));
        real phi3 = boundarySdf.sample(pt + vec_type(0.0f,  0.5f * h.y,  0.5f * h.z));
        real frac = fractionInside(phi0, phi1, phi2, phi3);
        real weight = math::clamp(1.0f - frac, 0.0f, 1.0f);

        if (weight < 0.01f && weight > 0.0f)
          weight = 0.01f;

        this->_weights[0][id] = weight;
      }
    );

    // v
    forEachIndex<3>(
      this->_weights[1].size(),
      [&](const Index3& index) {
        auto id = this->_weights[1].id(index);
        auto pt = vPos(index);
        real phi0 = boundarySdf.sample(pt + vec_type(-0.5f * h.x, 0.0f, -0.5f * h.z));
        real phi1 = boundarySdf.sample(pt + vec_type(-0.5f * h.x, 0.0f,  0.5f * h.z));
        real phi2 = boundarySdf.sample(pt + vec_type( 0.5f * h.x, 0.0f, -0.5f * h.z));
        real phi3 = boundarySdf.sample(pt + vec_type( 0.5f * h.x, 0.0f,  0.5f * h.z));
        real frac = fractionInside(phi0, phi1, phi2, phi3);
        real weight = math::clamp(1.0f - frac, 0.0f, 1.0f);

        if (weight < 0.01f && weight > 0.0f)
          weight = 0.01f;

        this->_weights[1][id] = weight;
      }
    );

    // w
    forEachIndex<3>(
      this->_weights[2].size(),
      [&](const Index3& index) {
        auto id = this->_weights[2].id(index);
        auto pt = uPos(index);
        real phi0 = boundarySdf.sample(pt + vec_type(-0.5f * h.x, -0.5f * h.y, 0.0f));
        real phi1 = boundarySdf.sample(pt + vec_type(-0.5f * h.x,  0.5f * h.y, 0.0f));
        real phi2 = boundarySdf.sample(pt + vec_type( 0.5f * h.x, -0.5f * h.y, 0.0f));
        real phi3 = boundarySdf.sample(pt + vec_type( 0.5f * h.x,  0.5f * h.y, 0.0f));
        real frac = fractionInside(phi0, phi1, phi2, phi3);
        real weight = math::clamp(1.0f - frac, 0.0f, 1.0f);

        if (weight < 0.01f && weight > 0.0f)
          weight = 0.01f;

        this->_weights[2][id] = weight;
      }
    );
  }
}; // GridFractionalSinglePhasePressureSolver<3, real>

} // end namespace cg

#endif // __GridFractionalSinglePhasePressureSolver_h
