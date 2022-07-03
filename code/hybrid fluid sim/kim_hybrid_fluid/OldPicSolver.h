#ifndef __OldPicSolver_h
#define __OldPicSolver_h

#include "PhysicsAnimation.h"
#include "FaceCenteredGrid.h"
#include "LinearArraySampler2.h"
#include "geometry/Grid2.h"
#include "geometry/Index2.h"
#include "geometry/ParticleSystem.h"
#include "VolumeParticleEmitter.h"
#include "Box.h"
#include "GridFractionalSinglePhasePressureSolver.h"
#include "GridFractionalBoundaryConditionSolver.h"
#include "GridBackwardEulerDiffusionSolver.h"
#include "CellCenteredScalarGrid.h"
#include "PointGridHashSearcher.h"
#include "geometry/KNNHelper.h"
#include "utils/Stopwatch.h"

namespace cg
{

#ifdef _DEBUG
#define debug(...) fprintf(stdout, __VA_ARGS__)
#else
#define debug(...) 
#endif

// no futuro PicSolverBase deve receber um grid nos argumentos de template
// esse grid deverá implementar algumas funcoes que o Pic usa pra funcionar
// atualmente está fixo com Grid Regular
template<size_t D, typename real, typename Allocator> class OldPicSolver;

template<typename real>
class OldPicSolver<2, real, ArrayAllocator> : public PhysicsAnimation
{
public:
  using vec_type = Vector2<real>;
  using FlipParticleSystem = ParticleSystem<2, real, ArrayAllocator, vec_type>;
  using grid_type = FaceCenteredGrid2<real>;
  using Searcher = PointGridHashSearcher<2, real, FlipParticleSystem>;

  enum
  {
    POSITION,
    VELOCITY,
  };

  OldPicSolver(const Index2& size, const vec_type& gridSpacing, const vec_type& origin, size_t capacity) :
    _particleSystem(capacity),
    _gravity(0.0f, -9.8f),
    _grid(new grid_type(size, gridSpacing, origin)),
    _searcher(new Searcher(Index2{64LL}, 2.4f * gridSpacing.max() / sqrt(2.0f))),
    _signedDistanceField(size, gridSpacing, origin, math::Limits<real>::inf())
  {
    // do nothing
  }

  virtual ~OldPicSolver()
  {
    // do nothing
  }
    // TODO: Constructor that takes grid info

  const FlipParticleSystem& particleSystem() const { return _particleSystem; };

  FlipParticleSystem& particleSystem() { return _particleSystem; };

  const grid_type& grid() const { return *_grid; }

  const vec_type& gravity() const { return _gravity; };

  void setGravity(const vec_type& gravity) { _gravity = gravity; };

  void setGrid(grid_type* newGrid) { _grid.operator=(newGrid); }

  real cfl(double timeInterval) const;

  real maxCfl() const { return _maxCfl; }

  void setMaxCfl(real cfl) { _maxCfl = math::max(cfl, math::Limits<real>::eps()); }

  const auto particleEmitter() const { return _particleEmitter; }

  void setParticleEmitter(ParticleEmitter<FlipParticleSystem>* newEmitter);

  real viscosityCoefficient() const { return _viscosityCoefficient; }

  void setViscosityCoefficient(real viscosity) { _viscosityCoefficient = math::max(viscosity, real(0.0f)); }

protected:
  vec_type _gravity;
  real _viscosityCoefficient = 0.0f;
  real _maxCfl = 5.0f;
  FlipParticleSystem _particleSystem;
  Reference<ParticleEmitter<FlipParticleSystem>> _particleEmitter;
  Reference<grid_type> _grid;
  Reference<PointGridHashSearcher<2, real, FlipParticleSystem>> _searcher;
  GridData<2, char> _uMarkers;
  GridData<2, char> _vMarkers;
  GridBackwardEulerDiffusionSolver<2, real, false> _diffusionSolver;
  GridFractionalSinglePhasePressureSolver<2, real> _pressureSolver;
  GridFractionalBoundaryConditionSolver<2, real> _boundaryConditionSolver;
  CellCenteredScalarGrid<2, real> _signedDistanceField;
  Reference<Collider<2, real>> _collider;

  /// Base class virtual functions

  void onAdvanceTimeStep(double timeInterval) override;

  void initialize() override
  {
    debug("Initializing...\n");
    updateParticleEmitter(0.0);
  }

  size_t numberOfSubTimeSteps(double timeInterval) const override;

  void applyBoundaryCondition();

  void extrapolateVelocityToAir();

  void extrapolateIntoCollider(CellCenteredScalarGrid<2, real>& grid);

  /// Pic Pipeline functions

  void computeExternalForces(double timeInterval);

  void computeViscosity(double timeInterval);

  void computePressure(double timeInterval);

  virtual void computeAdvection(double timeInterval);

  virtual void moveParticles(double timeInterval);

  void computeGravity(double timeInterval);

  virtual void transferFromParticlesToGrids();

  virtual void transferFromGridsToParticles();

private:
  void buildSignedDistanceField();

  void updateParticleEmitter(double timeInterval);

}; // OldPicSolver

template<typename real>
inline void
OldPicSolver<2, real, ArrayAllocator>::setParticleEmitter(ParticleEmitter<FlipParticleSystem>* newEmitter)
{
  _particleEmitter = newEmitter;
  _particleEmitter->setTarget(_particleSystem);
}

template<typename real>
using PicSolver2 = OldPicSolver<2, real, ArrayAllocator>;

using PicSolver2f = PicSolver2<float>;

template<typename real>
inline real
OldPicSolver<2, real, ArrayAllocator>::cfl(double timeInterval) const
{
  real maxVel = 0.0f;
  forEachIndex<2>(_grid->size(), [&](const Index2& index) {
    auto v = _grid->valueAtCellCenter(index) + timeInterval * _gravity;
    maxVel = math::max(maxVel, v.max());
  });

  return real(maxVel * timeInterval / _grid->gridSpacing().min());
}

template<typename real>
inline void cg::OldPicSolver<2, real, ArrayAllocator>::transferFromParticlesToGrids()
{
  using vec_type = Vector2<real>;

  const auto& spacing = this->_grid->gridSpacing();
  const auto& origin = this->_grid->origin();

  LinearArraySampler2<real> uSampler(
    (this->_grid->data<0>()),
    spacing,
    origin
  );
  LinearArraySampler2<real> vSampler(
    (this->_grid->data<1>()),
    spacing,
    origin
  );

  auto uSize = this->_grid->iSize<0>();
  auto vSize = this->_grid->iSize<1>();

  // fill velocity with zero
  this->_grid->fill(vec_type::null());
  GridData<2, real> uWeight;
  GridData<2, real> vWeight;
  uWeight.resize(uSize);
  vWeight.resize(vSize);
  
  _uMarkers.resize(uSize);
  _vMarkers.resize(vSize);

  fill(uWeight, (real)0.0f);
  fill(_uMarkers, (char)0);
  fill(vWeight, (real)0.0f);
  fill(_vMarkers, (char)0);

  auto numberOfParticles = this->_particleSystem.size();
  for (size_t i = 0; i < numberOfParticles; ++i)
  {
    std::array<Index2, 4> indices;
    std::array<real, 4> weights;

    uSampler.getCoordinatesAndWeights(this->_particleSystem[i], indices, weights);
    for (int j = 0; j < 4; ++j)
    {
      _grid->velocityAt<0>(indices[j])
        += _particleSystem.get<VELOCITY>(i).x * weights[j];
      uWeight[uWeight.id(indices[j])] += weights[j];
      _uMarkers[_uMarkers.id(indices[j])] = 1;
    }

    vSampler.getCoordinatesAndWeights(this->_particleSystem[i], indices, weights);
    for (int j = 0; j < 4; ++j)
    {
      _grid->velocityAt<1>(indices[j])
        += _particleSystem.get<VELOCITY>(i).y * weights[j];
      vWeight[vWeight.id(indices[j])] += weights[j];
      _vMarkers[_vMarkers.id(indices[j])] = 1;
    }
  }

  /*printf("[*] V-Grid velocities AFTER transfer\n");
  for (int64_t y = vSize.y - 1; y >= 0; --y)
  {
    for (int64_t x = 0; x < vSize.x; ++x)
    {
      printf("[%f]", _grid->velocityAt<1>(Index2{ x, y }));
    }
    printf("\n");
  }*/

  for (int64_t i = 0; i < uSize.prod(); ++i)
    if (uWeight[i] > 0.0f)
      _grid->velocityAt<0>(i) /= uWeight[i];

  for (int64_t i = 0; i < vSize.prod(); ++i)
    if (vWeight[i] > 0.0f)
      _grid->velocityAt<1>(i) /= vWeight[i];
}

template<typename real>
inline void cg::OldPicSolver<2, real, ArrayAllocator>::transferFromGridsToParticles()
{
  auto numberOfParticles = _particleSystem.size();
  auto& particles = _particleSystem;
  for (size_t i = 0; i < numberOfParticles; ++i)
  {
    particles.get<VELOCITY>(i) = _grid->sample(particles[i]);
  }
}

template<typename real>
void OldPicSolver<2, real, ArrayAllocator>::computeExternalForces(double timeInterval)
{
  //printf("->Computando forças externas...\n");
  computeGravity(timeInterval);
}

template<typename real>
void OldPicSolver<2, real, ArrayAllocator>::computeGravity(double timeInterval)
{
  if (this->_gravity.squaredNorm() > math::Limits<real>::eps())
  {
    auto uSize = this->_grid->iSize<0>();
    auto vSize = this->_grid->iSize<1>();

    if (!math::isZero(this->_gravity.x))
    {
      for (int64_t i = 0; i < uSize.prod(); ++i)
        _grid->velocityAt<0>(i) += timeInterval * _gravity.x;
    }

    if (!math::isZero(this->_gravity.y))
    {
      for (int64_t i = 0; i < vSize.prod(); ++i)
        _grid->velocityAt<1>(i) += timeInterval * _gravity.y;
    }

    /*printf("[*] V-Grid velocities AFTER gravity applied\n");
    for (int64_t y = vSize.y - 1; y >= 0; --y)
    {
      for (int64_t x = 0; x < vSize.x; ++x)
      {
        printf("[%f]", _grid->velocityAt<1>(Index2{ x, y }));
      }
      printf("\n");
    }*/

    applyBoundaryCondition();
  }
}

template<typename real>
void OldPicSolver<2, real, ArrayAllocator>::computeViscosity(double timeInterval)
{
  Stopwatch s;
  if (math::isPositive(_viscosityCoefficient))
  {
    s.start();
    _diffusionSolver.solve(
      _grid,
      _viscosityCoefficient,
      timeInterval,
      _grid,
      *_boundaryConditionSolver.colliderSdf(),
      _signedDistanceField
    );
    debug("[INFO] Viscosity solver took %lld ms\n", s.time());

    applyBoundaryCondition();
  }
}

template<typename real>
void OldPicSolver<2, real, ArrayAllocator>::computePressure(double timeInterval)
{
  Stopwatch s;
  s.start();
  _pressureSolver.solve(
    _grid,
    timeInterval,
    _grid,
    *_boundaryConditionSolver.colliderSdf(),
    _signedDistanceField,
    *_boundaryConditionSolver.colliderVelocityField()
  );
  debug("[INFO] Pressure solver took %lld ms\n", s.time());

  applyBoundaryCondition();
}

template<typename real>
void OldPicSolver<2, real, ArrayAllocator>::computeAdvection(double timeInterval)
{
  Stopwatch s;
  s.start();
  extrapolateVelocityToAir();
  debug("[INFO] ExtrapolateVelocityToAir took %lld ms\n", s.lap());

  applyBoundaryCondition();

  transferFromGridsToParticles();
  debug("[INFO] TransferFromGridsToParticles took %lld ms\n", s.lap());

  moveParticles(timeInterval);
  debug("[INFO] MoveParticles took %lld ms\n", s.lap());
}

template<typename real>
inline void OldPicSolver<2, real, ArrayAllocator>::moveParticles(double timeInterval)
{
  auto numberOfParticles = _particleSystem.size();
  
  for (size_t i = 0; i < numberOfParticles; ++i)
  {
    vec_type pt0{ _particleSystem[i] };
    vec_type pt1{ pt0 };
    vec_type vel{ _particleSystem.get<VELOCITY>(i) };
    const auto& bounds = _grid->bounds();
    vec_type boundsMin{ bounds.min() };
    vec_type boundsMax{ bounds.max() };

    // TODO: adaptive time-stepping
    unsigned int numSubSteps
      = static_cast<unsigned int>(math::max<real>(_maxCfl, 1.0f));
    real dt = static_cast<real>(timeInterval / numSubSteps);
    for (unsigned t = 0; t < numSubSteps; ++t)
    {
      vec_type vel0{ _grid->sample(pt0) };

      // mid-point rule
      vec_type midPt{ pt0 + 0.5 * dt * vel0 };
      vec_type midVel{ _grid->sample(midPt) };
      pt1 = pt0 + dt * midVel;

      pt0 = pt1;
    }

    // assuming that closed domain flag is set to ALL directions
    if (pt1.x - boundsMin.x <= static_cast<real>(0))
    {
      pt1.x = boundsMin.x;
      vel.x = 0.0;
    }
    if (math::isPositive(pt1.x - boundsMax.x))
    {
      pt1.x = boundsMax.x;
      vel.x = 0.0;
    }
    if (pt1.y - boundsMin.y <= static_cast<real>(0))
    {
      pt1.y = boundsMin.y;
      vel.y = 0.0;
    }
    if (math::isPositive(pt1.y - boundsMax.y))
    {
      pt1.y = boundsMax.y;
      vel.y = 0.0;
    }

    _particleSystem.set(i, pt1, vel);
  }

  auto col = _boundaryConditionSolver.collider();
  if (col != nullptr)
  {
    for (size_t i = 0; i < numberOfParticles; ++i)
    {
      col->resolveCollision(
        0.0f,
        0.0f,
        _particleSystem[i],
        _particleSystem.get<VELOCITY>(i)
      );
    }
  }
}

template<typename real>
void OldPicSolver<2, real, ArrayAllocator>::onAdvanceTimeStep(double timeInterval)
{
  Stopwatch s;
  s.start();

  if (_collider != nullptr)
    _collider->update(this->currentTime(), timeInterval);

  _boundaryConditionSolver.updateCollider(
    _collider,
    _grid->size(),
    _grid->gridSpacing(),
    _grid->origin()
  );

  applyBoundaryCondition();
  
  updateParticleEmitter(timeInterval);

  s.lap();
  transferFromParticlesToGrids();
  debug("Transfer from particles to grids: %lld ms\n", s.lap());

  buildSignedDistanceField();

  extrapolateVelocityToAir();

  applyBoundaryCondition();
  
  // main pipeline
  computeExternalForces(timeInterval);

  computeViscosity(timeInterval);
  
  computePressure(timeInterval);

  computeAdvection(timeInterval);
}

template<typename real>
inline size_t
OldPicSolver<2, real, ArrayAllocator>::numberOfSubTimeSteps(double timeInterval) const
{
  auto _cfl = cfl(timeInterval);
  return static_cast<size_t>(math::max<real>(std::ceil(_cfl / _maxCfl), 1.0f));
}

template<typename real>
inline void OldPicSolver<2, real, ArrayAllocator>::buildSignedDistanceField()
{
  auto maxH = _signedDistanceField.cellSize().max();
  auto radius = 1.2f * maxH / sqrt(2.0f);
  auto sdfBandRadius = 2.0f * radius;

  auto sdfSize = _signedDistanceField.dataSize();
  _searcher->build(_particleSystem);
  forEachIndex<2>(sdfSize, [&](const Index2& index) {
    auto pt = _signedDistanceField.dataPosition(index);
    auto minDist = sdfBandRadius;
    _searcher->forEachNearbyPoint(pt, sdfBandRadius, [&](size_t, const vec_type& x) {
      auto d = (pt - x).length();
      minDist = math::min(minDist, d);
    });
    _signedDistanceField[index] = minDist - radius;
  });

  extrapolateIntoCollider(_signedDistanceField);
}

template<typename real>
inline void
OldPicSolver<2, real, ArrayAllocator>::applyBoundaryCondition()
{
  unsigned int depth = static_cast<unsigned int>(std::ceil(_maxCfl));
  _boundaryConditionSolver.constrainVelocity(_grid, depth);
}

template<typename real>
inline void
OldPicSolver<2, real, ArrayAllocator>::extrapolateVelocityToAir()
{
  unsigned int depth = static_cast<unsigned int>(std::ceil(_maxCfl));
  extrapolateToRegion(*_grid->data<0>(), _uMarkers, depth, *_grid->data<0>());
  extrapolateToRegion(*_grid->data<1>(), _vMarkers, depth, *_grid->data<1>());
}

template<typename real>
inline void
OldPicSolver<2, real, ArrayAllocator>::extrapolateIntoCollider(CellCenteredScalarGrid<2, real>& grid)
{
  GridData<2, char> marker;
  marker.resize(grid.dataSize());

  for (size_t i = 0; i < grid.length(); ++i)
  {
    auto index = grid.index(i);
    auto colliderSdfQuery = _boundaryConditionSolver.colliderSdf()->sample(
      grid.dataPosition(index));
    if (isInsideSdf(colliderSdfQuery))
      marker[i] = 0;
    else
      marker[i] = 1;
  }
  
  unsigned int depth = static_cast<unsigned int>(std::ceil(_maxCfl));
  extrapolateToRegion(grid, marker, depth, grid);
}

template<typename real>
inline void
OldPicSolver<2, real, ArrayAllocator>::updateParticleEmitter(double timeInterval)
{
  if (_particleEmitter != nullptr)
    // TODO CurrentTimeInSeconds
    _particleEmitter->update(0.0, timeInterval);
}

} // end namespace cg

#endif // __OldPicSolver_h
