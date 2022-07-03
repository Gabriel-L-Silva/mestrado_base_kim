#ifndef __PicSolver_h
#define __PicSolver_h

#include "geometry/ParticleSystem.h"
#include "GridFluidSolver.h"
#include "PointGridHashSearcher.h"
#include "ParticleEmitter.h"

namespace cg
{

#ifdef _DEBUG
#define debug(...) fprintf(stdout, __VA_ARGS__)
#else
#define debug(...) 
#endif

template<size_t D, typename real, typename ArrayAllocator>
class PicSolver : public GridFluidSolver<D, real>
{
public:
  using Base = GridFluidSolver<D, real>;
  using vec = Vector<real, D>;
  using PicParticleSystem = ParticleSystem<D, real, ArrayAllocator, vec>;
  using Searcher = PointGridHashSearcher<D, real, PicParticleSystem>;
  template <typename T> using Ref = Reference<T>;

  PicSolver(const Index<D>& size, const vec& spacing, const vec& origin, size_t capacity = 100000)
    : Base{size, spacing, origin}, _particleSystem(capacity)
  {
    _searcher = new Searcher(Index<D>(64LL), 2.4f * spacing.max() / sqrt(2.0f));
    _signedDistanceField = new CellCenteredScalarGrid<D, real>(size, spacing, origin, math::Limits<real>::inf());
  }

  virtual ~PicSolver()
  {
    // do nothing
  }

  auto signedDistanceField() const { return _signedDistanceField; }

  const auto& particleSystem() const { return _particleSystem; }

  const auto particleEmitter() const { return _particleEmitter; }

  void setParticleEmitter(ParticleEmitter<PicParticleSystem>* emitter)
  {
    _particleEmitter = emitter;
    _particleEmitter->setTarget(_particleSystem);
  }

protected:
  std::array<GridData<D, char>, D> _markers;
  PicParticleSystem _particleSystem;

  void initialize() override;

  void onBeginAdvanceTimeStep(double timeInterval) override;

  void computeAdvection(double timeInterval) override;

  ScalarField<D, real>* fluidSdf() const override;

  // Transfers velocity field from particles to grids.
  virtual void transferFromParticlesToGrids();

  // Transfers velocity field from grids to particles.
  virtual void transferFromGridsToParticles();

  // Moves particles.
  virtual void moveParticles(double timeIntervalInSeconds);

private:
  Ref<ParticleEmitter<PicParticleSystem>> _particleEmitter;
  Ref<Searcher> _searcher;
  Ref<CellCenteredScalarGrid<D, real>> _signedDistanceField;

  void extrapolateVelocityToAir();

  void buildSignedDistanceField();

  void updateParticleEmitter(double timeInterval);

}; // PicSolver<D, real, ArrayAllocator>

template<size_t D, typename real, typename ArrayAllocator>
inline void
PicSolver<D, real, ArrayAllocator>::initialize()
{
  GridFluidSolver<D, real>::initialize();

  updateParticleEmitter(0.0);
}

template<size_t D, typename real, typename ArrayAllocator>
inline void
PicSolver<D, real, ArrayAllocator>::onBeginAdvanceTimeStep(double timeInterval)
{
  updateParticleEmitter(timeInterval);

  transferFromParticlesToGrids();

  buildSignedDistanceField();

  extrapolateVelocityToAir();

  this->applyBoundaryCondition();
}

template<size_t D, typename real, typename ArrayAllocator>
inline void
PicSolver<D, real, ArrayAllocator>::computeAdvection(double timeInterval)
{
  Stopwatch s;
  s.start();
  extrapolateVelocityToAir();
  debug("[INFO] ExtrapolateVelocityToAir took %lld ms\n", s.lap());

  this->applyBoundaryCondition();

  s.lap();
  transferFromGridsToParticles();
  debug("[INFO] TransferFromGridsToParticles took %lld ms\n", s.lap());

  moveParticles(timeInterval);
  debug("[INFO] MoveParticles took %lld ms\n", s.lap());
}

template<size_t D, typename real, typename ArrayAllocator>
inline ScalarField<D, real>*
PicSolver<D, real, ArrayAllocator>::fluidSdf() const
{
  return signedDistanceField().get();
}

template<size_t D, typename real, typename ArrayAllocator>
inline void
PicSolver<D, real, ArrayAllocator>::transferFromParticlesToGrids()
{
  if constexpr (D == 2)
  {
    auto vel = this->velocity();
    const auto& spacing = this->gridSpacing();
    const auto& origin = this->gridOrigin();

    LinearArraySampler2<real> uSampler(
      (this->velocity()->data<0>()),
      spacing,
      vel->iOrigin<0>()
    );
    LinearArraySampler2<real> vSampler(
      (this->velocity()->data<1>()),
      spacing,
      vel->iOrigin<1>()
    );

    auto uSize = this->velocity()->iSize<0>();
    auto vSize = this->velocity()->iSize<1>();

    // fill velocity with zero
    this->velocity()->fill(vec::null());
    GridData<2, real> uWeight;
    GridData<2, real> vWeight;
    uWeight.resize(uSize);
    vWeight.resize(vSize);

    _markers[0].resize(uSize);
    _markers[1].resize(vSize);

    fill(uWeight, (real)0.0f);
    fill(_markers[0], (char)0);
    fill(vWeight, (real)0.0f);
    fill(_markers[1], (char)0);

    auto numberOfParticles = this->_particleSystem.size();
    for (size_t i = 0; i < numberOfParticles; ++i)
    {
      std::array<Index2, 4> indices;
      std::array<real, 4> weights;

      uSampler.getCoordinatesAndWeights(this->_particleSystem[i], indices, weights);
      for (int j = 0; j < 4; ++j)
      {
        vel->velocityAt<0>(indices[j])
          += _particleSystem.get<1>(i).x * weights[j];
        uWeight[uWeight.id(indices[j])] += weights[j];
        _markers[0][_markers[0].id(indices[j])] = 1;
      }

      vSampler.getCoordinatesAndWeights(this->_particleSystem[i], indices, weights);
      for (int j = 0; j < 4; ++j)
      {
        vel->velocityAt<1>(indices[j])
          += _particleSystem.get<1>(i).y * weights[j];
        vWeight[vWeight.id(indices[j])] += weights[j];
        _markers[1][_markers[1].id(indices[j])] = 1;
      }
    }

    for (int64_t i = 0; i < uSize.prod(); ++i)
      if (uWeight[i] > 0.0f)
        vel->velocityAt<0>(i) /= uWeight[i];

    for (int64_t i = 0; i < vSize.prod(); ++i)
      if (vWeight[i] > 0.0f)
        vel->velocityAt<1>(i) /= vWeight[i];
  }
}

template<size_t D, typename real, typename ArrayAllocator>
inline void
PicSolver<D, real, ArrayAllocator>::transferFromGridsToParticles()
{
  auto numberOfParticles = _particleSystem.size();
  auto& particles = _particleSystem;
  auto vel = this->velocity();
  for (size_t i = 0; i < numberOfParticles; ++i)
  {
    particles.get<1>(i) = vel->sample(particles[i]);
  }
}

template<size_t D, typename real, typename ArrayAllocator>
inline void
PicSolver<D, real, ArrayAllocator>::moveParticles(double timeInterval)
{
  auto numberOfParticles = _particleSystem.size();
  auto velocity = this->velocity();

  for (size_t i = 0; i < numberOfParticles; ++i)
  {
    vec pt0{ _particleSystem[i] };
    vec pt1{ pt0 };
    vec vel{ _particleSystem.get<1>(i) };
    const auto& bounds = velocity->bounds();
    vec boundsMin{ bounds.min() };
    vec boundsMax{ bounds.max() };

    // Adaptive time-stepping
    unsigned int numSubSteps
      = static_cast<unsigned int>(math::max<real>(this->maxCfl(), 1.0f));
    real dt = static_cast<real>(timeInterval / numSubSteps);
    for (unsigned t = 0; t < numSubSteps; ++t)
    {
      vec vel0{ velocity->sample(pt0) };

      // mid-point rule
      vec midPt{ pt0 + 0.5 * dt * vel0 };
      vec midVel{ velocity->sample(midPt) };
      pt1 = pt0 + dt * midVel;

      pt0 = pt1;
    }

    // review this!!!!!!
    auto flag = this->closedDomainBoundaryFlag();
    if ((flag & constants::directionLeft) && (pt1.x <= boundsMin.x))
    {
      pt1.x = boundsMin.x;
      vel.x = 0.0f;
    }
    if ((flag & constants::directionRight) && (pt1.x >= boundsMax.x))
    {
      pt1.x = boundsMax.x;
      vel.x = 0.0f;
    }
    if ((flag & constants::directionDown) && (pt1.y <= boundsMin.y))
    {
      pt1.y = boundsMin.y;
      vel.y = 0.0f;
    }
    if ((flag & constants::directionUp) && (pt1.y >= boundsMax.y))
    {
      pt1.y = boundsMax.y;
      vel.y = 0.0f;
    }

    if constexpr (D == 3)
    {
      if ((flag & constants::directionBack) && (pt1.z <= boundsMin.z))
      {
        pt1.z = boundsMin.z;
        vel.z = 0.0f;
      }
      if ((flag & constants::directionFront) && (pt1.z >= boundsMax.z))
      {
        pt1.y = boundsMax.y;
        vel.y = 0.0f;
      }
    }

    _particleSystem.set(i, pt1, vel);
  }

  auto col = this->collider();
  if (col != nullptr)
  {
    for (size_t i = 0; i < numberOfParticles; ++i)
    {
      col->resolveCollision(
        0.0f,
        0.0f,
        _particleSystem[i],
        _particleSystem.get<1>(i)
      );
    }
  }
}

template<size_t D, typename real, typename ArrayAllocator>
inline void
PicSolver<D, real, ArrayAllocator>::extrapolateVelocityToAir()
{
  auto vel = this->velocity();

  auto depth = static_cast<unsigned int>(std::ceil(this->maxCfl()));
  extrapolateToRegion(*vel->data<0>(), _markers[0], depth, *vel->data<0>());
  extrapolateToRegion(*vel->data<1>(), _markers[1], depth, *vel->data<1>());
  if constexpr (D == 3)
    extrapolateToRegion(*vel->data<2>(), _markers[2], depth, *vel->data<2>());
}

template<size_t D, typename real, typename ArrayAllocator>
inline void
PicSolver<D, real, ArrayAllocator>::buildSignedDistanceField()
{
  auto maxH = _signedDistanceField->cellSize().max();
  auto radius = 1.2f * maxH / sqrt(2.0f);
  auto sdfBandRadius = 2.0f * radius;

  auto sdfSize = _signedDistanceField->dataSize();
  _searcher->build(_particleSystem);
  forEachIndex<D>(sdfSize, [&](const Index<D>& index) {
    auto pt = _signedDistanceField->dataPosition(index);
    auto minDist = sdfBandRadius;
    _searcher->forEachNearbyPoint(pt, sdfBandRadius, [&](size_t, const vec& x) {
      auto d = (pt - x).length();
      minDist = math::min(minDist, d);
      });
    (*_signedDistanceField)[index] = minDist - radius;
    });

  this->extrapolateIntoCollider(*_signedDistanceField);
}

template<size_t D, typename real, typename ArrayAllocator>
inline void
PicSolver<D, real, ArrayAllocator>::updateParticleEmitter(double timeInterval)
{
  if (_particleEmitter != nullptr)
    _particleEmitter->update(this->currentTime(), timeInterval);
}

} // end namespace cg

#endif // __PicSolver_h
