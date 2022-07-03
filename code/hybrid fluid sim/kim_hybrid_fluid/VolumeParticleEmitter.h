#ifndef __VolumeParticleEmitter_h
#define __VolumeParticleEmitter_h

#include <vector>
#include <algorithm>
#include "ParticleEmitter.h"
#include "TrianglePointGenerator.h"
#include "math/Surface.h"

namespace cg
{

// Simple Particle Emitter, lacking a lot of features, as it is now emits
// all particles inside the surface
/**
* Template class for a D-dimensional volumetric particle emitter.
* 
* This class emits particles from volumetric geometry.
* 
* \tparam D Defines the number of dimensions.
* \tparam real A floating point type.
* \tparam PointArray An array of points, or particles.
* 
* \todo Implement 3-D template specialization.
*/
template <size_t D, typename real, typename PointArray>
class VolumeParticleEmitter;

/**
* Template class for a 2-D volumetric particle emitter.
* 
* This class specializes the template class VolumeParticleEmitter to a 2-D
* case.
* 
* \tparam real A floating point type.
* \tparam PointArray An array of points, or particles.
* 
* \todo Allow random seed to be provided.
* \todo Use ImplicitSurface instead of Surface
* \todo Implement case where emitter is not one shot.
*/
template <typename real, typename PointArray>
class VolumeParticleEmitter<2, real, PointArray> final : public ParticleEmitter<PointArray>
{
public:
  using Base = ParticleEmitter<PointArray>; ///< Base class alias.
  using vec_type = Vector<real, 2>; ///< Vector type alias.
  using bounds_type = Bounds<real, 2>; ///< Bounds type alias.

  /**
  * Constructs an emitter that generates particles inside a given surface
  * that defines the volumetric geometry. These particles are placed in
  * \p particles. By default uses TrianglePointGenerator to produce particles.
  * 
  * \param[in, out] particles     The object to put particles into.
  * \param[in]  surface           The surface.
  * \param[in]  bounds            The max region.
  * \param[in]  spacing           The spacing between particles.
  * \param[in]  initialVel        The initial velocity of new particles.
  * \param[in]  linearVel         The linear velocity of the emitter.
  * \param[in]  angularVel        The angular velocity of the emitter.
  * \param[in]  maxParticles      The max number of particles to be emitted.
  * \param[in]  isOneShot         True if emitter should emit all particles
  *                               at once and then be disabled.
  * \param[in]  allowOverlapping  True if particles can overlap.
  */
  VolumeParticleEmitter(
    PointArray& particles,
    math::Surface<2, real>* surface,
    const bounds_type& bounds,
    real spacing,
    const vec_type& initialVel = vec_type::null(),
    const vec_type& linearVel = vec_type::null(),
    real angularVel = 0.0f,
    size_t maxParticles = math::Limits<size_t>::inf(),
    bool isOneShot = true,
    bool allowOverlapping = false
  ):
    Base(particles),
    _surface(surface),
    _bounds(bounds),
    _spacing(spacing),
    _initialVel(initialVel),
    _linearVel(linearVel),
    _angularVel(angularVel),
    _maxNumberOfParticles(maxParticles),
    _isOneShot(isOneShot),
    _allowOverlapping(allowOverlapping)
  {
    _pointsGen = new TrianglePointGenerator<real>();
  }

  /** Destructor */
  ~VolumeParticleEmitter()
  {
    delete _pointsGen;
  }

  /** Sets a new point generator. */
  void setPointGenerator(PointGenerator<2, real>* pointsGen) { _pointsGen = pointsGen; }

  /** Returns the surface. */
  const math::Surface<2, real>* surface() const { return _surface.get(); }

  /** Sets the surface. */
  void setSurface(math::Surface<2, real>* surface) { _surface = surface; }

  /** Returns the bounds in which particles can be generated. */
  const auto& maxRegion() const { return _bounds; }

  /** Sets the bounds in which particles can be generated. */
  void setMaxRegion(const bounds_type& bounds) { _bounds = bounds; }

  /** Returns \c true if particles can overlap, \c false otherwise. */
  bool allowOverlapping() const { return _allowOverlapping; }

  /** \p allowOverlapping defines if particles can overlap or not. */
  void setAllowOverlapping(bool allowOverlapping) { _allowOverlapping = allowOverlapping; }

  /** Returns \c true if emitter is one shot, \c false otherwise. */
  bool isOneShot() const { return _isOneShot; }

  /** Sets if emitter is one shot. */
  void setIsOneShot(bool oneShot) { _isOneShot = oneShot; }

  /** Returns the maximum number of particles that can be generated. */
  auto maxNumberOfParticles() const { return _maxNumberOfParticles; }

  /** Sets the maximum number of particles that can be generated. */
  void setMaxNumberOfParticles(size_t max) { _maxNumberOfParticles = math::max<size_t>(0, max); }

  /** Returns the spacing. */
  real spacing() const { return _spacing; }

  /** Sets the spacing. */
  void setSpacing(real newSpacing) { _spacing = math::abs(newSpacing); }

  /** Returns the initial velocity of the particles. */
  auto initialVelocity() const { return _initialVel; }

  /** Sets the initial velocity of the particles. */
  void setInitialVelocity(const vec_type& vel) { _initialVel = vel; }

  /** Returns the linear velocity of the emitter. */
  auto linearVelocity() const { return _linearVel; }

  /** Sets the linear velocity of the emitter. */
  void setLinearVelocity(const vec_type& vel) { _linearVel = vel; }

  /** Returns the angular velocity of the emitter. */
  real angularVelocity() const { return _angularVel; }

  /** Sets the angular velocity of the emitter. */
  void setAngularVelocity(real vel) { _angularVel = vel; }

private:
  /** Surface reference. */
  Reference<math::Surface<2, real>> _surface; // should be replaced by implicit surface
  /** Emitter bounds. */
  bounds_type _bounds;
  /** Particle spacing. */
  real _spacing;
  /** Particle initial velocity. */
  vec_type _initialVel;
  /** Emitter linear velocity. */
  vec_type _linearVel;
  /** Emitter angular velocity. */
  real _angularVel = 0.0f;
  /** Point generator. */
  PointGenerator<2, real>* _pointsGen;

  /** Max number of particles. */
  size_t _maxNumberOfParticles = math::Limits<size_t>::inf();
  /** Number of emitter particles. */
  size_t _numberOfEmittedParticles = 0;

  bool _isOneShot = true;
  bool _allowOverlapping = false;

  /**
  * Emits particles to the particle system.
  * 
  * \param[in] currentTimeInSeconds     Current simulation time.
  * \param[in] timeIntervalInSeconds    The time-step interval.
  */
  void onUpdate(double currentTimeInSeconds, double timeIntervalInSeconds) override;

  /**
  * \brief    Emits particles filling their information into
  *           \p positions and \p velocities.
  * 
  * Uses the point generator instance to produce particles within bounds.
  * The velocities calculated by VolumeParticleEmitter::velocityAt.
  */
  void emit(std::vector<vec_type>& positions, std::vector<vec_type>& velocities);

  /**
  * Calculates the velocity of a particle based on the emitter's properties.
  * 
  * \param[in] point Particle position.
  * \return Particle velocity.
  */
  vec_type velocityAt(const vec_type& point) const;

}; // VolumeParticleEmitter

template<typename real, typename PointArray>
inline void
VolumeParticleEmitter<2, real, PointArray>::onUpdate(double currentTimeInSeconds, double timeIntervalInSeconds)
{
  auto& particles = this->target();

  if (!this->isEnabled()) return;

  std::vector<vec_type> newPositions, newVelocities;

  emit(newPositions, newVelocities);

  for (size_t i = 0; i < newPositions.size(); ++i)
  {
    if (!particles.add(newPositions[i], newVelocities[i]))
    {
#ifdef _DEBUG
      printf("Early stopping in VolumeParticleEmitter, particle system capacity reached.\n");
#endif
      break;
    }
  }

  if (_isOneShot)
    this->setIsEnabled(false);
}

template<typename real, typename PointArray>
inline void
VolumeParticleEmitter<2, real, PointArray>::emit(std::vector<vec_type>& positions, std::vector<vec_type>& velocities)
{
  if (_surface == nullptr)
    return;

  vec_type min = _bounds.min(), max = _bounds.max();
  if (_surface->isBounded())
  {
    auto bb = _surface->bounds();
    const auto& _min = bb.min();
    const auto& _max = bb.max();
    max = vec_type{
      math::max(max.x, _max.x),
      math::max(max.y, _max.y)
    };
    min = vec_type{
      math::min(min.x, _min.x),
      math::min(min.y, _min.y)
    };
  }
  bounds_type region{ min, max };

  size_t numNewParticles = 0;

  if (_allowOverlapping || _isOneShot)
  {
    _pointsGen->forEachPoint(region, _spacing, [&](const vec_type& point) {
      if (_surface->isInside(point))
      {
        if (_numberOfEmittedParticles < _maxNumberOfParticles)
        {
          positions.push_back(point);
          ++_numberOfEmittedParticles;
          ++numNewParticles;
        }
        else
        {
          return false;
        }
      }
      return true;
    });
  }
  else
  {
    // TODO
  }

  velocities.resize(positions.size());
  for (size_t i = 0; i < positions.size(); ++i)
    velocities[i] = velocityAt(positions[i]);
}

template<typename real, typename PointArray>
inline Vector<real, 2>
VolumeParticleEmitter<2, real, PointArray>::velocityAt(const vec_type& point) const
{
  auto r = point - _surface->transform.transform(vec_type::null());
  return _linearVel + _angularVel * vec_type{ -r.y, r.x } + _initialVel;
}

} // end namespace cg

#endif // __VolumeParticleEmitter_h
