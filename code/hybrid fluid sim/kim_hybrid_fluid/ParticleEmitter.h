#ifndef __ParticleEmitter_h
#define __ParticleEmitter_h

#include <functional>
#include "geometry/ParticleSystem.h"

namespace cg
{

/**
* Abstract base class for particle emitter.
*/
template <typename PointArray>
class ParticleEmitter: public SharedObject
{
public:
  /**
  * Callback function type for update calls.
  * 
  * This type of callback function will take the emitter reference, current
  * time and time interval in seconds.
  */
  using OnBeginUpdateCallback = std::function<void(ParticleEmitter<PointArray>&, double, double)>;

  /**
  * \brief  Constructs a particle emitter that will generate particles into
  *         \p particles.
  */
  ParticleEmitter(PointArray& particles)
    : _particles(std::ref(particles))
  {
    // do nothing
  }

  /** Default destructor. */
  virtual ~ParticleEmitter()
  {
    // do nothing
  }

  /** Updates the emitter state from \p currentTimeInSeconds to the following
  * time-step.
  */
  void update(double currentTimeInSeconds, double timeIntervalInSeconds);

  /** Returns the target particle system to emit. */
  PointArray& target() const { return _particles.get(); }

  /** Sets the target particle system to emit. */
  void setTarget(PointArray& particles);

  /** Returns \c true if the emitter is enabled, \c false otherwise. */
  bool isEnabled() const { return _isEnabled; }

  /** Sets true/false to enable/disable the emitter. */
  void setIsEnabled(bool enabled) { _isEnabled = enabled; }

  /**
  * \brief    Sets the callback function to be called when
  *           ParticleEmitter::update method is invoked.
  * 
  * The callback function takes current simulation time in seconds.
  * Use this callback to track any motion or state changes related to this
  * emitter.
  * 
  * \param[in] callback The callback function.
  */
  void setOnBeginUpdateCallback(const OnBeginUpdateCallback& callback) { _onBeginUpdateCallback = callback; }

protected:
  /** Called when ParticleEmitter::setTarget is executed. */
  virtual void onSetTarget(PointArray& particles);

  /** Called when ParticleEmitter::update is executed. */
  virtual void onUpdate(double currentTimeInSeconds, double timeIntervalInSeconds) = 0;

private:
  bool _isEnabled = true;
  /** Particle system reference. */
  std::reference_wrapper<PointArray> _particles;
  /** Update callback. */
  OnBeginUpdateCallback _onBeginUpdateCallback;
};

template<typename PointArray>
inline void
ParticleEmitter<PointArray>::update(double currentTimeInSeconds, double timeIntervalInSeconds)
{
  if (_onBeginUpdateCallback)
    _onBeginUpdateCallback(*this, currentTimeInSeconds, timeIntervalInSeconds);
  onUpdate(currentTimeInSeconds, timeIntervalInSeconds);
}

template<typename PointArray>
inline void
ParticleEmitter<PointArray>::setTarget(PointArray& particles)
{
  _particles = particles;

  onSetTarget(particles);
}

template<typename PointArray>
inline void
ParticleEmitter<PointArray>::onSetTarget(PointArray& particles)
{
  // do nothing
}

} // end namespace cg

#endif // __ParticleEmitter_h
