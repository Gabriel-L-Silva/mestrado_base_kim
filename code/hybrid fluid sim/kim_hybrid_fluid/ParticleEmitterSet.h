#ifndef __ParticleEmitterSet_h
#define __ParticleEmitterSet_h

#include <vector>
#include "core/SharedObject.h"
#include "ParticleEmitter.h"

namespace cg
{

template <typename PointArray>
class ParticleEmitterSet final : public ParticleEmitter<PointArray>
{
public:
  using Base = ParticleEmitter<PointArray>;
  using EmitterVector = std::vector<Reference<ParticleEmitter<PointArray>>>;

  ParticleEmitterSet(PointArray& particles)
    : Base(particles)
  {
    // do nothing
  }

  explicit ParticleEmitterSet(const EmitterVector& emitters, PointArray& particles)
    : _emitters(emitters), Base(particles)
  {
    // do nothing
  }

  virtual ~ParticleEmitterSet()
  {
    // do nothing
  }

  void addEmitter(ParticleEmitter<PointArray>* emitter)
  {
    _emitters.push_back(emitter);
  }

private:
  EmitterVector _emitters;

  void onSetTarget(PointArray& particles) override;

  void onUpdate(double currentTimeInSeconds, double timeIntervalInSeconds) override;

}; // ParticleEmitterSet<PointArray>

template<typename PointArray>
inline void
ParticleEmitterSet<PointArray>::onSetTarget(PointArray& particles)
{
  for (auto& emitter : _emitters)
  {
    emitter->setTarget(particles);
  }
}

template<typename PointArray>
inline void
ParticleEmitterSet<PointArray>::onUpdate(double currentTimeInSeconds, double timeIntervalInSeconds)
{
  if (!this->isEnabled())
    return;

  for (auto& emitter : _emitters)
  {
    emitter->update(currentTimeInSeconds, timeIntervalInSeconds);
  }
}

} // end namespace cg

#endif // __ParticleEmitterSet_h
