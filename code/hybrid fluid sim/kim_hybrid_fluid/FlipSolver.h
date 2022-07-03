#ifndef __FlipSolver_h
#define __FlipSolver_h

#include "PicSolver.h"

namespace cg
{

/**
* Fluid-Implicit-Particle (FLIP) template implementation.
* 
* This class implements 2-D and 3-D Fluid-Implicit-Particle (FLIP) method by inheriting PicSolver.
* As it is a hybrid method, the solver uses particles to keep track of fluid particles.
* This implementation is based on Kim Doyub's Jet Fluid Simulation Library.
* 
* \tparam D Defines the number of dimensions
* \tparam real A floating point type
* \tparam ArrayAllocator Class responsible for allocating data
*/
template <size_t D, typename real, typename ArrayAllocator>
class FlipSolver : public PicSolver<D, real, ArrayAllocator>
{
public:
  using Base = PicSolver<D, real, ArrayAllocator>;
  using Base::Base;
  using vec = Vector<real, D>;

  /** Default destructor. */
  virtual ~FlipSolver()
  {
    // do nothing
  }

  /**
  * Returns PIC blending factor.
  */
  real picBlendingFactor() const
  {
    return _picBlendingFactor;
  }

  /**
  * Sets the PIC blending factor.
  * 
  * This function sets the PIC blending factor which mixes FLIP and PIC
  * results when transferring velocity from grids to particles in order to
  * reduce the noise. The factor can be a value between 0 and 1, where 0
  * means no blending and 1 means full PIC. Default is 0.
  * 
  * \param[in] factor Factor to be used.
  */
  void setPicBlendingFactor(real factor)
  {
    _picBlendingFactor = math::clamp<real>(factor, 0.0f, 1.0f);
  }

protected:
  /**
  * Transfers velocity field from particles to grids.
  * 
  * This method transfers the velocity using PicSolver's
  * transferFromParticlesToGrids and then stores the initial velocity in _delta.
  */
  void transferFromParticlesToGrids() override;

  /**
  * Transfers velocity field from grids to particles.
  * 
  * This method computes and then interpolates the velocity delta to update the
  * particle velocity.
  */
  void transferFromGridsToParticles() override;

private:
  /**
  * PIC blending factor, mixes results from PIC and FLIP.
  * 
  * This factor is a value between 0 and 1, where 0 means no blending and 1 means full PIC.
  * Default is 0.
  * 
  * \see setPicBlendingFactor
  */
  real _picBlendingFactor = 0.0f;

  /**
  * Array of GridData, stores the deltas used at transfer step.
  * 
  * This array has D instances of GridData, where each one stores the delta velocity for a dimension of the problem.
  */
  std::array<GridData<D, real>, D> _delta;

}; // FlipSolver<D, real, ArrayAllocator>

template<size_t D, typename real, typename ArrayAllocator>
inline void
FlipSolver<D, real, ArrayAllocator>::transferFromParticlesToGrids()
{
  Base::transferFromParticlesToGrids();
  auto vel = this->velocity();

  std::array<size_t, D> sizes;
  sizes[0] = vel->iSize<0>().prod();
  sizes[1] = vel->iSize<1>().prod();
  if constexpr (D == 3)
    sizes[2] = vel->iSize<2>().prod();
  _delta[0].resize(vel->iSize<0>());
  _delta[1].resize(vel->iSize<1>());
  if constexpr (D == 3)
    _delta[2].resize(vel->iSize<2>());

  for (size_t i = 0; i < sizes[0]; ++i)
    _delta[0][i] = vel->velocityAt<0>((int64_t)i);
  for (size_t i = 0; i < sizes[1]; ++i)
    _delta[1][i] = vel->velocityAt<1>((int64_t)i);

  if constexpr (D == 3)
  {
    for (size_t i = 0; i < sizes[2]; ++i)
      _delta[2][i] = vel->velocityAt<2>((int64_t)i);
  }

}

template<size_t D, typename real, typename ArrayAllocator>
inline void
FlipSolver<D, real, ArrayAllocator>::transferFromGridsToParticles()
{
  auto vel = this->velocity();
  auto numberOfParticles = this->particleSystem().size();

  std::array<size_t, D> sizes;
  sizes[0] = vel->iSize<0>().prod();
  sizes[1] = vel->iSize<1>().prod();
  if constexpr (D == 3)
    sizes[2] = vel->iSize<2>().prod();

  // Compute delta
  for (size_t i = 0; i < sizes[0]; ++i)
    _delta[0][i] = vel->velocityAt<0>((int64_t)i) - _delta[0][i];
  for (size_t i = 0; i < sizes[1]; ++i)
    _delta[1][i] = vel->velocityAt<1>((int64_t)i) - _delta[1][i];

  if constexpr (D == 3)
  {
    for (size_t i = 0; i < sizes[2]; ++i)
      _delta[2][i] = vel->velocityAt<2>((int64_t)i) - _delta[2][i];
  }
  // GAMBIARRA, REMOVER O MAIS CEDO POSSIVEL
  std::array<LinearArraySampler<real, D>, D> samplers;
  samplers[0] = LinearArraySampler<real, D>(nullptr, vel->gridSpacing(), vel->iOrigin<0>());
  samplers[0].setGridData(&_delta[0]);
  samplers[1] = LinearArraySampler<real, D>(nullptr, vel->gridSpacing(), vel->iOrigin<1>());
  samplers[1].setGridData(&_delta[1]);
  if constexpr (D == 3)
  {
    samplers[2] = LinearArraySampler<real, D>(nullptr, vel->gridSpacing(), vel->iOrigin<2>());
    samplers[2].setGridData(&_delta[2]);
  }
  
  auto sampler = [&](const vec& x) {
    vec r{ real(0.0f) };
    for (int i = 0; i < D; ++i)
      r[i] = samplers[i](x);
    return r;
  };

  auto& ps = this->_particleSystem;
  for (size_t i = 0; i < numberOfParticles; ++i)
  {
    auto flipVel = ps.get<1>(i) + sampler(ps[i]);
    if (math::isPositive(_picBlendingFactor))
    {
      auto picVel = vel->sample(ps[i]);
      flipVel = lerp(flipVel, picVel, _picBlendingFactor);
    }
    ps.get<1>(i) = flipVel;
  }
}

} // end namespace cg

#endif // __FlipSolver_h

#pragma once
