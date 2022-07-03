#ifndef __GridFluidSolver_h
#define __GridFluidSolver_h

#include "utils/Stopwatch.h"
#include "math/Vector3.h"
#include "PhysicsAnimation.h"
#include "FaceCenteredGrid.h"
#include "GridBackwardEulerDiffusionSolver.h"
#include "GridFractionalSinglePhasePressureSolver.h"
#include "GridFractionalBoundaryConditionSolver.h"
#include "Collider.h"
#include "Constants.h"


namespace cg
{

  template <size_t D, typename real>
  class GridFluidSolver : public PhysicsAnimation
  {
  public:
    using vec = Vector<real, D>;
    template <typename T> using Ref = Reference<T>;

    GridFluidSolver(const Index<D>& size, const vec& spacing, const vec& origin)
    {
      _velocity = new FaceCenteredGrid<D, real>(size, spacing, origin);

      // use Adaptive SubTimeStepping
      this->setIsUsingFixedSubTimeSteps(false);
    }

    virtual ~GridFluidSolver()
    {
      // do nothing
    }

    // Returns gravity for this solver.
    const auto& gravity() const { return _gravity; };

    // Sets the gravity for this solver.
    void setGravity(const vec& gravity) { _gravity = gravity; };

    // Returns the fluid viscosity.
    const auto& viscosityCoefficient() const { return _viscosityCoefficient; }

    // Sets the viscosity coefficient. Non-positive input will clamped to zero.
    void setViscosityCoefficient(real viscosity) { _viscosityCoefficient = math::max<real>(viscosity, 0.0f); }

    // Returns the CFL number from current velocity field for given time interval.
    real cfl(double timeInterval) const;

    // Returns the max allowed CFL number.
    real maxCfl() const { return _maxCfl; }

    // Sets the max allowed CFL number.
    void setMaxCfl(real newCfl) { _maxCfl = math::max(newCfl, math::Limits<real>::eps()); }

    // Returns the closed domain boundary flag.
    int closedDomainBoundaryFlag() const { return _closedDomainBoundaryFlag; }

    // Sets the closed domain boundary flag.
    void setClosedDomainBoundaryFlag(int flag);

    // Returns grid size.
    const auto& size() const { return _velocity->size(); }

    // Returns grid cell spacing.
    const auto& gridSpacing() const { return _velocity->gridSpacing(); }

    // Returns grid origin.
    const auto& gridOrigin() const { return _velocity->origin(); }

    // Returns the velocity field, represented by a FaceCenteredGrid
    const auto& velocity() const { return _velocity; }

    const auto& collider() const { return _collider; }

    void setCollider(Collider<D, real>* collider);

    /*const auto& emitter() const;
     TODO
    void setEmitte(const GridEmitter* emitter);*/

  protected:
    // PhysicsAnimation virtual functions
    void initialize() override;

    void onAdvanceTimeStep(double timeInterval) override;

    size_t numberOfSubTimeSteps(double timeInterval) const override;

    // Called at the beginning of a time-step.
    virtual void onBeginAdvanceTimeStep(double timeInterval);

    // Called at the end of a time-step.
    virtual void onEndAdvanceTimeStep(double timeInterval);

    virtual void computeExternalForces(double timeInterval);

    virtual void computeViscosity(double timeInterval);

    virtual void computePressure(double timeInterval);

    virtual void computeAdvection(double timeInterval);

    void computeGravity(double timeInterval);

    virtual ScalarField<D, real>* fluidSdf() const;

    void applyBoundaryCondition();

    void extrapolateIntoCollider(CellCenteredScalarGrid<D, real>& grid);

    ScalarField<D, real>* colliderSdf() const;

    VectorField<D, real>* colliderVelocityField() const;

  private:
    vec _gravity{ real(0.0f), real(-9.8f) };
    real _viscosityCoefficient{ 0.0f };
    real _maxCfl{ 5.0f };
    int _closedDomainBoundaryFlag{ constants::directionAll };

    Ref<FaceCenteredGrid<D, real>> _velocity;
    Ref<Collider<D, real>> _collider;
    // grid Emitter TODO

    // Solvers
    GridBackwardEulerDiffusionSolver<D, real, false> _diffusionSolver;
    GridFractionalSinglePhasePressureSolver<D, real> _pressureSolver;
    GridFractionalBoundaryConditionSolver<D, real> _boundaryConditionSolver;
    // advectionSolver TODO

    void beginAdvanceTimeStep(double timeInterval);

    void endAdvanceTimeStep(double timeInterval);

    void updateCollider(double timeInterval);

    void updateEmitter(double timeInterval);

  }; // GridFluidSolver<D, real>

  template<size_t D, typename real>
  inline real
    GridFluidSolver<D, real>::cfl(double timeInterval) const
  {
    real maxVel = 0.0f;
    forEachIndex<D>(_velocity->size(), [&](const Index<D>& index) {
      auto v = _velocity->valueAtCellCenter(index) + timeInterval * _gravity;
      maxVel = math::max(maxVel, v.max());
      });

    return real(maxVel * timeInterval / _velocity->gridSpacing().min());
  }

  template<size_t D, typename real>
  inline void
    GridFluidSolver<D, real>::setClosedDomainBoundaryFlag(int flag)
  {
    _closedDomainBoundaryFlag = flag;
    _boundaryConditionSolver.setClosedDomainBoundaryFlag(flag);
  }

  template<size_t D, typename real>
  inline void
    GridFluidSolver<D, real>::setCollider(Collider<D, real>* collider)
  {
    _collider = collider;
  }

  template<size_t D, typename real>
  inline void
    GridFluidSolver<D, real>::initialize()
  {
    updateCollider(0.0);

    updateEmitter(0.0);
  }

  template<size_t D, typename real>
  inline void
    GridFluidSolver<D, real>::onAdvanceTimeStep(double timeInterval)
  {
#ifdef _DEBUG
    // asserting min grid size
    assert(_velocity->size().min() > 0);
#endif // _DEBUG

    beginAdvanceTimeStep(timeInterval);

    computeExternalForces(timeInterval);

    computeViscosity(timeInterval);

    computePressure(timeInterval);

    computeAdvection(timeInterval);

    endAdvanceTimeStep(timeInterval);
  }

  template<size_t D, typename real>
  inline size_t
    GridFluidSolver<D, real>::numberOfSubTimeSteps(double timeInterval) const
  {
    auto _cfl = cfl(timeInterval);
    return static_cast<size_t>(math::max<real>(std::ceil(_cfl / _maxCfl), 1.0f));
  }

  template<size_t D, typename real>
  inline void
    GridFluidSolver<D, real>::onBeginAdvanceTimeStep(double timeInterval)
  {
    // do nothing
  }

  template<size_t D, typename real>
  inline void
    GridFluidSolver<D, real>::onEndAdvanceTimeStep(double timeInterval)
  {
    // do nothing
  }

  template<size_t D, typename real>
  inline void
    GridFluidSolver<D, real>::computeExternalForces(double timeInterval)
  {
    computeGravity(timeInterval);
  }

  template<size_t D, typename real>
  inline void
    GridFluidSolver<D, real>::computeViscosity(double timeInterval)
  {
    Stopwatch s;
    if (math::isPositive(_viscosityCoefficient))
    {
      s.start();
      _diffusionSolver.solve(
        _velocity,
        _viscosityCoefficient,
        timeInterval,
        _velocity,
        *colliderSdf(),
        *fluidSdf()
      );
      //debug("[INFO] Viscosity solver took %lld ms\n", s.time());

      applyBoundaryCondition();
    }
  }

  template<size_t D, typename real>
  inline void
    GridFluidSolver<D, real>::computePressure(double timeInterval)
  {
    Stopwatch s;
    s.start();
    _pressureSolver.solve(
      _velocity,
      timeInterval,
      _velocity,
      *colliderSdf(),
      *fluidSdf(),
      *colliderVelocityField()
    );
    //debug("[INFO] Pressure solver took %lld ms\n", s.time());

    applyBoundaryCondition();
  }

  template<size_t D, typename real>
  inline void
    GridFluidSolver<D, real>::computeAdvection(double timeInterval)
  {
    // TODO
  }

  template<size_t D, typename real>
  inline void
    GridFluidSolver<D, real>::computeGravity(double timeInterval)
  {
    if (this->_gravity.squaredNorm() > math::Limits<real>::eps())
    {
      std::array<int, D> sizes;
      sizes[0] = _velocity->iSize<0>().prod();
      sizes[1] = _velocity->iSize<1>().prod();
      if constexpr (D == 3)
        sizes[2] = _velocity->iSize<2>().prod();

      if (!math::isZero(_gravity.x))
      {
        for (int i = 0; i < sizes[0]; ++i)
          _velocity->velocityAt<0>(i) += timeInterval * _gravity.x;
      }

      if (!math::isZero(_gravity.y))
      {
        for (int i = 0; i < sizes[1]; ++i)
          _velocity->velocityAt<1>(i) += timeInterval * _gravity.y;
      }

      if constexpr (D == 3)
      {
        if (!math::isZero(_gravity.z))
        {
          for (int i = 0; i < sizes[2]; ++i)
            _velocity->velocityAt<2>(i) += timeInterval * _gravity.z;
        }
      }

      applyBoundaryCondition();
    }
  }

  template<size_t D, typename real>
  inline ScalarField<D, real>*
    GridFluidSolver<D, real>::fluidSdf() const
  {
    return new ConstantScalarField<D, real>(math::Limits<real>::inf());
  }

  template<size_t D, typename real>
  inline void
    GridFluidSolver<D, real>::applyBoundaryCondition()
  {
    auto depth = static_cast<unsigned int>(std::ceil(_maxCfl));
    _boundaryConditionSolver.constrainVelocity(_velocity, depth);
  }

  template<size_t D, typename real>
  inline void
    GridFluidSolver<D, real>::extrapolateIntoCollider(CellCenteredScalarGrid<D, real>& grid)
  {
    GridData<D, char> marker;
    marker.resize(grid.dataSize());

    for (size_t i = 0; i < grid.length(); ++i)
    {
      auto index = grid.index(i);
      auto colliderSdfQuery = colliderSdf()->sample(grid.dataPosition(index));
      if (isInsideSdf(colliderSdfQuery))
        marker[i] = 0;
      else
        marker[i] = 1;
    }

    auto depth = static_cast<unsigned int>(std::ceil(_maxCfl));
    extrapolateToRegion(grid, marker, depth, grid);
  }

  template<size_t D, typename real>
  inline ScalarField<D, real>*
    GridFluidSolver<D, real>::colliderSdf() const
  {
    return _boundaryConditionSolver.colliderSdf();
  }

  template<size_t D, typename real>
  inline VectorField<D, real>*
    GridFluidSolver<D, real>::colliderVelocityField() const
  {
    return _boundaryConditionSolver.colliderVelocityField();
  }

  template<size_t D, typename real>
  inline void
    GridFluidSolver<D, real>::beginAdvanceTimeStep(double timeInterval)
  {
    updateCollider(timeInterval);

    updateEmitter(timeInterval);

    _boundaryConditionSolver.updateCollider(
      _collider,
      _velocity->size(),
      _velocity->gridSpacing(),
      _velocity->origin()
    );

    applyBoundaryCondition();

    // Invoke callback
    onBeginAdvanceTimeStep(timeInterval);
  }

  template<size_t D, typename real>
  inline void
    GridFluidSolver<D, real>::endAdvanceTimeStep(double timeInterval)
  {
    // Invoke callback
    onEndAdvanceTimeStep(timeInterval);
  }

  template<size_t D, typename real>
  inline void
    GridFluidSolver<D, real>::updateCollider(double timeInterval)
  {
    if (_collider != nullptr)
      _collider->update(this->currentTime(), timeInterval);
  }

  template<size_t D, typename real>
  inline void
    GridFluidSolver<D, real>::updateEmitter(double timeInterval)
  {
    // TODO
  }

} // end namespace cg

#endif // __GridFluidSolver_h
