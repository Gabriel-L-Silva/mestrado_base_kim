#ifndef __GridSolver_h
#define __GridSolver_h

#include "utils/Stopwatch.h"
#include "math/Vector3.h"
#include "PhysicsAnimation.h"
#include "FaceCenteredGrid.h"
#include "CellCenteredScalarGrid.h"
#include "GridBackwardEulerDiffusionSolver.h"
#include "GridFractionalSinglePhasePressureSolver.h"
#include "GridFractionalBoundaryConditionSolver.h"
#include "Collider.h"
#include "Constants.h"


namespace cg
{

  template <size_t D, typename real>
  class GridSolver : public PhysicsAnimation
  {
  public:
    using vec = Vector<real, D>;
    template <typename T> using Ref = Reference<T>;

    GridSolver(const Index<D>& size, const vec& spacing, const vec& origin)
    {
      _velocity = new FaceCenteredGrid<D, real>(size+2, spacing, origin-spacing);
      _density = new CellCenteredScalarGrid<D, real>(size+2, spacing, origin-spacing);
      _temperature = new CellCenteredScalarGrid<D, real>(size + 2, spacing, origin - spacing);
      _solverSize = Index2{ size.x,size.y };
      // use Adaptive SubTimeStepping
      this->setIsUsingFixedSubTimeSteps(false);
    }

    virtual ~GridSolver()
    {
      // do nothing
    }

    enum class AdvectType
    {
      Velocity = 1,
      Density = 2,
    };

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
    const auto& size() const { return _solverSize; }

    const auto& minTemp() const { return _minTemp; }

    const auto& maxTemp() const { return _maxTemp; }

    // Returns grid cell spacing.
    const auto& gridSpacing() const { return _velocity->gridSpacing(); }

    // Returns grid origin.
    const auto& gridOrigin() const { return _velocity->origin(); }

    // Returns the velocity field, represented by a FaceCenteredGrid
    const auto& velocity() const { return _velocity; }

    // Returns the density field, represented by a CellCenteredScalarGrid
    const auto& density() const { return _density; }

    // Returns the temperature field, represented by a CellCenteredScalarGrid
    const auto& temperature() const { return _temperature; }

    const auto& collider() const { return _collider; }

    void setCollider(Collider<D, real>* collider);
    
    void updateTempMinMax();

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

    void divergent(CellCenteredScalarGrid<D, real>& div, Ref<FaceCenteredGrid<D, real>> input);

    virtual void computePressure(double timeInterval);

    void advectDensity(Index2 idx, double timeInterval);

    void computeAdvection(double timeInterval, AdvectType type);

    void advectVelocity(Index2 idx, double timeInterval);

    void computeBuoyancy(double timeInterval);

    void computeSource(double timeInterval);

    void computeGravity(double timeInterval);
    
    void densityStep(double timeInterval);
    
    void velocityStep(double timeInterval);

    virtual ScalarField<D, real>* fluidSdf() const;

    void applyBoundaryCondition();

    void applyBoundaryCondition(CellCenteredScalarGrid<D, real>& grid);

    void extrapolateIntoCollider(CellCenteredScalarGrid<D, real>& grid);

    void poisson_solver(CellCenteredScalarGrid<D, real>&, CellCenteredScalarGrid<D, real>&, double);
    
    ScalarField<D, real>* colliderSdf() const;

    VectorField<D, real>* colliderVelocityField() const;
  private:
    Index2 _solverSize{ 0 , 0};
    vec _gravity{ real(0.0f), real(-9.8f) };
    real _viscosityCoefficient{ 0.0f };
    real _densityBuoyancyFactor{ 0.005f };
    real _temperatureBuoyancyFactor{ 0.1f };
    real _maxCfl{ 5.0f };
    int _maxIterPoisson{ 10 };
    int _closedDomainBoundaryFlag{ constants::directionAll };
    float _minTemp{ math::Limits<real>::inf() }, _maxTemp{ 0.0f };

    Ref<FaceCenteredGrid<D, real>> _velocity;
    Ref<CellCenteredScalarGrid<D, real>> _density;
    Ref<CellCenteredScalarGrid<D, real>> _temperature;
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

  }; // GridSolver<D, real>

  template<size_t D, typename real>
  inline real
    GridSolver<D, real>::cfl(double timeInterval) const
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
    GridSolver<D, real>::setClosedDomainBoundaryFlag(int flag)
  {
    _closedDomainBoundaryFlag = flag;
    _boundaryConditionSolver.setClosedDomainBoundaryFlag(flag);
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::setCollider(Collider<D, real>* collider)
  {
    _collider = collider;
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::initialize()
  {
    updateCollider(0.0);

    updateEmitter(0.0);
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::densityStep(double timeInterval)
  {
    computeSource(timeInterval);

    computeViscosity(timeInterval);//k=0, passthrough

    computeAdvection(timeInterval, AdvectType::Density);
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::velocityStep(double timeInterval)
  {
#ifdef _DEBUG
    // asserting min grid size
    assert(_velocity->size().min() > 0);
#endif // _DEBUG

    beginAdvanceTimeStep(timeInterval);

    computeExternalForces(timeInterval);

    computeViscosity(timeInterval);//k=0, passthrough

    computePressure(timeInterval);

    computeAdvection(timeInterval, AdvectType::Velocity);

    computePressure(timeInterval);
    
    endAdvanceTimeStep(timeInterval);
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::onAdvanceTimeStep(double timeInterval)
  {
#ifdef _DEBUG
    // asserting min grid size
    assert(_velocity->size().min() > 0);
#endif // _DEBUG

    beginAdvanceTimeStep(timeInterval);

    velocityStep(timeInterval);

    densityStep(timeInterval);

    endAdvanceTimeStep(timeInterval);
  }

  template<size_t D, typename real>
  inline size_t
    GridSolver<D, real>::numberOfSubTimeSteps(double timeInterval) const
  {
    auto _cfl = cfl(timeInterval);
    return static_cast<size_t>(math::max<real>(std::ceil(_cfl / _maxCfl), 1.0f));
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::onBeginAdvanceTimeStep(double timeInterval)
  {
    // do nothing
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::onEndAdvanceTimeStep(double timeInterval)
  {
    // do nothing
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::computeExternalForces(double timeInterval)
  {
    computeBuoyancy(timeInterval);
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::computeViscosity(double timeInterval)
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
  void
  GridSolver<D, real>::divergent(CellCenteredScalarGrid<D, real>& div, Ref<FaceCenteredGrid<D, real>> input)
  {
    auto N = size();
    auto invX = 1.0f / N.x;
    auto invY = 1.0f / N.y;
    for (auto y = 1; y < N.y-1; y++)
    {
      for (auto x = 1; x < N.x-1; x++)
      {
        auto index = Index2(x, y);
        auto right = input->velocityAt<0>(index + Index2(1, 0));
        auto left = input->velocityAt<0>(index);
        auto up = input->velocityAt<1>(index + Index2(0, 1));
        auto down = input->velocityAt<1>(index);
        /*debug("%.2f\n", right);
        debug("%.2f\n", left);
        debug("%.2f\n", up);
        debug("%.2f\n", down);*/
        
        div[index] = invX * (right - left) + invY * (up - down);
      }
    }
    applyBoundaryCondition(div);
  }


  template<size_t D, typename real>
  void 
  GridSolver<D, real>::poisson_solver(CellCenteredScalarGrid<D, real>& div, CellCenteredScalarGrid<D, real>& x0, double tol)
  {
    auto N = size();
    auto invX = 1.0f / N.x;
    auto invY = 1.0f / N.y;
    auto old = new CellCenteredScalarGrid<D, real>(_density->size(), _density->cellSize(), _density->dataOrigin());
    for (int i = 0; i < _maxIterPoisson; i++) {
      forEachIndex<D>(x0.size(), [&](const Index<D>& index) {
        (*old)[index] = x0[index];
        });
      double accum = 0.0;
      for (auto y = 1; y < N.y-1; y++)
      {
        for (auto x = 1; x < N.x-1; x++)
        {
          auto index = Index2(x, y);
          /*debug("%.2f\n", x0[index + Index2(1, 0)]);
          debug("%.2f\n", x0[index + Index2(-1, 0)]);
          debug("%.2f\n", x0[index + Index2(0, 1)]);
          debug("%.2f\n", x0[index + Index2(0, -1)]);
          debug("%.2f\n", div[index]);*/
          x0[index] = 0.25 * (x0[index + Index2(1, 0)] + x0[index + Index2(-1, 0)] + x0[index + Index2(0, 1)] + x0[index + Index2(0, -1)] - div[index]);
          accum += abs((*old)[index] - x0[index]);
        }
      }
      
      if (accum < tol)
        break;
    }
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::computePressure(double timeInterval)
  {
    Stopwatch s;
    s.start();
    /*_pressureSolver.solve(
      _velocity,
      timeInterval,
      _velocity,
      *colliderSdf(),
      *fluidSdf(),
      *colliderVelocityField()
    );*/
    auto div = new CellCenteredScalarGrid<D, real>(_density->size(), _density->cellSize(), _density->dataOrigin());
    auto x0 = new CellCenteredScalarGrid<D, real>(_density->size(), _density->cellSize(), _density->dataOrigin());
    
    divergent(*div, _velocity);
    poisson_solver(*div, *x0, 10e-5);
    auto N = size();
    auto maxVel = 0.f;
    for (auto y = 1; y < N.y - 1; y++)
    {
      for (auto x = 1; x < N.x - 1; x++)
      {
        auto index = Index2(x, y);
        auto grad = x0->gradientAt(index);
        auto vel = vec{ _velocity->velocityAt<0>(index) , _velocity->velocityAt<1>(index) };
        if (vel.squaredNorm() > maxVel)
          maxVel = vel.squaredNorm();

        _velocity->velocityAt<0>(index) -= grad.x;
        _velocity->velocityAt<1>(index) -= grad.y;
      }
    }
    //debug("[INFO] Pressure solver took %lld ms\n", s.time());
    //debug("max vel: %.3f\n", maxVel);

    applyBoundaryCondition();
  }

  template<size_t D, typename real>
  inline void GridSolver<D, real>::advectVelocity(Index2 idx, double timeInterval)
  {
    auto vel = vec(_velocity->velocityAt<0>(idx), _velocity->velocityAt<1>(idx));
    auto bounds = _velocity->bounds();
    auto cellSize = gridSpacing();
    auto minX = bounds.min().x + cellSize.x;
    auto maxX = bounds.max().x - cellSize.x;
    auto minY = bounds.min().y + cellSize.y;
    auto maxY = bounds.max().y - cellSize.y;
    
    auto uPos = _velocity->positionInSpace<0>();

    auto newUPos = uPos(idx) - vel * timeInterval;
    auto midUPt = uPos(idx) - vel * 0.5 * timeInterval;

    //Stop backtracing at cell face
    newUPos.x = newUPos.x < minX ? minX : newUPos.x;
    newUPos.x = newUPos.x > maxX ? maxX : newUPos.x;
    newUPos.y = newUPos.y < minY ? minY : newUPos.y;
    newUPos.y = newUPos.y > maxY ? maxY : newUPos.y;

    _velocity->velocityAt<0>(idx) = _velocity->sample(midUPt).x;

    auto vPos = _velocity->positionInSpace<1>();

    auto newVPos = vPos(idx) - vel * timeInterval;
    auto midVPt = vPos(idx) - vel * 0.5f * timeInterval;

    //Stop backtracing at cell face
    newVPos.x = newVPos.x < minX ? minX : newVPos.x;
    newVPos.x = newVPos.x > maxX ? maxX : newVPos.x;
    newVPos.y = newVPos.y < minY ? minY : newVPos.y;
    newVPos.y = newVPos.y > maxY ? maxY : newVPos.y;

    _velocity->velocityAt<1>(idx) = _velocity->sample(midVPt).y;
  }

  template<size_t D, typename real>
  inline void GridSolver<D, real>::advectDensity(Index2 idx, double timeInterval)
  {
    auto vel = vec(_velocity->velocityAt<0>(idx), _velocity->velocityAt<1>(idx));
    auto pos = _density->dataPosition(idx);
    auto bounds = _density->bounds();
    auto cellSize = _density->cellSize();

    //if (idx == Index2{ 64,64 })
      //debug("aq");
    auto minX = bounds.min().x + cellSize.x;
    auto maxX = bounds.max().x - cellSize.x;
    auto minY = bounds.min().y + cellSize.y;
    auto maxY = bounds.max().y - cellSize.y;

    //Stop backtracing at cell face
    auto midPtVel = _velocity->sample(pos - vel * 0.5f * timeInterval);
    auto newPos = pos - midPtVel * timeInterval;
    newPos.x = newPos.x < minX ? minX : newPos.x;
    newPos.x = newPos.x > maxX ? maxX : newPos.x;
    newPos.y = newPos.y < minY ? minY : newPos.y;
    newPos.y = newPos.y > maxY ? maxY : newPos.y;

    (*_density)[idx] = _density->sample(newPos - gridSpacing() * .5f);
    (*_temperature)[idx] = _temperature->sample(newPos - gridSpacing() * .5f);
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::computeAdvection(double timeInterval, AdvectType type)
  {
    auto N = size().x;
    Index2 index;
    if(type == AdvectType::Density)
      for (index.y = 1; index.y <= N; ++index.y)
        for (index.x = 1; index.x <= N; ++index.x)
          advectDensity(index, timeInterval);
    else if (type == AdvectType::Velocity)
      for (index.y = 1; index.y <= N; ++index.y)
        for (index.x = 1; index.x <= N; ++index.x)
          advectVelocity(index, timeInterval);

    applyBoundaryCondition();
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::updateTempMinMax()
  {
    forEachIndex<D>(_temperature->size(), [&](const Index<D>& index) {
      auto t = (*_temperature)[index];
      if (t < _minTemp)
        _minTemp = t;
      if (t > _maxTemp)
        _maxTemp = t;
      });
  }
  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::computeBuoyancy(double timeInterval)
  {
    auto vPos = _velocity->positionInSpace<1>();
    auto N = size().x;
    auto Tamb = 0.0;
    forEachIndex<D>(_temperature->size(), [&](const Index<D>& index) {
      auto t = (*_temperature)[index];
      if (t < _minTemp)
        _minTemp = t;
      if (t > _maxTemp)
        _maxTemp = t;
      Tamb += t;
      });
    Tamb /= N * N;

    Index2 index;
    for (index.y = 1; index.y <= N; ++index.y)
    {
      for (index.x = 1; index.x <= N; ++index.x)
      {
        auto pos = vPos(index);
        _velocity->velocityAt<1>(index) += timeInterval * (_densityBuoyancyFactor * _density->sample(pos) +
          + _temperatureBuoyancyFactor * (_temperature->sample(pos) - Tamb));
      }
    }

    applyBoundaryCondition();
  }

  template<size_t D, typename real>
  inline void GridSolver<D, real>::computeSource(double timeInterval)
  {    
    applyBoundaryCondition();
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::computeGravity(double timeInterval)
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
    GridSolver<D, real>::fluidSdf() const
  {
    return new ConstantScalarField<D, real>(math::Limits<real>::inf());
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::applyBoundaryCondition()
  {
    /*auto depth = static_cast<unsigned int>(std::ceil(_maxCfl));
    _boundaryConditionSolver.constrainVelocity(_velocity, depth);*/
    auto N = size().x;
    for (int i = 1; i <= N; i++)
    {
      _velocity->velocityAt<0>(Index2(1, i)) = 0;
      _velocity->velocityAt<0>(Index2(N+1, i)) = 0;

      _velocity->velocityAt<1>(Index2(i, 1)) = 0;
      _velocity->velocityAt<1>(Index2(i, N+1)) = 0;

      (*_density)[Index2(0, i)] = (*_density)[Index2(1, i)];
      (*_density)[Index2(N+1, i)] = (*_density)[Index2(N, i)];
      (*_density)[Index2(i, 0)] = (*_density)[Index2(i, 1)];
      (*_density)[Index2(i, N+1)] = (*_density)[Index2(i, N)];
    }
    (*_density)[Index2(0, 0)] = .5f* ((*_density)[Index2(1, 0)]+ (*_density)[Index2(0, 1)]);
    (*_density)[Index2(0, N+1)] = .5f * ((*_density)[Index2(1, N+1)] + (*_density)[Index2(0, N)]);
    (*_density)[Index2(N+1, 0)] = .5f * ((*_density)[Index2(N, 0)] + (*_density)[Index2(N+1, 1)]);
    (*_density)[Index2(N+1, N+1)] = .5f * ((*_density)[Index2(N, N+1)] + (*_density)[Index2(N+1, N)]);
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::applyBoundaryCondition(CellCenteredScalarGrid<D, real>& grid)
  {
    /*auto depth = static_cast<unsigned int>(std::ceil(_maxCfl));
    _boundaryConditionSolver.constrainVelocity(_velocity, depth);*/
    auto N = size().x;
    for (int i = 1; i <= N; i++)
    {
      grid[Index2(0, i)] = grid[Index2(1, i)];
      grid[Index2(N + 1, i)] = grid[Index2(N, i)];
      grid[Index2(i, 0)] = grid[Index2(i, 1)];
      grid[Index2(i, N + 1)] = grid[Index2(i, N)];
    }
    grid[Index2(0, 0)] = .5f * (grid[Index2(1, 0)] + grid[Index2(0, 1)]);
    grid[Index2(0, N + 1)] = .5f * (grid[Index2(1, N + 1)] + grid[Index2(0, N)]);
    grid[Index2(N + 1, 0)] = .5f * (grid[Index2(N, 0)] + grid[Index2(N + 1, 1)]);
    grid[Index2(N + 1, N + 1)] = .5f * (grid[Index2(N, N + 1)] + grid[Index2(N + 1, N)]);
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::extrapolateIntoCollider(CellCenteredScalarGrid<D, real>& grid)
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
    GridSolver<D, real>::colliderSdf() const
  {
    return _boundaryConditionSolver.colliderSdf();
  }

  template<size_t D, typename real>
  inline VectorField<D, real>*
    GridSolver<D, real>::colliderVelocityField() const
  {
    return _boundaryConditionSolver.colliderVelocityField();
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::beginAdvanceTimeStep(double timeInterval)
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
    GridSolver<D, real>::endAdvanceTimeStep(double timeInterval)
  {
    // Invoke callback
    onEndAdvanceTimeStep(timeInterval);
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::updateCollider(double timeInterval)
  {
    if (_collider != nullptr)
      _collider->update(this->currentTime(), timeInterval);
  }

  template<size_t D, typename real>
  inline void
    GridSolver<D, real>::updateEmitter(double timeInterval)
  {
    // TODO
  }

} // end namespace cg

#endif // __GridSolver_h
