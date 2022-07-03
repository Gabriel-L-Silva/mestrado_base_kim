#ifndef __GridFractionalBoundaryConditionSolver_h
#define __GridFractionalBoundaryConditionSolver_h

#include "GridBoundaryConditionSolver.h"
#include "CellCenteredScalarGrid.h"
#include "CustomVectorField.h"
#include "GridUtils.h"

namespace cg
{

template <size_t D, typename real>
class GridFractionalBoundaryConditionSolver : public GridBoundaryConditionSolver<D, real>
{
public:
  using vec_type = Vector<real, D>;
  using id_type = typename Index<D>::base_type;
  using Base = GridBoundaryConditionSolver<D, real>;

  virtual ~GridFractionalBoundaryConditionSolver()
  {
    if (_colliderVel)
      delete _colliderVel;
  }

  void constrainVelocity(Reference<FaceCenteredGrid<D, real>> grid, unsigned extrapolationDepth = 5) override;

  ScalarField<D, real>* colliderSdf() const override
  {
    return _colliderSdf.get();
  }

  // DEFINE COLLIDERVEL TODO!!
  VectorField<D, real>* colliderVelocityField() const override
  {
    return _colliderVel;
  }

protected:
  void onColliderUpdated(const Index<D>& size, const vec_type& spacing, const vec_type& origin) override;

private:
  Reference<CellCenteredScalarGrid<D, real>> _colliderSdf;
  CustomVectorField<D, real>* _colliderVel;

}; // GridFractionalBoundaryConditionSolver

template<size_t D, typename real>
inline void
GridFractionalBoundaryConditionSolver<D, real>::constrainVelocity(Reference<FaceCenteredGrid<D, real>> grid, unsigned extrapolationDepth)
{
  auto size = grid->size();

  if (_colliderSdf == nullptr || _colliderSdf->size() != size)
  {
    this->updateCollider(
      this->collider(), size, grid->gridSpacing(), grid->origin());
  }

  // preparing data
  std::array<
    std::function<vec_type(const Index<D>&)>, D> positions;
  positions[0] = grid->positionInSpace<0>();
  positions[1] = grid->positionInSpace<1>();
  if constexpr (D == 3)
    positions[2] = grid->positionInSpace<2>();

  std::array<GridData<D, real>, D> temps;
  temps[0].resize(grid->iSize<0>());
  temps[1].resize(grid->iSize<1>());
  if constexpr (D == 3)
    temps[2].resize(grid->iSize<2>());

  std::array<GridData<D, char>, D> markers;
  markers[0].resize(grid->iSize<0>());
  markers[1].resize(grid->iSize<1>());
  if constexpr (D == 3)
    markers[2].resize(grid->iSize<2>());

  // initialize data
  for (size_t i = 0; i < D; ++i)
  {
    fill(temps[i], static_cast<real>(0.0f));
    fill(markers[i], (char)1);
  }

  auto h = grid->gridSpacing();

  grid->forEachIndex<0>([&](const Index<D>& index) {
    auto pt = positions[0](index);
    auto c = vec_type(static_cast<real>(0.0f));
    c.x = h.x * 0.5f;
    auto phi0 = _colliderSdf->sample(pt - c);
    auto phi1 = _colliderSdf->sample(pt + c);
    auto frac = fractionInsideSdf(phi0, phi1);
    frac = 1.0f - math::clamp<real>(frac, 0.0f, 1.0f);

    auto id = markers[0].id(index);
    if (frac > 0.0f)
      markers[0][id] = 1;
    else
    {
      auto colliderVel = this->collider()->velocityAt(pt);
      grid->velocityAt<0>(index) = colliderVel.x;
      markers[0][id] = 0;
    }
  });

  grid->forEachIndex<1>([&](const Index<D>& index) {
    auto pt = positions[1](index);
    auto c = vec_type(static_cast <real>(0.0f));
    c.y = h.y * 0.5f;
    auto phi0 = _colliderSdf->sample(pt - c);
    auto phi1 = _colliderSdf->sample(pt + c);
    auto frac = fractionInsideSdf(phi0, phi1);
    frac = 1.0f - math::clamp<real>(frac, 0.0f, 1.0f);

    auto id = markers[1].id(index);
    if (frac > 0.0f)
      markers[1][id] = 1;
    else
    {
      auto colliderVel = this->collider()->velocityAt(pt);
      grid->velocityAt<1>(index) = colliderVel.y;
      markers[1][id] = 0;
    }
  });

  if constexpr (D == 3)
  {
    grid->forEachIndex<2>([&](const Index<D>& index) {
      auto pt = positions[2](index);
      auto c = vec_type(static_cast <real>(0.0f));
      c.z = h.z * 0.5f;
      auto phi0 = _colliderSdf->sample(pt - c);
      auto phi1 = _colliderSdf->sample(pt + c);
      auto frac = fractionInsideSdf(phi0, phi1);
      frac = 1.0f - math::clamp<real>(frac, 0.0f, 1.0f);

      auto id = markers[2].id(index);
      if (frac > 0.0f)
        markers[2][id] = 1;
      else
      {
        auto colliderVel = this->collider()->velocityAt(pt);
        grid->velocityAt<2>(index) = colliderVel.z;
        markers[2][id] = 0;
      }
    });
  }
  // Free-slip: Extrapolate fluid velocity into the collider
  extrapolateToRegion(*grid->data<0>(), markers[0], extrapolationDepth, *grid->data<0>());
  extrapolateToRegion(*grid->data<1>(), markers[1], extrapolationDepth, *grid->data<1>());
  if constexpr (D == 3)
    extrapolateToRegion(*grid->data<2>(), markers[2], extrapolationDepth, *grid->data<2>());

  // No-flux: project the extrapolated velocity to the collider's surface
  grid->forEachIndex<0>([&](const Index<D>& index) {
    auto pt = positions[0](index);
    auto id = temps[0].id(index);
    if (isInsideSdf(_colliderSdf->sample(pt)))
    {
      auto colliderVel = this->collider()->velocityAt(pt);
      auto vel = grid->sample(pt);
      auto g = _colliderSdf->gradient(pt);
      if (g.squaredNorm() > 0.0f)
      {
        auto velr = vel - colliderVel;
        auto velt = projectAndApplyFriction(
          velr, g.normalize(), this->collider()->frictionCoefficient());
        auto velp = velt + colliderVel;
        temps[0][id] = velp.x;
      }
      else
        temps[0][id] = colliderVel.x;
    }
    else
      temps[0][id] = grid->velocityAt<0>(index);
    grid->velocityAt<0>(index) = temps[0][id];
  });

  grid->forEachIndex<1>([&](const Index<D>& index) {
    auto pt = positions[1](index);
    auto id = temps[1].id(index);
    if (isInsideSdf(_colliderSdf->sample(pt)))
    {
      auto colliderVel = this->collider()->velocityAt(pt);
      auto vel = grid->sample(pt);
      auto g = _colliderSdf->gradient(pt);
      if (g.squaredNorm() > 0.0f)
      {
        auto velr = vel - colliderVel;
        auto velt = projectAndApplyFriction(
          velr, g.normalize(), this->collider()->frictionCoefficient());
        auto velp = velt + colliderVel;
        temps[1][id] = velp.y;
      }
      else
        temps[1][id] = colliderVel.y;
    }
    else
      temps[1][id] = grid->velocityAt<1>(index);
    grid->velocityAt<1>(index) = temps[1][id];
  });

  if constexpr(D == 3)
    grid->forEachIndex<2>([&](const Index<D>& index) {
      auto pt = positions[2](index);
      auto id = temps[2].id(index);
      if (isInsideSdf(_colliderSdf->sample(pt)))
      {
        auto colliderVel = this->collider()->velocityAt(pt);
        auto vel = grid->sample(pt);
        auto g = _colliderSdf->gradient(pt);
        if (g.squaredNorm() > 0.0f)
        {
          auto velr = vel - colliderVel;
          auto velt = projectAndApplyFriction(
            velr, g.normalize(), this->collider()->frictionCoefficient());
          auto velp = velt + colliderVel;
          temps[2][id] = velp.z;
        }
        else
          temps[2][id] = colliderVel.z;
      }
      else
        temps[2][id] = grid->velocityAt<2>(index);
      grid->velocityAt<2>(index) = temps[2][id];
    });

  // No-flux: Project velocity on the domain boundary if closed
  auto flag = this->closedDomainBoundaryFlag();
  if (flag & constants::directionLeft)
  {
    if constexpr (D == 2)
      for (id_type j = 0; j < grid->iSize<0>().y; ++j)
        grid->velocityAt<0>(Index2{ 0, j }) = 0;
    else
      for (id_type k = 0; k < grid->iSize<0>().z; ++k)
        for (id_type j = 0; j < grid->iSize<0>().y; ++j)
          grid->velocityAt<0>(Index3{ 0, j, k }) = 0;
  }

  if (flag & constants::directionRight)
  {
    if constexpr (D == 2)
      for (id_type j = 0; j < grid->iSize<0>().y; ++j)
        grid->velocityAt<0>(Index2{ grid->iSize<0>().x - 1, j }) = 0;
    else
      for (id_type k = 0; k < grid->iSize<0>().z; ++k)
        for (id_type j = 0; j < grid->iSize<0>().y; ++j)
          grid->velocityAt<0>(Index3{ grid->iSize<0>().x - 1, j, k }) = 0;
  }

  if (flag & constants::directionDown)
  {
    if constexpr (D == 2)
      for (id_type i = 0; i < grid->iSize<1>().x; ++i)
        grid->velocityAt<1>(Index2{ i, 0 }) = 0;
    else
      for (id_type k = 0; k < grid->iSize<1>().z; ++k)
        for (id_type i = 0; i < grid->iSize<1>().x; ++i)
          grid->velocityAt<1>(Index3{ i, 0, k }) = 0;
  }

  if (flag & constants::directionUp)
  {
    if constexpr (D == 2)
      for (id_type i = 0; i < grid->iSize<1>().x; ++i)
        grid->velocityAt<1>(Index2{ i, grid->iSize<1>().y - 1 }) = 0;
    else
      for (id_type k = 0; k < grid->iSize<1>().z; ++k)
        for (id_type i = 0; i < grid->iSize<1>().x; ++i)
          grid->velocityAt<1>(Index3{ i, grid->iSize<1>().y - 1, k }) = 0;
  }

  if constexpr (D == 3)
  {
    if (flag & constants::directionBack)
    {
      for (id_type j = 0; j < grid->iSize<2>().y; ++j)
        for (id_type i = 0; i < grid->iSize<2>().x; ++i)
          grid->velocityAt<2>(Index3{ i, j, 0 }) = 0;
    }
    
    if (flag & constants::directionFront)
    {
      for (id_type j = 0; j < grid->iSize<2>().y; ++j)
        for (id_type i = 0; i < grid->iSize<2>().x; ++i)
          grid->velocityAt<2>(Index3{ i, j, grid->iSize<2>().z - 1 }) = 0;
    }
  }
}

template<size_t D, typename real>
inline void
GridFractionalBoundaryConditionSolver<D, real>::onColliderUpdated(const Index<D>& size, const vec_type& spacing, const vec_type& origin)
{
  if (_colliderSdf == nullptr)
    _colliderSdf = new CellCenteredScalarGrid<D, real>(size, spacing, origin);

  if (this->collider() != nullptr)
  {
    // TODO
    auto surface = this->collider()->surface();
  }
  else
  {
    for (auto& v : *_colliderSdf)
      v = math::Limits<real>::inf();

    if (_colliderVel != nullptr)
      delete _colliderVel;
    _colliderVel = new CustomVectorField<D, real>(
      [](const vec_type&) {
        return vec_type::null();
      },
      this->gridSpacing().x
    );
  }
}

} // end namespace cg

#endif // __GridFractionalBoundaryConditionSolver_h
