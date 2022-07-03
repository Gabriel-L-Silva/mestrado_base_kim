#ifndef __GridBoundaryConditionSolver_h
#define __GridBoundaryConditionSolver_h

#include "FaceCenteredGrid.h"
#include "ScalarField.h"
#include "Collider.h"
#include "Constants.h"

namespace cg
{

template <size_t D, typename real>
class GridBoundaryConditionSolver : public SharedObject
{
public:
  using vec_type = typename FaceCenteredGrid<D, real>::vec_type;

  virtual ~GridBoundaryConditionSolver()
  {
    // do nothing
  }

  const auto collider() const
  {
    return _collider;
  }

  void updateCollider(
    const Reference<Collider<D, real>> newCollider,
    const Index<D>& size,
    const vec_type& spacing,
    const vec_type& origin
  )
  {
    _collider = newCollider;
    _gridSize = size;
    _gridSpacing = spacing;
    _gridOrigin = origin;

    onColliderUpdated(size, spacing, origin);
  }

  // Returns the closed domain boundary flag.
  int closedDomainBoundaryFlag() const
  {
    return _closedDomainBoundaryFlag;
  }

  void setClosedDomainBoundaryFlag(int flag)
  {
    _closedDomainBoundaryFlag = flag;
  }

  //! Sets the closed domain boundary flag.
  //void setClosedDomainBoundaryFlag(int flag);

  virtual void constrainVelocity(Reference<FaceCenteredGrid<D, real>> grid, unsigned extrapolationDepth = 5) = 0;

  virtual ScalarField<D, real>* colliderSdf() const = 0;

  virtual VectorField<D, real>* colliderVelocityField() const = 0;

protected:
  virtual void onColliderUpdated(const Index<D>& size, const vec_type& spacing, const vec_type& origin) = 0;

  const Index<D>& gridSize() const
  {
    return _gridSize;
  }

  const vec_type& gridSpacing() const
  {
    return _gridSpacing;
  }

  const vec_type& gridOrigin() const
  {
    return _gridOrigin;
  }

private:
  Reference<Collider<D, real>> _collider;
  Index<D> _gridSize;
  vec_type _gridSpacing;
  vec_type _gridOrigin;
  int _closedDomainBoundaryFlag = cg::constants::directionAll;

}; // GridBoundaryConditionSolver

} // end namespace cg

#endif // __GridBoundaryConditionSolver_h
