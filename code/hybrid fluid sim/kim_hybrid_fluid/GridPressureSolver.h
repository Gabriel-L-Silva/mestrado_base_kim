#ifndef __GridPressureSolver_h
#define __GridPressureSolver_h

#include "FaceCenteredGrid.h"
#include "ConstantScalarField.h"
#include "ConstantVectorField.h"

namespace cg
{

/// <summary>
/// Abstract template base class for grid-based pressure solver.
/// </summary>
/// <typeparam name="real">Floating point type</typeparam>
template <size_t D, typename real>
class GridPressureSolver : public SharedObject
{
public:
  template <typename T> using Ref = Reference<T>;
  GridPressureSolver()
  {
    // do nothing
  }

  virtual ~GridPressureSolver()
  {
    // do nothing
  }

  virtual void solve(
    const Ref<FaceCenteredGrid<D, real>> input,
    double timeInterval,
    Ref<FaceCenteredGrid<D, real>> output,
    const ScalarField<D, real>& boundarySdf = ConstantScalarField<D, real>(math::Limits<real>::inf()),
    const ScalarField<D, real>& fluidSdf = ConstantScalarField<D, real>(-math::Limits<real>::inf()),
    const VectorField<D, real>& boundaryVelocity = ConstantVectorField<D, real>(vec_type{ real(0.0f) })) = 0;

};

} // end namespace cg

#endif // __GridPressureSolver_h

