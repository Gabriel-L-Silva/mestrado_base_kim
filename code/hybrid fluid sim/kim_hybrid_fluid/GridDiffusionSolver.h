#ifndef __GridDiffusionSolver_h
#define __GridDiffusionSolver_h

#include "FaceCenteredGrid.h"
#include "ScalarGrid.h"
#include "ConstantScalarField.h"

namespace cg
{
template <size_t D, typename real>
class GridDiffusionSolver : public SharedObject
{
public:
  using GridReference = Reference<FaceCenteredGrid<D, real>>;
  GridDiffusionSolver()
  {
    // do nothing
  }

  virtual ~GridDiffusionSolver()
  {
    // do nothing
  }


  virtual void solve(
    const GridReference source,
    real diffusionCoefficient,
    double timeInterval,
    GridReference dest,
    const ScalarField<D, real>& boundarySdf
    = ConstantScalarField<D, real>(real(0.0)),
    const ScalarField<D, real>& fluidSdf
    = ConstantScalarField<D, real>(real(0.0))
  ) = 0;

}; // GridDiffusionSolver<D, real>

} // end namespace cg

#endif // __GridDiffusionSolver_h
