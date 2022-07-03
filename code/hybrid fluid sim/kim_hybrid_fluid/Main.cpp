#include "graphics/Application.h"
#include "GLSimulationWindow.h"
#include <cstdlib>
#include <inttypes.h>

using namespace cg;

int
main(int argc, char** argv)
{
  return cg::Application{ new GLSimulationWindow<float>("SimulationWindow", 721, 720) }.run(argc, argv);
  /*Index2 size{ 3, 3 };
  auto backwardEuler = GridBackwardEulerDiffusionSolver<2, float, false>();
  Reference<FaceCenteredGrid<2, float>> grid = new FaceCenteredGrid<2, float>(size, vec2f{ 1.0f }, vec2f::null());
  Reference<CellCenteredScalarGrid<2, float>> fluidSdf = new CellCenteredScalarGrid<2, float>(size, vec2f{ 1.0f }, vec2f::null());
  for (auto& d : *fluidSdf)
  {
    d = -0.5f;
  }

  grid->fill(vec2f{ 10.0f });
  grid->velocityAt<0>(0) = 0.0f;
  grid->velocityAt<0>(3) = -2.0f;
  grid->velocityAt<0>(6) = 8.0f;
  grid->velocityAt<0>(7) = -8.0f;


  backwardEuler.solve(grid, 1.0f, 1.0 / 100.0, grid, ConstantScalarField<2, float>(10000.0f), *fluidSdf);*/

  return 0;
}