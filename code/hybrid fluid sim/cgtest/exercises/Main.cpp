#include "graphics/Application.h"
#include "ExTriangle.h"
#include "ExRandomTriangles.h"
#include "ExRGBCube.h"
#include "ExTexture.h"

template <typename T, typename ...Args>
inline int
run(int argc, char** argv, Args... args)
{
  return cg::Application{new T{1280, 720, args...}}.run(argc, argv);
}

int
main(int argc, char** argv)
{
  //return run<ExTriangle>(argc, argv);
  //return run<ExRandomTriangles>(argc, argv, 1000);
  return run<ExRGBCube>(argc, argv);
  //return run<ExTexture>(argc, argv);
}