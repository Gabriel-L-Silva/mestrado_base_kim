#include "graphics/Application.h"
#include "FramebufferTest.h"
#include "GLTestSuite.h"
#include "LineTest.h"
#include "QuadTest.h"
#include "TextureTest.h"
#include "TriangleTest.h"

inline auto
testSuite()
{
  return new GLTestSuite({new LineTest,
    new TriangleTest,
    new QuadTest,
    new TextureTest,
    new MultiTextureTest,
    new FramebufferTest});
}

int
main(int argc, char** argv)
{
  return cg::Application{testSuite()}.run(argc, argv);
}
