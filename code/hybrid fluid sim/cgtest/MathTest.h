#include "geometry/Point3.h"
#include "graphics/Color.h"
#include "math/Matrix4x4.h"

inline void
testVec3f()
{
  using namespace cg;

  puts("**vec3f test**");

  vec3f a{1, 2, 3};
  vec3f b{vec4f{4, 5, 6, 1}};

  a.print("a=");
  b.print("b=");
  (a + b).print("a+b=");
  (a - b).print("a-b=");
  (a * b).print("a*b=");
  (a * 2).print("a*2=");
  (2 * a).print("2*a=");
  (-a).print("-a=");
  a.inverse().print("1/a=");
  printf("length(a)=%g\n", a.length());
  printf("a.b=%g\n", a.dot(b));
  (a.cross(b)).print("axb=");
}

inline void
testVec4f()
{
  using namespace cg;

  puts("**vec4f test**");

  vec4f a{1, 2, 3};
  vec4f b{vec3f{4, 5, 6}};

  a.print("a=");
  b.print("b=");
  (a + b).print("a+b=");
  (a - b).print("a-b=");
  (a * b).print("a*b=");
  (a * 2).print("a*2=");
  (2 * a).print("2*a=");
  (-a).print("-a=");
  a.inverse().print("1/a=");
  printf("length(a)=%g\n", a.length());
  printf("a.b=%g\n", a.dot(b));
}

inline void
testQuatf()
{
  using namespace cg;

  puts("**quatf test**");

  quatf a{45, vec3f::up()};
  quatf b{quatf::eulerAngles(30, 45, 60)};

  a.print("a=");
  b.print("b=");
  (a + b).print("a+b=");
  (a - b).print("a-b=");
  (a * b).print("a*b=");
  (a * 2).print("a*2=");
  (2 * a).print("2*a=");
  (-a).print("-a=");
  (~a).print("~a=");
  a.inverse().print("1/a=");
  printf("length(a)=%g\n", a.length());
  (a * vec3f::up()).print("a*y=");
}

inline void
testMat3f()
{
  using namespace cg;

  puts("**mat3f test**");

  mat3f a{quatf{45, vec3f::up()}};
  mat3f b;

  a.inverse(b);
  a.print("a=");
  b.print("b=");
  (a * b).print("a*b=");
  (a * 2).print("a*2=");
  (2 * a).print("2*a=");
  (a * vec3f::up()).print("a*y=");
}

inline void
testMat4f()
{
  using namespace cg;

  puts("**mat4f test**");

  mat4f a{mat4f::TRS(vec3f::null(), quatf{45, vec3f::up()}, vec3f{1.f})};
  mat4f c{mat4f::rotation(quatf::identity(), vec3f::null())};
  mat4f b;

  a.inverse(b);
  a.print("a=");
  b.print("b=");
  c.print("c=");
  (a * b).print("a*b=");
  (a * 2).print("a*2=");
  (2 * a).print("2*a=");
  (c.transform3x4(vec3f::up())).print("a*y=");
  mat4f::lookAt(vec3f{0, 0, 1}, vec3f::null(), vec3f::up()).print("lookAt=");
  mat4f::perspective(60, 1, 0.1f, 100).print("perspective=");
  mat4f::frustum(-1, 1, -1, 1, 0.1f, 100).print("frustum=");
  mat4f::ortho(-1, 1, -1, 1, 0.1f, 100).print("ortho=");
}

inline void
testColor()
{
  using namespace cg;

  puts("**color test**");

  Color a{38, 80, 120};
  Color b{vec4f{0.5f, 0.5f, 0.5f}};

  printf("a=rgb(%g,%g,%g)\n", a.x, a.y, a.z);
  a.print("a=");
  b.print("b=");
  (a * b).print("a*b=");
  (a * 2).print("a*2=");
  (2 * a).print("2*a=");
}

inline void
testPoint2f()
{
  using namespace cg;

  puts("**Point2f test**");

  Point2f a{1, 2};
  vec2f v{3, 4};
  Point2f b{v};

  a.print("a=");
  b.print("b=");
  (a + b).print("a+b=");
  (a + v).print("a+v=");
  (a - b).print("a-b=");
  (a - v).print("a-v=");
  (a * 2).print("a*2=");
  (2 * a).print("2*a=");
  (-a).print("-a=");
  (a = v).print("a=v=");
  printf(a == b ? "a==b\n" : "a!=b\n");
}

inline void
testPoint3f()
{
  using namespace cg;

  puts("**Point3f test**");

  Point3f a{1, 2, 3};
  vec3f v{4, 5, 6};
  Point3f b{v};

  a.print("a=");
  b.print("b=");
  (a + b).print("a+b=");
  (a + v).print("a+v=");
  (a - b).print("a-b=");
  (a - v).print("a-v=");
  (a * 2).print("a*2=");
  (2 * a).print("2*a=");
  (-a).print("-a=");
  (a = v).print("a=v=");
  printf(a == b ? "a==b\n" : "a!=b\n");
}

int
testMath()
{
  /*
  testVec3f();
  printf("\n");
  testVec4f();
  printf("\n");
  testQuatf();
  printf("\n");
  testMat3f();
  printf("\n");
  testMat4f();
  printf("\n");
  testColor();
  printf("\n");
  */
  testPoint2f();
  printf("\n");
  testPoint3f();
  puts("\nPress any key to exit...");
  getchar();
  return 0;
}
