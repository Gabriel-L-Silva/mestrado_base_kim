#ifndef __Utils_h
#define __Utils_h

#include "math/Vector2.h"

inline float
frand()
{
  return (float)((double)rand() / RAND_MAX);
}

inline float
frand(float min, float max)
{
  return min + frand() * (max - min);
}

inline cg::vec2f
vrand()
{
  return {frand(), frand()};
}

inline cg::vec2f
vrand(float min, float max)
{
  return {frand(min, max), frand(min, max)};
}

#endif // __Utils_h
