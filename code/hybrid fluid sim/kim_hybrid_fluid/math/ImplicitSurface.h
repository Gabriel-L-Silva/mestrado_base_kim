#ifndef __ImplicitSurface_h
#define __ImplicitSurface_h

#include "Surface.h"

namespace cg::math
{

template <size_t D, typename real>
class ImplicitSurface : public Surface<D, real>
{
public:
  using Base = Surface<D, real>;

  ImplicitSurface(const Transform<D>& t)
    : Base(t)
  {
    // do nothing
  }

  virtual ~ImplicitSurface()
  {

  }

  // TODO!!!!!!!!!!!!!!
}; // ImplicitSurface

} // end namespace cg::math

#endif // __ImplicitSurface_h