#ifndef __VectorField_h
#define __VectorField_h

#include "math/Vector4.h"

namespace cg
{

template <size_t D, typename real>
class VectorField
{
public:
  ASSERT_REAL(real, "*VectorField: real must be float or double");
  static_assert(
    D == 2 || D == 3, "Not implemented - D should be either 2 or 3."
    );

  using vec_type = Vector<real, D>;

  VectorField()
  {
    // do nothing
  }

  virtual ~VectorField()
  {
    // do nothing
  }

  // Retorna o valor interpolado para a posicao x
  virtual vec_type sample(const vec_type& x) const = 0;

  // Retorna o divergente para a posicao x
  virtual real divergence(const vec_type& x) const
  {
    return static_cast<real>(0.0);
  }

  // Retorna o rotacional para posicao x
  virtual vec_type curl(const vec_type& x) const
  {
    return vec_type::null();
  }
}; // VectorField

} // end namespace cg

#endif // __VectorField_h

