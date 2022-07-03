#ifndef __ScalarField_h
#define __ScalarField_h

#include "math/Vector4.h"

namespace cg
{

template <size_t D, typename real>
class ScalarField
{
public:
  ASSERT_REAL(real, "*ScalarField: real must be float or double");
  static_assert(
    D == 2 || D == 3, "Not implemented - D should be either 2 or 3."
  );

  using vec_type = Vector<real, D>;

  ScalarField()
  {
    // do nothing
  }

  virtual ~ScalarField()
  {
    // do nothing
  }

  // Retorna o valor interpolado para a posicao x
  virtual real sample(const vec_type& x) const = 0;

  // Retorna o vetor gradiente para a posicao x
  virtual vec_type gradient(const vec_type& x) const
  {
    return vec_type();
  }

  // Retorna o Laplaciano para posicao x
  virtual real laplacian(const vec_type& x) const
  {
    return static_cast<real>(0.0);
  }
};

}  // end namespace cg

#endif  // __ScalarField_h
