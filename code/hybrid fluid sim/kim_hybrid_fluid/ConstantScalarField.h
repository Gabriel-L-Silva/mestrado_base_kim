#ifndef __ConstantScalarField_h
#define __ConstantScalarField_h

#include "ScalarField.h"

namespace cg
{

/**
* D-dimensional constant scalar field.
* 
* \tparam D Defines the number of dimensions.
* \tparam real A floating point type.
*/
template <size_t D, typename real>
class ConstantScalarField final : public ScalarField<D, real>
{
public:
  using vec_type = Vector<real, D>; ///< Vector type alias.

  /** Constructs a constant scalar field with given \p value. */
  explicit ConstantScalarField(real value):
    _value(value)
  {
    // do nothing
  }

  /** Returns the sampled value at position \p x. */
  real sample(const vec_type& x) const
  {
    return _value;
  }

private:
  real _value = static_cast<real>(0.0);
}; // ConstantScalarField

} // end namespace cg

#endif // __ConstantScalarField_h
