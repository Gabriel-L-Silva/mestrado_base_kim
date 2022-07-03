#ifndef __ConstantVectorField_h
#define __ConstantVectorField_h

#include "VectorField.h"

namespace cg
{

/**
* D-dimensional constant vector field.
* 
* \tparam D Defines the number of dimensions.
* \tparam real A floating point type.
*/
template <size_t D, typename real>
class ConstantVectorField: public VectorField<D, real>
{
public:
  using vec_type = Vector<real, D>; ///< Vector type alias.

  /** Constructs a constant vector field with given \p value. */
  explicit ConstantVectorField(const vec_type& value):
    _value(value)
  {
    // do nothing
  }

  /** Returns the sampled value at position \p x. */
  vec_type sample(const vec_type& x) const
  {
    return _value;
  }

private:
  /** Constant value for this vector field. */
  vec_type _value;
}; // ConstantVectorField

} // end namespace cg

#endif // __ConstantVectorField_h
