#ifndef __CustomVectorField_h
#define __CustomVectorField_h

#include <functional>
#include "VectorField.h"

namespace cg
{

/**
* Custom D-dimensional vector field template class.
* 
* This class implements a vector field with a user-provided function.
* 
* \tparam D Defines the number of dimensions.
* \tparam real A floating point type.
*/
template <size_t D, typename real>
class CustomVectorField final: public VectorField<D, real>
{
public:
  using vec_type = Vector<real, D>; ///< Vector type alias.

  /**
  * Constructs a field with given function.
  * 
  * This constructor creates a field with user-provided function.
  * To compute derivatives, such as gradient and Laplacian, finite
  * differencing is used. Thus, the differencing resolution also can be
  * provided as the last parameter.
  */
  CustomVectorField(
    const std::function<vec_type(const vec_type&)>& func,
    real resolution = 1e-3
  ) : _customFunc(func),
    _resolution(resolution)
  {
    // do nothing
  }

  /**
  * Constructs a field with given field and gradient function.
  * 
  * This constructor creates a field with user-provided field and gradient
  * function. To compute Laplacian, finite differencing is used. Thus, the
  * differencing resolution also can be provided as the last parameter.
  */
  CustomVectorField(
    const std::function<vec_type(const vec_type&)>& func,
    const std::function<real(const vec_type&)>& divergenceFunc,
    real resolution = 1e-3
  ) : _customFunc(func),
    _divergenceFunc(divergenceFunc),
    _resolution(resolution)
  {
    // do nothing
  }


  /** Constructs a field with given field, gradient, and Laplacian function. */
  CustomVectorField(
    const std::function<vec_type(const vec_type&)>& func,
    const std::function<real(const vec_type&)>& divergenceFunc,
    const std::function<vec_type(const vec_type&)>& curlFunc
  ) : _customFunc(func),
    _divergenceFunc(divergenceFunc),
    _curlFunc(curlFunc)
  {
    // do nothing
  }

  /** Returns the sampled value at given position \p x. */
  vec_type sample(const vec_type& x) const override
  {
    return _customFunc(x);
  }

  /** Returns the divergence at given position \p x. */
  real divergence(const vec_type& x) const override
  {
    if (_divergenceFunc)
      return _divergenceFunc(x);
    else
    {
      real acc = 0.0f;
      for (size_t i = 0; i < D; ++i)
      {
        auto _x = vec_type::null();
        _x[int(i)] = 0.5f * _resolution;
        acc += (_customFunc(x + _x) - _customFunc(x - _x))[int(i)];
      }
      return acc / _resolution;
    }
  }

  /**
  * Returns the curl at given position \p x.
  * 
  * \todo Code the case where no curl function was provided.
  */
  vec_type curl(const vec_type& x) const override
  {
    if (_curlFunc)
      return _curlFunc(x);
    else
    {
      // TODO!!
      return vec_type::null();
    }
  }

private:
  /** Field function. */
  std::function<vec_type(const vec_type&)> _customFunc;
  /** Divergence function. */
  std::function<real(const vec_type&)> _divergenceFunc;
  /** Curl function. */
  std::function<vec_type(const vec_type&)> _curlFunc;
  /** Finite differencing resolution. */
  real _resolution = 1e-3;

}; // CustomVectorField

} // end namespace cg

#endif // __CustomVectorField_h
