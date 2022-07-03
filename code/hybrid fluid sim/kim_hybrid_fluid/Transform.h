#ifndef __Transform_h
#define __Transform_h

#include "math/Matrix4x4.h"

namespace cg
{ // begin namespace cg

template <size_t D, typename real>
class Transform;

template <typename real>
class Transform<2, real> final
{
public:
  using vec2 = Vector2<real>;

  Transform():
    _position{real(0.0f)},
    _angle{real(0.0f)},
    _cosAngle{real(1.0f)},
    _sinAngle{real(0.0f)}
  {
    // do nothing
  }
  const vec2& position() const
  {
    return _position;
  }

  const real& rotation() const
  {
    return _angle;
  }

  const real& eulerAngles() const
  {
    return _angle;
  }

  const vec2& scale() const
  {
    return _scale;
  }

  /// Sets the local position of this transform.
  void setPosition(const vec2& position)
  {
    _position = position;
  }

  /// Sets the local rotation of this transform.
  void setRotation(real angle)
  {
    setEulerAngles(angle);
  }

  /// Sets the local Euler angles (in degrees) of this transform.
  void setEulerAngles(real angle)
  {
    _angle = angle;
    _cosAngle = std::cos(angle);
    _sinAngle = std::sin(angle);
  }

  /// Sets the local scale of this transform.
  void setScale(const vec2& scale)
  {
    _scale = scale;
  }

  /// Sets the local uniform scale of this transform.
  void setScale(real scale)
  {
    setScale(vec2{ scale });
  }

  /// Returns the direction of the world Y axis of this transform.
  vec2 up() const
  {
    return vec2{ -_sinAngle, _cosAngle };
  }

  /// Returns the direction of the world Z axis of this transform.
  vec2 right() const
  {
    return vec2{ _cosAngle, _sinAngle };
  }

  /// Translates this transform.
  void translate(const vec2& t)
  {
    setPosition(_position + t);
  }

  /// Rotates this transform.
  void rotate(const float& angle)
  {
    // TODO
  }

  /// Transforms \c p from local space to world space.
  vec2 transform(const vec2& p) const
  {
    return vec2{
      _cosAngle * p.x - _sinAngle * p.y + _position.x,
      _sinAngle * p.x + _cosAngle * p.y + _position.y
    };
  }

  /// Transforms \c p from world space to local space.
  vec2 inverseTransform(const vec2& p) const
  {
    const auto t = p - _position;
    return vec2{
       _cosAngle * t.x + _sinAngle * t.y,
      -_sinAngle * t.x + _cosAngle * t.y
    };
  }

  /// Transforms \c v from local space to world space.
  vec2 transformVector(const vec2& v) const
  {
    return vec2{
      _cosAngle * v.x - _sinAngle * v.y,
      _sinAngle * v.x + _cosAngle * v.y
    };
  }

  /// Transforms \c v from world space to local space.
  vec2 inverseTransformVector(const vec2& v) const
  {
    return vec2{
       _cosAngle * v.x + _sinAngle * v.y,
      -_sinAngle * v.x + _cosAngle * v.y
    };
  }

  /// Transforms \c d from world space to local space.
  vec2 transformDirection(const vec2& d) const
  {
    return inverseTransformVector(d);
  }

  /// Sets this transform as an identity transform.
  void reset();

private:
  vec2 _position;
  vec2 _scale;
  real _angle;
  real _cosAngle = 1.0f;
  real _sinAngle = 0.0f;
}; // Transform<2, real>

template<typename real>
inline void
Transform<2, real>::reset()
{
  _position = vec2::null();
  _scale = vec2{ real(1.0f) };
  _angle = _sinAngle = real(0.0f);
  _cosAngle = real(1.0f);
}

template <typename real>
class Transform<3, real> final
{
public:
  using vec = Vector3<real>;
  using quat = Quaternion<real>;
  using mat4 = Matrix4x4<real>;

  /// Constructs an identity transform.
  Transform():
    _position{real(0.0f)},
    _rotation{quat::identity()},
    _scale{real(1.0f)},
    _matrix{real(1.0f)}
  {
    _inverseMatrix = _matrix;
  }

  /// Returns the world position of this transform.
  const vec& position() const
  {
    return _position;
  }

  /// Returns the world rotation of this transform.
  const quat& rotation() const
  {
    return _rotation;
  }

  /// Returns the world Euler angles (in degrees) of this transform.
  const vec& eulerAngles() const
  {
    return _rotation.eulerAngles();
  }

  /// Returns the local scale of this transform.
  const vec& scale() const
  {
    return _scale;
  }

  /// Sets the local position of this transform.
  void setPosition(const vec& position)
  {
    _position = position;
    update();
  }

  /// Sets the local rotation of this transform.
  void setRotation(const quat& rotation)
  {
    _rotation = rotation;
    update();
  }

  /// Sets the local Euler angles (in degrees) of this transform.
  void setEulerAngles(const vec& angles)
  {
    _rotation = quat::eulerAngles(angles);
    update();
  }

  /// Sets the local scale of this transform.
  void setScale(const vec& scale)
  {
    _scale = scale;
    update();
  }

  /// Sets the local uniform scale of this transform.
  void setScale(real scale)
  {
    setScale(vec{ scale });
  }

  /// Returns the direction of the world Z axis of this transform.
  vec forward() const
  {
    return _rotation * vec{ 0, 0, 1 };
  }

  /// Returns the direction of the world Y axis of this transform.
  vec up() const
  {
    return _rotation * vec::up();
  }

  /// Returns the direction of the world Z axis of this transform.
  vec right() const
  {
    return _rotation * vec{ 1, 0, 0 };
  }

  /// Translates this transform.
  void translate(const vec& t)
  {
    setPosition(_position + t);
  }

  /// Rotates this transform.
  void rotate(const vec& angles)
  {
    rotate(quat::eulerAngles(angles));
  }

  /// Rotates this transform.
  void rotate(const vec& axis, real angle)
  {
    rotate(quat{ angle, axis });
  }

  /// Returns the local to world _matrix of this transform.
  const mat4& localToWorldMatrix() const
  {
    return _matrix;
  }

  /// Returns the world to local _matrix of this transform.
  const mat4& worldToLocalMatrix() const
  {
    return _inverseMatrix;
  }

  /// Transforms \c p from local space to world space.
  vec transform(const vec& p) const
  {
    return _matrix.transform3x4(p);
  }

  /// Transforms \c p from world space to local space.
  vec inverseTransform(const vec& p) const
  {
    return _inverseMatrix.transform3x4(p);
  }

  /// Transforms \c v from local space to world space.
  vec transformVector(const vec& v) const
  {
    return _matrix.transformVector(v);
  }

  /// Transforms \c v from world space to local space.
  vec inverseTransformVector(const vec& v) const
  {
    return _inverseMatrix.transformVector(v);
  }

  /// Transforms \c d from world space to local space.
  vec transformDirection(const vec& d) const
  {
    return _rotation.rotate(d);
  }

  /// Sets this transform as an identity transform.
  void reset();

private:
  vec _position;
  quat _rotation;
  vec _scale;
  mat4 _matrix;
  mat4 _inverseMatrix;

  mat4 localMatrix() const;

  mat4 inverseLocalMatrix() const;

  void rotate(const quat&);
  void update();

}; // Transform<3, real>

template <typename real>
inline Matrix4x4<real>
Transform<3, real>::localMatrix() const
{
  return mat4::TRS(_position, _rotation, _scale);
}

template <typename real>
inline Matrix4x4<real>
Transform<3, real>::inverseLocalMatrix() const
{
  Matrix3x3<real> r{ _rotation };

  r[0] *= math::inverse(_scale[0]);
  r[1] *= math::inverse(_scale[1]);
  r[2] *= math::inverse(_scale[2]);

  mat4 m;

  m[0].set(r[0][0], r[1][0], r[2][0]);
  m[1].set(r[0][1], r[1][1], r[2][1]);
  m[2].set(r[0][2], r[1][2], r[2][2]);
  m[3][0] = -(r[0].dot(_position));
  m[3][1] = -(r[1].dot(_position));
  m[3][2] = -(r[2].dot(_position));
  m[3][3] = 1.0f;
  return m;
}

template <typename real>
inline void
Transform<3, real>::rotate(const quat& q)
{
  setRotation(_rotation * (_rotation.inverse() * q * _rotation));
}

template <typename real>
inline void
Transform<3, real>::reset()
{
  _position = vec{ 0.0f };
  _rotation = quat::identity();
  _scale = vec{ 1.0f };
  _matrix = _inverseMatrix = mat4{ 1.0f };
}

template <typename real>
inline void
Transform<3, real>::update()
{
  _matrix = localMatrix();
  _inverseMatrix = inverseLocalMatrix();
}

} // end namespace cg

#endif // __Transform_h
