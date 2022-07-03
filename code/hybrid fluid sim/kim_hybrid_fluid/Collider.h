#ifndef __Collider_h
#define __Collider_h

#include "math/Surface.h"
#include <functional>
#include <cassert>

namespace cg
{

/**
* Abstract template class for generic collider.
* 
* This class defines basic interfaces for colliders. Most of the
* functionalities are implemented at this template class, except the member
* method Collider<D, real>::velocityAt. The subclasses should provide a
* Surface<D, real> instance to define the collider surface using the
* Collider<D, real>::setSurface.
* 
* \tparam D Defines the number of dimensions.
* \tparam real A floating point type.
*/
template <size_t D, typename real>
class Collider : public SharedObject
{
public:
	using vec_type = Vector<real, D>; ///< Vector type alias.
	/**
	* Callback function type for update calls.
	*
	* The callbacks of this type take the collider pointer, current time and the
	* time interval in seconds.
	*/
	using OnBeginUpdateCallback = std::function<void(Collider<D, real>*, real, real)>;

	/** Default constructor. */
	Collider()
	{
		// do nothing
	}


	/** Default destructor. */
	virtual ~Collider()
	{
		// do nothing
	}

	/** Returns the velocity of the collider at point \p p. */
	virtual vec_type velocityAt(const vec_type& p) const = 0;


	/**
	* Resolves collision for given point.
	* 
	* \see Fluid Engine Development by Doyub Kim.
	* 
	* \param[in] radius Radius of the colliding point.
	* \param[in] restitutionCoefficient Defines the restitution effect.
	* \param[in, out] position Input and output position of the point.
	* \param[in, out] velocity Input and output velocity of the point.
	*/
	void resolveCollision(
		real radius,
		real restitutionCoefficient,
		vec_type& position,
		vec_type& velocity
	);

	/** Returns the friction coefficient. */
	real frictionCoefficient() const
	{
		return _frictionCoefficient;
	}

	/**
	* Sets the friction coefficient.
	* 
	* This method assigns the friction coefficient to the collider. Any
	* negative inputs will be set to zero.
	* 
	* \param[in] frictionCoefficient New friction coefficient.
	*/
	void setFrictionCoefficient(real frictionCoefficient)
	{
		_frictionCoefficient = math::max(frictionCoefficient, (real)0.0f);
	}

	/** Returns a Reference to this collider's surface. */
	auto surface() const
	{
		return _surface;
	}

	/** Updates the collider state. */
	void update(double currentTime, double timeInterval)
	{
#ifdef _DEBUG
		assert(_surface.get());
#endif // _DEBUG
		if (_onUpdateCallback)
			_onUpdateCallback(this, currentTime, timeInterval);
	}

	/**
	* \brief Sets the callback function that should be called when
	* Collider<D, real>::update is invoked.
	* 
	* \see OnBeginUpdateCallback
	* 
	* \param[in] callback The callback function.
	*/
	void setOnBeginUpdateCallback(const OnBeginUpdateCallback& callback)
	{
		_onUpdateCallback = callback;
	}

protected:
	/** Internal query result structure. */
	struct ColliderQueryResult final {
		real distance;
		vec_type point;
		vec_type normal;
		vec_type velocity;
	};

	/** Assigns the surface instance for this collider. */
	void setSurface(const Reference<math::Surface<D, real>> surface)
	{
		_surface = surface;
	}

	/** Computes the closest point's information. */
	void getClosestPoint(const Reference<math::Surface<D, real>> surface, const vec_type& queryPoint, ColliderQueryResult& result) const;

	/** Returns \c true if given point is inside the surface. */
	bool isPenetrating(const ColliderQueryResult& colliderPoint, const vec_type& position, real radius);

private:
	/** Reference to surface. */
	Reference<math::Surface<D, real>> _surface;
	/** Collider friction coefficient. */
	real _frictionCoefficient = 0.0f;
	/** Update callback. */
	OnBeginUpdateCallback _onUpdateCallback;

}; // Collider

template<size_t D, typename real>
inline void
Collider<D, real>::resolveCollision(real radius, real restitutionCoefficient, vec_type& position, vec_type& velocity)
{
#ifdef _DEBUG
	assert(_surface.get());
#endif // _DEBUG

	ColliderQueryResult colliderPoint;
	getClosestPoint(_surface, position, colliderPoint);

	if (isPenetrating(colliderPoint, position, radius))
	{
		auto targetNormal = colliderPoint.normal;
		auto targetPoint = colliderPoint.point + radius * targetNormal;
		auto colliderVelAtTargetPoint = colliderPoint.velocity;

		auto relativeVel = velocity - colliderVelAtTargetPoint;
		auto normalDotRelativeVel = targetNormal.dot(relativeVel);
		auto relativeVelN = normalDotRelativeVel * targetNormal;
		auto relativeVelT = relativeVel - relativeVelN;

		if (normalDotRelativeVel < 0.0f)
		{
			auto deltaRelativeVelN = (-restitutionCoefficient - 1.0f) * relativeVelN;
			relativeVelN = relativeVelN * (-restitutionCoefficient);

			if (relativeVelT.squaredNorm() > 0.0f)
			{
				auto frictionScale = math::max<real>(1.0f - _frictionCoefficient * deltaRelativeVelN.length() / relativeVelT.length(), 0.0f);
				relativeVelT *= frictionScale;
			}

			velocity = relativeVelN + relativeVelT + colliderVelAtTargetPoint;
		}
		position = targetPoint;
	}
}

template<size_t D, typename real>
inline void
Collider<D, real>::getClosestPoint(const Reference<math::Surface<D, real>> surface, const vec_type& queryPoint, ColliderQueryResult& result) const
{
	result.distance = surface->closestDistance(queryPoint);
	result.point = surface->closestPoint(queryPoint);
	result.normal = surface->closestNormal(queryPoint);
	result.velocity = velocityAt(queryPoint);
}

template<size_t D, typename real>
inline bool
Collider<D, real>::isPenetrating(const ColliderQueryResult& colliderPoint, const vec_type& position, real radius)
{
	// If the new candidate position of the particle is inside
	// the volume defined by the surface OR the new distance to the surface is
	// less than the particle's radius, this particle is in colliding state.
	return _surface->isInside(position) || math::isNegative(colliderPoint.distance - radius);
}

} // end namespace cg

#endif // __Collider_h
