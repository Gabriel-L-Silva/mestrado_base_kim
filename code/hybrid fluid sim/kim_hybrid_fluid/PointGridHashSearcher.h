#ifndef __PointGridHashSearcher_h
#define __PointGridHashSearcher_h

#include <vector>
#include <functional>
#include "geometry/ParticleSystem.h"
#include "geometry/Index3.h"

namespace cg
{

template <size_t D, typename real, typename PointArray>
class PointGridHashSearcher final: public SharedObject
{
public:
  using vec_type = Vector<real, D>;
  using id_type = typename Index<D>::base_type;

  PointGridHashSearcher(const Index<D>& size, real gridSpacing)
    : _size(size), _gridSpacing(gridSpacing)
  {
    // TODO: clamp size to positive values only
    // do nothing
  }

  void build(const PointArray& points)
  {
    _buckets.clear();
    _points.clear();

    if (points.size() == 0)
      return;

    // allocate memory
    _buckets.resize(_size.prod());
    _points.resize(points.size());
    
    // put points into buckets
    for (size_t i = 0; i < points.size(); ++i)
    {
      _points[i] = points[i];
      size_t key = getHashKeyFromPosition(points[i]);
      _buckets[key].push_back(i);
    }
  }

  /// <summary>
  /// Invokes the callback for each nearby point around the origin within given radius
  /// </summary>
  /// <param name="o">: origin</param>
  /// <param name="radius">: radius to search</param>
  /// <param name="callback">: callback function</param>
  void forEachNearbyPoint(const vec_type& o, real radius, const std::function<void(size_t, const vec_type&)>& callback) const
  {
    if (_buckets.empty())
      return;

    constexpr size_t n = 1 << D;
    size_t nearbyKeys[n];
    getNearbyKeys(o, nearbyKeys);

    const double queryRadiusSquared = radius * radius;

    for (size_t i = 0; i < n; ++i)
    {
      const auto& bucket = _buckets[nearbyKeys[i]];
      size_t pointsInBucket = bucket.size();

      for (size_t j = 0; j < pointsInBucket; ++j)
      {
        auto pointIndex = bucket[j];
        auto rSquared = (_points[pointIndex] - o).squaredNorm();
        if (rSquared <= queryRadiusSquared)
        {
          callback(pointIndex, _points[pointIndex]);
        }
      }
    }
  }

private:
  real _gridSpacing = 1.0f;
  Index<D> _size = Index<D>(1LL);
  std::vector<vec_type> _points;
  std::vector<std::vector<id_type>> _buckets;

  Index<D> getBucketIndex(const vec_type& position) const
  {
    auto invSpacing = math::inverse(_gridSpacing);
    return Index<D>(position * invSpacing);
  }

  size_t getHashKeyFromBucketIndex(const Index<D>& index) const
  {
    auto wrappedIndex = index;
    for (size_t i = 0; i < D; ++i)
    {
      wrappedIndex[i] = index[(int)i] % _size[i];
      if (wrappedIndex[i] < 0)
        wrappedIndex[i] += _size[i];
    }

    if constexpr (D == 3)
      return static_cast<size_t>((wrappedIndex.z * _size.y + wrappedIndex.y) * _size.x + wrappedIndex.x);
    else
      return static_cast<size_t>(wrappedIndex.y * _size.x + wrappedIndex.x);
  }

  size_t getHashKeyFromPosition(const vec_type& position) const
  {
    auto index = getBucketIndex(position);
    return getHashKeyFromBucketIndex(index);
  }

  void getNearbyKeys(const vec_type& position, size_t* nearbyKeys) const
  {
    constexpr size_t n = 1 << D;
    Index<D> originIndex = getBucketIndex(position);
    Index<D> nearbyBucketIndices[n];

    for (size_t i = 0; i < n; ++i)
    {
      nearbyBucketIndices[i] = originIndex;
    }

    if ((originIndex.x + 0.5f) * _gridSpacing <= position.x)
    {
      if constexpr (D == 2)
      {
        nearbyBucketIndices[2].x += 1;
        nearbyBucketIndices[3].x += 1;
      }
      else
      {
        nearbyBucketIndices[4].x += 1;
        nearbyBucketIndices[5].x += 1;
        nearbyBucketIndices[6].x += 1;
        nearbyBucketIndices[7].x += 1;
      }
    }
    else
    {
      if constexpr (D == 2)
      {
        nearbyBucketIndices[2].x -= 1;
        nearbyBucketIndices[3].x -= 1;
      }
      else
      {
        nearbyBucketIndices[4].x -= 1;
        nearbyBucketIndices[5].x -= 1;
        nearbyBucketIndices[6].x -= 1;
        nearbyBucketIndices[7].x -= 1;
      }
    }

    if ((originIndex.y + 0.5f) * _gridSpacing <= position.y)
    {
      if constexpr (D == 2)
      {
        nearbyBucketIndices[1].y += 1;
        nearbyBucketIndices[3].y += 1;
      }
      else
      {
        nearbyBucketIndices[2].y += 1;
        nearbyBucketIndices[3].y += 1;
        nearbyBucketIndices[6].y += 1;
        nearbyBucketIndices[7].y += 1;
      }
    }
    else
    {
      if constexpr (D == 2)
      {
        nearbyBucketIndices[1].y -= 1;
        nearbyBucketIndices[3].y -= 1;
      }
      else
      {
        nearbyBucketIndices[2].y -= 1;
        nearbyBucketIndices[3].y -= 1;
        nearbyBucketIndices[6].y -= 1;
        nearbyBucketIndices[7].y -= 1;
      }
    }
    if constexpr (D == 3)
    {
      if ((originIndex.z + 0.5f) * _gridSpacing <= position.z)
      {
        nearbyBucketIndices[1].z += 1;
        nearbyBucketIndices[3].z += 1;
        nearbyBucketIndices[5].z += 1;
        nearbyBucketIndices[7].z += 1;
      }
      else
      {
        nearbyBucketIndices[1].z -= 1;
        nearbyBucketIndices[3].z -= 1;
        nearbyBucketIndices[5].z -= 1;
        nearbyBucketIndices[7].z -= 1;
      }
    }
    
    for (size_t i = 0; i < n; ++i)
    {
      nearbyKeys[i] = getHashKeyFromBucketIndex(nearbyBucketIndices[i]);
    }
  }
};

} // end namespace cg

#endif // __PointGridHashSearcher_h
