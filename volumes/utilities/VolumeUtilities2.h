/*
 * volume_utilities2.h
 */

#ifndef VOLUME_UTILITIES2_H_
#define VOLUME_UTILITIES2_H_

#include "base/Vector3D.h"
#include "base/Global.h"
#include "volumes/PlacedBox.h"
#include "volumes/LogicalVolume.h"
#include <cstdio>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {
namespace volumeUtilities {

/**
 * @brief Returns a random point, based on a sampling rectangular volume.
 * @details Mostly used for benchmarks and navigation tests
 * @param size is a Vector3D containing the rectangular dimensions of the sampling volume
 * @param scale an optional scale factor (default is 1)
 * @return a random output point
 */
template <typename RngEngine>
VECGEOM_FORCE_INLINE
Vector3D<Precision> SamplePoint(Vector3D<Precision> const &size, RngEngine &rngengine, const Precision scale = 1)
{
  std::uniform_real_distribution<double> dist(0, 2.);
  const Vector3D<Precision> ret(scale * (1. - dist(rngengine)) * size[0], scale * (1. - dist(rngengine)) * size[1],
                                scale * (1. - dist(rngengine)) * size[2]);
  return ret;
}

// *** The following functions allow to give an external generator
// *** which should make these functions usable in parallel

/**
 * @brief Fills the volume with 3D points which are _not_ contained in
 *    any daughters of the input mother volume.
 * @details Requires a proper bounding box from the input volume.
 *    Point coordinates are local to input mother volume.
 * @param volume is the input mother volume containing all output points.
 * @param points is the output container, provided by the caller.
 */
template <typename RandomEngine, typename TrackContainer>
VECGEOM_FORCE_INLINE
void FillUncontainedPoints(VPlacedVolume const &volume, RandomEngine &rngengine, TrackContainer &points)
{
  static double lastUncontCap = 0.0;
  double uncontainedCapacity  = UncontainedCapacity(volume);
  if (uncontainedCapacity != lastUncontCap) {
    printf("Uncontained capacity for %s: %g units\n", volume.GetLabel().c_str(), uncontainedCapacity);
    lastUncontCap = uncontainedCapacity;
  }
  if (uncontainedCapacity <= 1000 * kTolerance) {
    std::cout << "\nVolUtil: FillUncontPts: ERROR: Volume provided <" << volume.GetLabel()
              << "> does not have uncontained capacity!  Aborting.\n";
    assert(false);
  }

  const int size = points.capacity();
  points.resize(points.capacity());

  Vector3D<Precision> lower, upper, offset;
  volume.Extent(lower, upper);
  offset                        = 0.5 * (upper + lower);
  const Vector3D<Precision> dim = 0.5 * (upper - lower);

  int tries = 0;
  for (int i = 0; i < size; ++i) {
    bool contained;
    Vector3D<Precision> point;
    tries = 0;
    do {
      // ensure that point is contained in mother volume
      do {
        ++tries;
        if (tries % 1000000 == 0) {
          printf("%s line %i: Warning: %i tries to find uncontained points... volume=%s.  Please check.\n", __FILE__,
                 __LINE__, tries, volume.GetLabel().c_str());
        }

        point = offset + SamplePoint(dim, rngengine);
      } while (!volume.UnplacedContains(point));
      points.set(i, point);

      contained = false;
      int kk    = 0;
      for (Vector<Daughter>::const_iterator j = volume.GetDaughters().cbegin(), jEnd = volume.GetDaughters().cend();
           j != jEnd; ++j, ++kk) {
        if ((*j)->Contains(points[i])) {
          contained = true;
          break;
        }
      }
    } while (contained);
  }
}

template <typename RandomEngine, typename TrackContainer>
VECGEOM_FORCE_INLINE
void FillUncontainedPoints(LogicalVolume const &volume, RandomEngine &rngengine, TrackContainer &points)
{
  VPlacedVolume const *const placed = volume.Place();
  FillUncontainedPoints(*placed, rngengine, points);
  delete placed;
}

} // end namespace volumeUtilities
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif /* VOLUME_UTILITIES2_H_ */
