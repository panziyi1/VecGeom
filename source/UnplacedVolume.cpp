#include "VecGeom/volumes/UnplacedVolume.h"
#include "VecGeom/volumes/PlacedVolume.h"
#include "VecGeom/base/SOA3D.h"
#include "VecGeom/volumes/utilities/VolumeUtilities.h"

#include "VecGeom/volumes/UnplacedBox.h"
#include "VecGeom/volumes/UnplacedTrd.h"
#include "VecGeom/volumes/UnplacedTube.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

using Vec3D  = Vector3D<Precision>;
using Real_v = vecgeom::VectorBackend::Real_v;

// generic implementation for SamplePointOnSurface
Vector3D<Precision> VUnplacedVolume::SamplePointOnSurface() const
{
  Vector3D<Precision> surfacepoint;
  SOA3D<Precision> points(1);
  volumeUtilities::FillRandomPoints(*this, points);

  Vector3D<Precision> dir = volumeUtilities::SampleDirection();
  surfacepoint            = points[0] + DistanceToOut(points[0], dir) * dir;

  // assert( Inside(surfacepoint) == vecgeom::kSurface );
  return surfacepoint;
}

// trivial implementations for the interface functions
// (since we are moving to these interfaces only gradually)

// ---------------- Normal ---------------------------------------------------------------------

VECCORE_ATT_HOST_DEVICE
bool VUnplacedVolume::Normal(Vector3D<Precision> const &p, Vector3D<Precision> &normal) const
{
#ifndef VECCORE_CUDA
  throw std::runtime_error("unimplemented function called");
#endif
  return false;
}

// ----------------- Extent --------------------------------------------------------------------
VECCORE_ATT_HOST_DEVICE
void VUnplacedVolume::Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const
{
#ifndef VECCORE_CUDA
  throw std::runtime_error("unimplemented function called");
#endif
}

// estimating the surface area by sampling
// based on the six-point method of G4
Precision VUnplacedVolume::EstimateSurfaceArea(int nStat) const
{
  static const Precision s2 = 1. / Sqrt(2.);
  static const Precision s3 = 1. / Sqrt(3.);

  // Predefined directions
  //
  static const Vec3D directions[64] = {
      Vec3D(0, 0, 0),      Vec3D(-1, 0, 0),    // (  ,  ,  ) ( -,  ,  )
      Vec3D(1, 0, 0),      Vec3D(-1, 0, 0),    // ( +,  ,  ) (-+,  ,  )
      Vec3D(0, -1, 0),     Vec3D(-s2, -s2, 0), // (  , -,  ) ( -, -,  )
      Vec3D(s2, -s2, 0),   Vec3D(0, -1, 0),    // ( +, -,  ) (-+, -,  )

      Vec3D(0, 1, 0),      Vec3D(-s2, s2, 0), // (  , +,  ) ( -, +,  )
      Vec3D(s2, s2, 0),    Vec3D(0, 1, 0),    // ( +, +,  ) (-+, +,  )
      Vec3D(0, -1, 0),     Vec3D(-1, 0, 0),   // (  ,-+,  ) ( -,-+,  )
      Vec3D(1, 0, 0),      Vec3D(-1, 0, 0),   // ( +,-+,  ) (-+,-+,  )

      Vec3D(0, 0, -1),     Vec3D(-s2, 0, -s2),   // (  ,  , -) ( -,  , -)
      Vec3D(s2, 0, -s2),   Vec3D(0, 0, -1),      // ( +,  , -) (-+,  , -)
      Vec3D(0, -s2, -s2),  Vec3D(-s3, -s3, -s3), // (  , -, -) ( -, -, -)
      Vec3D(s3, -s3, -s3), Vec3D(0, -s2, -s2),   // ( +, -, -) (-+, -, -)

      Vec3D(0, s2, -s2),   Vec3D(-s3, s3, -s3), // (  , +, -) ( -, +, -)
      Vec3D(s3, s3, -s3),  Vec3D(0, s2, -s2),   // ( +, +, -) (-+, +, -)
      Vec3D(0, 0, -1),     Vec3D(-s2, 0, -s2),  // (  ,-+, -) ( -,-+, -)
      Vec3D(s2, 0, -s2),   Vec3D(0, 0, -1),     // ( +,-+, -) (-+,-+, -)

      Vec3D(0, 0, 1),      Vec3D(-s2, 0, s2),   // (  ,  , +) ( -,  , +)
      Vec3D(s2, 0, s2),    Vec3D(0, 0, 1),      // ( +,  , +) (-+,  , +)
      Vec3D(0, -s2, s2),   Vec3D(-s3, -s3, s3), // (  , -, +) ( -, -, +)
      Vec3D(s3, -s3, s3),  Vec3D(0, -s2, s2),   // ( +, -, +) (-+, -, +)

      Vec3D(0, s2, s2),    Vec3D(-s3, s3, s3), // (  , +, +) ( -, +, +)
      Vec3D(s3, s3, s3),   Vec3D(0, s2, s2),   // ( +, +, +) (-+, +, +)
      Vec3D(0, 0, 1),      Vec3D(-s2, 0, s2),  // (  ,-+, +) ( -,-+, +)
      Vec3D(s2, 0, s2),    Vec3D(0, 0, 1),     // ( +,-+, +) (-+,-+, +)

      Vec3D(0, 0, -1),     Vec3D(-1, 0, 0),    // (  ,  ,-+) ( -,  ,-+)
      Vec3D(1, 0, 0),      Vec3D(-1, 0, 0),    // ( +,  ,-+) (-+,  ,-+)
      Vec3D(0, -1, 0),     Vec3D(-s2, -s2, 0), // (  , -,-+) ( -, -,-+)
      Vec3D(s2, -s2, 0),   Vec3D(0, -1, 0),    // ( +, -,-+) (-+, -,-+)

      Vec3D(0, 1, 0),      Vec3D(-s2, s2, 0), // (  , +,-+) ( -, +,-+)
      Vec3D(s2, s2, 0),    Vec3D(0, 1, 0),    // ( +, +,-+) (-+, +,-+)
      Vec3D(0, -1, 0),     Vec3D(-1, 0, 0),   // (  ,-+,-+) ( -,-+,-+)
      Vec3D(1, 0, 0),      Vec3D(-1, 0, 0),   // ( +,-+,-+) (-+,-+,-+)
  };

  // Get bounding box
  //
  Vec3D bmin, bmax;
  this->Extent(bmin, bmax);
  Vec3D bdim = bmax - bmin;

  // Define statistics and shell thickness
  //
  int npoints      = (nStat < 1000) ? 1000 : nStat;
  Precision coeff  = 0.5 / Cbrt(Precision(npoints));
  Precision eps    = coeff * bdim.Min(); // half thickness
  Precision twoeps = 2. * eps;
  Precision del    = 1.8 * eps; // six-point offset - should be more than sqrt(3.)

  // Enlarge bounding box by eps
  //
  bmin -= Vec3D(eps);
  bdim += Vec3D(twoeps);

  // Calculate surface area
  //
  int icount = 0;
  for (int i = 0; i < npoints; ++i) {
    Precision px = bmin.x() + bdim.x() * RNG::Instance().uniform();
    Precision py = bmin.y() + bdim.y() * RNG::Instance().uniform();
    Precision pz = bmin.z() + bdim.z() * RNG::Instance().uniform();
    Vec3D p(px, py, pz);
    EnumInside in  = this->Inside(p);
    Precision dist = 0;
    if (in == EInside::kInside) {
      if (this->SafetyToOut(p) >= eps) continue;
      int icase = 0;
      if (this->Inside(Vec3D(px - del, py, pz)) != EInside::kInside) icase += 1;
      if (this->Inside(Vec3D(px + del, py, pz)) != EInside::kInside) icase += 2;
      if (this->Inside(Vec3D(px, py - del, pz)) != EInside::kInside) icase += 4;
      if (this->Inside(Vec3D(px, py + del, pz)) != EInside::kInside) icase += 8;
      if (this->Inside(Vec3D(px, py, pz - del)) != EInside::kInside) icase += 16;
      if (this->Inside(Vec3D(px, py, pz + del)) != EInside::kInside) icase += 32;
      if (icase == 0) continue;
      Vec3D v = directions[icase];
      dist    = this->DistanceToOut(p, v);
      Vec3D n;
      this->Normal(p + v * dist, n);
      dist *= v.Dot(n);
    } else if (in == EInside::kOutside) {
      if (this->SafetyToIn(p) >= eps) continue;
      int icase = 0;
      if (this->Inside(Vec3D(px - del, py, pz)) != EInside::kOutside) icase += 1;
      if (this->Inside(Vec3D(px + del, py, pz)) != EInside::kOutside) icase += 2;
      if (this->Inside(Vec3D(px, py - del, pz)) != EInside::kOutside) icase += 4;
      if (this->Inside(Vec3D(px, py + del, pz)) != EInside::kOutside) icase += 8;
      if (this->Inside(Vec3D(px, py, pz - del)) != EInside::kOutside) icase += 16;
      if (this->Inside(Vec3D(px, py, pz + del)) != EInside::kOutside) icase += 32;
      if (icase == 0) continue;
      Vec3D v = directions[icase];
      dist    = this->DistanceToIn(p, v);
      if (dist == kInfLength) continue;
      Vec3D n;
      this->Normal(p + v * dist, n);
      dist *= -(v.Dot(n));
    }
    if (dist < eps) icount++;
  }
  return bdim.x() * bdim.y() * bdim.z() * icount / npoints / twoeps;
}

// estimating the cubic volume by sampling
// based on the method of G4
Precision VUnplacedVolume::EstimateCapacity(int nStat) const
{
  Precision epsilon = 1E-4;

  // limits
  if (nStat < 100) nStat = 100;

  Vector3D<Precision> lower, upper, offset;
  this->Extent(lower, upper);
  offset                        = 0.5 * (upper + lower);
  const Vector3D<Precision> dim = 0.5 * (upper - lower);

  int insidecounter = 0;
  for (int i = 0; i < nStat; i++) {
    auto p = offset + volumeUtilities::SamplePoint(dim);
    if (this->Contains(p)) insidecounter++;
  }
  return 8. * (dim[0] + epsilon) * (dim[1] + epsilon) * (dim[2] + epsilon) * insidecounter / nStat;
}

std::ostream &operator<<(std::ostream &os, VUnplacedVolume const &vol)
{
  vol.Print(os);
  return os;
}

#define DISPATCH_TO_UNPLACED_VOLUME_KIND(kind, class, func, args) \
  case VolumeTypes::kind: \
    return static_cast<const class *>(this)->class::func args

#define DISPATCH_TO_UNPLACED_VOLUME(func, args) \
  switch (this->fType) { \
  DISPATCH_TO_UNPLACED_VOLUME_KIND(kBox, UnplacedBox, func, args); \
  DISPATCH_TO_UNPLACED_VOLUME_KIND(kTrd, UnplacedTrd, func, args); \
  DISPATCH_TO_UNPLACED_VOLUME_KIND(kTube, UnplacedTube, func, args); \
  default: break; \
  }
  // return this->func args

bool VUnplacedVolume::Contains(Vector3D<Precision> const &pos) const
{
  DISPATCH_TO_UNPLACED_VOLUME(Contains, (pos));
  return 0;
}

EnumInside VUnplacedVolume::Inside(Vector3D<Precision> const &pos) const
{
  DISPATCH_TO_UNPLACED_VOLUME(Inside, (pos));
  return kOutside;
}

Precision VUnplacedVolume::DistanceToOut(Vector3D<Precision> const &pos, Vector3D<Precision> const &dir,
                                         Precision step_max) const
{
  DISPATCH_TO_UNPLACED_VOLUME(DistanceToOut, (pos, dir, step_max));
  return step_max;
}

Precision VUnplacedVolume::SafetyToOut(Vector3D<Precision> const &pos) const
{
  DISPATCH_TO_UNPLACED_VOLUME(SafetyToOut, (pos));
  return 0;
}

Precision VUnplacedVolume::DistanceToIn(Vector3D<Precision> const &position, Vector3D<Precision> const &direction,
                                        const Precision step_max) const
{
  DISPATCH_TO_UNPLACED_VOLUME(DistanceToIn, (position, direction, step_max));
  return step_max;
}

Precision VUnplacedVolume::SafetyToIn(Vector3D<Precision> const &pos) const
{
  DISPATCH_TO_UNPLACED_VOLUME(SafetyToIn, (pos));
  return 0;
}

#ifndef VECCORE_CUDA

VPlacedVolume *VUnplacedVolume::PlaceVolume(char const *const label, LogicalVolume const *const volume,
                                            Transformation3D const *const transformation,
                                            VPlacedVolume *const placement) const
{
  VPlacedVolume *const placed = PlaceVolume(volume, transformation, placement);
  placed->set_label(label);
  return placed;
}

#endif
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom
