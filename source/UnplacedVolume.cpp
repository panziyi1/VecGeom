#include "volumes/UnplacedVolume.h"
#include "volumes/PlacedVolume.h"

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

using Real_v = vecgeom::VectorBackend::Real_v;

// trivial implementations for the interface functions
// (since we are moving to these interfaces only gradually)

// ---------------- DistanceToOut functions -----------------------------------------------------

// the USolid/GEANT4-like interface for DistanceToOut (returning also exiting normal)
VECCORE_ATT_HOST_DEVICE
Precision VUnplacedVolume::DistanceToOut(Vector3D<Precision> const &p, Vector3D<Precision> const &d,
                                         Vector3D<Precision> &normal, bool &convex, Precision step_max) const
{
#ifndef VECCORE_CUDA
  throw std::runtime_error("unimplemented function called");
#endif
  return -1.;
}

// an explicit SIMD interface
VECCORE_ATT_HOST_DEVICE
Real_v VUnplacedVolume::DistanceToOutVec(Vector3D<Real_v> const &p, Vector3D<Real_v> const &d,
                                         Real_v const &step_max) const
{
#ifndef VECCORE_CUDA
  throw std::runtime_error("unimplemented function called");
#endif
  return Real_v(-1.);
}

// the container/basket interface (possibly to be deprecated)
void VUnplacedVolume::DistanceToOut(SOA3D<Precision> const &points, SOA3D<Precision> const &directions,
                                    Precision const *const step_max, Precision *const output) const
{
#ifndef VECCORE_CUDA
  throw std::runtime_error("unimplemented function called");
#endif
}

// ---------------- SafetyToOut functions -----------------------------------------------------
// an explicit SIMD interface
Real_v VUnplacedVolume::SafetyToOutVec(Vector3D<Real_v> const &p) const
{
#ifndef VECCORE_CUDA
  throw std::runtime_error("unimplemented function called");
#endif
  return Real_v(-1.);
}

// the container/basket interface (possibly to be deprecated)
void VUnplacedVolume::SafetyToOut(SOA3D<Precision> const &points, Precision *const output) const
{
#ifndef VECCORE_CUDA
  throw std::runtime_error("unimplemented function called");
#endif
}

// ---------------- DistanceToIn functions -----------------------------------------------------
VECCORE_ATT_HOST_DEVICE
Real_v VUnplacedVolume::DistanceToInVec(Vector3D<Real_v> const &position, Vector3D<Real_v> const &direction,
                                        const Real_v &step_max) const
{
#ifndef VECCORE_CUDA
  throw std::runtime_error("unimplemented function called");
#endif
  return Real_v(-1.);
}

// ---------------- SafetyToIn functions -------------------------------------------------------
// explicit SIMD interface
Real_v VUnplacedVolume::SafetyToInVec(Vector3D<Real_v> const &p) const
{
#ifndef VECCORE_CUDA
  throw std::runtime_error("unimplemented function called");
#endif
  return Real_v(-1.);
}

// ---------------- Normal ---------------------------------------------------------------------

VECCORE_ATT_HOST_DEVICE
bool VUnplacedVolume::Normal(Vector3D<Precision> const &p, Vector3D<Precision> &normal) const
{
#ifndef VECCORE_CUDA
  throw std::runtime_error("unimplemented function called");
#endif
  return false;
}

// ---------------- SamplePointOnSurface ----------------------------------------------------------
Vector3D<Precision> VUnplacedVolume::SamplePointOnSurface() const
{
  throw std::runtime_error("unimplemented function called");
  return Vector3D<Precision>();
}

// ----------------- Extent --------------------------------------------------------------------
VECCORE_ATT_HOST_DEVICE
void VUnplacedVolume::Extent(Vector3D<Precision> &aMin, Vector3D<Precision> &aMax) const
{
#ifndef VECCORE_CUDA
  throw std::runtime_error("unimplemented function called");
#endif
}

std::ostream &operator<<(std::ostream &os, VUnplacedVolume const &vol)
{
  vol.Print(os);
  return os;
}

#ifndef VECCORE_CUDA

VPlacedVolume *VUnplacedVolume::PlaceVolume(LogicalVolume const *const volume,
                                            Transformation3D const *const transformation,
                                            VPlacedVolume *const placement) const
{

  const TranslationCode trans_code = transformation->GenerateTranslationCode();
  const RotationCode rot_code      = transformation->GenerateRotationCode();

  return SpecializedVolume(volume, transformation, trans_code, rot_code, placement);
}

VPlacedVolume *VUnplacedVolume::PlaceVolume(char const *const label, LogicalVolume const *const volume,
                                            Transformation3D const *const transformation,
                                            VPlacedVolume *const placement) const
{
  VPlacedVolume *const placed = PlaceVolume(volume, transformation, placement);
  placed->set_label(label);
  return placed;
}

#endif
}
} // End global namespace
