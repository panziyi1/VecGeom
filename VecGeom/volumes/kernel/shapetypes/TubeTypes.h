/**
 * @file TubeTypes.h
 * @author Georgios Bitzes (georgios.bitzes@cern.ch)
 *
 * Contains all possible tube types
 **/

#ifndef VECGEOM_VOLUMES_KERNEL_SHAPETYPES_TUBETYPES_H_
#define VECGEOM_VOLUMES_KERNEL_SHAPETYPES_TUBETYPES_H_

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_NS_CONV(TubeTypes, struct, UniversalTube, UniversalTube);

inline namespace VECGEOM_IMPL_NAMESPACE {
namespace TubeTypes {

#define DEFINE_TUBE_TYPE(name)                      \
  struct name {                                     \
    VECCORE_ATT_HOST_DEVICE                         \
    static char const *toString() { return #name; } \
  }

// A tube that encompasses all cases - not specialized and will do extra checks at runtime
DEFINE_TUBE_TYPE(UniversalTube);

#undef DEFINE_TUBE_TYPE

// Mapping of tube types to certain characteristics
enum ETreatmentType { kYes = 0, kNo, kUnknown };

// asking for phi treatment
template <typename T>
struct NeedsPhiTreatment {
  static const ETreatmentType value = kYes;
};

template <>
struct NeedsPhiTreatment<UniversalTube> {
  static const ETreatmentType value = kUnknown;
};

template <typename T, typename UnplacedTube>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
bool checkPhiTreatment(const UnplacedTube &tube)
{
  if (NeedsPhiTreatment<T>::value != kUnknown)
    return NeedsPhiTreatment<T>::value == kYes;
  else
    return tube.fDphi < vecgeom::kTwoPi;
}

// asking for rmin treatment
template <typename T>
struct NeedsRminTreatment {
  static const ETreatmentType value = kUnknown; // unknown will trigger real rmin check
};

template <>
struct NeedsRminTreatment<UniversalTube> {
  static const ETreatmentType value = kUnknown;
};

template <typename T, typename UnplacedTube>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
bool checkRminTreatment(const UnplacedTube &tube)
{
  if (NeedsRminTreatment<T>::value != kUnknown)
    return NeedsRminTreatment<T>::value == kYes;
  else
    return tube.rmin() > 0;
}

// sector size
enum EAngleType { kNoAngle = 0, kSmallerThanPi, kOnePi, kBiggerThanPi, kUnknownAngle };

template <typename T>
struct SectorType {
  static const EAngleType value = kNoAngle;
};

template <>
struct SectorType<UniversalTube> {
  static const EAngleType value = kUnknownAngle;
};

} // End namespace TubeTypes
} // End VECGEOM_IMPL_NAMESPACE
} // End global namespace vecgeom

#endif // VECGEOM_VOLUMES_KERNEL_SHAPETYPES_TUBETYPES_H_
