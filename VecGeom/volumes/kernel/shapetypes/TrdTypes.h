// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// This file defines specialized types for the Trd shape
/// @file volumes/kernel/shapetypes/TrdTypes.h
/// @author Georgios Bitzes

#ifndef VECGEOM_VOLUMES_KERNEL_SHAPETYPES_TRDTYPES_H_
#define VECGEOM_VOLUMES_KERNEL_SHAPETYPES_TRDTYPES_H_

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_NS_CONV(TrdTypes, struct, UniversalTrd, UniversalTrd);

inline namespace VECGEOM_IMPL_NAMESPACE {
namespace TrdTypes {

#define DEFINE_TRD_TYPE(name)                       \
  struct name {                                     \
    VECCORE_ATT_HOST_DEVICE                         \
    static char const *toString() { return #name; } \
  }

// A Trd that includes all cases but does runtime checks
DEFINE_TRD_TYPE(UniversalTrd);

#undef DEFINE_TRD_TYPE

enum ETreatmentType { kYes = 0, kNo, kUnknown };

template <typename T>
struct HasVaryingY {
  static const ETreatmentType value = kUnknown;
};

template <typename T>
VECGEOM_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
bool checkVaryingY(const TrdStruct<Precision> &trd)
{
  if (HasVaryingY<T>::value != kUnknown)
    return HasVaryingY<T>::value == kYes;
  else
    return trd.fDY1 != trd.fDY2;
}

} // namespace TrdTypes
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_KERNEL_SHAPETYPES_TRDTYPES_H_
