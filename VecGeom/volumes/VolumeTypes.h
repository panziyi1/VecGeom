// This file is part of VecGeom and is distributed under the
// conditions in the file LICENSE.txt in the top directory.
// For the full list of authors see CONTRIBUTORS.txt and `git log`.

/// \brief A list of all volume types supported by VecGeom.
/// Used for dynamic dispatch without virtual function calls.
/// \file VecGeom/volumes/VolumeTypes.h
/// \author created by Jonas Hahnfeld

#ifndef VECGEOM_VOLUMES_VOLUMETYPES_H_
#define VECGEOM_VOLUMES_VOLUMETYPES_H_

namespace vecgeom {
enum VolumeTypes {
  kUnknown = 0,

  kBox        = 1,
  kCone       = 2,
  kCutTube    = 3,
  kOrb        = 4,
  kPolycone   = 5,
  kPolyhedron = 6,
  kSphere     = 7,
  kTorus2     = 8,
  kTrapezoid  = 9,
  kTrd        = 10,
  kTube       = 11,

  kBooleanIntersection = 12,
  kBooleanUnion        = 13,
  kBooleanSubtraction  = 14,
};
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_VOLUMETYPES_H_
