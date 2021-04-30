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

  kBUnion          = 1,
  kBSubtraction    = 2,
  kBIntersection   = 3,
  kBox             = 4,
  kCoaxialCones    = 5,
  kCone            = 6,
  kCutTube         = 7,
  kEllipsoid       = 8,
  kEllipticalCone  = 9,
  kEllipticalTube  = 10,
  kExtruded        = 11,
  kGenericPolycone = 12,
  kGenTrap         = 13,
  kHype            = 14,
  kMultiUnion      = 15,
  kOrb             = 16,
  kParaboloid      = 17,
  kParallelepiped  = 18,
  kPolycone        = 19,
  kPolyhedron      = 20,
  kScaledShape     = 21,
  kSExtru          = 22,
  kSphere          = 23,
  kTessellated     = 24,
  kTet             = 25,
  kTorus2          = 26,
  kTrapezoid       = 27,
  kTrd             = 28,
  kTube            = 29,
};
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_VOLUMETYPES_H_
