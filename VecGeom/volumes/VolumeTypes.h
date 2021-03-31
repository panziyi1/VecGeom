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

  kBoolean         = 1,
  kBox             = 2,
  kCoaxialCones    = 3,
  kCone            = 4,
  kCutTube         = 5,
  kEllipsoid       = 6,
  kEllipticalCone  = 7,
  kEllipticalTube  = 8,
  kExtruded        = 9,
  kGenericPolycone = 10,
  kGenTrap         = 11,
  kHype            = 12,
  kMultiUnion      = 13,
  kOrb             = 14,
  kParaboloid      = 15,
  kParallelepiped  = 16,
  kPolycone        = 17,
  kPolyhedron      = 18,
  kScaledShape     = 19,
  kSExtru          = 20,
  kSphere          = 21,
  kTessellated     = 22,
  kTet             = 23,
  kTorus2          = 24,
  kTrapezoid       = 25,
  kTrd             = 26,
  kTube            = 27,
};
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_VOLUMETYPES_H_
