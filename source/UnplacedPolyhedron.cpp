/// \file UnplacedPolyhedron.cpp
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#include "volumes/UnplacedPolyhedron.h"

#include "volumes/PlacedPolyhedron.h"
#include "volumes/SpecializedPolyhedron.h"

#include <cmath>

namespace VECGEOM_NAMESPACE {

#ifdef VECGEOM_STD_CXX11
UnplacedPolyhedron::UnplacedPolyhedron(
    const int sideCount,
    const int zPlaneCount,
    Precision zPlanes[],
    Precision rMin[],
    Precision rMax[])
    : fSideCount(sideCount), fHasInnerRadii(false),
      fZBounds{zPlanes[0]-kTolerance, zPlanes[zPlaneCount-1]+kTolerance},
      fEndCaps(2), fEndCapsOuterRadii{0, 0}, fSegments(zPlaneCount-1),
      fZPlanes(zPlaneCount), fRMin(zPlaneCount), fRMax(zPlaneCount),
      fPhiSections(sideCount+1), fBoundingTube(0, 1, 1, 0, 360) {

  typedef Vector3D<Precision> Vec_t;

  // Sanity check of input parameters
  Assert(zPlaneCount > 1, "Need at least two z-planes to construct polyhedron"
         " segments.\n");
  Assert(fSideCount > 2, "Need at least three sides to construct polyhedron"
         " segments.\n");

  copy(zPlanes, zPlanes+zPlaneCount, &fZPlanes[0]);
  copy(rMin, rMin+zPlaneCount, &fRMin[0]);
  copy(rMax, rMax+zPlaneCount, &fRMax[0]);
  // Initialize segments
  for (int i = 0; i < zPlaneCount-1; ++i) {
    fSegments[i].outer = Quadrilaterals(fSideCount);
    Assert(zPlanes[i] <= zPlanes[i+1], "Polyhedron Z-planes must be "
           "monotonically increasing.\n");
    if (rMin[i] > kTolerance || rMin[i+1] > kTolerance) {
      fHasInnerRadii = true;
      fSegments[i].hasInnerRadius = true;
      fSegments[i].inner = Quadrilaterals(fSideCount);
    } else {
      fSegments[i].hasInnerRadius = false;
    }
  }

  fEndCaps.Set(0, Vector3D<Precision>(0, 0, -1),
                  Vector3D<Precision>(0, 0, zPlanes[0]));
  fEndCaps.Set(1, Vector3D<Precision>(0, 0, 1),
                  Vector3D<Precision>(0, 0, zPlanes[zPlaneCount-1]));

  // Compute the cylindrical coordinate phi along which the corners are placed
  Precision deltaPhi = kTwoPi / sideCount;
  Precision *vertixPhi = new Precision[sideCount];
  for (int i = 0; i < sideCount; ++i) {
    vertixPhi[i] = i*deltaPhi;
    fPhiSections.set(i,
        Vec_t::FromCylindrical(1., vertixPhi[i], 0).Normalized().FixZeroes());
  }
  fPhiSections.set(sideCount, fPhiSections[0]); // Close the loop

  // Specified radii are to the sides, not to the corners. Change these values,
  // as corners and not sides are used to build the structure
  Precision cosHalfDeltaPhi = cos(0.5*deltaPhi);
  Precision innerRadius = kInfinity, outerRadius = -kInfinity;
  for (int i = 0; i < zPlaneCount; ++i) {
    // Use distance to side for minimizing inner radius
    if (rMin[i] < innerRadius) innerRadius = rMin[i];
    rMin[i] /= cosHalfDeltaPhi;
    rMax[i] /= cosHalfDeltaPhi;
    Assert(rMin[i] >= 0 && rMax[i] > kTolerance, "Invalid radius provided to "
           "unplaced polyhedron constructor.");
    // Use distance to corner for minimizing outer radius
    if (rMax[i] > outerRadius) outerRadius = rMax[i];
  }
  fEndCapsOuterRadii[0] = rMax[0];
  fEndCapsOuterRadii[1] = rMax[zPlaneCount-1];
  // Create bounding tube with biggest outer radius and smallest inner radius
  Precision boundingTubeZ = zPlanes[zPlaneCount-1] - zPlanes[0] + 2.*kTolerance;
  fBoundingTube = UnplacedTube(innerRadius - kTolerance,
                               outerRadius + kTolerance, 0.5*boundingTubeZ,
                               0, kTwoPi);
  fBoundingTubeOffset = Abs(zPlanes[0] + 0.5*boundingTubeZ);

  // Ease indexing into twodimensional vertix array
  auto GetVertix = [&sideCount] (int plane, int corner) {
    return plane*(sideCount+1) + corner;
  };

  // Precompute all vertices to ensure that there are no numerical cracks in the
  // surface.
  Vec_t *outerVertices = 0, *innerVertices = 0;
  outerVertices = new Vec_t[zPlaneCount*(sideCount+1)];
  if (fHasInnerRadii) innerVertices = new Vec_t[zPlaneCount*(sideCount+1)];
  for (int i = 0; i < zPlaneCount; ++i) {
    for (int j = 0; j < sideCount; ++j) {
      int index = GetVertix(i, j);
      outerVertices[index] =
          Vec_t::FromCylindrical(rMax[i], vertixPhi[j], zPlanes[i]).FixZeroes();
      if (fHasInnerRadii) {
        innerVertices[index] =
            Vec_t::FromCylindrical(
                rMin[i], vertixPhi[j], zPlanes[i]).FixZeroes();
      }
    }
    // Close the loop for easy indexing when building the sides
    outerVertices[GetVertix(i, sideCount)] = outerVertices[GetVertix(i, 0)];
    if (fHasInnerRadii) {
      innerVertices[GetVertix(i, sideCount)] = innerVertices[GetVertix(i, 0)];
    }
  }
  delete[] vertixPhi;

  // Build segments by drawing quadrilaterals between vertices
  for (int i = 0; i < zPlaneCount-1; ++i) {
    Vec_t corner0, corner1, corner2, corner3;
    // Sides of outer shell
    for (int j = 0; j < sideCount; ++j) {
      corner0 = outerVertices[GetVertix(i, j)];
      corner1 = outerVertices[GetVertix(i, j+1)];
      corner2 = outerVertices[GetVertix(i+1, j+1)];
      corner3 = outerVertices[GetVertix(i+1, j)];
      fSegments[i].outer.Set(j, corner0, corner1, corner2, corner3);
    }
    // Make sure the normals are pointing away from the Z-axis
    bool tiltingUp = rMax[i+1] - rMax[i] < 0;
    fSegments[i].outer.FixNormalSign(2, tiltingUp);
    // Sides of inner shell (if needed)
    if (fSegments[i].hasInnerRadius) {
      for (int j = 0; j < sideCount; ++j) {
        corner0 = innerVertices[GetVertix(i, j)];
        corner1 = innerVertices[GetVertix(i, j+1)];
        corner2 = innerVertices[GetVertix(i+1, j+1)];
        corner3 = innerVertices[GetVertix(i+1, j)];
        fSegments[i].inner.Set(j, corner0, corner1, corner2, corner3);
      }
      // Make sure the normals are pointing away from the Z-axis
      bool tiltingUp = rMin[i+1] - rMin[i] < 0;
      fSegments[i].inner.FixNormalSign(2, tiltingUp);
    }
  }

  delete[] outerVertices;
  if (fHasInnerRadii) delete[] innerVertices;
}
#endif

VECGEOM_CUDA_HEADER_DEVICE
VPlacedVolume* UnplacedPolyhedron::SpecializedVolume(
    LogicalVolume const *const volume,
    Transformation3D const *const transformation,
    const TranslationCode trans_code, const RotationCode rot_code,
#ifdef VECGEOM_NVCC
    const int id,
#endif
    VPlacedVolume *const placement) const {

  UnplacedPolyhedron const *unplaced =
      static_cast<UnplacedPolyhedron const *>(volume->unplaced_volume());

  bool hasInner = unplaced->HasInnerRadii();

#ifndef VECGEOM_NVCC
  #define POLYHEDRON_CREATE_SPECIALIZATION(INNER) \
  if (hasInner == INNER) { \
    if (placement) { \
      return new(placement) \
             SpecializedPolyhedron<INNER>(volume, transformation); \
    } else { \
      return new SpecializedPolyhedron<INNER>(volume, transformation); \
    } \
  }
#else
  #define POLYHEDRON_CREATE_SPECIALIZATION(INNER) \
  if (hasInner == INNER) { \
    if (placement) { \
      return new(placement) \
             SpecializedPolyhedron<INNER>(volume, transformation, id); \
    } else { \
      return new \
             SpecializedPolyhedron<INNER>(volume, transformation, id); \
    } \
  }
#endif

  POLYHEDRON_CREATE_SPECIALIZATION(true);
  POLYHEDRON_CREATE_SPECIALIZATION(false);

#ifndef VECGEOM_NVCC
  if (placement) {
    return new(placement)
           SpecializedPolyhedron<true>(volume, transformation);
  } else {
    return new SpecializedPolyhedron<true>(volume, transformation);
  }
#else
  if (placement) {
    return new(placement)
           SpecializedPolyhedron<true>(volume, transformation, id);
  } else {
    return new SpecializedPolyhedron<true>(volume, transformation, id);
  }
#endif

  #undef POLYHEDRON_CREATE_SPECIALIZATION
}

VECGEOM_CUDA_HEADER_BOTH
void UnplacedPolyhedron::Print() const {
  printf("UnplacedPolyhedron {%i sides, %i segments, %s}",
         fSideCount, fSegments.size(),
         (fHasInnerRadii) ? "has inner radii" : "no inner radii");
}

void UnplacedPolyhedron::Print(std::ostream &os) const {
  os << "UnplacedPolyhedron {" << fSideCount << " sides, " << fSegments.size()
     << " segments, "
     << ((fHasInnerRadii) ? "has inner radii" : "no inner radii") << "}";
}

} // End global namespace
