/// \file QuadrilateralFacet.h
/// \author Mihaela Gheata (mihaela.gheata@cern.ch)

#ifndef VECGEOM_VOLUMES_QUADRILATERALCFACET_H_
#define VECGEOM_VOLUMES_QUADRILATERALCFACET_H_

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE(struct, QuadrilateralCFacet, typename);

inline namespace VECGEOM_IMPL_NAMESPACE {
//______________________________________________________________________________
// Basic structure of indices to 4 vertices making a convex quadrilateral.
// The vertices making the quadrilateral have to be given in anti-clockwise
// order looking from the outsider of the solid where it belongs.
//______________________________________________________________________________
template <typename T = double>
struct QuadrilateralCFacet {
  Vector3D<T> fVertices[4]; ///< vertices of the triangle
  Vector3D<T> fCenter;      ///< Center of the triangle
  int fIndices[4];          ///< indices for 3 distinct vertices
  int fNeighbors[4];        ///< indices to triangle neighbors
  T fSurfaceArea = 0;       ///< surface area
  Vector3D<T> fNormal;      ///< normal vector pointing outside
  T fDistance;              ///< distance between the origin and the triangle plane
  Vector3D<T> fSideVectors[4]; ///< side vectors perpendicular to edges

  VECCORE_ATT_HOST_DEVICE
  QuadrilateralCFacet() { fNeighbors.reserve(4); }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  bool SetVertices(Vector3D<T> const &vtx0, Vector3D<T> const &vtx1, Vector3D<T> const &vtx2, Vector3D<T> const &vtx3,
                   int ind0 = 0, int ind1 = 0, int ind2 = 0, int ind3 = 0)
  {
    fVertices[0] = vtx0;
    fVertices[1] = vtx1;
    fVertices[2] = vtx2;
    fVertices[3] = vtx3;
    fIndices[0]  = ind0;
    fIndices[1]  = ind1;
    fIndices[2]  = ind2;
    fIndices[3]  = ind3;
    // Check validity
    // We accept degeneration of 2 vertices only (triangle)
    int ismall = -1;
    int nsmall = 0;
    for (int i=0; i<4; ++i) {
      int j = (i+1)%4;
      if ((fVertices[j] - fVertices[i]).Mag() < kTolerance) {
        ismall = i;
        nsmall++;
      }
    }
    if (nsmall > 1) {
      // TO DO: add more verbosity
      std::cout << "Quadrilateral degenerated: Length of sides of facet are too small." << std::endl;
      return false;
    }
    // Compute normal using non-zero segments
    Vector3D<T> e1 = fVertices[(ismall+2)%4] - fVertices[(ismall+1)%4];
    Vector3D<T> e2 = fVertices[(ismall+3)%4] - fVertices[(ismall+2)%4];
  
    fNormal = e1.Cross(e2).Unit();
    // Compute side vectors
    for (int i        = 0; i < 4; i++)
      fSideVectors[i] = fNormal.Cross(fVertices[(i + 1) % 4] - fVertices[i]).Normalized();
    if (nsmall > 0)
      fSideVectors[ismall] = fSideVectors[(ismall+1)%4];
    // Distace to facet
    fDistance    = -fNormal.Dot(vtx0);
    fSurfaceArea = 0.;
    for (int i=0; i<3; ++i) {
      Vector3D<T> e1 = fVertices[(i+1)%4] - fVertices[i];
      Vector3D<T> e2 = fVertices[(i+2)%4] - fVertices[(i+1)%4];
      fSurfaceArea += 0.5 * (e1.Cross(e2)).Mag();

    }
    if (fSurfaceArea < kTolerance * kTolerance) {
      // TO DO: add more verbosity
      std::cout << "Flat triangle." << std::endl;
      return false;
    }
    // Center of the quadrilateral
    fCenter = (vtx0 + vtx1 + vtx2 + vtx3) / 4.;
    return true;
  }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  bool IsSurrounded() const { return fNeighbors.size() > 3; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  int IsNeighbor(QuadrilateralCFacet const &other)
  {
    // Check if a segment is common
    int ncommon = 0;
    for (int ind1 = 0; ind1 < 4; ++ind1) {
      for (int ind2 = 0; ind2 < 4; ++ind2) {
        if (fIndices[ind1] == other.fIndices[ind2]) ncommon++;
      }
    }
    return ncommon;
  }

  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  bool Contains(Vector3D<T> const &point) const
  {
    // Check id point within the triangle plane is inside the triangle.
    bool inside = true;
    for (size_t i = 0; i < 4; ++i) {
      T saf = (point - fVertices[i]).Dot(fSideVectors[i]);
      inside &= saf > -kTolerance;
    }
    return inside;
  }

  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  T DistPlane(Vector3D<T> const &point) const
  {
    // Returns distance from point to plane. This is positive if the point is on
    // the outside halfspace, negative otherwise.
    return (point.Dot(fNormal) + fDistance);
  }

  VECCORE_ATT_HOST_DEVICE
  T DistanceToIn(Vector3D<T> const &point, Vector3D<T> const &direction) const
  {
    T ndd      = NonZero(direction.Dot(fNormal));
    T saf      = DistPlane(point);
    bool valid = ndd < 0. && saf > -kTolerance;
    if (!valid) return InfinityLength<T>();
    T distance = -saf / ndd;
    // Propagate the point with the distance to the plane.
    Vector3D<T> point_prop = point + distance * direction;
    // Check if propagated points hit the triangle
    if (!Contains(point_prop)) return InfinityLength<T>();
    return distance;
  }

  VECCORE_ATT_HOST_DEVICE
  Precision DistanceToOut(Vector3D<T> const &point, Vector3D<T> const &direction) const
  {
    T ndd      = NonZero(direction.Dot(fNormal));
    T saf      = DistPlane(point);
    bool valid = ndd > 0. && saf < kTolerance;
    if (!valid) return InfinityLength<T>();
    T distance = -saf / ndd;
    // Propagate the point with the distance to the plane.
    Vector3D<T> point_prop = point + distance * direction;
    // Check if propagated points hit the triangle
    if (!Contains(point_prop)) return InfinityLength<T>();
    return distance;
  }

  template <bool ToIn>
  VECCORE_ATT_HOST_DEVICE
  T SafetySq(Vector3D<T> const &point) const
  {
    T safety = DistPlane(point);
    // Find the projection of the point on each plane
    Vector3D<T> intersection = point - safety * fNormal;
    bool withinBound         = Contains(intersection);
    if (ToIn)
      withinBound &= safety > -kTolerance;
    else
      withinBound &= safety < kTolerance;
    safety *= safety;
    if (withinBound) return safety;

    Vector3D<T> safety_outbound = InfinityLength<T>();
    for (int ivert = 0; ivert < 4; ++ivert) {
      safety_outbound[ivert] =
          DistanceToLineSegmentSquared<kScalar>(fVertices[ivert], fVertices[(ivert + 1) % 4], point);
    }
    return (safety_outbound.Min());
  }
};

std::ostream &operator<<(std::ostream &os, QuadrilateralCFacet<double> const &facet);

} // end VECGEOM_IMPL_NAMESPACE
} // end namespace vecgeom

#endif
