/// \file QuadrilateralFacet.h
/// \author Mihaela Gheata (mihaela.gheata@cern.ch)

#ifndef VECGEOM_VOLUMES_TILE_H_
#define VECGEOM_VOLUMES_TILE_H_

namespace vecgeom {

VECGEOM_DEVICE_DECLARE_CONV_TEMPLATE_2t(struct, Tile, size_t, typename);

inline namespace VECGEOM_IMPL_NAMESPACE {
//______________________________________________________________________________
// Basic facet tile structure having NVERT vertices making a convex polygon.
// The vertices making the tile have to be given in anti-clockwise
// order looking from the outsider of the solid where it belongs.
//______________________________________________________________________________
template <size_t NVERT, typename T = double>
struct Tile {
  int fNvert = 0;                  ///< the tile is fully defined after adding the last vertex
  Vector3D<T> fVertices[NVERT];    ///< vertices of the tile
  Vector3D<T> fCenter;             ///< Center of the tile
  int fIndices[NVERT];             ///< indices for 3 distinct vertices
  int fNeighbors[NVERT];           ///< indices to triangle neighbors
  T fSurfaceArea = 0;              ///< surface area
  Vector3D<T> fNormal;             ///< normal vector pointing outside
  T fDistance;                     ///< distance between the origin and the triangle plane
  Vector3D<T> fSideVectors[NVERT]; ///< side vectors perpendicular to edges

  VECCORE_ATT_HOST_DEVICE
  Tile() { fNeighbors.reserve(NVERT); }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  bool AddVertex(Vector3D<T> const &vtx, int ind)
  {
    fVertices[fNvert] = vtx;
    fIndices[fNvert] = ind;
    fNvert++;
    if ( fNvert < NVERT ) return true;
    // Check validity
    // Get number of different vertices
    size_t nvert = NVERT;
    for (size_t i = 0; i < NVERT; ++i) {
      const Vector3D<T> vi = fVertices[(i+1)%NVERT] - fVertices[i];
      if ((vi.Mag2() < kTolerance) {
        nvert--;
      }
    }
    
    if (nvert < 3) {
      std::cout << "Tile degenerated: Length of sides of facet are too small." << std::endl;
      return false;
    }

    // Compute normal using non-zero segments

    bool degenerated = true;
    for (size_t i = 0; i < NVERT - 1; ++i) {
      Vector3D<T> e1 = fVertices[i+1] - fVertices[i];
      if (e1.Mag2() < kTolerance) continue;
      for (size_t j = i + 1; j < NVERT, ++j) {
        Vector3D<T> e2 = fVertices[(j + 1) % NVERT] - fVertices[j];
        if (e2.Mag2() < kTolerance) continue;
        fNormal = e1.Cross(e2);
        // e1 and e2 may be colinear
        if (fNormal.Mag2() < kTolerance) continue;
        fNormal.Normalize();
        degenerated = false;
      }
    }

    if (degenerated) {
      std::cout << "Tile degenerated: Length of sides of facet are too small." << std::endl;
      return false;
    }

    // Compute side vectors
    for (size_t i = 0; i < NVERT; ++i) {
      Vector3D<T> e1 = fVertices[(i+1) % NVERT] - fVertices[i];
      if (e1.Mag2() < kTolerance) continue;
      fSideVectors[i] = fNormal.Cross(e1).Normalized();
      fDistance    = -fNormal.Dot(fVertices[i]);
      for (size_t j = i + 1; j < i + NVERT, ++j) {
        Vector3D<T> e2 = fVertices[(j+1) % NVERT] - fVertices[j % NVERT];
        if (e2.Mag2() < kTolerance) 
          fSideVectors[j % NVERT] = fSideVectors[(j-1) % NVERT];
        else
          fSideVectors[j % NVERT] = fNormal.Cross(e2).Normalized();
      }
      break;
    }

    // Compute surface area
    fSurfaceArea = 0.;
    for (int i = 0; i < NVERT; ++i) {
      Vector3D<T> e1 = fVertices[(i+1)%4] - fVertices[i];
      Vector3D<T> e2 = fVertices[(i+2)%4] - fVertices[(i+1)%4];
      fSurfaceArea += 0.5 * (e1.Cross(e2)).Mag();
    }
    assert(fSurfaceArea < kTolerance * kTolerance);

    // Center of the tile
    for (int i = 0; i < NVERT; ++i)
      fCenter += fVertices[i];
    fCenter /= NVERT;
    return true;    
  }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  bool IsSurrounded() const { return fNeighbors.size() >= NVERT; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  int IsNeighbor(Tile<NVERT, T> const &other)
  {
    // Check if a segment is common
    int ncommon = 0;
    for (int ind1 = 0; ind1 < NVERT; ++ind1) {
      for (int ind2 = 0; ind2 < NVERT; ++ind2) {
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
    for (size_t i = 0; i < NVERT; ++i) {
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
    for (int ivert = 0; ivert < NVERT; ++ivert) {
      safety_outbound[ivert] =
          DistanceToLineSegmentSquared<kScalar>(fVertices[ivert], fVertices[(ivert + 1) % NVERT], point);
    }
    return (safety_outbound.Min());
  }
};

std::ostream &operator<<(std::ostream &os, Tile<NVERT, double> const &facet);

} // end VECGEOM_IMPL_NAMESPACE
} // end namespace vecgeom

#endif // VECGEOM_VOLUMES_TILE_H_
