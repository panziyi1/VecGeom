/// \file TessellatedStruct.h
/// \author Mihaela Gheata (mihaela.gheata@cern.ch)

#ifndef VECGEOM_VOLUMES_TESSELLATEDSTRUCT_H_
#define VECGEOM_VOLUMES_TESSELLATEDSTRUCT_H_

#include <VecCore/VecCore>

#include "base/Global.h"
#include "base/Vector3D.h"

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

constexpr size_t kVecSize = vecCore::VectorSize<vecgeom::VectorBackend::Real_v>();

// Basic structure of indices to 3 vertices making a triangle.
// The vertices making the triangle have to be given in anti-clockwise
// order looking from the outsider of the solid where it belongs.
// This helper structure is used temporarily by the tesselated solid in the
// creation and clustering phase.

template <typename T = double>
struct TriangleFacet {
  Vector3D<T> fVertices[3];      ///< vertices of the triangle
  Vector3D<int> fIndices;        ///< indices for 3 distinct vertices
  Vector<int> fNeighbors;        ///< indices to triangle neighbors
  T fSurfaceArea = 0;            ///< surface area
  Vector3D<T> fNormal;           ///< normal vector pointing outside
  T fDistance;                   ///< distance between the origin and the triangle plane

  VECGEOM_CUDA_HEADER_BOTH
  TriangleFacet()
  {
    fVertices.reserve(3);
    fNeighbors.reserve(3);
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  bool SetVertices(Vector3D<T> const &vtx0, Vector3D<T> const &vtx1, Vector3D<T> const &vtx2, int ind0, int ind1,
                   int ind2)
  {
    fVertices.push_back(vtx0);
    fVertices.push_back(vtx1);
    fVertices.push_back(vtx2);
    fInsices.Set(ind0, ind1, ind2);
    // Check validity
    Vector3D<T> e1 = vtx1 - vtx0;
    Vector3D<T> e2 = vtx2 - vtx0;
    double eMag1   = e1.Mag();
    double eMag2   = e2.Mag();
    double eMag3   = (e2 - e1).Mag();

    if (eMag1 <= kTolerance || eMag2 <= kTolerance || eMag3 <= kTolerance) {
      // TO DO: add more verbosity
      std::cout << "Length of sides of facet are too small." << std::endl;
      return false;
    }
    // Compute normal
    fNormal      = e1.Cross(e2).Unit();
    fSurfaceArea = 0.5 * (e1.Cross(e2)).Mag();
    if (fSurfaceArea < kTolerance * kTolerance) {
      // TO DO: add more verbosity
      std::cout << "Flat triangle." << std::endl;
      return false;
    }
    // Any more fields to fill to be added here
  }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  bool IsSurrounded() const { return fNeighbors.size() > 2; }

  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_FORCE_INLINE
  bool IsNeighbor(TriangleFacet const &other)
  {
    // Check if a segment is common
    int ncommon = 0;
    int i1      = 0;
    int i2      = 0;
    for (int ind1 = 0; ind1 < 3; ++ind1) {
      for (int ind2 = 0; ind2 < 3; ++ind2) {
        if ((fIndices[ind1] == other.fIndices[ind2]) {
          ncommon++;
          i1 = ind1;
          i2 = ind2;
        }
      }
    }
    // In case we detect a single common vertex, it is still possible that the
    // facets are neighbors
    // if (ncommon == 1) DetectCollision(other)
    return (ncommon == 2);
  }
};

// Structure used for vectorizing queries on groups of triangles
template <typename Real_v>
struct TessellatedCluster {
  using T = vecCore::ScalarType<Real_v>::Type;

  //____________________________________________________________________________
  Vector3D<Real_v> fNormals;        ///< Normals to facet components
  Real_v fDistances;                ///< Distances from origin to facets
  Vector3D<Real_v> fSideVectors[3]; ///< Side vectors of the triangular facets
  Real_v fVertices[3];              ///< Vertices stored in SIMD format
  
  Vector3D<T> fMinExtent;           ///< Minimum extent
  Vector3D<T> fMaxExtent;           ///< Maximum extent
  //____________________________________________________________________________

  VECGEOM_CUDA_HEADER_BOTH
  TesselatedCluster() {}

  /** @brief Fill the components 'i' of the cluster with facet data
    * @param index Triangle index, equivalent to SIMD lane index
    * @param facet Triangle facet data
    */
  VECGEOM_CUDA_HEADER_BOTH
  AddFacet(size_t index, TriangleFacet<T> const &facet)
  {
    // Fill the facet normal by accessing individual SIMD lanes
    assert(index <= kVecSize);
    fNormals.x()[index] = facet.fNormal.x();
    // Fill the distance to the plane
    fDistances[index] = facet.fDistance;
    // Compute side vectors and fill them using the store operation per SIMD lane
    for (size_t ivert = 0; ivert < 3; ++ivert) {
      Vector3D<T> c0 = facet.fVertices[ivert];
      Vector3D<T> c1 = facet.fVertices[(ivert+1) % 3];
      Vector3D<T> sideVector = facet.fNormal.Cross(c1 - c0).Normalized();
      fSideVectors[ivert].x()[index] = sideVector.x();
      fSideVectors[ivert].y()[index] = sideVector.y();
      fSideVectors[ivert].z()[index] = sideVector.z();
      fVertices[ivert].x()[index] = c0.x();
      fVertices[ivert].y()[index] = c0.y();
      fVertices[ivert].z()[index] = c0.z();
    }
  }

  // === Navigation functionality === //
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  void InsideCluster(Vector3D<Real_v> const &point,
                     typename vecCore::Mask<Real_v> inside) const
  {
    // Check if the points are inside some of the triangles. The points are assumed
    // to be already propagated on the triangle planes.
    using Bool_v = vecCore::Mask<Real_v>;
    
    inside = Bool_v(true);
    for (size_t i = 0; i < 3; ++i) {
      Real_v saf = (point - fVertices[i]).Dot(fSideVectors[i]);
      inside &= saf < Real_v(kTolerance);
    }
  }
  
  VECGEOM_FORCE_INLINE
  VECGEOM_CUDA_HEADER_BOTH
  Real_v DistPlanes(Vector3D<Real_v> const &point) const
  {
    // Returns distance from point to plane. This is positive if the point is on
    // the outside halfspace, negative otherwise.
    return (point.Dot(fNormals) + fDistances);
  }

  VECGEOM_CUDA_HEADER_BOTH
  T DistanceToIn(Vector3D<T> const &point, Vector3D<T> const &direction,
                    T const &/*stepMax*/, T &distance, int &isurf)
  {
    using Bool_v = vecCore::Mask<Real_v>;

    Real_v distance = InfinityLength<Real_v>();
    isurf = -1;
    Real_v pointv(point);
    Real_v dirv(direction);
    Real_v ndd   = NonZero(dirv.Dot(fNormals));
    Real_v saf = DistPlanes(pointv);
    Bool_v valid = ndd < Real_v(0.) && saf > Real_v(-kTolerance);
    if ( vecCore::EarlyReturnAllowed() && vecCore::MaskEmpty(valid))
      return InfinityLength<T>();

    vecCore__MaskedAssignFunc(distance, valid, -saf / ndd);
    // Since we can make no assumptions on convexity, we need to actually check
    // which surface is actually crossed. First propagate the point with the
    // distance to each plane.
    pointv += distance * dirv;
    // Check if propagated points hit the triangles
    Bool_v hit;
    InsideCluster(pointv, hit);
    valid &= hit;
    // Now we need to return the minimum distance for the hit facets
    if ( vecCore::MaskEmpty(valid) ) return;
    T distmin = InfinityLength<T>();
    for (int i = 0;  i < kVecSize; ++i) {
      if (valid[i] && distance[i] < distmin) {
        distmin = distance[i];
        isurf = i;
      }
    }
    return distmin;
  }
};

template <typename T = double>
class TessellatedStruct {

private:
  bool fSolidClosed;      ///< Closure of the solid
  int fNsurrounded = 0;   ///< Number of surrounded facets
  T fCubicVolume   = 0;   ///< cubic volume
  Vector3D<T> fMinExtent; ///< Minimum extent
  Vector3D<T> fMaxExtent; ///< Maximum extent

  // Here we have a pointer to the aligned bbox structure
  // ABBoxanager *fABBoxManager;

  Vector<Vector3D<T>> fVertices;                ///< Vector of unique vertices
  Vector<TriangleFacet> fFacets;                ///< Vector of triangular facets
  Vector<TessellatedCluster<T, kNV>> fClusters; ///< Vector of facet clusters

protected:
  VECGEOM_CUDA_HEADER_BOTH
  int AddVertex(Vector3D<T> const &vtx)
  {
    // This method needs to check if the vertex is duplicated. The
    // index of the vertex is returned.

    // The trivial loop makes AddVertex a N squared problem...
    // *** TO DO: uniform grid store/search, ...
    constexpr Precision tolerancesq = kTolerance * kTolerance;
    int ivert                       = 0;
    for (auto vertex : fVertices) {
      if ((vtx - vertex).Mag2() < tolerancesq) return ivert;
      ivert++;
    }
    fVertices.push_back(vtx);
    return ivert;
  }

  VECGEOM_CUDA_HEADER_BOTH
  void FindNeighbors(TriangleFacet &facet, int ifacet)
  {
    // Loop non-closed facets.
    int nfacets = fFacets.size();
    for (int icrt = 0; icrt < nfacets; ++icrt) {
      if (fFacets[icrt].IsSurrounded()) continue;
      if (facet.IsNeighbor(fFacets[icrt])) {
        facet.fNeighbors.push_back(icrt);
        fFacets[icrt].fNeighbors.push_back(ifacet);
        if (fFacets[icrt].IsSurrounded()) fNsurrounded++;
        if (facet.IsSurrounded()) {
          fNsurrounded++ return;
        }
      }
    }
  }

public:
  VECGEOM_CUDA_HEADER_BOTH
  TesselatedStruct() { fClusterSize = vecCore::VectorSize<vecgeom::VectorBackend::Real_v>(); }

  /* @brief Methods for adding a new facet
   * @detailed The method akes 4 parameters to define the three fVertices:
   *      1) UFacetvertexType = "ABSOLUTE": in this case Pt0, vt1 and vt2 are
   *         the 3 fVertices in anti-clockwise order looking from the outsider.
   *      2) UFacetvertexType = "RELATIVE": in this case the first vertex is Pt0,
   *         the second vertex is Pt0+vt1 and the third vertex is Pt0+vt2, all
   *         in anti-clockwise order when looking from the outsider.
   */
  VECGEOM_CUDA_HEADER_BOTH
  bool AddTriangularFacet(Vector3D<T> const &vt0, Vector3D<T> const &vt1, Vector3D<T> const &vt2, bool absolute = true)
  {
    TriangleFacet facet;
    int ind0, ind1, ind2;
    ind0 = AddVertex(vt0);
    if (absolute) {
      ind1 = AddVertex(vt1);
      ind2 = AddVertex(vt2);
    } else {
      ind1 = AddVertex(vt0 + vt1);
      ind2 = AddVertex(vt0 + vt1 + vt2);
    }
    if (!facet.AddVertices(fVertices[ind0], fVertices[ind1], fVertices[ind2], ind0, ind1, ind2)) return false;
    FindNeighbors(facet);
    fFacets.push_back(facet);
    if (fNsurrounded == fFacets.size()) fClosed = true;
  }

  VECGEOM_CUDA_HEADER_BOTH
  bool AddQuadrilateralFacet(Vector3D<T> const &vt0, Vector3D<T> const &vt1, Vector3D<T> const &vt2,
                             Vector3D<T> const &vt3, bool absolute = true)
  {
    // We should check the quadrilateral convexity to correctly define the
    // triangle facets
    // CheckConvexity()vt0, vt1, vt2, vt3, absolute);
    TriangleFacet facet1, facet2;
    int ind0, ind1, ind2, ind22;
    ind0 = AddVertex(vt0);
    if (absolute) {
      ind1  = AddVertex(vt1);
      ind2  = AddVertex(vt2);
      ind22 = AddVertex(vt3);
    } else {
      ind1  = AddVertex(vt0 + vt1);
      ind2  = AddVertex(vt0 + vt1 + vt2);
      ind22 = AddVertex(vt0 + vt1 + vt2 + vt3);
    }
    bool added1 = facet1.AddVertices(fVertices[ind0], fVertices[ind1], fVertices[ind2], ind0, ind1, ind2);
    bool added2 = facet2.AddVertices(fVertices[ind0], fVertices[ind2], fVertices[ind22], ind0, ind2, ind22);

    // The two facets are neighbors
    int ifacet1 = fFacets.size();
    if (added1) facet2.fNeighbors.push_back(ifacet1);
    if (added2) facet1.fNeighbors.push_back(ifacet1 + 1);
    if (added1) {
      FindNeighbors(facet1);
      fFacets.push_back(facet1);
    }
    if (added2) {
      FindNeighbors(facet2);
      fFacets.push_back(facet2);
    }
    if (fNsurrounded == fFacets.size()) fClosed = true;
    return (added1 || added2);
  }
};
}
} // end

#endif