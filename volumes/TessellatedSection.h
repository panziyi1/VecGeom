/// \file TessellatedSection.h
/// \author Mihaela Gheata (mihaela.gheata@cern.ch)

#ifndef VECGEOM_VOLUMES_TESSELLATEDSECTION_H_
#define VECGEOM_VOLUMES_TESSELLATEDSECTION_H_

#include "volumes/TessellatedCluster.h"

#include "management/HybridManager2.h"
#include "navigation/HybridNavigator2.h"
#include "management/ABBoxManager.h"

namespace vecgeom {

inline namespace VECGEOM_IMPL_NAMESPACE {

// Navigation helper for adjacent quadrilateral facets forming a closed
// section in between 2 Z planes. A base class defines the scalar navigation
// interfaces.

template <typename T>
class TessellatedSectionBase {

  using Real_v = vecgeom::VectorBackend::Real_v;

public:
  VECCORE_ATT_HOST_DEVICE
  virtual Inside_t Inside(Vector3D<Real_v> const &point) const = 0;

  VECCORE_ATT_HOST_DEVICE
  virtual bool Contains(Vector3D<Real_v> const &point) const = 0;

  VECCORE_ATT_HOST_DEVICE
  virtual T DistanceToIn(Vector3D<Real_v> const &point, Vector3D<Real_v> const &direction, int &isurf) const = 0;

  VECCORE_ATT_HOST_DEVICE
  virtual T DistanceToOut(Vector3D<Real_v> const &point, Vector3D<Real_v> const &direction, int &isurf) const = 0;

  VECCORE_ATT_HOST_DEVICE
  virtual double SafetyToInSq(Vector3D<Real_v> const &point, int &isurf) const = 0;

  VECCORE_ATT_HOST_DEVICE
  virtual double SafetyToOutSq(Vector3D<Real_v> const &point, int &isurf) const = 0;

  VECCORE_ATT_HOST_DEVICE
  virtual void Normal(Vector3D<T> const &point, Vector3D<T> &normal, bool &valid) const = 0;
};

// Derived templated tessellated section on convexity
// making a flat surface at same z.

template <typename T, bool Convex = true>
class TessellatedSection /* : public TessellatedSectionBase<T> */ {
  // Here we should be able to use vecgeom::Vector
  template <typename U>
  using vector_t  = std::vector<U>;
  using Real_v    = vecgeom::VectorBackend::Real_v;
  using Facet_t   = QuadrilateralFacet<T>;
  using Cluster_t = TessellatedCluster<4, Real_v>;

private:
  int fNfacets;           ///< Number of triangle facets on the section
  T fZmin        = 0;     ///< Minimum Z
  T fZmax        = 0;     ///< Maximum Z
  T fCubicVolume = 0;     ///< Cubic volume
  T fSurfaceArea = 0;     ///< Surface area
  bool fConvex   = false; ///< Convexity
  Vector3D<T> fMinExtent; ///< Minimum extent
  Vector3D<T> fMaxExtent; ///< Maximum extent
  Vector3D<T> fTestDir;   ///< Test direction for generalized Contains/Inside

  vector_t<Vector3D<T>> fVertices; ///< Vector of unique vertices
  vector_t<Facet_t *> fFacets;     ///< Vector of quadrilateral convex facets
  vector_t<Cluster_t *> fClusters; ///< Vector of facet clusters

protected:
  VECCORE_ATT_HOST_DEVICE
  void AddFacet(Facet_t *facet)
  {
    // Method adding a facet to the structure. The vertices are added to the
    // list of all vertices (including duplications) and the extent is re-adjusted.
    fFacets.push_back(facet);
    // Adjust extent
    using vecCore::math::Min;
    using vecCore::math::Max;
    T xmin        = Min(Min(facet->fVertices[0].x(), facet->fVertices[1].x()),
                 Min(facet->fVertices[2].x(), facet->fVertices[3].x()));
    T ymin        = Min(Min(facet->fVertices[0].y(), facet->fVertices[1].y()),
                 Min(facet->fVertices[2].y(), facet->fVertices[3].y()));
    fMinExtent[0] = Min(fMinExtent[0], xmin);
    fMinExtent[1] = Min(fMinExtent[1], ymin);
    fMinExtent[2] = fZmin;
    T xmax        = Max(Max(facet->fVertices[0].x(), facet->fVertices[1].x()),
                 Max(facet->fVertices[2].x(), facet->fVertices[3].x()));
    T ymax        = Max(Min(facet->fVertices[0].y(), facet->fVertices[1].y()),
                 Max(facet->fVertices[2].y(), facet->fVertices[3].y()));
    fMaxExtent[0] = Max(fMaxExtent[0], xmax);
    fMaxExtent[1] = Max(fMaxExtent[1], ymax);
    fMaxExtent[2] = fZmax;
    // Check if we can create a Tessellated cluster
    int nfacets = fFacets.size();
    assert(nfacets <= fNfacets && "Cannot add extra facets to section");
    if (nfacets % kVecSize == 0 || nfacets == fNfacets) {
      int istart         = nfacets - (nfacets - 1) % kVecSize - 1;
      int i              = 0;
      Cluster_t *cluster = new Cluster_t();
      for (; istart < nfacets; ++istart) {
        cluster->AddFacet(i++, fFacets[istart], istart);
      }
      // The last cluster may not be yet full: fill with last facet
      for (; i < kVecSize; ++i)
        cluster->AddFacet(i, facet, nfacets - 1);
      fClusters.push_back(cluster);
    }
    if (nfacets == fNfacets) {
      if (Convex) {
        assert(CalculateConvexity() == Convex);
      }
    }
  }

  /** &brief Calculate convexity of the section with respect to itself */
  VECCORE_ATT_HOST_DEVICE
  bool CalculateConvexity()
  {
    int nconvex = 0;
    for (size_t i = 0; i < fNfacets; ++i) {
      Facet_t *facet = fFacets[i];
      bool convex    = true;
      for (size_t j = 0; j < fNfacets; ++j) {
        if (j == i) continue;
        for (size_t ivert = 0; ivert < 4; ++ivert) {
          convex &= facet->DistPlane(fFacets[j]->fVertices[ivert]) < kTolerance;
        }
        if (!convex) continue;
      }
      facet->fConvex = convex;
      if (convex) nconvex++;
    }
    for (auto cluster : fClusters)
      cluster->CalculateConvexity();
    if (nconvex == fNfacets) return true;
    return false;
  }

public:
  VECCORE_ATT_HOST_DEVICE
  TessellatedSection(int nfacets, T zmin, T zmax) : fNfacets(nfacets), fZmin(zmin), fZmax(zmax)
  {
    assert(zmax > zmin && "zmin is greater than zmax");
    double testangle = 1.e-3;
    fTestDir.Set(vecCore::math::Cos(testangle), vecCore::math::Sin(testangle), 0.);
    fMinExtent.Set(InfinityLength<T>());
    fMaxExtent.Set(-InfinityLength<T>());
  }

  VECCORE_ATT_HOST_DEVICE
  bool AddQuadrilateralFacet(Vector3D<T> const &vt0, Vector3D<T> const &vt1, Vector3D<T> const &vt2,
                             Vector3D<T> const &vt3, bool absolute = true)
  {
    // Quadrilateral facet, normal pointing outside
    Facet_t *facet = new Facet_t;
    if (absolute) {
      if (!facet->SetVertices(vt0, vt1, vt2, vt3)) {
        delete facet;
        return false;
      }
      AddFacet(facet);
    } else {
      if (!facet->SetVertices(vt0, vt0 + vt1, vt0 + vt1 + vt2, vt0 + vt1 + vt2 + vt3)) {
        delete facet;
        return false;
      }
      AddFacet(facet);
    }
    return true;
  }

  /** @brief Fast check using the extent if the point is outside */
  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  bool IsOutside(Vector3D<T> const &point)
  {
    return ((point - fMinExtent).Min() < -kTolerance || (point - fMaxExtent).Max() > kTolerance);
  }

  VECCORE_ATT_HOST_DEVICE
  Inside_t Inside(Vector3D<Real_v> const &point) const
  {
    // All lanes of point contain the same scalar point
    using Bool_v = vecCore::Mask<Real_v>;

    // Assuming the fast check on extent was already done using the scalar point
    size_t nclusters = fClusters.size();
    if (Convex) {
      // Convex polygone on top/bottom
      Real_v distPlanes;
      Bool_v inside(true), outside(false);
      for (size_t i = 0; i < nclusters; ++i) {
        distPlanes = fClusters[i]->DistPlanes(point);
        outside |= distPlanes > Real_v(kTolerance);
        //        if (!vecCore::MaskEmpty(outside)) return kOutside;
        inside &= distPlanes < -kTolerance;
      }
      if (!vecCore::MaskEmpty(outside)) return kOutside;
      if (vecCore::MaskFull(inside)) return kInside;
      return kSurface;
    }

    // The general case.
    T distanceToIn  = InfinityLength<T>();
    T distanceToOut = InfinityLength<T>();
    int isurfToIn   = -1;
    int isurfToOut  = -1;

    T clusterToIn, clusterToOut;
    int icrtToIn, icrtToOut;
    for (size_t i = 0; i < nclusters; ++i) {
      fClusters[i]->DistanceToCluster(point, fTestDir, clusterToIn, clusterToOut, icrtToIn, icrtToOut);

      // Update distanceToIn/Out
      if (icrtToIn >= 0 && clusterToIn < distanceToIn) {
        distanceToIn = clusterToIn;
        isurfToIn    = icrtToIn;
      }

      if (icrtToOut >= 0 && clusterToOut < distanceToOut) {
        distanceToOut = clusterToOut;
        isurfToOut    = icrtToOut;
      }
    }
    if (isurfToOut < 0) return kOutside;
    if (isurfToIn >= 0 && distanceToIn < distanceToOut &&
        distanceToIn * fTestDir.Dot(fFacets[isurfToIn]->fNormal) < -kTolerance)
      return kOutside;

    if (distanceToOut < 0 || distanceToOut * fTestDir.Dot(fFacets[isurfToOut]->fNormal) < kTolerance) return kSurface;

    if (isurfToIn < 0 || distanceToOut < distanceToIn) return kInside;

    if (distanceToIn < 0 || distanceToIn * fTestDir.Dot(fFacets[isurfToIn]->fNormal) > -kTolerance) return kSurface;
  }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  size_t GetNfacets() const { return fFacets.size(); }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  size_t GetNclusters() const { return fClusters.size(); }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Cluster_t const &GetCluster(size_t i) const { return *fClusters[i]; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  Facet_t const &GetFacet(size_t i) const { return *fFacets[i]; }

  VECCORE_ATT_HOST_DEVICE
  bool Contains(Vector3D<Real_v> const &point) const
  {
    using Bool_v = vecCore::Mask<Real_v>;

    // Check if point is in the bounding box
    // if ((point - fMinExtent).Min() < 0. || (point - fMaxExtent).Max() > 0.) return kOutside;

    size_t nclusters = fClusters.size();
    if (Convex) {
      // Convex polygone on top/bottom
      Bool_v outside(false);
      for (size_t i = 0; i < nclusters; ++i) {
        Real_v distPlanes = fClusters[i]->DistPlanes(point);
        outside |= distPlanes > Real_v(0);
        if (!vecCore::MaskEmpty(outside)) return kOutside;
      }
      return kInside;
    }

    // The general case.
    T distanceToIn  = InfinityLength<T>();
    T distanceToOut = InfinityLength<T>();
    int isurfToIn   = -1;
    int isurfToOut  = -1;

    T clusterToIn, clusterToOut;
    int icrtToIn, icrtToOut;
    for (size_t i = 0; i < nclusters; ++i) {
      fClusters[i]->DistanceToCluster(point, fTestDir, clusterToIn, clusterToOut, icrtToIn, icrtToOut);

      // Update distanceToIn/Out
      if (icrtToIn >= 0 && clusterToIn < distanceToIn) {
        distanceToIn = clusterToIn;
        isurfToIn    = icrtToIn;
      }

      if (icrtToOut >= 0 && clusterToOut < distanceToOut) {
        distanceToOut = clusterToOut;
        isurfToOut    = icrtToOut;
      }
    }
    if (isurfToOut < 0) return kOutside;
    if (isurfToIn >= 0 && distanceToIn < distanceToOut && distanceToIn * fTestDir.Dot(fFacets[isurfToIn]->fNormal) < 0)
      return kOutside;

    return kInside;
  }

  VECCORE_ATT_HOST_DEVICE
  T DistanceToInConvex(Vector3D<T> const &point, Vector3D<T> const &direction, T stepmax) const
  {
    // Compute distance to segment from outside point.
    T dz = 0.5 * (fZmax - fZmin);
    T pz = point.z() - 0.5 * (fZmax + fZmin);
    if ((vecCore::math::Abs(pz) - dz) > -kTolerance && pz * direction.z() >= 0) return InfinityLength<T>();
    const T invz = -1. / NonZero(direction.z());
    const T ddz  = (invz < 0) ? dz : -dz;
    T distance   = (pz + ddz) * invz;
    T limit      = vecCore::math::Min((pz - ddz) * invz, stepmax);

    const int nclusters = fClusters.size();
    for (int i = 0; i < nclusters; ++i) {
      // bool hitcluster =
      fClusters[i]->DistanceToInConvex(point, direction, distance, limit);
    }
    return distance;
  }

  VECCORE_ATT_HOST_DEVICE
  T DistanceToOut(Vector3D<Real_v> const &point, Vector3D<Real_v> const &direction, int &isurf) const
  {
    // Compute distance to segment from point inside, returning also the crossed
    // facet.
    isurf               = -1;
    T distance          = InfinityLength<T>();
    T stepmax           = InfinityLength<T>();
    const int nclusters = fClusters.size();
    for (int i = 0; i < nclusters; ++i) {
      int isurfcrt = -1;
      T distcrt;
      fClusters[i]->DistanceToOut(point, direction, stepmax, distcrt, isurfcrt);
      if (distcrt < distance) {
        distance = distcrt;
        isurf    = isurfcrt;
      }
    }
    return distance;
  }

  VECCORE_ATT_HOST_DEVICE
  T SafetyToInSq(Vector3D<Real_v> const &point, int &isurf) const
  {
    // Compute safety squared to segment from point outside, returning also the crossed
    // facet.
    isurf               = -1;
    T safetysq          = InfinityLength<T>();
    const int nclusters = fClusters.size();
    for (int i = 0; i < nclusters; ++i) {
      int isurfcrt     = -1;
      const T safsqcrt = fClusters[i]->template SafetySq<true>(point, isurfcrt);
      if (safsqcrt < safetysq) {
        safetysq = safsqcrt;
        isurf    = isurfcrt;
      }
    }
    return safetysq;
  }

  VECCORE_ATT_HOST_DEVICE
  T SafetyToOutSq(Vector3D<Real_v> const &point, int &isurf) const
  {
    // Compute safety squared to segment from point inside, returning also the crossed
    // facet.
    isurf               = -1;
    T safetysq          = InfinityLength<T>();
    const int nclusters = fClusters.size();
    for (int i = 0; i < nclusters; ++i) {
      int isurfcrt     = -1;
      const T safsqcrt = fClusters[i]->template SafetySq<false>(point, isurfcrt);
      if (safsqcrt < safetysq) {
        safetysq = safsqcrt;
        isurf    = isurfcrt;
      }
    }
    return safetysq;
  }

  VECCORE_ATT_HOST_DEVICE
  void Normal(Vector3D<T> const & /*point*/, Vector3D<T> & /*normal*/, bool & /*valid*/) const
  {
    // Compute normal to segment surface in given point near surface.
  }
};

std::ostream &operator<<(std::ostream &os, TessellatedSection<double, true> const &ts);

} // namespace VECGEOM_IMPL_NAMESPACE
} // end namespace vecgeom

#endif
