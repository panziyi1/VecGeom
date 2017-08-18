#ifndef VECGEOM_ENABLE_CUDA

#include <VecCore/VecCore>

#include "test/benchmark/ArgParser.h"
#include "volumes/utilities/VolumeUtilities.h"
#include "base/Stopwatch.h"
#include "volumes/TessellatedStruct.h"
#include "volumes/kernel/TessellatedImplementation.h"
#include "test/core/TessellatedOrb.h"

#ifdef VECGEOM_ROOT
#include "utilities/Visualizer.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TH1.h"
#include "TGraph.h"
#include "TFile.h"
#include "volumes/Box.h"
#endif

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

/* Simple test for the scalability of creation of the tessellated structure.
   An orb is split into ngrid theta and phi regions; each cell is represented
   as a quadrilateral. The solid will contain 2*(ngrid-1)*ngrid triangle facets */

using namespace vecgeom;
using Real_v = vecgeom::VectorBackend::Real_v;

void RandomDirection(Vector3D<double> &direction)
{
  double phi    = RNG::Instance().uniform(0., 2. * kPi);
  double theta  = std::acos(1. - 2. * RNG::Instance().uniform(0, 1));
  direction.x() = std::sin(theta) * std::cos(phi);
  direction.y() = std::sin(theta) * std::sin(phi);
  direction.z() = std::cos(theta);
}

void RandomPointInBBox(Vector3D<double> &point, TessellatedStruct<double> const &tsl)
{
  Vector3D<double> rnd(RNG::Instance().uniform(0, 1), RNG::Instance().uniform(0, 1), RNG::Instance().uniform(0, 1));
  point = tsl.fMinExtent + rnd * (tsl.fMaxExtent - tsl.fMinExtent);
}

#ifdef VECGEOM_ROOT
void AddFacetToVisualizer(TriangleFacet<double> const *facet, Visualizer &visualizer)
{
  TPolyLine3D pl(3);
  pl.SetLineColor(kBlue);
  for (int i = 0; i < 3; i++)
    pl.SetNextPoint(facet->fVertices[i].x(), facet->fVertices[i].y(), facet->fVertices[i].z());
  visualizer.AddLine(pl);
}

void DrawCluster(TessellatedStruct<double> const &tsl, int icluster, Visualizer &visualizer, bool boxonly = false)
{
  // Draw only segments of the facets which are not shared within the cluster
  TPolyLine3D pl(2);
  pl.SetLineColor(kBlue);
  if (boxonly) {
    Vector3D<double> minext = tsl.fClusters[icluster]->fMinExtent;
    Vector3D<double> maxext = tsl.fClusters[icluster]->fMaxExtent;
    Vector3D<double> dext   = maxext - minext;
    pl.SetPoint(0, minext.x(), minext.y(), minext.z());
    pl.SetPoint(1, minext.x() + dext.x(), minext.y(), minext.z());
    visualizer.AddLine(pl);
    pl.SetPoint(0, minext.x() + dext.x(), minext.y(), minext.z());
    pl.SetPoint(1, minext.x() + dext.x(), minext.y() + dext.y(), minext.z());
    visualizer.AddLine(pl);
    pl.SetPoint(0, minext.x() + dext.x(), minext.y() + dext.y(), minext.z());
    pl.SetPoint(1, minext.x(), minext.y() + dext.y(), minext.z());
    visualizer.AddLine(pl);
    pl.SetPoint(0, minext.x(), minext.y() + dext.y(), minext.z());
    pl.SetPoint(1, minext.x(), minext.y(), minext.z());
    visualizer.AddLine(pl);

    pl.SetPoint(0, minext.x(), minext.y(), minext.z());
    pl.SetPoint(1, minext.x(), minext.y(), minext.z() + dext.z());
    visualizer.AddLine(pl);
    pl.SetPoint(0, minext.x() + dext.x(), minext.y(), minext.z());
    pl.SetPoint(1, minext.x() + dext.x(), minext.y(), minext.z() + dext.z());
    visualizer.AddLine(pl);
    pl.SetPoint(0, minext.x() + dext.x(), minext.y() + dext.y(), minext.z());
    pl.SetPoint(1, minext.x() + dext.x(), minext.y() + dext.y(), minext.z() + dext.z());
    visualizer.AddLine(pl);
    pl.SetPoint(0, minext.x(), minext.y() + dext.y(), minext.z());
    pl.SetPoint(1, minext.x(), minext.y() + dext.y(), minext.z() + dext.z());
    visualizer.AddLine(pl);

    pl.SetPoint(0, minext.x(), minext.y(), minext.z() + dext.z());
    pl.SetPoint(1, minext.x() + dext.x(), minext.y(), minext.z() + dext.z());
    visualizer.AddLine(pl);
    pl.SetPoint(0, minext.x() + dext.x(), minext.y(), minext.z() + dext.z());
    pl.SetPoint(1, minext.x() + dext.x(), minext.y() + dext.y(), minext.z() + dext.z());
    visualizer.AddLine(pl);
    pl.SetPoint(0, minext.x() + dext.x(), minext.y() + dext.y(), minext.z() + dext.z());
    pl.SetPoint(1, minext.x(), minext.y() + dext.y(), minext.z() + dext.z());
    visualizer.AddLine(pl);
    pl.SetPoint(0, minext.x(), minext.y() + dext.y(), minext.z() + dext.z());
    pl.SetPoint(1, minext.x(), minext.y(), minext.z() + dext.z());
    visualizer.AddLine(pl);
    return;
  }

  unsigned nfacets = 0;
  unsigned ifacet  = 0;
  unsigned iother  = 0;
  TriangleFacet<double> *facets[kVecSize];
  while (ifacet < kVecSize) {
    bool add = true;
    for (unsigned i = 0; i < nfacets; ++i) {
      if (tsl.fClusters[icluster]->fFacets[ifacet] == facets[i]) {
        ifacet++;
        add = false;
        break;
      }
    }
    if (add) facets[nfacets++] = tsl.fClusters[icluster]->fFacets[ifacet++];
  }
  // Loop facets
  ifacet = 0;
  int ivert[2];
  while (ifacet < nfacets) {
    // loop segments
    for (int iseg = 0; iseg < 3; iseg++) {
      bool shared = false;
      ivert[0]    = facets[ifacet]->fIndices[iseg];
      ivert[1]    = facets[ifacet]->fIndices[(iseg + 1) % 3];
      // loop remaining facets
      for (iother = 0; iother < nfacets; iother++) {
        if (iother == ifacet) continue;
        // check if the other facet has the 2 vertices
        if (facets[iother]->fIndices[0] != ivert[0] && facets[iother]->fIndices[1] != ivert[0] &&
            facets[iother]->fIndices[2] != ivert[0])
          continue;
        if (facets[iother]->fIndices[0] != ivert[1] && facets[iother]->fIndices[1] != ivert[1] &&
            facets[iother]->fIndices[2] != ivert[1])
          continue;
        // The line is shared
        shared = true;
        break;
      }
      if (shared) continue;
      // Add the line segment to the visualizer
      pl.SetPoint(0, facets[ifacet]->fVertices[iseg].x(), facets[ifacet]->fVertices[iseg].y(),
                  facets[ifacet]->fVertices[iseg].z());
      pl.SetPoint(1, facets[ifacet]->fVertices[(iseg + 1) % 3].x(), facets[ifacet]->fVertices[(iseg + 1) % 3].y(),
                  facets[ifacet]->fVertices[(iseg + 1) % 3].z());
      visualizer.AddLine(pl);
    }
    ifacet++;
  }
}
#endif // VECGEOM_ROOT
#endif // VECGEOM_CUDA

int main(int argc, char *argv[])
{
#ifndef VECGEOM_ENABLE_CUDA
  using namespace vecgeom;
  //  using Real_v = typename VectorBackend::Real_v;

  OPTION_INT(ngrid, 100);
  OPTION_INT(npoints, 10000);
#ifdef VECGEOM_ROOT
  OPTION_INT(vis, 0);
  OPTION_INT(scalability, 0);
  constexpr double r = 10.;
  int ngrid1         = 10;
  int i              = 0;
  const double sqrt2 = vecCore::math::Sqrt(2.);
#endif
  Vector3D<double> start(0, 0, 0);
  Vector3D<double> point;

  Vector3D<double> *dirs = new Vector3D<double>[npoints];
  for (int i = 0; i < npoints; ++i)
    RandomDirection(dirs[i]);

#ifdef VECGEOM_ROOT
  TGraph *gtime = nullptr;
  if (scalability) {
    gtime = new TGraph(14);
    while (ngrid1 < 1000) {
      SimpleTessellated *stsl   = new SimpleTessellated("test_VecGeomTessellated");
      UnplacedTessellated *tsl1 = (UnplacedTessellated *)stsl->GetUnplacedVolume();
      int nfacets1              = TessellatedOrb(r, ngrid1, *tsl1);
      // Close the solid
      Stopwatch timer;
      timer.Start();
      tsl1->Close();
      double tbuild = timer.Stop();
      gtime->SetPoint(i++, nfacets1, tbuild);
      // Check Distance performance
      timer.Start();
      for (int i = 0; i < npoints; ++i)
        stsl->DistanceToIn(start, dirs[i]);
      double trun = timer.Stop();
      printf("n=%d ngrid=%d nfacets=%d  build time=%g run time=%g\n", i, ngrid1, nfacets1, tbuild, trun);
      delete tsl1;
      ngrid1 = sqrt2 * double(ngrid1);
    }
  }
#endif
  SimpleTessellated *stsl1             = new SimpleTessellated("test_VecGeomTessellated");
  UnplacedTessellated *utsl            = (UnplacedTessellated *)stsl1->GetUnplacedVolume();
  TessellatedStruct<double> const &tsl = utsl->GetStruct();
  TessellatedOrb(r, ngrid, *utsl);
  utsl->Close();
  std::cout << "=== Tessellated solid statistics: nfacets = " << tsl.fFacets.size()
            << "  nclusters = " << tsl.fClusters.size() << "  kVecSize = " << kVecSize << std::endl;
  std::cout << "    cluster distribution: ";
  for (unsigned i = 1; i <= kVecSize; ++i) {
    std::cout << i << ": " << tsl.fNcldist[i] << " | ";
  }
  std::cout << "\n";

// Visualize the facets
#ifdef VECGEOM_ROOT
  // Analyze clusters
  int nblobs, nfacets;
  int nfacetstot    = 0;
  TH1F *hdispersion = new TH1F("hdispersion", "Cluster dispersion", 100, 0., 10.);
  TH1I *hblobs      = new TH1I("hblobs", "Blobs in clusters", 8, 0, 8);
  for (unsigned icluster = 0; icluster < tsl.fClusters.size(); ++icluster) {
    double dispersion = tsl.fClusters[icluster]->ComputeSparsity(nblobs, nfacets);
    nfacetstot += nfacets;
    hdispersion->Fill(dispersion);
    hblobs->Fill(nblobs);
  }
  printf("Number of facets = %d/%lu\n", nfacetstot, tsl.fFacets.size());
  TFile *file = TFile::Open("dispersion.root", "RECREATE");
  hdispersion->Write();
  hblobs->Write();
  if (gtime) gtime->Write();
  file->Write();

  if (vis) {
    Visualizer visualizer;
    // Visualize bounding box
    Vector3D<double> deltas = 0.5 * (tsl.fMaxExtent - tsl.fMinExtent);
    Vector3D<double> origin = 0.5 * (tsl.fMaxExtent + tsl.fMinExtent);
    SimpleBox box("bbox", deltas.x(), deltas.y(), deltas.z());
    visualizer.AddVolume(box, Transformation3D(origin.x(), origin.y(), origin.z()));

    // Visualize facets

    for (auto facet : tsl.fFacets)
      AddFacetToVisualizer(facet, visualizer);

    // Visualize clusters
    //    for (unsigned icluster = 0; icluster < tsl.fClusters.size(); ++icluster)
    //      DrawCluster(tsl, icluster, visualizer, false);

    TPolyMarker3D pm(npoints);
    pm.SetMarkerColor(kRed);
    pm.SetMarkerStyle(7);
    // Test contains function

    for (int i = 0; i < npoints; ++i) {
      RandomPointInBBox(point, tsl);
      if (0) {
        // Visualize a specific point/direction
        point.Set(-8, 8, 0);
        Vector3D<double> direction(-0.74608321159322855, -0.28587882094198169, -0.60135941093123035);
        pm.SetNextPoint(point[0], point[1], point[2]);
        TPolyLine3D pl(2);
        pl.SetLineColor(kRed);
        pl.SetNextPoint(point[0], point[1], point[2]);
        point += direction * 25;
        pm.SetNextPoint(point[0], point[1], point[2]);
        pl.SetNextPoint(point[0], point[1], point[2]);
        visualizer.AddLine(pl);
        break;
      }
      bool contains;
      TessellatedImplementation::Contains<double, bool>(tsl, point, contains);
      if (contains) pm.SetNextPoint(point[0], point[1], point[2]);
    }

    // Test distance to out from origin

    /*
        for (int i = 0; i < npoints; ++i) {
          tsl.DistanceToSolid<false>(start, dirs[i], InfinityLength<double>(), distance, ifacet);
          point = start + distance * dirs[i];
          pm.SetNextPoint(point[0], point[1], point[2]);
        }
    */
    delete[] dirs;

    visualizer.AddPoints(pm);
    visualizer.Show();
#endif
    return 0;
  }

#endif

  return 0;
}