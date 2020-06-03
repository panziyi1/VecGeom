/*
 * RaytraceBenchmark.cpp
 *
 *  Created on: May 8, 2020
 *      Author: andrei.gheata@cern.ch
 */

#include <iomanip>
#include <VecGeom/management/GeoManager.h>
#include <VecGeom/management/NavIndexTable.h>
#include <VecGeom/navigation/NavigationState.h>
#include <VecGeom/base/Stopwatch.h>
#include "ArgParser.h"

#ifdef VECGEOM_GDML
#include "Frontend.h"
#endif

using namespace vecgeom;

void TestNavIndexCPU(vecgeom::cxx::VPlacedVolume const* const world, int maxdepth, int npasses);

#ifdef VECGEOM_ENABLE_CUDA
void TestNavIndexGPU(vecgeom::cxx::VPlacedVolume const* const world, int maxdepth, int npasses);
#endif

namespace visitorcxx {

class GlobalToLocalVisitor {
private:
  int fType = 0; ///< visit type: 0=no computation 1=master-to-local NavStatePath 2=master-to-local NavStateIndex
  Vector3D<Precision> fGlobalPoint; ///< point coordinates in the global frame
  Vector3D<Precision> fLocal;       ///< make sure we store the result somewhere

public:
  GlobalToLocalVisitor() {}

  VECCORE_ATT_HOST_DEVICE
  void SetType(int type) { fType = type; }

  VECCORE_ATT_HOST_DEVICE
  void SetGlobalPoint(Precision x, Precision y, Precision z) { fGlobalPoint.Set(x, y, z); }

  VECCORE_ATT_HOST_DEVICE
  void apply(NavStatePath *state, NavIndex_t nav_ind)
  {
    switch (fType) {
    case 0:
      break;
    case 1:
      fLocal = state->GlobalToLocal(fGlobalPoint);
      break;
    case 2:
      fLocal = NavStateIndex::GlobalToLocalImpl(nav_ind, fGlobalPoint);
      break;
    }
  }
};

/// Traverses the geometry tree keeping track of the state context (volume path or navigation state)
/// and applies the injected Visitor
VECCORE_ATT_HOST_DEVICE
template <typename Visitor>
void visitAllPlacedVolumesPassNavIndex(VPlacedVolume const *currentvolume, Visitor *visitor, NavStatePath *state,
                                       NavIndex_t nav_ind)
{
  if (currentvolume != NULL) {
    state->Push(currentvolume);
    visitor->apply(state, nav_ind);
    for (auto daughter : currentvolume->GetDaughters()) {
      auto nav_ind_d = NavStateIndex::PushImpl(nav_ind, daughter);
      visitAllPlacedVolumesPassNavIndex(daughter, visitor, state, nav_ind_d);
    }
    state->Pop();
  }
}

} // namespace visitorcxx

void TestNavIndexCPU(vecgeom::cxx::VPlacedVolume const* const world, int maxdepth, int npasses)
{
  // Check performance
  using namespace visitorcxx;

  Stopwatch timer;
  NavStatePath *state = NavStatePath::MakeInstance(maxdepth);
  state->Clear();
  GlobalToLocalVisitor visitor;

  NavIndex_t nav_ind_top = 1; // The navigation index corresponding to the world

  // empty pass to subtract the overhead
  visitor.SetType(0);
  if (npasses > 1) {
    // Warming phase
    for (auto i = 0; i < 5; ++i)
      visitAllPlacedVolumesPassNavIndex(world, &visitor, state, nav_ind_top);
  }

  timer.Start();
  for (auto i = 0; i < npasses; ++i)
    visitAllPlacedVolumesPassNavIndex(world, &visitor, state, nav_ind_top);
  auto tbaseline = timer.Stop();

  // Compute the time for NavStatePath::GlobalToLocal
  visitor.SetType(1);
  if (npasses > 1) {
    // Warming phase
    for (auto i = 0; i < 5; ++i)
      visitAllPlacedVolumesPassNavIndex(world, &visitor, state, nav_ind_top);
  }
  timer.Start();
  for (auto i = 0; i < npasses; ++i)
    visitAllPlacedVolumesPassNavIndex(world, &visitor, state, nav_ind_top);
  auto tnavstate = timer.Stop();
  std::cout << "NavStatePath::GlobalToLocal took: " << tnavstate - tbaseline << " sec.\n";

  // Compute the time for NavStateIndex::GlobalToLocal
  visitor.SetType(2);
  if (npasses > 1) {
    // Warming phase
    for (auto i = 0; i < 5; ++i)
      visitAllPlacedVolumesPassNavIndex(world, &visitor, state, nav_ind_top);
  }
  timer.Start();
  for (auto i = 0; i < npasses; ++i)
    visitAllPlacedVolumesPassNavIndex(world, &visitor, state, nav_ind_top);
  auto tnavindex = timer.Stop();
  std::cout << "NavStateIndex::GlobalToLocal took: " << tnavindex - tbaseline << " sec.\n";

  std::cout << "Speedup per GlobalToLocal averaged over all states: " << std::setprecision(3)
            << (tnavstate - tbaseline) / (tnavindex - tbaseline) << "\n";

  NavStatePath::ReleaseInstance(state);
}

int main(int argc, char *argv[])
{
#ifndef VECGEOM_GDML
  std::cout << "### VecGeom must be compiled with GDML support to run this.\n";
  return 1;
#endif

  OPTION_STRING(gdml_name, "default.gdml");
  OPTION_INT(max_depth, 0);
  OPTION_INT(validate, 0);
  OPTION_INT(npasses, 10);
  OPTION_INT(on_gpu, 0);

  Stopwatch timer;
// Try to open the input file
#ifdef VECGEOM_GDML
  timer.Start();
  bool load = vgdml::Frontend::Load(gdml_name.c_str(), false);
  if (!load) return 2;
  auto tload = timer.Stop();
  std::cout << "Loading " << gdml_name << " took " << tload << " sec.\n";
#endif

  auto world = GeoManager::Instance().GetWorld();
  if (!world) return 3;
  int maxdepth = GeoManager::Instance().getMaxDepth();
  std::cout << "geometry max depth: " << maxdepth
            << "  total physical nodes: " << GeoManager::Instance().GetTotalNodeCount() << "\n";

  timer.Start();
  bool success = GeoManager::Instance().MakeNavIndexTable(max_depth, validate);
  auto tbuild  = timer.Stop();
  std::cout << "Building the navigation index table took: " << tbuild << " sec.\n";

  float frac_build = 100. * tbuild / tload;
  std::cout << "Navigation table build time as fraction of the GDML load time: " << std::setprecision(2) << frac_build
            << " %\n";
  
  if (on_gpu) {
  #ifdef VECGEOM_ENABLE_CUDA
    TestNavIndexGPU(GeoManager::Instance().GetWorld(), maxdepth, npasses);
  #else
    std::cout << "=== Cannot run the test on GPU since VecGeom CUDA support not compiled.\n";
    return 1;
  #endif
  } else {
    TestNavIndexCPU(GeoManager::Instance().GetWorld(), maxdepth, npasses);
  }

  if (!success) return 1;
  return 0;
}
