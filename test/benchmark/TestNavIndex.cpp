/*
 * RaytraceBenchmark.cpp
 *
 *  Created on: May 8, 2020
 *      Author: andrei.gheata@cern.ch
 */

#include <VecGeom/management/GeoManager.h>
#include <VecGeom/management/NavIndexTable.h>
#include <VecGeom/navigation/NavigationState.h>
#include <VecGeom/base/Stopwatch.h>
#include "ArgParser.h"

#ifdef VECGEOM_GDML
#include "Frontend.h"
#endif

using namespace vecgeom;

class GlobalToLocalVisitor {
private:
  NavIndexTable *fTable = nullptr; ///< pointer to navigation index table
  int fType = 0; ///< visit type: 0=no computation 1=master-to-local NavigationState 2=master-to-local NavIndex_t
  Vector3D<Precision> fGlobalPoint; ///< point coordinates in the global frame
  Vector3D<Precision> fLocal;       ///< make sure we store the result somewhere

public:
  GlobalToLocalVisitor(NavIndexTable *table) { fTable = NavIndexTable::Instance(); }

  void SetType(int type) { fType = type; }
  void SetGlobalPoint(Precision x, Precision y, Precision z) { fGlobalPoint.Set(x, y, z); }

  void apply(NavigationState *state, NavIndex_t nav_ind)
  {
    switch (fType) {
    case 0:
      break;
    case 1:
      fLocal = state->GlobalToLocal(fGlobalPoint);
      break;
    case 2:
      fLocal = fTable->GlobalToLocal(nav_ind, fGlobalPoint);
      break;
    }
  }
};

/// Traverses the geometry tree keeping track of the state context (volume path or navigation state)
/// and applies the injected Visitor
template <typename Visitor>
void visitAllPlacedVolumesPassNavIndex(VPlacedVolume const *currentvolume, Visitor *visitor, NavigationState *state,
                                       NavIndex_t nav_ind)
{
  if (currentvolume != NULL) {
    state->Push(currentvolume);
    visitor->apply(state, nav_ind);
    int size = currentvolume->GetDaughters().size();
    for (int i = 0; i < size; ++i) {
      auto nav_ind_d = NavIndexTable::Instance()->Push(nav_ind, i);
      visitAllPlacedVolumesPassNavIndex(currentvolume->GetDaughters().operator[](i), visitor, state, nav_ind_d);
    }
    state->Pop();
  }
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
  std::cout << "geometry max depth: " << GeoManager::Instance().getMaxDepth()
            << "  total physical nodes: " << GeoManager::Instance().GetTotalNodeCount() << "\n";

  timer.Start();
  bool success = GeoManager::Instance().MakeNavIndexTable(max_depth, validate);
  auto tbuild  = timer.Stop();
  std::cout << "Building the navigation index table took: " << tbuild << " sec.\n";

  float frac_build = 100. * tbuild / tload;

  // Check performance
  NavigationState *state = NavigationState::MakeInstance(GeoManager::Instance().getMaxDepth());
  state->Clear();
  auto visitor = new GlobalToLocalVisitor(NavIndexTable::Instance());

  NavIndex_t nav_ind_top = 1; // The navigation index corresponding to the world

  // empty pass to subtract the overhead
  visitor->SetType(0);
  timer.Start();
  for (auto i = 0; i < npasses; ++i)
    visitAllPlacedVolumesPassNavIndex(GeoManager::Instance().GetWorld(), visitor, state, nav_ind_top);
  auto tbaseline = timer.Stop();

  // Compute the time for NavigationState::GlobalToLocal
  visitor->SetType(1);
  timer.Start();
  for (auto i = 0; i < npasses; ++i)
    visitAllPlacedVolumesPassNavIndex(GeoManager::Instance().GetWorld(), visitor, state, nav_ind_top);
  auto tnavstate = timer.Stop();
  std::cout << "NavigationState::GlobalToLocal took: " << tnavstate - tbaseline << " sec.\n";

  // Compute the time for NavIndexTable::GlobalToLocal
  visitor->SetType(2);
  timer.Start();
  for (auto i = 0; i < npasses; ++i)
    visitAllPlacedVolumesPassNavIndex(GeoManager::Instance().GetWorld(), visitor, state, nav_ind_top);
  auto tnavindex = timer.Stop();
  std::cout << "NavIndexTable::GlobalToLocal took: " << tnavindex - tbaseline << " sec.\n";

  std::cout << "Navigation table build time as fraction of the GDML load time: " << std::setprecision(2) << frac_build
            << " %\n";
  std::cout << "Speedup per GlobalToLocal averaged over all states: " << std::setprecision(3)
            << (tnavstate - tbaseline) / (tnavindex - tbaseline) << "\n";

  NavigationState::ReleaseInstance(state);

  if (!success) return 1;
  return 0;
}
