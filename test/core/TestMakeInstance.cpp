#include "VecGeom/volumes/UnplacedBox.h"
#include "VecGeom/volumes/UnplacedTube.h"
#include "VecGeom/volumes/UnplacedCone.h"
#include "VecGeom/volumes/UnplacedOrb.h"
#include "VecGeom/volumes/LogicalVolume.h"
#include "VecGeom/management/GeoManager.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>

using namespace vecgeom;

// a unit test checking the factory mechanism to produce specialized unplaced
// volumes
int main()
{
  // BOX: IS TRIVIAL
  auto ubox = GeoManager::MakeInstance<UnplacedBox>(1., 1., 2.);
  assert(ubox != nullptr);
  assert(dynamic_cast<UnplacedBox *>(ubox));
  // let me try to make a specialized placed box
  Transformation3D placement(0, 0, 0);
  LogicalVolume lv("mybox", ubox);
  auto pv = lv.Place(&placement);
  assert(pv->Contains(Vector3D<double>(0, 0, 0)));

  // ORB: IS TRIVIAL
  auto uorb = GeoManager::MakeInstance<UnplacedOrb>(1.);
  assert(uorb != nullptr);
  assert(dynamic_cast<UnplacedOrb *>(uorb));

  // CHECK THE TUBE CASES
  {
    // an ordinary tube without inner radius
    auto utube = GeoManager::MakeInstance<UnplacedTube>(0., 1., 1., 0., 2. * M_PI);
    assert(utube != nullptr);
    assert(dynamic_cast<UnplacedTube *>(utube));

    // let me try to make a specialized placed hollow tube
    Transformation3D placement(0, 0, 0);
    LogicalVolume lv("mytube", utube);
    auto pv = lv.Place(&placement);
    auto c  = pv->Contains(Vector3D<double>(0, 0, 0));
    assert(c);
  }

  {
    // an ordinary hollow tube
    auto utube = GeoManager::MakeInstance<UnplacedTube>(0.5, 1., 1., 0., 2. * M_PI);
    assert(utube != nullptr);
    assert(dynamic_cast<UnplacedTube *>(utube));
  }

  // CHECK THE CONE CASES
  {
    // an ordinary cone without inner radii
    auto ucone = GeoManager::MakeInstance<UnplacedCone>(0., 1., 0., 1., 2., 0., kTwoPi);
    assert(ucone != nullptr);
    assert(dynamic_cast<UnplacedCone *>(ucone));

    // let me try to make a specialized placed hollow cone
    Transformation3D placement(0, 0, 0);
    LogicalVolume lv("mycone", ucone);
    auto pv = lv.Place(&placement);
    auto c  = pv->Contains(Vector3D<double>(0, 0, 0));
    assert(c);
  }

  {
    // an ordinary hollow cone
    auto ucone = GeoManager::MakeInstance<UnplacedCone>(0.5, 1., 0.4, 1., 1.8, 0., kTwoPi);
    assert(ucone != nullptr);
    assert(dynamic_cast<UnplacedCone *>(ucone));
  }

  {
    // a hollow cone with a smaller than PI sector
    auto ucone = GeoManager::MakeInstance<UnplacedCone>(0.5, 1., 0.4, 1., 1.8, 0., kPi / 3.);
    assert(ucone != nullptr);
    assert(dynamic_cast<UnplacedCone *>(ucone));
  }

  std::cout << "test passed \n";
  return 0;
}
