#include "volumes/UnplacedAssembly.h"
#include "volumes/PlacedAssembly.h"
#include "navigation/NavigationState.h"
#include "volumes/Box.h"

// make sure that assert are included even in Release mode
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

using namespace vecgeom;

int main()
{

  // make a simple assembly of 2 boxes
  UnplacedBox *b    = new UnplacedBox(10., 10., 10.);
  LogicalVolume *lb = new LogicalVolume("boxlv", b);

  UnplacedAssembly *ass = new UnplacedAssembly();
  assert(ass->GetLogicalVolume() == nullptr);

  // this statement makes lv an assembly
  LogicalVolume *lv = new LogicalVolume("assemblylv", ass);
  assert(ass->GetLogicalVolume() == lv);

  // make simple assembly out of 2 placed boxes
  ass->AddVolume(lb->Place(new Transformation3D(-20., 0., 0.)));
  ass->AddVolume(lb->Place(new Transformation3D(20., 0., 0.)));
  assert(ass->GetNVolumes() == 2);

  // check that the bounding box is initialized
  assert(ass->GetLowerCorner().x() > -kInfinity);

  VPlacedVolume *const pv = lv->Place();

  // check the assembly property
  assert(pv->GetUnplacedVolume()->IsAssembly());
  assert(ass->IsAssembly());
  assert(!b->IsAssembly());

  GeoManager::Instance().SetWorld(pv);
  GeoManager::Instance().CloseGeometry();

  // verify correct conversion of Unplaced to Placed type
  PlacedAssembly const *pa = dynamic_cast<PlacedAssembly const *>(pv);
  assert(pa != nullptr);

  // some checks on Contains, Safety and DistanceToIn
  {
    Vector3D<Precision> p(0., 0., 0.);
    Vector3D<Precision> p2(20., 0., 0.);
    Vector3D<Precision> lp(0., 0., 0.);
    NavigationState *state = NavigationState::MakeInstance(10);
    std::cerr << pa->Contains(p, lp, *state) << "\n";
    assert(!pa->Contains(p, lp, *state));
    std::cerr << pa->Contains(p2, lp, *state) << "\n";
    assert(pa->Contains(p2, lp, *state));
    std::cerr << pa->Contains(p) << "\n";
    assert(!pa->Contains(p));
    std::cerr << pa->Contains(p2) << "\n";
    assert(pa->Contains(p2));

    assert(pa->SafetyToIn(Vector3D<Precision>(-10, 0, 0)) == 0.);
    assert(pa->SafetyToIn(Vector3D<Precision>(0, 0, 0)) == 10.);

    assert(pa->DistanceToIn(Vector3D<Precision>(-40, 0, 0), Vector3D<Precision>(1., 0, 0)) == 10);
    assert(pa->DistanceToIn(Vector3D<Precision>(0, -40, 0), Vector3D<Precision>(0, 1, 0)) == kInfinity);
    assert(pa->DistanceToIn(Vector3D<Precision>(0, 0, 0), Vector3D<Precision>(1., 0, 0)) == 10);
    assert(pa->DistanceToIn(Vector3D<Precision>(0, 0, 0), Vector3D<Precision>(-1., 0, 0)) == 10);
  }

  if (pv->GetUnplacedVolume()->IsAssembly()) {
    Vector3D<Precision> p(20., 0., 0.);
    Vector3D<Precision> lp(0., 0., 0.);
    NavigationState *state = NavigationState::MakeInstance(10);
    state->Clear();
    static_cast<PlacedAssembly const *>(pv)->Contains(p, lp, *state);
    state->Print();
  }

  // check Capacity and Surface Area
  assert(ass->Capacity() == 2. * b->Capacity());
  assert(ass->SurfaceArea() == 2. * b->SurfaceArea());
  assert(((PlacedAssembly *)pa)->Capacity() == 2. * b->Capacity());
  assert(((PlacedAssembly *)pa)->SurfaceArea() == 2. * b->SurfaceArea());

  // check Extent
  Vector3D<Precision> emin;
  Vector3D<Precision> emax;
  ass->Extent(emin, emax);
  assert(emin.x() <= -30);
  assert(emin.y() <= -10);
  assert(emin.z() <= -10);
  assert(emax.x() >= 30);
  assert(emax.y() >= 10);
  assert(emax.z() >= 10);

  assert(emin == ass->GetLowerCorner());
  assert(emax == ass->GetUpperCorner());

  // check Normal
  // TBD

  // check some properties of BoxImplementation
  BoxStruct<Precision> bs(kInfinity, kInfinity, kInfinity);
  Vector3D<Precision> p(0., 0., 0.);
  Vector3D<Precision> d(1., 0., 0.);
  Precision dist;
  BoxImplementation::DistanceToIn(bs, p, d, kInfinity, dist);
  assert(dist < kInfinity);
  bool cont;
  BoxImplementation::Contains(bs, p, cont);
  assert(cont);
  std::cerr << "dist " << dist << "\n";

  Vector3D<Precision> corners[2];
  corners[0].Set(-kInfinity, -kInfinity, -kInfinity);
  corners[1].Set(kInfinity, kInfinity, kInfinity);
  assert(BoxImplementation::Intersect(corners, p, d, 0, kInfinity));

  return 0;
}