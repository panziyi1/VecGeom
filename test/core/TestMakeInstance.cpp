#include "volumes/UnplacedBox.h"
#include "volumes/UnplacedTube.h"
#include "management/GeoManager.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>

using namespace vecgeom;

int main()
{
  auto ubox = GeoManager::MakeInstance<UnplacedBox>(1., 1., 2.);
  assert(ubox != nullptr);
  assert(dynamic_cast<UnplacedBox *>(ubox));

  auto utube = GeoManager::MakeInstance<UnplacedTube>(0., 1., 1., 0., 2. * M_PI);
  assert(utube != nullptr);
  assert(dynamic_cast<SUnplacedTube<TubeTypes::NonHollowTube> *>(utube));

  std::cout << "test passed \n";
  return 0;
}
