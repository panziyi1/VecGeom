#include "RootPersistencyComponentsTest.h"

using namespace std;
using namespace vecgeom;

void usage(std::string name)
{
  cerr << "Usage: " << name << " <test>\n"
       << "Available tests:\n"
       << "\t1 -> perform all tests\n"
       << "\t2 -> UnplacedBox test\n"
       << "\t3 -> LogicalVolume test\n"
       << "\t4 -> std:vector test" << endl;
}

int main(int argc, char *argv[])
{
  if (argc != 2) {
    usage(argv[0]);
    return 0;
  }

  switch (atoi(argv[1])) {
  case 1:
    bool u, l, v;
    u = unplaced_test();
    l = logical_test();
    v = vector_test();
    return (u && l && v) ? 0 : 2;
  case 2:
    return unplaced_test() ? 0 : 2;
  case 3:
    return logical_test() ? 0 : 2;
  case 4:
    return vector_test() ? 0 : 2;
  default:
    usage(argv[0]);
  }
  return 1;
}
