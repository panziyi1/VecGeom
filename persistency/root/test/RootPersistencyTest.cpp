#include "RootPersistencyTest.h"

using namespace std;
using namespace vecgeom;

void usage()
{
  cout << "Usage:" << endl;
  cout << "\t./RootPersistencyTest test" << endl << endl;
  cout << "available tests:" << endl;
  cout << "u -> UnplacedBox test" << endl;
  cout << "l -> LogicalVolume test" << endl;
  cout << "v -> vecgeom::Vector test" << endl;
}

int main(int argc, char *argv[])
{
  if (argc != 2) {
    usage();
    return 0;
  }

  switch (argv[1][0]) {
  case 'u':
    unplaced_test();
    break;
  case 'l':
    logical_test();
    break;
  case 'v':
    vector_test();
    break;
  default:
    usage();
  }
  return 0;
}