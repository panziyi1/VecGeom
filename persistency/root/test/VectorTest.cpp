#include "RootPersistencyComponentsTest.h"

using namespace std;
using namespace vecgeom;
#include <base/Vector.h>
#include <base/Array.h>
#include <base/SOA3D.h>
#include <volumes/Planes.h>
#include <volumes/Quadrilaterals.h>

bool vector_test()
{
  cout << "///// Running vector_test /////" << endl << endl;

  RootPersistencyProxy(); // calling the proxy

//  auto vec = new vecgeom::Array<double>(3);
  auto arr = new vecgeom::Array<double>(3);
  (*arr)[0] = 1.24;
  (*arr)[1] = 2.14;
  (*arr)[2] = 70;
  std::cout << "arr: " << (*arr)[0] << ", " << (*arr)[1] << ", " << (*arr)[2] << "\n";
  
  auto vec = new vecgeom::Vector<double>(3);
  vec->push_back(1.24);
  vec->push_back(2.14);
  vec->push_back(70);  
  std::cout << "vec: " << (*vec)[0] << ", " << (*vec)[1] << ", " << (*vec)[2] << "\n";

  auto soa = new vecgeom::SOA3D<double>(10);
  std::cout << "SOA3D:\n";
  for (auto i = 0; i < 10; ++i) {
    soa->set(i, i, i+1, i+2);
    std::cout << "   " << i << ": " << soa->x(i) << ", " << soa->y(i) << ", " << soa->z(i) << "\n";
  }

  auto planes = new vecgeom::Planes(4);
  planes->Set(0, {-1, 0, 0}, 0.5);
  planes->Set(1, {0, -1, 0}, 0.5);
  planes->Set(2, {1, 0, 0}, 0.5);
  planes->Set(3, {0, 1, 0}, 0.5);
  std::cout << "Planes: " << *planes << std::endl;

  auto quads = new vecgeom::Quadrilaterals(4);
  quads->Set(0, {-1,-1,-1}, {-1,-1,1}, {-1,1,1}, {-1,1,-1}, 0.5);
  quads->Set(1, {-1,-1,-1}, {1,-1,-1}, {1,1,-11}, {-1,-1,1}, 0.5);
  quads->Set(2, {1,-1,-1}, {1,1,1}, {1,1,1}, {1,1,-1}, 0.5);
  quads->Set(3, {-1,1,-1}, {1,1,-1}, {1,1,1}, {-1,1,1}, 0.5);
  std::cout << *quads << std::endl;

  cout << "writing on vec.root...\n\n" << endl;
  TFile fo("vec.root", "RECREATE");
  //gDebug = 10;
  fo.WriteObject(arr, "arr_saved");
  fo.WriteObject(vec, "vec_saved");
  fo.WriteObject(soa, "soa_saved");
  fo.WriteObject(planes, "planes_saved");
  fo.WriteObject(quads, "quads_saved");
  fo.Close();

  cout << "reading from vec.root" << endl << endl;
  TFile fi("vec.root");

  vecgeom::Array<double> *rarr;
  vecgeom::Vector<double> *rvec;
  vecgeom::SOA3D<double> *rsoa;
  vecgeom::Planes *rplanes;
  vecgeom::Quads *rquads;

  // gDebug = 3;

  fi.GetObject("arr_saved", rarr);
  std::cout << "rvec: " << (*rarr)[0] << ", " << (*rarr)[1] << ", " << (*rarr)[2] << "\n";
  fi.GetObject("vec_saved", rvec);
  std::cout << "rvec: " << (*rvec)[0] << ", " << (*rvec)[1] << ", " << (*rvec)[2] << "\n";
  fi.GetObject("soa_saved", rsoa);
  std::cout << "SOA3D:\n";
  for (auto i = 0; i < 10; ++i)
    std::cout << "   " << i << ": " << rsoa->x(i) << ", " << rsoa->y(i) << ", " << rsoa->z(i) << "\n";
  fi.GetObject("planes_saved", rplanes);
  std::cout << "Planes: " << *rplanes << std::endl;
  fi.GetObject("quads_saved", rquads);
  std::cout << "Quads: " << *rquads << std::endl;
  return true;
  /*

  // testing
  bool all_test_ok = true;
  bool test_ok;

  // [1]
  cout << "[1] comparing vector<double> size\n\n"
       << ">> before\n"
       << "vector<double>.size(): " << vec->size() << "\n\n"
       << ">> after\n"
       << "vector<double>.size(): " << rvec->size() << "\n-----------------" << endl;
  test_ok = (vec->size() == rvec->size());
  if (test_ok)
    cout << "test passed\n\n" << endl;
  else {
    cout << "! test not passed\n\n" << endl;
    all_test_ok = false;
  }

  // [2]
  cout << "[2] comparing elements of vector<double>\n"
       << ">> before\n"
       << "vector<double>: ";

  for (auto i : *vec)
    cout << i << ",";

  cout << "\n\n"
       << ">> after\n"
       << "vector<double>: ";

  for (auto i : *rvec)
    cout << i << ",";

  cout << "\n-----------------" << endl;
  test_ok = true;
  for (auto i = vec->begin(), j = rvec->begin(), end = vec->end(), rend = rvec->end(); i != end && j != rend; i++, j++)
    if (*i != *j) test_ok = false;

  if (test_ok)
    cout << "test passed\n\n" << endl;
  else {
    cout << "! test not passed\n\n" << endl;
    all_test_ok = false;
  }
  return all_test_ok;
  */
}