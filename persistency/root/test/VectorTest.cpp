#include "RootPersistencyComponentsTest.h"

using namespace std;
using namespace vecgeom;
#include <vector>

bool vector_test()
{
  cout << "///// Running vector_test /////" << endl << endl;

  RootPersistencyProxy(); // calling the proxy

  auto a = new double(9.548);

  auto vec  = new std::vector<double>();
  auto vecp = new std::vector<double *>();
  vec->push_back(1.24);
  vec->push_back(2.14);
  vec->push_back(70);
  vecp->push_back(a);

  cout << "writing on vec.root" << endl;
  TFile fo("vec.root", "RECREATE");
  fo.WriteObject(vec, "vec_saved");
  fo.WriteObject(vecp, "vecp_saved");
  fo.Close();

  cout << "reading from vec.root\n\n" << endl << endl;
  TFile fi("vec.root");

  std::vector<double> *rvec;
  std::vector<double *> *rvecp;

  fi.GetObject("vec_saved", rvec);
  fi.GetObject("vecp_saved", rvecp);

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

  // [3]
  cout << "[3] comparing vector<double*> size\n\n"
       << ">> before\n"
       << "vector<double*>.size(): " << vecp->size() << "\n\n"
       << ">> after\n"
       << "vector<double*>.size(): " << rvecp->size() << "\n-----------------" << endl;
  test_ok = (vecp->size() == rvecp->size());

  if (test_ok)
    cout << "test passed\n\n" << endl;
  else {
    cout << "! test not passed\n\n" << endl;
    all_test_ok = false;
  }

  // [4]
  cout << "[4] comparing elements of vector<double*>\n"
       << ">> before\n"
       << "vector<double*>: ";

  for (auto i : *vecp)
    cout << *i << ",";

  cout << "\n\n"
       << ">> after\n"
       << "vector<double*>: ";

  for (auto i : *rvecp)
    cout << *i << ",";

  cout << "\n-----------------" << endl;
  test_ok = true;
  for (auto i = vecp->begin(), j = rvecp->begin(), end = vecp->end(), rend = rvecp->end(); i != end && j != rend;
       i++, j++)
    if (*i != *j) test_ok = false;

  if (test_ok)
    cout << "test passed\n\n" << endl;
  else {
    cout << "! test not passed\n\n" << endl;
    all_test_ok = false;
  }

  return all_test_ok;
}