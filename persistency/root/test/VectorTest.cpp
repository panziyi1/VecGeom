#include "RootPersistencyTest.h"

using namespace std;
using namespace vecgeom;

void vector_test()
{
  cout << "Running vector_test" << endl << endl;

  RootPersistencyProxy(); // calling the proxy

  double a = 9.548;

  auto vec  = new Vector<double>();
  auto vecp = new Vector<double *>();
  vec->push_back(1.24);
  vec->push_back(2.14);
  vec->push_back(70);
  vecp->push_back(&a);

  cout << "writing on vec.root" << endl << endl;
  TFile fo("vec.root", "RECREATE");
  fo.WriteObject(vec, "vec_saved");
  fo.WriteObject(vecp, "vecp_saved");
  fo.Close();

  cout << "reading from vec.root" << endl << endl;
  TFile fi("vec.root");

  Vector<double> *rvec;
  Vector<double *> *rvecp;

  fi.GetObject("vec_saved", rvec);
  fi.GetObject("vecp_saved", rvecp);

  cout << "before Vector<double>:" << vec->size() << endl;
  cout << "after Vector<double>:" << rvec->size() << endl;
  for (auto i : *rvec)
    cout << i << ", ";

  cout << endl << "before Vector<double*>:" << vecp->size() << endl;
  cout << "after Vector<double*>:" << rvecp->size() << endl;
  for (auto i : *rvecp)
    cout << *i << "(" << i << ")"
         << ", ";

  cout << endl;
}