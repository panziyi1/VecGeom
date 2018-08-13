#include "RootPersistencyComponentsTest.h"

using namespace std;
using namespace vecgeom;

bool unplaced_test()
{
  cout << "///// Running unplaced_test /////" << endl << endl;

  auto box  = new UnplacedBox(.05, .06, .07);
  auto par  = new UnplacedParaboloid(.08, .09, .10);
  auto pal  = new UnplacedParallelepiped(.11, .12, .13, 14, 15, 16);
  auto box2 = new UnplacedBox(.17, .18, .19);

  cout << "writing on unplaced.root" << endl << endl;

  TFile fo("unplaced.root", "RECREATE");

  fo.WriteObject(par, "par_saved");
  fo.WriteObject(box, "box_saved");
  fo.WriteObject(pal, "pal_saved");
  fo.WriteObject(box2, "box2_saved");

  fo.Close();

  // reading back

  cout << "reading from unplaced.root" << endl;
  TFile fi("unplaced.root");

  UnplacedBox *rbox;
  UnplacedParaboloid *rpar;
  UnplacedParallelepiped *rpal;
  UnplacedBox *rbox2;

  fi.GetObject("box_saved", rbox);
  fi.GetObject("par_saved", rpar);
  fi.GetObject("pal_saved", rpal);
  fi.GetObject("box2_saved", rbox2);
  // testing
  bool all_test_ok = true;
  bool test_ok;

  // [1]
  cout << "[1] comparing UnplacedBox box\n\n"
       << ">> before\n",
      box->Print();
  cout << "\n"
       << ">> after\n";
  rbox->Print();

  test_ok = (box->dimensions() == rbox->dimensions());
  cout << "\n\n";
  if (test_ok)
    cout << "test passed\n\n" << endl;
  else {
    cout << "! test not passed\n\n" << endl;
    all_test_ok = false;
  }

  // [2]
  cout << "[2] comparing UnplacedParaboloid\n\n"
       << ">> before\n",
      par->Print();
  cout << "\n"
       << ">> after\n";
  rpar->Print();
  cout << "\n\n";
  test_ok = (par->GetRlo() == rpar->GetRlo() && par->GetRhi() == rpar->GetRhi() && par->GetDz() == rpar->GetDz() &&
             par->GetA() == rpar->GetA() && par->GetB() == rpar->GetB());

  if (test_ok)
    cout << "test passed\n\n" << endl;
  else {
    cout << "! test not passed\n\n" << endl;
    all_test_ok = false;
  }
  // [3]
  cout << "[3] comparing UnplacedParallelepiped\n\n"
       << ">> before\n",
      pal->Print();
  cout << "\n"
       << ">> after\n";
  rpal->Print();
  cout << "\n\n";
  test_ok = (pal->GetX() == rpal->GetX() && pal->GetY() == rpal->GetY() && pal->GetZ() == rpal->GetZ() &&
             pal->GetTanAlpha() == rpal->GetTanAlpha() && pal->GetTanThetaCosPhi() == rpal->GetTanThetaCosPhi() &&
             pal->GetTanThetaSinPhi() == rpal->GetTanThetaSinPhi());

  if (test_ok)
    cout << "test passed\n\n" << endl;
  else {
    cout << "! test not passed\n\n" << endl;
    all_test_ok = false;
  }

  // [4]
  cout << "[4] comparing UnplacedBox box2\n\n"
       << ">> before\n",
      box2->Print();
  cout << "\n"
       << ">> after\n";
  rbox2->Print();

  test_ok = (box2->dimensions() == rbox2->dimensions());
  cout << "\n\n";
  if (test_ok)
    cout << "test passed\n\n" << endl;
  else {
    cout << "! test not passed\n\n" << endl;
    all_test_ok = false;
  }
  return all_test_ok;
}