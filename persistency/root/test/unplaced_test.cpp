#include "RootPersistencyTest.h"

using namespace std;
using namespace vecgeom;

void unplaced_test()
{
  cout << "Running unplaced_test" << endl << endl;

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
  TFile fi("geo.root");

  UnplacedBox *rbox;
  UnplacedParaboloid *rpar;
  UnplacedParallelepiped *rpal;
  UnplacedBox *rbox2;

  fi.GetObject("box_saved", rbox);
  fi.GetObject("par_saved", rpar);
  fi.GetObject("pal_saved", rpal);
  fi.GetObject("box2_saved", rbox2);

  // printing
  box->Print();
  cout << endl;
  rbox->Print();
  cout << endl << endl;

  par->Print();
  cout << endl;
  rpar->Print();
  cout << endl << endl;

  pal->Print();
  cout << endl;
  rpal->Print();
  cout << endl << endl;

  box2->Print();
  cout << endl;
  rbox2->Print();
  cout << endl << endl;
}