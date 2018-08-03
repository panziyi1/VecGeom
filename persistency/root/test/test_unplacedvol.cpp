#include "test_unplacedvol.hh"

using namespace vecgeom;
using namespace std;

int main(int argc, char *argv[])
{
  if (argc > 1 && strcmp(argv[1], "w") == 0) { // if w -> write the file

    TFile fo("geo.root", "RECREATE");

    // example unplaced volumes:

    cout << "writing on geo.root:";

    auto box = new UnplacedBox(.1, .2, .1);
    cout << "\n\t";
    box->Print();
    fo.WriteObject(box, "box_saved");

    auto par = new UnplacedParaboloid(.1, .2, .1);
    cout << "\n\t";
    par->Print();
    fo.WriteObject(par, "par_saved");

    auto pal = new UnplacedParallelepiped();
    cout << "\n\t";
    pal->Print();
    fo.WriteObject(pal, "pal_saved");

    auto trd = new UnplacedTrd(3, 4, 5, 6, 7);
    cout << "\n\t";
    trd->Print();
    fo.WriteObject(trd, "trd_saved");

    auto gen_trp = new UnplacedGenTrap();
    cout << "\n\t";
    gen_trp->Print();
    fo.WriteObject(gen_trp, "gen_trp_saved");

    auto trapez = new UnplacedTrapezoid();
    cout << "\n\t";
    trapez->Print();
    fo.WriteObject(trapez, "trapez_saved");

//    auto hype = new UnplacedHype(0.2, .3, .4, .4, .01);
//    cout << "\n\t";
//    hype->Print();
//    fo.WriteObject(hype, "hype_saved");

    auto orb = new UnplacedOrb(.215);
    cout << "\n\t";
    orb->Print();
    fo.WriteObject(orb, "orb_saved");

    auto cuttube = new UnplacedCutTube();
    cout << "\n\t";
    cuttube->Print();
    fo.WriteObject(cuttube, "cuttube_saved");

    auto sphere = new UnplacedSphere(0.5, 1, 2, 0.4, 0.2, 0.1);
    cout << "\n\t";
    sphere->Print();
    fo.WriteObject(sphere, "sphere_saved");

  } else { // otherwise read it

    TFile fi("geo.root");

    cout << "reading from geo.root:";

    UnplacedBox *box;
    fi.GetObject("box_saved", box);
    cout << "\n\t";
    box->Print();

    UnplacedParaboloid *par;
    fi.GetObject("par_saved", par);
    cout << "\n\t";
    par->Print();

    UnplacedParallelepiped *pal;
    fi.GetObject("pal_saved", pal);
    cout << "\n\t";
    pal->Print();

    UnplacedTrd *trd;
    fi.GetObject("trd_saved", trd);
    cout << "\n\t";
    trd->Print();

    UnplacedGenTrap *gen_trp;
    fi.GetObject("gen_trp_saved", gen_trp);
    cout << "\n\t";
    gen_trp->Print();

    UnplacedTrapezoid *trapez;
    fi.GetObject("trapez_saved", trapez);
    cout << "\n\t";
    trapez->Print();

//    UnplacedHype *hype;
//    fi.GetObject("hype_saved", hype);
//    cout << "\n\t";
//    hype->Print();

    UnplacedOrb *orb;
    fi.GetObject("orb_saved", orb);
    cout << "\n\t";
    orb->Print();

    UnplacedCutTube *cuttube;
    fi.GetObject("cuttube_saved", cuttube);
    cout << "\n\t";
    cuttube->Print();

    UnplacedSphere *sphere;
    fi.GetObject("sphere_saved", sphere);
    cout << "\n\t";
    sphere->Print();
    cout << "\n\t";
  }

  return 0;
}
