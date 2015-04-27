//
//
// TestBox
//             Ensure asserts are compiled in

#undef NDEBUG
#include "base/Global.h"
#include "base/Vector3D.h"
#include "volumes/Box.h"
#include "volumes/Hype.h"
#include "ApproxEqual.h"
#ifdef VECGEOM_USOLIDS
#include "UBox.hh"
#include "G4Hype.hh"
#include "UVector3.hh"
#endif

#include <cassert>
#include <cmath>
#include <iomanip> 

#define PI 3.14159265358979323846


template <class Hype_t, class Vec_t = vecgeom::Vector3D<vecgeom::Precision> >
bool TestHype() {
    
    
    Hype_t test("Solid VecGeomHype #test",10,15, PI/4,PI/4, 50); //5,6,0.523598775598298,0.523598775598298,10);
    Vec_t testOutPoint(18,0,0);
	double Dist = test.SafetyFromOutside(testOutPoint);
    std::cout<<"SafetyFromOutside : "<<Dist<<std::endl;
	
	Vec_t testInPoint(12,0,0);
	Dist = test.SafetyFromInside(testInPoint);
    std::cout<<"SafetyFromInside : "<<Dist<<std::endl;

    std::cout<<"Capacity is : "<<test.Capacity()<<std::endl;
    std::cout<<"SurfaceArea is : "<<test.SurfaceArea()<<std::endl;

	Vec_t testPoint(18,0,0);
	Vec_t normal;
	test.Normal(testPoint,normal);
	std::cout<<"Normal Calculated from VecGeom is  : "<<normal<<std::endl;
	
    return true;
}

int main() {
    
//#ifdef VECGEOM_USOLIDS
//  assert(TestHype<UHype>());
//  std::cout << "UHype passed\n";
//#endif
  std::cout<<"-------------------------------------------------------------------------------------------------"<<std::endl;
  //std::cout<<"*************************************************************************************************"<<std::endl;
  assert(TestHype<vecgeom::SimpleHype>());
  //vecgeom::SimpleHype test("Solid VecGeomHype #test",5, 30, 6, 30, 10); //5,6,0.523598775598298,0.523598775598298,10);


  std::cout << "VecGeomHype passed\n";
  return 0;
}

