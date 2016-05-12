#include "base/Vector3D.h"
#include "backend/Backend.h"
#include "volumes/Tube.h"
 
using namespace vecgeom;
 
int main()
{
    //using Backend = kVc;
    using Backend = kUmeSimd;
    using Real_v = Backend::precision_v;
 
    SOA3D<double> dirs(32);
    SOA3D<double> points(32);
 
    for(auto i=0;i<32;++i) {
      points.set(i,static_cast<double>(i),static_cast<double>(i),static_cast<double>(i));
      dirs.set(i,static_cast<double>(1.0),static_cast<double>(0),static_cast<double>(0));
    }
 
   // void *daughter;
    double rmax = 5.;
    double dz = 5.;
    UnplacedTube *uTube = new UnplacedTube(0, rmax, dz, 0, kTwoPi);
    
    Real_v tot_dist=0.;
    
    for (int ipoints=0; ipoints < 10000000; ipoints++)
    {  
      for (auto i = 0; i < 32; i += Real_v::Size) {
         Real_v distance, distance1;
         Vector3D<Real_v> pos(Real_v(points.x() + i), Real_v(points.y() + i), Real_v(points.z() + i));
         Vector3D<Real_v> dir(Real_v(dirs.x() + i), Real_v(dirs.y() + i), Real_v(dirs.z() + i));
__asm__ __volatile__ ("TUBE_START:") ;
         TubeImplementation<translation::kIdentity, rotation::kIdentity, TubeTypes::UniversalTube>::DistanceToIn<Backend>(
              *static_cast<UnplacedTube const *>(uTube), Transformation3D::kIdentity, pos, dir, kInfinity, distance1);
__asm__ __volatile__ ("TUBE_END:"); 
         tot_dist += distance1;
       }
     }
     std::cout<<tot_dist[0]<<std::endl;
    return 0;
 
}
