/*
 * XRayBenchmark.cpp
 *
 * this benchmark performs an X-Ray scan of a (logical volume
 * in a) detector
 *
 * the benchmark stresses the distance functions of the volumes as well as
 * the basic higher level navigation functionality
 */

#include "management/GeoManager.h"
#include "volumes/LogicalVolume.h"
#include "volumes/PlacedBox.h"

#include "base/Global.h"
#include "base/Vector3D.h"
#include "base/Stopwatch.h"
#include "navigation/SimpleNavigator.h"
#include "base/Transformation3D.h"
#include "base/SOA3D.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <map>
#include <cassert>
#include <sstream>

#ifdef VECGEOM_OPENMP
 #include "omp.h"
#endif

#ifdef CALLGRIND
#include "base/callgrind.h"
#endif


#undef NDEBUG

#define VERBOSE false //true or false
#define WRITE_FILE_NAME "volumeImage.bmp" // output image name

using namespace vecgeom;

#pragma pack(push, 1)

typedef struct tFILE_HEADER
{
  unsigned short bfType;
  unsigned long bfSize;
  unsigned short bfReserved1;
  unsigned short bfReserved2;
  unsigned long bfOffBits;
} FILE_HEADER;

#pragma pack(pop)

typedef struct tINFO_HEADER
{
   unsigned long biSize;
   unsigned long biWidth;
   unsigned long biHeight;
   unsigned short biPlanes;
   unsigned short biBitCount;
   unsigned long biCompression;
   unsigned long biSizeImage;
   unsigned long biXPelsPerMeter;
   unsigned long biYPelsPerMeter;
   unsigned long biClrUsed;
   unsigned long biClrImportant;
} INFO_HEADER;

typedef struct tMY_BITMAP
{
  FILE_HEADER  bmpFileHeader;
  INFO_HEADER  bmpInfoHeader;
  unsigned char* bmpPalette;
  unsigned char* bmpRawData;
} MY_BITMAP;



// produce a bmp image out of pixel information given in volume_results
int make_bmp_header( );
int make_bmp(int const * image_result, char const *, int data_size_x, int data_size_y, bool linear = true);
int make_diff_bmp(int const * image1, int const * image2, char const *, int sizex, int sizey);


template<typename Nav_t = SimpleNavigator>
void XRayWithVecGeom(int axis,
                  Vector3D<Precision> origin,
                  Vector3D<Precision> bbox,
                  Vector3D<Precision> dir,
                  double axis1_start, double axis1_end,
                  double axis2_start, double axis2_end,
                  int data_size_x,
                  int data_size_y,
                  double pixel_axis,
                  int * image) {



if(VERBOSE){
    std::cout << "from [" << axis1_start << ";" << axis2_start
              << "] to [" << axis1_end   << ";" << axis2_end << "]\n";
    std::cout << "Xpixels " << data_size_x << " YPixels " << data_size_y << "\n";

    std::cout << pixel_axis << "\n";
}
    double pixel_width_1 = (axis1_end-axis1_start)/data_size_x;
    double pixel_width_2 = (axis2_end-axis2_start)/data_size_y;

    if(VERBOSE){
        std::cout << pixel_width_1 << "\n";
        std::cout << pixel_width_2 << "\n";
    }

    
  #ifdef VECGEOM_OPENMP 
   int Nthreads = omp_get_max_threads();
   NavigationState * newnavstateArray[Nthreads];
   NavigationState * curnavstateArray[Nthreads];
   for(size_t index=0;index<Nthreads;++index){
     std::cerr << "allocating for thread " << index << "\n";
     newnavstateArray[index]=NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );
     curnavstateArray[index]=NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );
   }
   #endif 

#pragma omp parallel for collapse(2) schedule(dynamic)
    for( int pixel_count_2 = 0; pixel_count_2 < data_size_y; ++pixel_count_2 ){
        for( int pixel_count_1 = 0; pixel_count_1 < data_size_x; ++pixel_count_1 )
        {



             #ifdef VECGEOM_OPENMP 
	         size_t threadid=omp_get_thread_num();
	  	 NavigationState * newnavstate = newnavstateArray[threadid];
	  	 NavigationState * curnavstate = curnavstateArray[threadid];
             #else	
     	         NavigationState * newnavstate = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );
		 NavigationState * curnavstate = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );
	     #endif
	  double axis2_count = axis2_start + pixel_count_2 * pixel_width_2 + 1E-6;
            double axis1_count = axis1_start + pixel_count_1 * pixel_width_1 + 1E-6;

            if(VERBOSE) {
                std::cout << "\n OutputPoint("<< axis1_count<< ", "<< axis2_count<< ")\n";
            }

            // set start point of XRay
            Vector3D<Precision> p;
            if( axis== 1 )
              p.Set( origin[0]-bbox[0], axis1_count, axis2_count);
            else if( axis== 2)
              p.Set( axis1_count, origin[1]-bbox[1], axis2_count);
            else if( axis== 3)
              p.Set( axis1_count, axis2_count, origin[2]-bbox[2]);

            SimpleNavigator nav;
            curnavstate->Clear();
            nav.LocatePoint( GeoManager::Instance().GetWorld(), p, *curnavstate, true );


            double distancetravelled=0.;
            int crossedvolumecount=0;

            if(VERBOSE) {
              std::cout << " StartPoint(" << p[0] << ", " << p[1] << ", " << p[2] << ")";
              std::cout << " Direction <" << dir[0] << ", " << dir[1] << ", " << dir[2] << ">"<< std::endl;
            }

            while( ! curnavstate->IsOutside() ) {
                double step = 0;
                newnavstate->Clear();
                Nav_t navigator;
                navigator.FindNextBoundaryAndStep( p,
                        dir,
                        *curnavstate,
                        *newnavstate,
                        vecgeom::kInfinity, step);

                //std::cout << "step " << step << "\n";
                distancetravelled+=step;

                // here we have to propagate particle ourselves and adjust navigation state
                p = p + dir*(step + 1E-6);


                newnavstate->CopyTo(curnavstate);

                // Increase passed_volume
                // TODO: correct counting of travel in "world" bounding box
                if(step>0) crossedvolumecount++;

             } // end while

             ///////////////////////////////////
             // Store the number of passed volume at 'volume_result'
             *(image+pixel_count_2*data_size_x+pixel_count_1) = crossedvolumecount;
             #ifndef VECGEOM_OPENMP 
	         NavigationState::ReleaseInstance( curnavstate );
	         NavigationState::ReleaseInstance( newnavstate );
             #endif  

      } // end inner loop
    } // end outer loop

  #ifdef VECGEOM_OPENMP 
    for(size_t index=0;index<Nthreads;++index){
      NavigationState::ReleaseInstance( curnavstateArray[index] );
      NavigationState::ReleaseInstance( newnavstateArray[index] );
    }
  #endif
} // end XRayWithVecGeom


// stressing the vector interface of navigator
void XRayWithVecGeom_VecNav(int axis,
                  Vector3D<Precision> origin,
                  Vector3D<Precision> bbox,
                  Vector3D<Precision> dir,
                  double axis1_start, double axis1_end,
                  double axis2_start, double axis2_end,
                  int data_size_x,
                  int data_size_y,
                  double pixel_axis,
                  int * image) {

if(VERBOSE){
    std::cout << "from [" << axis1_start << ";" << axis2_start
              << "] to [" << axis1_end   << ";" << axis2_end << "]\n";
    std::cout << "Xpixels " << data_size_x << " YPixels " << data_size_y << "\n";

    std::cout << pixel_axis << "\n";
}
    double pixel_width_1 = (axis1_end-axis1_start)/data_size_x;
    double pixel_width_2 = (axis2_end-axis2_start)/data_size_y;

    if(VERBOSE){
        std::cout << pixel_width_1 << "\n";
        std::cout << pixel_width_2 << "\n";
    }

    // we need N navstates ( where N should be a multiple of the SIMD width )
    unsigned int N = 512;
    NavigationState ** newnavstates = new NavigationState*[N];
    NavigationState ** curnavstates = new NavigationState*[N];
    for( unsigned int j=0;j<N;++j ){
        newnavstates[j] = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );
        curnavstates[j] = NavigationState::MakeInstance( GeoManager::Instance().getMaxDepth() );
    }

    SOA3D<Precision> points(N);
    SOA3D<Precision> dirs(N);
    SOA3D<Precision> workspaceforlocalpoints(N);
    SOA3D<Precision> workspaceforlocaldirs(N);

    // initialize dirs from dir
    for( unsigned int j=0; j<N; ++j )
        dirs.set(j, dir.x(), dir.y(),dir.z());

    double * steps    = new double[N];
    double * psteps   = new double[N];
    double * safeties = new double[N];
    int * nextnodeworkspace = new int[N]; // some workspace for the navigator; not important here
    // initialize physical steps to infinity
    for(unsigned int j=0;j<N;++j)
        psteps[j]=vecgeom::kInfinity;

    for( int pixel_count_2 = 0; pixel_count_2 < data_size_y; ++pixel_count_2 ){
        for( int pixel_count_1 = 0; pixel_count_1 < data_size_x; ++pixel_count_1 )
        {
            double axis1_count = axis1_start + pixel_count_1 * pixel_width_1 + 1E-6;
            double axis2_count = axis2_start + pixel_count_2 * pixel_width_2 + 1E-6;

            if(VERBOSE) {
                std::cout << "\n OutputPoint("<< axis1_count<< ", "<< axis2_count<< ")\n";
            }
            // std::cout << pixel_count_1 << " " << pixel_count_2 << "\n";
            // set start points of XRay; points should be in a SOA/AOS

            SimpleNavigator nav;
            // initialize points and locate them is serialized
            for( unsigned int j=0; j<N; ++j ){

                if( axis== 1 )
                  points.set( j, origin[0]-bbox[0], axis1_count, axis2_count );
                else if( axis== 2)
                  points.set( j, axis1_count, origin[1]-bbox[1], axis2_count );
                else if( axis== 3)
                  points.set( j, axis1_count, axis2_count, origin[2]-bbox[2] );

                curnavstates[j]->Clear();
                nav.LocatePoint( GeoManager::Instance().GetWorld(), points[j], *curnavstates[j], true );
            }
            double distancetravelled=0.;
            int crossedvolumecount=0;
            if(VERBOSE) {
              std::cout << " StartPoint(" << points[0].x() << ", " << points[1].y() << ", " <<
                      points[2].z() << ")";
              std::cout << " Direction <" << dirs[0].x() << ", " << dirs[1].y() << ", " << dirs[2].z() << ">"<< std::endl;
            }

            // we do the while loop only over the first "particle index"
            // the rest of the particles should follow exactly the same path
            while( ! curnavstates[0]->IsOutside() ) {
                nav.FindNextBoundaryAndStep(
                        points,
                        dirs,
                        workspaceforlocalpoints,
                        workspaceforlocaldirs,
                        curnavstates,
                        newnavstates,
                        psteps,
                        safeties,
                        steps,
                        nextnodeworkspace);

                std::cout << "step " << steps[0] << "\n";
                distancetravelled+=steps[0];

                // TODO: DO HERE AN ASSERTION THAT ALL STEPS AGREE

//
//              std::cout << "GOING FROM "
//                       << curnavstate->Top()->GetLabel() << "(";
//                        curnavstate->Top()->PrintType();
//                     std::cout << ") to ";
//
//                if( newnavstate->Top() ){
//                    std::cout << newnavstate->Top()->GetLabel() << "(";
//                    newnavstate->Top()->PrintType();
//                    std::cout << ")";
//                }
//                else
//                    std::cout << "outside ";
//
//                if ( curnavstate->Top() == newnavstate->Top() ) {
//                    std::cout << " CROSSING PROBLEM \n";
//                    curnavstate->Print();
//                    newnavstate->Print();
//                }
//
//                std::cout << "# step " << step << " crossed "<< crossedvolumecount << "\n";

                // here we have to propagate particle ourselves and adjust navigation state
                // propagate points
                for(unsigned int j=0;j<N;++j){
                    points.set(j, points[j] + dirs[j]*(steps[0] + 1E-6));
                    newnavstates[j]->CopyTo(curnavstates[j]);
                }

                // Increase passed_volume
                // TODO: correct counting of travel in "world" bounding box
                if(steps[0]>0) crossedvolumecount++;
             } // end while

             ///////////////////////////////////
             // Store the number of passed volume at 'volume_result'
             *(image+pixel_count_2*data_size_x+pixel_count_1) = crossedvolumecount;

      } // end inner loop
    } // end outer loop

    for( unsigned int j=0; j<N ; ++j ) {
        NavigationState::ReleaseInstance( curnavstates[j] );
        NavigationState::ReleaseInstance( newnavstates[j] );
    }
} // end XRayWithVecGeomVectorInterface



//////////////////////////////////
// main function
int main(int argc, char * argv[])
{
  int axis= 0;

  double axis1_start= 0.;
  double axis1_end= 0.;

  double axis2_start= 0.;
  double axis2_end= 0.;

  double pixel_width= 0;
  double pixel_axis= 1.;

  if( argc < 3 )
  {
    std::cerr<< std::endl;
    std::cerr<< "Need to give volumename, axis and number of axis"<< std::endl;
    std::cerr<< "USAGE : ./XRayBenchmark [VolumeName] [ViewDirection(Axis)]"
             << "[PixelWidth(OutputImageSize)]"
             << std::endl;
    std::cerr<< "  ex) ./XRayBenchmark BSCTrap y 95"<< std::endl;
    return 1;
  }

  GeoManager &gGeoManager = GeoManager::Instance();
  gGeoManager.LoadGeometryFromSharedLib("libCMSGeom.so");

  std::string testvolume( argv[1] );

  if( strcmp(argv[2], "x")==0 )
    axis= 1;
  else if( strcmp(argv[2], "y")==0 )
    axis= 2;
  else if( strcmp(argv[2], "z")==0 )
    axis= 3;
  else
  {
    std::cerr<< "Incorrect axis"<< std::endl<< std::endl;
    return 1;
  }

  std::stringstream imagenamebase;
  imagenamebase << "volumeImage_" << testvolume;
  if(axis==1) imagenamebase << "x"; 
  if(axis==2) imagenamebase << "y"; 
  if(axis==3) imagenamebase << "z"; 

  pixel_width= atof(argv[3]);
/*
  for(auto i= 3; i< argc; i++)
  {
    if( ! strcmp(argv[i], "--usolids") )
      usolids= true;
    if( ! strcmp(argv[i], "--vecgeom") )
      usolids= false;
    if( ! strcmp(argv[i], "--novoxel") )
      voxelize = false;
  }
*/

  std::cout<<argv[1]<<std::endl;
  VPlacedVolume * foundvolume = gGeoManager.FindPlacedVolume(argv[1]);
  LogicalVolume * foundlvolume = gGeoManager.FindLogicalVolume(argv[1]);
/*
  std::vector<vecgeom::VPlacedVolume *> v2; 
  gGeoManager.getAllPlacedVolumes( v2 );
  std::cout<<"Have placed volumes %ld\n"<<v2.size()<<std::endl;
 
  for (int i=0; i<v2.size();i++)
  {
     std::cout<<"Volume label "<< (v2[i])->GetLabel()<<std::endl; 
  }
*/
  if( foundvolume ) {
 
  //std::cout<<"Volume found"<<std::endl;

  double dx = 1.;
  double dy = 1.;
  double dz = 1.;

  Vector3D<Precision> minExtent,maxExtent;
  foundvolume->Extent(minExtent,maxExtent);
  Vector3D<Precision> delta = maxExtent -minExtent;
  dx = delta.x();
  dy = delta.y();
  dz = delta.z();
  
  if (VERBOSE) {
     std::cout<<dx<<" "<<dy<<" "<<dz<<std::endl;
     std::cout<<"minExtent"<<minExtent.x()<<" "<< minExtent.y()<<" "<<minExtent.z()<<std::endl;
     std::cout<<"maxExtent"<<maxExtent.x()<<" "<< maxExtent.y()<<" "<<maxExtent.z()<<std::endl;
  }

  //const Vector3D<Precision> vv = foundvolume->bounding_box()->dimensions();
  UnplacedBox worldUnplaced = UnplacedBox(delta);

  LogicalVolume newworld = LogicalVolume("newworld", &worldUnplaced);


  newworld.PlaceDaughter(foundvolume);

  VPlacedVolume *newworldPlaced = newworld.Place();

  gGeoManager.SetWorld(newworldPlaced);


  double origin[3]= {0.,};
  Vector3D<Precision> orig = (maxExtent + minExtent)/2.;
  origin[0] = orig.x();
  origin[1] = orig.y();
  origin[2] = orig.z();



    Vector3D<Precision> p;
    Vector3D<Precision> dir;
    
    if(axis== 1)
    {
      dir.Set(1., 0., 0.);
      //Transformation3D trans( 0, 0, 0, 5, 5, 5);
      //trans.Print();
     // dir = trans.TransformDirection( Vector3D<Precision> (1,0,0));

      axis1_start= origin[1]-dy;
      axis1_end= origin[1]+ dy;
      axis2_start= origin[2]-dz;
      axis2_end= origin[2]+ dz;
      pixel_axis= (dy*2)/pixel_width;
    }
    else if(axis== 2)
    {
      dir.Set(0., 1., 0.);
      //vecgeom::Transformation3D trans( 0, 0, 0, 5, 5, 5);
      //dir = trans.TransformDirection(dir);
      axis1_start= origin[0]-dx;
      axis1_end= origin[0]+ dx;
      axis2_start= origin[2]-dz;
      axis2_end= origin[2]+ dz;
      pixel_axis= (dx*2)/pixel_width;
    }
    else if(axis== 3)
    {
      dir.Set(0., 0., 1.);
      //vecgeom::Transformation3D trans( 0, 0, 0, 5, 5, 5);
      //dir = trans.TransformDirection(dir);
      axis1_start= origin[0] -dx;
      axis1_end= origin[0]+ dx;
      axis2_start= origin[1]-dy;
      axis2_end= origin[1]+ dy;
      pixel_axis= (dx*2)/pixel_width;
    }

    // init data for image
    int data_size_x= (axis1_end-axis1_start)/pixel_axis;
    int data_size_y= (axis2_end-axis2_start)/pixel_axis;

    int *volume_result_VecGeom= (int*) new int[data_size_y * data_size_x*3];


    gGeoManager.CloseGeometry();

    Stopwatch timer;

    timer.Start();
    //XRayWithVecGeom_VecNav( axis,

    XRayWithVecGeom<SimpleNavigator>( axis,
               Vector3D<Precision>(origin[0],origin[1],origin[2]),
               Vector3D<Precision>(dx,dy,dz),
               dir,
               axis1_start, axis1_end,
               axis2_start, axis2_end,
               data_size_x, data_size_y,
               pixel_axis,
               volume_result_VecGeom );
   timer.Stop();

    std::stringstream VecGeomimage;
    VecGeomimage << imagenamebase.str();
    VecGeomimage << "_VecGeom.bmp";
    make_bmp(volume_result_VecGeom, VecGeomimage.str().c_str(), data_size_x, data_size_y);


    std::cout << std::endl;
    std::cout << " VecGeom Elapsed time : "<< timer.Elapsed() << std::endl;
  }
  else {
    std::cout<<"volume not found"<<std::endl;
  }
  return 0;
}


void make_bmp_header( MY_BITMAP * pBitmap, unsigned char * bmpBuf, int sizex, int sizey )
{
  int width_4= (sizex+ 3)&~3;
  unsigned int len= 0;

  // bitmap file header
  pBitmap->bmpFileHeader.bfType=0x4d42;
  pBitmap->bmpFileHeader.bfSize=sizey* width_4* 3+ 54;
  pBitmap->bmpFileHeader.bfReserved1= 0;
  pBitmap->bmpFileHeader.bfReserved2= 0;
  pBitmap->bmpFileHeader.bfOffBits= 54;

  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfType, 2);
  len+= 2;
  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfSize, 4);
  len+= 4;
  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfReserved1, 2);
  len+= 2;
  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfReserved2, 2);
  len+= 2;
  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfOffBits, 4);
  len+= 4;

  // bitmap information header
  pBitmap->bmpInfoHeader.biSize= 40;
  pBitmap->bmpInfoHeader.biWidth= width_4;
  pBitmap->bmpInfoHeader.biHeight= sizey;
  pBitmap->bmpInfoHeader.biPlanes= 1;
  pBitmap->bmpInfoHeader.biBitCount= 24;
  pBitmap->bmpInfoHeader.biCompression= 0;
  pBitmap->bmpInfoHeader.biSizeImage= sizey* width_4* 3;
  pBitmap->bmpInfoHeader.biXPelsPerMeter= 0;
  pBitmap->bmpInfoHeader.biYPelsPerMeter= 0;
  pBitmap->bmpInfoHeader.biClrUsed= 0;
  pBitmap->bmpInfoHeader.biClrImportant=0;


  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biSize, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biWidth, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biHeight, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biPlanes, 2);
  len+= 2;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biBitCount, 2);
  len+= 2;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biCompression, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biSizeImage, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biXPelsPerMeter, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biYPelsPerMeter, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biClrUsed, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biClrImportant, 4);
  len+= 4;
}


int make_bmp(int const * volume_result, char const * name, int data_size_x, int data_size_y, bool linear )
{

  MY_BITMAP* pBitmap= new MY_BITMAP;
  FILE *pBitmapFile;
  int width_4= (data_size_x+ 3)&~3;
  unsigned char* bmpBuf;

  bmpBuf = (unsigned char*) new unsigned char[data_size_y* width_4* 3+ 54];
  printf("\n Write bitmap...\n");

  unsigned int len= 0;

  // bitmap file header
  pBitmap->bmpFileHeader.bfType=0x4d42;
  pBitmap->bmpFileHeader.bfSize=data_size_y* width_4* 3+ 54;
  pBitmap->bmpFileHeader.bfReserved1= 0;
  pBitmap->bmpFileHeader.bfReserved2= 0;
  pBitmap->bmpFileHeader.bfOffBits= 54;
  
  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfType, 2);
  len+= 2;
  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfSize, 4);
  len+= 4;
  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfReserved1, 2);
  len+= 2;
  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfReserved2, 2);
  len+= 2;
  memcpy(bmpBuf + len, &pBitmap->bmpFileHeader.bfOffBits, 4);
  len+= 4;

  // bitmap information header
  pBitmap->bmpInfoHeader.biSize= 40;
  pBitmap->bmpInfoHeader.biWidth= width_4;
  pBitmap->bmpInfoHeader.biHeight= data_size_y;
  pBitmap->bmpInfoHeader.biPlanes= 1;
  pBitmap->bmpInfoHeader.biBitCount= 24;
  pBitmap->bmpInfoHeader.biCompression= 0;
  pBitmap->bmpInfoHeader.biSizeImage= data_size_y* width_4* 3;
  pBitmap->bmpInfoHeader.biXPelsPerMeter= 0;
  pBitmap->bmpInfoHeader.biYPelsPerMeter= 0;
  pBitmap->bmpInfoHeader.biClrUsed= 0;
  pBitmap->bmpInfoHeader.biClrImportant=0;


  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biSize, 4); 
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biWidth, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biHeight, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biPlanes, 2);
  len+= 2;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biBitCount, 2);
  len+= 2;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biCompression, 4); 
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biSizeImage, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biXPelsPerMeter, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biYPelsPerMeter, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biClrUsed, 4);
  len+= 4;
  memcpy(bmpBuf+len, &pBitmap->bmpInfoHeader.biClrImportant, 4);
  len+= 4;

  // find out maxcount before doing the picture
  int maxcount = 0;
  int x=0,y=0,origin_x=0;
  while( y< data_size_y )
  {
     while( origin_x< data_size_x )
     {
       int value = *(volume_result+y*data_size_x+origin_x);
       maxcount = ( value > maxcount )? value : maxcount;

       x++;
       origin_x++;
     }
     y++;
     x = 0;
     origin_x = 0;
  }
//  maxcount = std::log(maxcount + 1);

  x= 0;
  y= 0;
  origin_x= 0;

  int padding= width_4- data_size_x;
  int padding_idx= padding;
  unsigned char *imgdata= (unsigned char*) new unsigned char[data_size_y*width_4*3];

  int totalcount = 0;

  while( y< data_size_y )
  {
    while( origin_x< data_size_x )
    {
      int value = *(volume_result+y*data_size_x+origin_x);
      totalcount += value;

      //*(imgdata+y*width_4*3+x*3+0)= (value *50) % 256;
      //*(imgdata+y*width_4*3+x*3+1)= (value *40) % 256;
      //*(imgdata+y*width_4*3+x*3+2)= (value *30) % 256;

      //*(imgdata+y*width_4*3+x*3+0)= (std::log(value)/(1.*maxcount)) * 256;
      //*(imgdata+y*width_4*3+x*3+1)= (std::log(value)/(1.2*maxcount)) * 256;
      //*(imgdata+y*width_4*3+x*3+2)= (std::log(value)/(1.4*maxcount)) * 256;
if( linear ){
      *(imgdata+y*width_4*3+x*3+0)= (value/(1.*maxcount)) * 256;
      *(imgdata+y*width_4*3+x*3+1)= (value/(1.*maxcount)) * 256;
      *(imgdata+y*width_4*3+x*3+2)= (value/(1.*maxcount)) * 256;
}
else {
          *(imgdata+y*width_4*3+x*3+0)= (log(value+1))/(1.*(log(1+maxcount))) * 256;
          *(imgdata+y*width_4*3+x*3+1)= (log(value+1))/(1.*(log(1+maxcount))) * 256;
          *(imgdata+y*width_4*3+x*3+2)= (log(value+1))/(1.*(log(1+maxcount))) * 256;
}
  x++;
  origin_x++;

  while( origin_x== data_size_x && padding_idx)
      {
// padding 4-byte at bitmap image
        *(imgdata+y*width_4*3+x*3+0)= 0;
        *(imgdata+y*width_4*3+x*3+1)= 0;
        *(imgdata+y*width_4*3+x*3+2)= 0;
        x++;
        padding_idx--;
      }
      padding_idx= padding;
    }
    y++;
    x= 0;
    origin_x= 0;
  }
  
  memcpy(bmpBuf + 54, imgdata, width_4* data_size_y* 3);

  pBitmapFile = fopen(name, "wb");
  fwrite(bmpBuf, sizeof(char), width_4*data_size_y*3+54, pBitmapFile);


  fclose(pBitmapFile);
  delete[] imgdata;
  delete[] bmpBuf;
  delete pBitmap;

  std::cout << " wrote image file " << name <<  "\n";
  std::cout << " total count " << totalcount << "\n";
  std::cout << " max count " << maxcount << "\n";
  return 0;
}



int make_diff_bmp(int const * image1, int const * image2, char const * name, int data_size_x, int data_size_y )
{

  MY_BITMAP* pBitmap= new MY_BITMAP;
  FILE *pBitmapFile;
  int width_4 = (data_size_x + 3)&~3;
  unsigned char* bmpBuf = (unsigned char*) new unsigned char[data_size_y* width_4* 3+ 54];

  // init buffer and write header
  make_bmp_header(pBitmap, bmpBuf, data_size_x, data_size_y);

  std::cerr << "AFTER HEADER " << bmpBuf << "\n";

  // TODO: verify the 2 images have same dimensions

  // find out maxcount before doing the picture
  int maxdiff = 0;
  int mindiff = 0;
  int x=0,y=0,origin_x=0;
  while( y< data_size_y )
  {
     while( origin_x< data_size_x )
     {
       int value = *(image1+y*data_size_x+origin_x) - *(image2+y*data_size_x + origin_x);
       maxdiff = ( value > maxdiff )? value : maxdiff;
       mindiff = ( value < mindiff )? value : mindiff;
       x++;
       origin_x++;
     }
     y++;
     x = 0;
     origin_x = 0;
  }


  x= 0;
  y= 0;
  origin_x= 0;

  int padding= width_4 - data_size_x;
  int padding_idx= padding;
  unsigned char *imgdata= (unsigned char*) new unsigned char[data_size_y*width_4*3];

  while( y< data_size_y )
  {
    while( origin_x< data_size_x )
    {
      int value = *(image1 + y*data_size_x + origin_x) - *(image2 + y*data_size_x + origin_x);

      if( value >=0 ){
          *(imgdata+y*width_4*3+x*3+0)= (value/(1.*maxdiff)) * 256;
          *(imgdata+y*width_4*3+x*3+1)= 0;// (value/(1.*maxcount)) * 256;
          *(imgdata+y*width_4*3+x*3+2)= 0;//(value/(1.*maxcount)) * 256;}
      }
      else
      {
          *(imgdata+y*width_4*3+x*3+0)= 0;
          *(imgdata+y*width_4*3+x*3+1)= 0;// (value/(1.*maxcount)) * 256;
          *(imgdata+y*width_4*3+x*3+2)= (value/(1.*mindiff)) * 256;//(value/(1.*maxcount)) * 256;}
      }
      x++;
      origin_x++;

      while( origin_x== data_size_x && padding_idx)
      {
        // padding 4-byte at bitmap image
        *(imgdata+y*width_4*3+x*3+0)= 0;
        *(imgdata+y*width_4*3+x*3+1)= 0;
        *(imgdata+y*width_4*3+x*3+2)= 0;
        x++;
        padding_idx--;
      }
      padding_idx= padding;
    }
    y++;
    x= 0;
    origin_x= 0;
  }

  memcpy(bmpBuf + 54, imgdata, width_4* data_size_y* 3);
  pBitmapFile = fopen(name, "wb");
  fwrite(bmpBuf, sizeof(char), width_4*data_size_y*3+54, pBitmapFile);

  fclose(pBitmapFile);
  delete[] imgdata;
  delete[] bmpBuf;
  delete pBitmap;

  std::cout << " wrote image file " << name <<  "\n";
  return 0;
}
