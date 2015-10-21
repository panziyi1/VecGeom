/// author: Andrei Gheata (andrei.gheata@cern.ch) Oct 14, 2015
#include <sstream>
#include <string>
#include "volumes/Tube.h"
#include "volumes/LogicalVolume.h"
#include "management/GeoManager.h"
#include "base/Stopwatch.h"
#include "SC15Navigators.h"
#include "navigation/SimpleNavigator.h"
#include "base/Vector3D.h"
#include "base/SOA3D.h"
#include <iostream>

#ifdef VECGEOM_OPENMP
 #include "omp.h"
#endif

using namespace vecgeom;

#define ALIGN_PADDING kAlignmentBoundary

#define SETINNERNAV(depth) \
  if(layer==depth) vol->SetUserExtensionPtr( (void*) InnerMostTubeNavigator<depth+1>::Instance() );

#define SETLAYERNAV(depth) \
  if(layer==depth) vol->SetUserExtensionPtr( (void*) LayerNavigator<depth+1>::Instance() );

#define SETWORLDNAV(depth) \
  if(layer==0) vol->SetUserExtensionPtr( (void*) WorldNavigator<true>::Instance() );

#define WRITE_FILE_NAME "volumeImage.bmp" // output image name

static VNavigator const* GetNavigator(LogicalVolume const*lvol){
  return (VNavigator const*)lvol->GetUserExtensionPtr();
}

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

constexpr unsigned int kMaxVec = 256;
constexpr unsigned int kMaxDepth = 20;

inline int round_up_align(int num) {
  int remainder = num % ALIGN_PADDING;
  if (remainder == 0) return num;
  return (num + ALIGN_PADDING - remainder);
}
inline int VolIndex(const NavigationState *state) { return (state->GetCurrentLevel() - 1); }
// Decode pixels into coordinates
inline void PixelToCoord(unsigned int ix, unsigned int iy, Precision start[3]) {
// to implement
}
inline unsigned int EncodePixel(unsigned int ix, unsigned int iy) { return ( (iy << 16) | ix); }
inline void DecodePixel(unsigned int const pixel, unsigned int &ix, unsigned int &iy) { 
  ix = pixel & 0xFFFF; iy = pixel >> 16;
}  
void make_bmp(int const *, char const *, int, int, bool linear = true);

// This structure handles a pixel window and handles inline conversions between pixel and world
// coordinates.
struct Window
{
  int              pixel_width;    // Number of pixels on each axis
  int              axis_;          // Shooting axis along which the window is seen (1,2 or 3)
  int              nslices;        // Number of sub-window slices
  Vector3D<Precision> dir_;        // Shooting direction
  Precision        offset;         // Offset on the shooting axis
  Precision        axis1_start;    // Start coordinate on first axis
  Precision        axis1_end;      // End coordinate on first axis
  Precision        axis2_start;    // Start coordinate on second axis
  Precision        axis2_end;      // End coordinate on second axis
  
  int              data_size_x;    // Size of image in x (pixels)
  int              data_size_y;    // Size of image in y (pixels)
  Precision        pixel_width_1;  // Width of one pixel on first axis
  Precision        pixel_width_2;  // Width of one pixel on second axis
  int              dslice;         // Width in pixels of a slice window
  
  Window(VPlacedVolume const *vol, int axis, int width, int nsw) : pixel_width(width), axis_(axis-1), nslices(nsw), dir_(),
      offset(0), axis1_start(0), axis1_end(0), axis2_start(0), axis2_end(0), 
      data_size_x(0), data_size_y(0), pixel_width_1(0), pixel_width_2(0), dslice(0) {
    Vector3D<Precision> minExtent,maxExtent;
    vol->Extent(minExtent,maxExtent);
  
    Precision dx = (maxExtent - minExtent).x(); 
    Precision dy = (maxExtent - minExtent).y(); 
    Precision dz = (maxExtent - minExtent).z(); 
  
    Vector3D<Precision> orig = (maxExtent + minExtent)/2.;

    Precision pixel_axis= 1.;
    switch (axis_) {
      case (0):
        dir_.Set(1., 0., 0.);
	     offset = orig.x()-dx;
        axis1_start= orig.y()-dy;
        axis1_end= orig.y() + dy;
        axis2_start= orig.z()-dz;
        axis2_end= orig.z()+ dz;
        pixel_axis= (dy*2)/pixel_width;
        break;
      case (1):
        dir_.Set(0., 1., 0.);
	     offset = orig.y()-dy;
        axis1_start= orig.z()-dz;
        axis1_end= orig.z()+ dz;
        axis2_start= orig.x()-dx;
        axis2_end= orig.x()+ dx;
        pixel_axis= (dz*2)/pixel_width;
        break;
      case (2):
        dir_.Set(0., 0., 1.);
	     offset = orig.z()-dz;
        axis1_start= orig.x() -dx;
        axis1_end= orig.x()+ dx;
        axis2_start= orig.y()-dy;
        axis2_end= orig.y()+ dy;
        pixel_axis= (dx*2)/pixel_width;
        break;
      default:
        std::cout << "Unknown axis" << std::endl;  
    }
    // init data for image
    data_size_x= (axis1_end-axis1_start)/pixel_axis;
    data_size_y= (axis2_end-axis2_start)/pixel_axis;
    pixel_width_1 = (axis1_end-axis1_start)/data_size_x;
    pixel_width_2 = (axis2_end-axis2_start)/data_size_y;    
    assert(data_size_x > nslices);
    dslice = data_size_x/nslices;
    std::cout << "pixel_width= " << pixel_width << std::endl;
    std::cout << "direction= " << dir_ << std::endl;
    std::cout << "offset= " << offset << std::endl;
    std::cout << "data_size_x= " << data_size_x << "  data_size_y= " << data_size_y << std::endl;
    std::cout << "axis1_start= " << axis1_start << "  axis1_end= " << axis1_end << std::endl;
    std::cout << "axis2_start= " << axis2_start << "  axis2_end= " << axis2_end << std::endl;
    std::cout << "pixel_width_1= " << pixel_width_1 << " pixel_width_2= " << pixel_width_2 << std::endl;
  }
  
  void GetSubwindow(int islice, int &i0, int &ni) const {
    // Get sub-window coordinates and widths
    assert(islice < nslices);
    i0 = islice*dslice;
    ni = dslice;
    if ( i0+ni > data_size_x) ni = data_size_x - i0;
  }
  
  inline void GetCoordinates(int i, int j, Vector3D<Precision> &point) const {
    point[axis_] = offset;
    point[(axis_+1)%3] = axis1_start + i * pixel_width_1 + 1E-6;
    point[(axis_+2)%3] = axis2_start + j * pixel_width_2 + 1E-6;
  }
};

// This structure contains track data propagated after one step
struct TrackState
{
  unsigned int     bindex;
  unsigned int     pixel;
  unsigned int     nbound;
  Precision        xp;
  Precision        yp;
  Precision        zp;
  
  inline TrackState(unsigned int b, unsigned int pix, unsigned int nb, Precision x, Precision y, Precision z) {
    bindex = b; pixel = pix; nbound = nb; xp = x; yp = y; zp = z; }
};  

class Stack {
public:
  std::vector<TrackState> states_;
  Stack(int size) : states_() { states_.reserve(size); }
  
  inline bool empty() const { return states_.empty(); }
  inline void Push(unsigned int b, unsigned int pixel, unsigned int nbound, Precision xp, Precision yp, Precision zp) {
    states_.emplace_back(b, pixel, nbound, xp, yp, zp);
  }  
  inline TrackState &Pop() { TrackState &tmp = states_.back(); states_.pop_back(); return tmp;}
};
  
// Each stepper will create an array of baskets (one per volume)
class Stepper;

class Basket
{
public:
  int              size_;             // current size of the basket
  int              threshold_;        // threshold to start transport
  bool             ready_;            // signal that basket is ready for transport
  char            *fBuf;              // Buffer holding tracks data
  // Arrays of size kMaxVec aligned in buffer
  unsigned int    *pixel;             // pixel X/Y coordinates encoded in a single long
  unsigned int    *nbound;            // Number of crossed boundaries
  Precision       *xp;                // current coordinates for tracks in the basket
  Precision       *yp;                //
  Precision       *zp;                //
  Precision       *xdir;              // current directions - now all the same for all baskets
  Precision       *ydir;              //
  Precision       *zdir;              //
  NavStatePool    *navstates_;        // navigation states for tracks - now all the same for the same basket
  VNavigator const *specialnav_;      // Special navigator for this basket(volume)
  int             *pixmap_;           // Pixel map to be filled
  const Window    *window_;           // Shooting window
  
  Basket(int thres) : size_(0), threshold_(thres), ready_(false) {
     int size_max = round_up_align(kMaxVec); // multiple of ALIGN_PADDING - all sub-arrays automatically aligned
     int bufsize = size_max*(2*sizeof(unsigned int) + 6*sizeof(Precision));
     fBuf = (char *)_mm_malloc(bufsize, ALIGN_PADDING);
     memset(fBuf, 0, bufsize);
     char *buf = fBuf;
     pixel = (unsigned int *)buf;
     buf += size_max*sizeof(unsigned int);
     nbound = (unsigned int *)buf;
     buf += size_max*sizeof(unsigned int);
     xp = (Precision *)buf;
     buf += size_max*sizeof(Precision);
     yp = (Precision *)buf;
     buf += size_max*sizeof(Precision);
     zp = (Precision *)buf;
     buf += size_max*sizeof(Precision);
     xdir = (Precision *)buf;
     buf += size_max*sizeof(Precision);
     ydir = (Precision *)buf;
     buf += size_max*sizeof(Precision);
     zdir = (Precision *)buf;
     buf += size_max*sizeof(Precision);
     navstates_ = new NavStatePool(thres, GeoManager::Instance().getMaxDepth());
     specialnav_ = nullptr;
  }
  
  ~Basket() {
    _mm_free(fBuf);
    delete navstates_;
  }
  
  inline void FillPixmap(int ix, int iy, int nbound) {
    // Fill a pixel
    *(pixmap_+iy*window_->data_size_x+ix) = nbound;
  }

  void FlushTrack(int itr) {
    // Track fully transported - collect its "signal" in the output
    unsigned int px, py;
    DecodePixel(pixel[itr], px, py);
    FillPixmap(px, py, nbound[itr]+1);
  }

  bool AddTrack(TrackState const &state) {
    assert (!ready_);
    pixel[size_] = state.pixel;
    xp[size_] = state.xp;
    yp[size_] = state.yp;
    zp[size_] = state.zp;
    nbound[size_] = state.nbound;
    // No need to copy direction
    size_++;
    if (size_ == threshold_) ready_ = true;
    return ready_;
  }
  
  void SetDirection(const Vector3D<Precision> &dir) { 
    for (auto i=0; i<threshold_; ++i) {
      xdir[i] = dir.x(); ydir[i] = dir.y(), zdir[i] = dir.z();
    }
  }
  
  void SetNavStates(const NavigationState *state) {
    specialnav_ = GetNavigator(state->Top()->GetLogicalVolume());
    for (auto i=0; i<threshold_; ++i) state->CopyTo((*navstates_)[i]);
  }
  
  void Step(Precision *psteps, Precision *steps, NavStatePool *newnavstates, Stack &stack) {
    // Perform a single step and pushes the new states to the stack, flushing exiting tracks
//    VNavigator const *specialnav = GetNavigator((*navstates_)[0]->Top()->GetLogicalVolume());
    specialnav_->ComputeStepsAndPropagatedStates(SOA3D<Precision>(xp,yp,zp,size_), SOA3D<Precision>(xdir,ydir,zdir,size_), 
        psteps, *navstates_, *newnavstates, steps);

    // Push all tracks to stack using construct in place
    for(auto itr=0; itr<size_; ++itr) {
      NavigationState *state = (*newnavstates)[itr];
      if (state->IsOutside()) FlushTrack(itr);
      else stack.Push(VolIndex(state), pixel[itr], nbound[itr]+1, 
                 xp[itr] + xdir[itr]*steps[itr] + 1E-6,
                 yp[itr] + ydir[itr]*steps[itr] + 1E-6,
                 zp[itr] + zdir[itr]*steps[itr] + 1E-6);
    }
    // Reset the basket content
    ready_ = false;
    size_ = 0;
  }  
};

/// Each thread will create a stepper
class Stepper
{
public:
  int              vecsize_;          // Vector size (most baskets will get transported with this vector size)
  int              nvolumes_;         // number of geometry volumes
  Stack            stack_;            // Stack of transported track states (per step)
  Precision       *psteps_;           // proposed steps
  Precision       *steps_;            // current steps
  Basket          *baskets_[kMaxDepth];  // Array of baskets (one per volume)
  NavStatePool    *outstates_;        // navigation states for crossing tracks
  
  Stepper(int vsize, const Window* window, int* pixmap) : vecsize_(vsize), nvolumes_(0), stack_(vsize*kMaxDepth),
        psteps_(nullptr), steps_(nullptr) {
    psteps_ = (Precision*)_mm_malloc(round_up_align(kMaxVec)*sizeof(Precision), ALIGN_PADDING);
    steps_ = (Precision*)_mm_malloc(round_up_align(kMaxVec)*sizeof(Precision), ALIGN_PADDING);
    NavigationState *state = NavigationState::MakeInstance(GeoManager::Instance().getMaxDepth());
    VPlacedVolume const *node = GeoManager::Instance().GetWorld();
    nvolumes_ = 0;
    while (node) {
      // Update path to match current level
      state->Push(node);
      baskets_[nvolumes_] = new Basket(vsize);
      // Synchronize all states with state
      baskets_[nvolumes_]->SetNavStates(state);
      // Synchronize all directions with dir
      baskets_[nvolumes_]->SetDirection(window->dir_);
      baskets_[nvolumes_]->window_ = window;
      baskets_[nvolumes_]->pixmap_ = pixmap;
      nvolumes_++;
      // Go to next level
      LogicalVolume const* lvol = node->GetLogicalVolume();
      Vector<VPlacedVolume const*> const &daughters = lvol->GetDaughters();
      int nd = daughters.size();
      assert(nd <= 1); // yeah, for now...
      if (nd==0) {node = nullptr; continue;}
      node = daughters[0];      
    }
    // Setup arrays of proposed steps and returned ones
    for (auto i=0; i<vsize; ++i) {
      psteps_[i] = kInfinity;
      steps_[i] = 0;
    }
    // Setup output states
    outstates_ =  new NavStatePool(vsize, GeoManager::Instance().getMaxDepth());
  }
  
  ~Stepper() {
    _mm_free(psteps_);
    _mm_free(steps_);
    for (auto i=0; i<nvolumes_; ++i) delete baskets_[i];
    delete outstates_;
  }
  
  void BasketTransport(Basket *basket) {
    // Main iterative transport loop. Won't do tails.
    // Transport one step and flush basket
    basket->Step(psteps_, steps_, outstates_, stack_);
    while (!stack_.empty()) {
      TrackState const &state = stack_.Pop();
      bool ready = baskets_[state.bindex]->AddTrack(state);
      if (ready) BasketTransport(baskets_[state.bindex]);
    }
  }
  
  void TailsTransport() {
    // Transport partially filled baskets until all tracks exit the setup
    bool loop = true;
    while (loop) {
      loop = false;
      for (auto i=0; i<nvolumes_; ++i) {
        if (baskets_[i]->size_) {
	  loop = true;
	  BasketTransport(baskets_[i]);
	}
      }	
    }
  }
  

  void TransportTask(Window const &window, int islice, int vecsize, int* volume_result)
  {
    // Transport thread entry point
    int istart1, npix1;
    std::cout << "Thread " << islice << std::endl;
    window.GetSubwindow(islice, istart1, npix1);

    // The stepper variable should be using TLS
    Basket *top = baskets_[0];
    TrackState state(0, 0, 0, 0, 0, 0);
    Vector3D<Precision> startpos;
    for( int pixel_count_2 = 0; pixel_count_2 < window.data_size_y; ++pixel_count_2 ) {
      for (int pixel_count_1 = istart1; pixel_count_1<istart1+npix1; ++pixel_count_1) {
        window.GetCoordinates(pixel_count_1, pixel_count_2, startpos);
        // Encode pixel coordinates
        state.xp = startpos[0];
        state.yp = startpos[1];
        state.zp = startpos[2];
        state.pixel = EncodePixel(pixel_count_1, pixel_count_2);
        if (top->AddTrack(state)) BasketTransport(top);
      }  
    }
    TailsTransport();
  }

};


void AssignNavigatorToVolume( LogicalVolume *vol, int layer, int maxlayers ){
  assert( maxlayers <= 20 ); // we are only listing 10 template specializations up to depth 10 here
  if( layer == maxlayers - 1 ) // it is the last layer
  {
    SETINNERNAV(0);
    SETINNERNAV(1);
    SETINNERNAV(2);
    SETINNERNAV(3);
    SETINNERNAV(4);
    SETINNERNAV(5);
    SETINNERNAV(6);
    SETINNERNAV(7);
    SETINNERNAV(8);
    SETINNERNAV(9);
    SETINNERNAV(10);
    SETINNERNAV(11);
    SETINNERNAV(12);
    SETINNERNAV(13);
    SETINNERNAV(14);
    SETINNERNAV(15);
    SETINNERNAV(16);
    SETINNERNAV(17);
    SETINNERNAV(18);
    SETINNERNAV(19);
    SETINNERNAV(20);
  }
  if( layer < maxlayers - 1 ){
    SETLAYERNAV(0);
    SETLAYERNAV(1);
    SETLAYERNAV(2);
    SETLAYERNAV(3);
    SETLAYERNAV(4);
    SETLAYERNAV(5);
    SETLAYERNAV(6);
    SETLAYERNAV(7);
    SETLAYERNAV(8);
    SETLAYERNAV(9);
    SETLAYERNAV(10);
    SETLAYERNAV(11);
    SETLAYERNAV(12);
    SETLAYERNAV(13);
    SETLAYERNAV(14);
    SETLAYERNAV(15);
    SETLAYERNAV(16);
    SETLAYERNAV(17);
    SETLAYERNAV(18);
    SETLAYERNAV(18);
    SETLAYERNAV(20);
  }
}


VPlacedVolume *CreateSimpleTracker(unsigned int nlayers) {
  std::cout << "Creating SimpleTracker geometry having " << nlayers << " layers"
            << std::endl;
  assert(nlayers < kMaxDepth-1);          
  // World size
  const double world_size = 500.;

  // Top volume
  UnplacedBox *uTop = new UnplacedBox(world_size + 2, world_size + 2, world_size + 2);
  LogicalVolume *top = new LogicalVolume("world", uTop);
  top->SetUserExtensionPtr( (void *) WorldNavigator<true>::Instance() );
  
// Cylindrical layers
  double rmax, dz;
  std::string layerBase = "layer_";
  double deltaR = world_size/nlayers;

  LogicalVolume *mother = top;
  for (unsigned int layer=0; layer<nlayers; ++layer) {
    rmax = world_size - layer * deltaR;
    dz = world_size - 0.1*layer*deltaR;
    std::ostringstream layerName;
    layerName << layerBase << layer;
    UnplacedTube *uLayer = new UnplacedTube(0, rmax, dz, 0, kTwoPi);
    LogicalVolume *layerVol = new LogicalVolume(layerName.str().c_str(), uLayer);
    
    AssignNavigatorToVolume(layerVol, layer, nlayers);

    // Place in mother
    mother->PlaceDaughter(layerName.str().c_str(), layerVol, &Transformation3D::kIdentity);
    // change mother to current lvol (to make a real hierarchy)
    mother = layerVol;
  }  

  VPlacedVolume *world = top->Place();
  GeoManager::Instance().SetWorld(world);
  GeoManager::Instance().CloseGeometry();
  return world;
}


int ScalarNavigation(Vector3D<Precision> &p, Vector3D<Precision> const &dir, NavigationState *&curnavstate, NavigationState *&newnavstate) {

  int crossedvolumecount=0;
  while( ! curnavstate->IsOutside() ) {
    //
    //std::cout << "tracking in " << curnavstate->Top()->GetLogicalVolume()->GetName() << "\n";
      // get navigator object and move point
    VNavigator const *specialnav = GetNavigator(curnavstate->Top()->GetLogicalVolume());
    double step = specialnav->ComputeStepAndPropagatedState(p, dir, kInfinity, *curnavstate, *newnavstate);
    //std::cout << "step " << step << "\n";
    p = p + dir * (step + 1E-6);

    // pointer swap is enough
    auto *tmp = curnavstate;
    curnavstate = newnavstate;
    newnavstate = tmp;
    if (step>0.0) crossedvolumecount++;
  }

  // now test the vector progression (of coherent rays)
  // this is just testing the interface and makes sure that no trivial things go wrong
  return crossedvolumecount;
}

int ScalarNavigation_NonSpecialized(Vector3D<Precision> &p, Vector3D<Precision> const &dir, NavigationState *&curnavstate, NavigationState *&newnavstate) {
  SimpleNavigator nav;
  int crossedvolumecount=0;
  while( ! curnavstate->IsOutside() ) {
    //
    double step = 0;
    nav.FindNextBoundaryAndStep(p, dir, *curnavstate, *newnavstate, kInfinity, step);
    //std::cout << "step " << step << "\n";
    p = p + dir * (step + 1E-6);

    // pointer swap is enough
    auto *tmp = curnavstate;
    curnavstate = newnavstate;
    newnavstate = tmp;
    if (step>0.0) crossedvolumecount++;
  }

  // now test the vector progression (of coherent rays)
  // this is just testing the interface and makes sure that no trivial things go wrong
  return crossedvolumecount;
}



////////////////// VECTOR NAVIGATION

int VectorNavigation (SOA3D<Precision> &points, SOA3D<Precision> const &dirs, int np, NavStatePool *&curnavstates, NavStatePool *&newnavstates, double *psteps, double *steps) {
// we assume that that curnavstates are already initialized correctly

  int crossedvolumecount=0;
  while(!(*curnavstates)[0]->IsOutside()) {
      //
     // std::cout << "tracking in " << (*curnavstates)[0]->Top()->GetLogicalVolume()->GetName() << "\n";
        // get navigator object and move point
      VNavigator const *specialnav = GetNavigator((*curnavstates)[0]->Top()->GetLogicalVolume());
      specialnav->ComputeStepsAndPropagatedStates(points, dirs, psteps, *curnavstates, *newnavstates, steps);

      // transport + crosschecks + pointerswap
      // We should vectorize this transport !!

      for(auto i=0; i<np; ++i){
        //assert(steps[i]==steps[0]);
        assert((*newnavstates)[i]->Top()==(*newnavstates)[0]->Top());

        points.set(i, points[i] + dirs[i]*(steps[i] + 1E-6));

      }

     // SOA3D<Precision>::SOA_vector_transport(points,dirs,steps);

      // pointer swap to update navstates
      auto tmp=curnavstates;
      curnavstates = newnavstates;
      newnavstates = tmp;
      if (steps[0]>0.0) crossedvolumecount++;
  }

    return crossedvolumecount;
}

// a test function to verify correct functioning of the navigators
// for a simple test case
void TestScalarNavigation() {
  // setup point and direction in world
  Vector3D<Precision> p(-500.,0,0);
  Vector3D<Precision> dir(1.,0,0);
  // init navstates
  NavigationState * curnavstate = NavigationState::MakeInstance(GeoManager::Instance().getMaxDepth());
  NavigationState * newnavstate = NavigationState::MakeInstance(GeoManager::Instance().getMaxDepth());
  NavigationState * worldnavstate = NavigationState::MakeInstance(GeoManager::Instance().getMaxDepth());
  SimpleNavigator nav;
  nav.LocatePoint( GeoManager::Instance().GetWorld(), Vector3D<Precision>(-500,-500,-500), *worldnavstate, true );
  worldnavstate->CopyToFixedSize<NavigationState::SizeOf(1)>(curnavstate);

  std::cout<<"crossedvolumecount "<<ScalarNavigation(p,dir,curnavstate, newnavstate)<<std::endl;

  NavigationState::ReleaseInstance(curnavstate);
  NavigationState::ReleaseInstance(newnavstate);
}

void TestVectorNavigation() {
  // use vector navigation
  // choose 11 particles to also make sure we have tail treatment
  auto np=11;
  SOA3D<Precision> points(np);
  SOA3D<Precision> dirs(np);
  for(auto i=0;i<np;++i){
      points.set(i, -51,0,0);
      dirs.set(i, 1,0,0);
  }


  NavStatePool *curnavstates = new NavStatePool(np, GeoManager::Instance().getMaxDepth());
  NavStatePool *newnavstates = new NavStatePool(np, GeoManager::Instance().getMaxDepth());
  double *steps    = (double*) _mm_malloc(np*sizeof(double),64);
  double *psteps    = (double*) _mm_malloc(np*sizeof(double),64);
  SimpleNavigator nav;
  nav.LocatePoint( GeoManager::Instance().GetWorld(), points[0], *(*curnavstates)[0], true );
  for(auto i=1;i<np;++i){
      (*curnavstates)[0]->CopyTo((*curnavstates)[i]);
      assert(curnavstates->operator[](i)->Top() == curnavstates->operator[](0)->Top());
  }

  std::cout << "VecNav crossedvolumecount "
            << VectorNavigation(points, dirs, np, curnavstates, newnavstates, psteps, steps) << std::endl;

  _mm_free(steps);
  _mm_free(psteps);
  delete curnavstates;
  delete newnavstates;
}
// reproducing the pixel-by-pixel XRayBenchmark
// target: show speed gain from specialized navigators
// in comparision to current XRay-Benchmark
void XRayBenchmark(int axis, int pixel_width) {

  std::stringstream imagenamebase;
  imagenamebase << "simpleTrackerImage_";
  if(axis==1) imagenamebase << "x";
  if(axis==2) imagenamebase << "y";
  if(axis==3) imagenamebase << "z";

  VPlacedVolume *vol = GeoManager::Instance().FindPlacedVolume(2);
  Window window(vol, axis, pixel_width, 1);

  int *volume_result= (int*) new int[window.data_size_y * window.data_size_x*3];
    // init navstates
    NavigationState * curnavstate = NavigationState::MakeInstance(GeoManager::Instance().getMaxDepth());
    NavigationState * newnavstate = NavigationState::MakeInstance(GeoManager::Instance().getMaxDepth());

    NavigationState *worldnavstate = NavigationState::MakeInstance(GeoManager::Instance().getMaxDepth());
    worldnavstate->Push(GeoManager::Instance().GetWorld());
    assert( worldnavstate->Top() == GeoManager::Instance().GetWorld());

    Vector3D<Precision> cpoint;
    Stopwatch timer;
    timer.Start();

    for( int pixel_count_2 = 0; pixel_count_2 < window.data_size_y; ++pixel_count_2 ) {
       for( int pixel_count_1 = 0; pixel_count_1 < window.data_size_x; ++pixel_count_1 ) {
          window.GetCoordinates(pixel_count_1, pixel_count_2, cpoint);
          // Start is always in the top volume
          worldnavstate->CopyToFixedSize<NavigationState::SizeOf(1)>(curnavstate);
          *(volume_result+pixel_count_2*window.data_size_x+pixel_count_1) = ScalarNavigation(cpoint,window.dir_,curnavstate, newnavstate);
          // *(volume_result+pixel_count_2*data_size_x+pixel_count_1) = ScalarNavigation_NonSpecialized(p,dir,curnavstate, newnavstate);
       } // end inner loop
    } // end outer loop

   timer.Stop();
   std::cout << " XRay Elapsed time : "<< timer.Elapsed() << std::endl;

  std::stringstream VecGeomimage;
  VecGeomimage << imagenamebase.str();
  VecGeomimage << "_VecGeom.bmp";
  make_bmp(volume_result, VecGeomimage.str().c_str(), window.data_size_x, window.data_size_y);

  NavigationState::ReleaseInstance(curnavstate);
  NavigationState::ReleaseInstance(newnavstate);
}



void XRayBenchmarkVecNav(int axis, int pixel_width, int vecsize) {

  const int N=vecsize;
  std::stringstream imagenamebase;
  imagenamebase << "simpleTrackerimage_";
  if(axis==1) imagenamebase << "x";
  if(axis==2) imagenamebase << "y";
  if(axis==3) imagenamebase << "z";
  
  VPlacedVolume *vol = GeoManager::Instance().FindPlacedVolume(2);
  Window window(vol, axis, pixel_width, 1);

  int *volume_result= (int*) new int[window.data_size_y * window.data_size_x*3];
   


  SOA3D<Precision> dirs(N);
  for(auto i=0;i<N;++i){
      dirs.set(i, window.dir_.x(), window.dir_.y(),window.dir_.z());
  }

  NavStatePool *curnavstates = new NavStatePool(N, GeoManager::Instance().getMaxDepth());
  NavStatePool *newnavstates = new NavStatePool(N, GeoManager::Instance().getMaxDepth());
  double *steps    = (double*) _mm_malloc(N*sizeof(double),64);
  double *psteps    = (double*) _mm_malloc(N*sizeof(double),64);

  SOA3D<Precision> points(N);

  NavigationState *worldnavstate = NavigationState::MakeInstance(GeoManager::Instance().getMaxDepth());
  worldnavstate->Push(GeoManager::Instance().GetWorld());
  assert( worldnavstate->Top() == GeoManager::Instance().GetWorld());

  Vector3D<Precision> cpoint;
  Stopwatch timer;
  timer.Start();

    for( int pixel_count_2 = 0; pixel_count_2 < window.data_size_y; ++pixel_count_2 ) {
       for( int pixel_count_1 = 0; pixel_count_1 < window.data_size_x; ++pixel_count_1 ) {
          window.GetCoordinates(pixel_count_1, pixel_count_2, cpoint);
          for(auto i=0;i<N;++i){
             points.set(i, cpoint[0], cpoint[1], cpoint[2]);
          }
          // init initial nav state from the top volume
          for(auto i=0;i<N;++i){
              worldnavstate->CopyToFixedSize<NavigationState::SizeOf(1)>(curnavstates->operator[](i));
          }

          *(volume_result+pixel_count_2*window.data_size_x+pixel_count_1) = VectorNavigation(points,dirs,N,curnavstates,newnavstates,psteps,steps);
      } // end inner loop
   } // end outer loop 
   timer.Stop();
   std::cout << " XRayVecNav Elapsed time /Nrep : "<< timer.Elapsed()/N << std::endl;

    std::stringstream VecGeomimage;
    VecGeomimage << imagenamebase.str();
    VecGeomimage << "_Vector_VecGeom.bmp";
    make_bmp(volume_result, VecGeomimage.str().c_str(), window.data_size_x, window.data_size_y);

    _mm_free(steps);
    _mm_free(psteps);
    delete curnavstates;
    delete newnavstates;

}



void XRayBenchmarkBasketized(int axis, int pixel_width, int vecsize, int nthreads) {
  // to be filled in by Andrei  
  std::stringstream imagenamebase;
  imagenamebase << "simpleTrackerimage_";
  if(axis==1) imagenamebase << "x";
  if(axis==2) imagenamebase << "y";
  if(axis==3) imagenamebase << "z";

  VPlacedVolume *vol = GeoManager::Instance().FindPlacedVolume(2);
  Window window(vol, axis, pixel_width, nthreads);

  // Create the pixmap
  int *volume_result= (int*) new int[window.data_size_y * window.data_size_x*3];

  // Create steppers for each slice. This is not best for NUMA - there should be rather the thread code doing this
  // Another policy would be to have as many steppers as threads, making sure that the same thread picks the same stepper
  Stepper **steppers = new Stepper*[nthreads];
  for (auto i=0; i<nthreads; ++i) steppers[i] = new Stepper(vecsize, &window, volume_result);

  // This is the main parallelizable loop
  Stopwatch timer;
  timer.Start();

//#ifdef VECGEOM_OPENMP
  #pragma omp parallel
//#endif  
  {
//#ifdef VECGEOM_OPENMP
    #pragma omp for schedule(dynamic) nowait
//#endif
    for (auto i=0; i<nthreads; ++i) {
      steppers[i]->TransportTask(window, i, vecsize, volume_result);
    }    
  }
  timer.Stop();
  std::cout << " XRayBasketized Elapsed time : "<< timer.Elapsed() << std::endl;

  std::stringstream VecGeomimage;
  VecGeomimage << imagenamebase.str();
  VecGeomimage << "_basket_VecGeom.bmp";
  make_bmp(volume_result, VecGeomimage.str().c_str(), window.data_size_x, window.data_size_y);
  
  for (auto i=0; i<nthreads; ++i) delete steppers[i];
  delete [] steppers;
}

// produce a bmp image out of pixel information given in volume_results
void make_bmp_header(MY_BITMAP * pBitmap, unsigned char * bmpBuf, int sizex, int sizey)
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





void make_bmp(int const * volume_result, char const *name, int data_size_x, int data_size_y, bool linear)
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
}


int main(int argc, char* argv[]) {
  int nlayers = 10;
  if (argc > 1) {
    nlayers = std::atoi(argv[1]);
  }
  CreateSimpleTracker(nlayers);

  int axis= 1;

  if( strcmp(argv[2], "x")==0 )
    axis= 1;
  else if( strcmp(argv[2], "y")==0 )
    axis= 2;
  else if( strcmp(argv[2], "z")==0 )
    axis= 3;  

  int pixel_width = 500;
  if (argc > 3) pixel_width = atoi(argv[3]);
  int vecsize = 32;
  if (argc > 4) vecsize = atoi(argv[4]);
  
  int nthreads = 1;
  if (argc > 5) nthreads = atoi(argv[5]);
  
//  TestScalarNavigation();
//  TestVectorNavigation();
  std::cout << "Running with: img. width=" << pixel_width << " vecsize=" << vecsize << " nthreads=" << nthreads << " on axis " << argv[2] << std::endl;
  XRayBenchmark(axis, pixel_width);
  XRayBenchmarkVecNav(axis, pixel_width, vecsize);
  XRayBenchmarkBasketized(axis, pixel_width, vecsize, nthreads);
}
