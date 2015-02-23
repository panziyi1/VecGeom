//===--- GeantBasket.h - Geant-V --------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantBasket.h
 * @brief Implementation of baskets of tracks for Geant-V prototype
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_BASKET
#define GEANT_BASKET

#include "TObject.h"
#include "TGeoExtension.h"
#include "TGeoVolume.h"
#include "GeantTrack.h"
#include "priority_queue.h"
#include "mpmc_bounded_queue.h"

class TGeoVolume;
class GeantBasketMgr;

/**
 * @brief Class GeantBasket descripting basic operations with baskets
 * @details Basket of tracks in the same volume which are transported by a single thread
 */
class GeantBasket : public TObject {
public:

  /**
   * @enum EbasketFlags
   * @details Flags marking up different types of baskets
   */
  enum EBasketFlags {
    kMixed = BIT(14)        /** Basket mixing tracks from different volumes */
  };

protected:
  GeantBasketMgr *fManager; /** Manager for the basket */
  GeantTrack_v fTracksIn;   /** Vector of input tracks */
  GeantTrack_v fTracksOut;  /** Vector of output tracks */
//GeantHit_v        fHits;  /** Vector of produced hits */
#if __cplusplus >= 201103L
  std::atomic_int fAddingOp; /** Number of concurrent track adding operations */
#endif
  Int_t fThreshold; /** Current transport threshold */

private:

  /** @todo Still not implemented */
  GeantBasket(const GeantBasket &);

  /** @todo Still not implemented operator = */
  GeantBasket &operator=(const GeantBasket &);
public:
  
  /** @brief Default GeantBasket constructor */
  GeantBasket();
  
  /**
   * @brief GeantBasket standard constructor 
   * 
   * @param size Imitial size of input/output track arrays
   * @param mgr  Basket manager handling this basket
   */
  GeantBasket(Int_t size, GeantBasketMgr *mgr);
  
  /** @brief GeantBasket destructor */
  virtual ~GeantBasket();

  /**
   * @brief Add a scalar track to basket input.
   * @details Add concurrently track from generator or physics process
   * 
   * @param track Reference to track object
   */
  void AddTrack(GeantTrack &track);
  
  /**
   * @bref Add a track from vector container to basket input.
   * @details Add concurrently a track to basket input.
   * 
   * @param tracks Array of tracks to copy from.
   * @param itr Track id.
   */
  void AddTrack(GeantTrack_v &tracks, Int_t itr);
  
  /**
   * @brief Function to add multiple tracks to basket
   * @details Add multiple tracks from a track_v array
   * 
   * @param tracks Tracks from a track_v array 
   * @param istart Start track id 
   * @param iend End track id
   */
  void AddTracks(GeantTrack_v &tracks, Int_t istart, Int_t iend);
  
  /** @brief Virtual function for clearing the basket */
  virtual void Clear(Option_t *option = "");
  
  /**
   * @brief Check if a basket contains tracks in a given event range
   * 
   * @param evstart Start event id.
   * @param nevents Number of events (default 1)
   */
  Bool_t Contains(Int_t evstart, Int_t nevents = 1) const;
  
  /**
   * @brief Function returning the number of input tracks
   * @return Number of input tracks
   */
  Int_t GetNinput() const { return fTracksIn.GetNtracks(); }
  
  /**
   * @brief Function returning the number of output tracks
   * @return Number of output tracks
   */
  Int_t GetNoutput() const { return fTracksOut.GetNtracks(); }
  
  /**
   * @brief Function returning a reference to the vector of input tracks
   * @return Reference to input vector of tracks
   */
  GeantTrack_v &GetInputTracks() { return fTracksIn; }
  
  /**
   * @brief Function returning a reference to the vector of output tracks
   * @return Reference to output vector of tracks
   */
  GeantTrack_v &GetOutputTracks() { return fTracksOut; }
  
  /**
   * @brief Function returning the manager of the basket
   * @return Pointer to manager of the basket
   */
  GeantBasketMgr *GetBasketMgr() const { return fManager; }
  
  /**
   * @brief Function for defining basket transportability threshold
   * @return  Value of transportability threshold
   */
  Int_t GetThreshold() const { return fThreshold; }
  
  /**
   * @brief Function returning the volume for this basket
   * @return Pointer to associated logical volume
   */
  TGeoVolume *GetVolume() const;
  
  /**
   * @brief Function returning the mixed tracks property.
   * @return Boolean value if the tracks are mixed from several volumes
   */
  Bool_t IsMixed() const { return TObject::TestBit(kMixed); }
  
  /**
   * @brief Atomic snapshot of the number of concurrent track adding operations
   * @return Number of track adding operations
   */
  inline Bool_t IsAddingOp() const { return (fAddingOp.load()); }
  
  /**
   * @brief Atomic increment of the number of track adding operations
   * @return Number of concurrent track adding operations
   */
  inline Int_t LockAddingOp() { return ++fAddingOp; }
  
  /**
   * @brief Atomic decrement of the number of track adding operations
   * @return Number of concurrent track adding operations remaining
   */
  inline Int_t UnLockAddingOp() { return --fAddingOp; }
  
  /**
   * @brief Print the basket content
   */
  virtual void Print(Option_t *option = "") const;
  
  /**
   * @brief Print the parameters for a given track
   * 
   * @param itr Track id.
   * @param input Refer to input or output track (default input)
   */
  void PrintTrack(Int_t itr, Bool_t input = kTRUE) const;
  
  /** @brief Recycle this basket */
  void Recycle();
  
  /**
   * @brief  Function that providing the size of this basket in bytes
   * @return Sum of sizes of all tracks (input and output) + data members
   */
  size_t Sizeof() const {
    return fTracksIn.Sizeof() + fTracksOut.Sizeof() + sizeof(TObject) + sizeof(GeantBasketMgr *) +
           sizeof(std::atomic_int);
  }

  /**
   * @brief Flag the basket to contain tracks that are mixed from different volumes
   * 
   * @param flag Boolean flag.
   */ 
  void SetMixed(Bool_t flag) { TObject::SetBit(kMixed, flag); }
  
  /**
   * @brief Function to change the transportability threshold for the basket
   * 
   * @param threshold New threshold value
   */
  void SetThreshold(Int_t threshold);

  ClassDef(GeantBasket, 1) // A basket containing tracks in the same geomety volume
};

class GeantScheduler;

/**
 * @brief Class managing all baskets for a given logical volume
 * @details Basket manager for a given volume. Holds a list of free baskets stored in a
 * concurrent queue, a current basket to be filled and a priority basket used
 * in priority mode.
 */
class GeantBasketMgr : public TGeoExtension {
protected:
  GeantScheduler *fScheduler; /** Scheduler for this basket */
  TGeoVolume *fVolume;        /** Volume for which applies */
  Int_t fNumber;              /** Number matching the volume index */
  Int_t fBcap;                /** Maximum capacity of baskets held */
  Int_t fQcap;                /** Queue capacity */
#if __cplusplus >= 201103L
  std::atomic_int fThreshold; /** Adjustable transportability threshold */
  std::atomic_int fNbaskets;  /** Number of baskets for this volume */
  std::atomic_int fNused;     /** Number of baskets in use */
  typedef std::atomic<GeantBasket *> atomic_basket;
  atomic_basket fCBasket;  /** Current basket being filled */
  atomic_basket fPBasket;  /** Current priority basket being filled */
  std::atomic_flag fLock;  /** Atomic lock for stealing current basket */
  std::atomic_flag fQLock; /** Atomic lock for increasing queue size */
#endif
  mpmc_bounded_queue<GeantBasket *> *fBaskets;   /** Queue of available baskets */
  Geant::priority_queue<GeantBasket *> *fFeeder; /** Feeder queue to which baskets get injected */
  TMutex fMutex;                                 /** Mutex for this basket manager */
private:

  /** @todo Still not implemented */
  GeantBasketMgr(const GeantBasketMgr &);

  /** @todo Still not implemented operator = */          
  GeantBasketMgr &operator=(const GeantBasketMgr &);
#if __cplusplus >= 201103L

  /**
   * @brief Attempt to steal the current filled basket and replace it
   *  
   * @param current Current atomic basket to be replaced
   * @return Released basket pointer if operation succeeded, 0 if not
   */
  GeantBasket *StealAndReplace(atomic_basket &current);

  /**
   * @brief The caller thread steals temporarily the basket to mark a track addition
   * 
   * @param current Current atomic basket to be pinned to a thread for track adding
   * @return Current basket being pinned
   */
  GeantBasket *StealAndPin(atomic_basket &current);

  /**
   * @brief Attempt to steal a global basket matching the content
   * 
   * @param global Global atomic basket
   * @param content Content of GeantBasket 
   * @return Flag marking the success/failure of the steal operation
   */
  Bool_t StealMatching(atomic_basket &global, GeantBasket *content);
#endif

  /**
   * @brief Function for increasing the current queue size 
   * 
   * @param newsize New size of queue  
   */
  Bool_t IncreaseQueueSize(Int_t newsize);

public:

  /**
   * @brief Get next free basket 
   * @details Get pointer to next free basket
   */
  GeantBasket *GetNextBasket();

public:

  /** @brief GeantBasketMgr dummy constructor */
  GeantBasketMgr()
      : fScheduler(0), fVolume(0), fNumber(0), fBcap(0), fQcap(0), fThreshold(0), fNbaskets(0),
        fNused(0), fCBasket(0), fPBasket(0), fLock(), fQLock(), fBaskets(0), fFeeder(0), fMutex() {}

  /** @brief GeantBasketMgr normal constructor 
   *
   * @param sch Scheduler dealing with this basket manager
   * @param vol Volume associated with this
   * @param number Number for the basket manager
  */
  GeantBasketMgr(GeantScheduler *sch, TGeoVolume *vol, Int_t number);
  
  /** @brief Destructor of GeantBasketMgr */
  virtual ~GeantBasketMgr();
  
  /**
   * @brief Grab function
   * @details Interface of TGeoExtension for getting a reference to this from TGeoVolume
   * @return Pointer to the base class
   */
  virtual TGeoExtension *Grab() { return this; }

  /** 
   * @brief Release function
   * @details Interface of TGeoExtension to signal releasing ownership of this from TGeoVolume
   */
  virtual void Release() const {}

  /**
   * @brief Function adding a track to basket up to the basket threshold
   * 
   * @param track  Track that should be added to basket
   * @param priority Set priority (by default kFALSE)
   */
  Int_t AddTrack(GeantTrack &track, Bool_t priority = kFALSE);

  /**
   * @brief Function adding an indexed track to basket up to the basket threshold
   * 
   * @param trackv Array of tracks containing the track to be added
   * @param itr Track id
   * @param priority Set priority (by default kFALSE)
   */
  Int_t AddTrack(GeantTrack_v &trackv, Int_t itr, Bool_t priority = kFALSE);

  /**
   * @brief Thread local garbage collection of tracks from prioritized events
   *
   * @param evmin Minimum event index
   * @param evmax Maximum event index
   * @param gc Garbage collector basket
   */
  Int_t GarbageCollectEvents(Int_t evmin, Int_t evmax, GeantBasketMgr *gc);

  /**
   * @brief Garbage collection of prioritized tracks in an event range
   * 
   * @param evmin Minimum event index
   * @param evmax Maximum event index
   */
  Int_t CollectPrioritizedTracks(Int_t evmin, Int_t evmax);

  /**
   * @brief Function cleaning a number of free baskets
   * 
   * @param ntoclean Number of baskets to be cleaned
   */
  void CleanBaskets(Int_t ntoclean);

  /** @brief Function flushing to the work queue only priority baskets */
  Int_t FlushPriorityBasket();

  /** @brief Function doing full collection to work queue of non-empty baskets*/
  Int_t GarbageCollect();

  /**
   * @brief Function that set capacity of baskets
   * 
   * @param capacity Capacity of baskets to be set
   */
  void SetBcap(Int_t capacity) { fBcap = capacity; }

  /**
   * @brief Function that returns the capacity of baskets
   * @return Maximum capacity of baskets held
   */
  Int_t GetBcap() const { return fBcap; }

#if __cplusplus >= 201103L

   /**
    * @brief Snapshot of the number of baskets
    * @return Number of baskets
    */
  Int_t GetNbaskets() const { return fNbaskets.load(); }

  /**
   * @brief Snapshot of the number of baskets in use
   * @return number of baskets in use
   */
  Int_t GetNused() const { return fNused.load(); }

  /**
   * @brief Snapshot of the current basket threshold 
   * @return Threshold value
   */
  Int_t GetThreshold() const { return fThreshold.load(); }

  /**
   * @brief Function to set transportability threshold
   * 
   * @param thr Threshold that should be set
   */
  void SetThreshold(Int_t thr) { fThreshold.store(thr); }

  /**
   * @brief Function that will load the current basket
   * @return Load current basket
   */
  GeantBasket *GetCBasket() const { return fCBasket.load(); }

  /**
   * @brief Function that will load the priority basket
   * @return Load priority basket 
   */
  GeantBasket *GetPBasket() const { return fPBasket.load(); }

  /**
   * @brief Function to set current basket
   * 
   * @param basket Basket to be set as current
   */
  void SetCBasket(GeantBasket *basket) { fCBasket.store(basket); }

  /**
   * @brief Function to set the priority basket 
   * 
   * @param basket Basket that will be set as priority basket
   */
  void SetPBasket(GeantBasket *basket) { fPBasket.store(basket); }
#endif

  /**
   * @brief Function that returns the scheduler
   * @return Scheduler for basket
   */
  GeantScheduler *GetScheduler() const { return fScheduler; }

  /** @brief Function that returns the name of volume */
  const char *GetName() const { return (fVolume) ? fVolume->GetName() : ClassName(); }

  /**
   * @brief Function that returns the number assigned to basket
   * @return Number assigned to basket
   */
  Int_t GetNumber() const { return fNumber; }

  /**
   * @brief Function that returns the associated volume pointer
   * @return Volume for which applies basket
   */
  TGeoVolume *GetVolume() const { return fVolume; }

  /** @brief Print the current basket */
  virtual void Print(Option_t *option = "") const;

  /**
   * @brief Recycles a given basket
   * 
   * @param b Pointer to current GeantBasket for recycling
   */
  void RecycleBasket(GeantBasket *b);

  /**
   * @brief Function setting the feeder work queue
   * 
   * @param queue priority_queue for GeantBasket
   */
  void SetFeederQueue(Geant::priority_queue<GeantBasket *> *queue) { fFeeder = queue; }

  /**
   * @brief Function returning the size of the basket being filled
   * @return Size of basket
   */
  size_t Sizeof() const {
    GeantBasket *c = GetCBasket();
    return (c) ? (sizeof(GeantBasketMgr) + (GetNbaskets()) * c->Sizeof()) : sizeof(GeantBasketMgr);
  }

  /** @brief Print size */
  void PrintSize() const;

  /**
   * @brief Function GetFeederQueue() 
   * @return Feeder queue to which baskets get injected
   */
  Geant::priority_queue<GeantBasket *> *GetFeederQueue() const { return fFeeder; }

  ClassDef(GeantBasketMgr, 0) // A path in geometry represented by the array of indices
};
#endif
