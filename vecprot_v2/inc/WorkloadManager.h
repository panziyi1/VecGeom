//===--- WorkloadManager.h - Geant-V ----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file WorkloadManager.h
 * @brief Definition of workload manager in Geant-V prototype  
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_WORKLOADMANAGER
#define GEANT_WORKLOADMANAGER

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "priority_queue.h"
#include "condition_locker.h"

#include "GeantTrack.h"

class GeantBasketMgr;
class GeantBasket;
class GeantScheduler;
class TaskBroker;


/**
 * @brief WorlloadManager class
 * @details Main work manager class. This creates and manages all the worker threads,
 * has pointers to the synchronization objects, but also to the currently
 * transported baskets.
 */
class WorkloadManager : public TObject {
public:

  /**
   * @brief Monitoring type
   */
enum EGeantMonitoringType {
  kMonQueue = 0,
  kMonMemory,
  kMonBasketsPerVol,
  kMonConcurrency,
  kMonTracksPerEvent
};

protected:
  Int_t fNthreads;         /** Number of managed threads */
  Int_t fNbaskets;         /** Total number of baskets */
  Int_t fBasketGeneration; /** Basket generation */
  Int_t fNbasketgen;       /** Number of baskets to transport in the current generation */
  Int_t fNidle;            /** Number of idle workers */
  Int_t fNminThreshold;    /** Minimum number of tracks in a basket to trigger transport */
  Int_t fNqueued;          /** Number of chunks queued */
  Int_t *fBtogo;           /** Array of baskets to be processed in the next generation */
  bool fStarted;         /** Start flag */
  bool fStopped;         /** Stop flag */
  Geant::priority_queue<GeantBasket *> *fFeederQ;      /** Queue of transportable baskets */
  Geant::priority_queue<GeantBasket *> *fTransportedQ; /** Queue of transported baskets */
  Geant::priority_queue<GeantBasket *> *fDoneQ;        /** Thread "all work done" queue */
  static WorkloadManager *fgInstance; /** Singleton instance */
  TList *fListThreads;                /** List of threads */
  bool fFlushed;                    /** Buffer flushed */
  bool fFilling;                    /** Worker queue is filling */
  Int_t  fMonQueue;                   /** Monitor the work queue */
  Int_t  fMonMemory;                  /** Monitor the memory */
  Int_t  fMonBasketsPerVol;           /** Monitor baskets per volume */
  Int_t  fMonConcurrency;             /** Monitor concurrency */
  Int_t  fMonTracksPerEvent;          /** Monitor tracks status per event */
  GeantScheduler *fScheduler;         /** Main basket scheduler */

  TaskBroker *fBroker; /** Pointer to the coprocessor broker, this could be made a collection. */
  Int_t *fWaiting;     /** ![fNthreads+1] Threads in waiting flag */
  condition_locker fSchLocker; /** Scheduler locker */
  condition_locker fGbcLocker; /** Garbage collector locker */

  /** 
   * @brief WorkloadManager parameterized constructor
   * 
   * @param  nthreads Number of threads foe workload manager
   */
  WorkloadManager(Int_t nthreads);

public:

  /** @brief WorkloadManager destructor */
  virtual ~WorkloadManager();

  /** @brief Create basket function */
  void CreateBaskets();

  /** @brief Function for feeder queue of transportable baskets */
  Geant::priority_queue<GeantBasket *> *FeederQueue() const { return fFeederQ; }
  
  /**
   * @brief Function that provides trasported queue for baskets
   */
  Geant::priority_queue<GeantBasket *> *TransportedQueue() const { return fTransportedQ; }
  
  /**
   * @brief Function that provides thread's "all work done" queue
   */
  Geant::priority_queue<GeantBasket *> *DoneQueue() const { return fDoneQ; }
  //   GeantObjectPool<VolumePath_t>
  //   rr_pool<VolumePath_t>
  //                      *NavStates() const   {return fNavStates;}
  
  /** @brief Function that returns number of managed threads */
  Int_t GetNthreads() const { return fNthreads; }

  /** @brief Function that returns total number of baskets */
  Int_t GetNbaskets() const { return fNbaskets; }

  /** @brief Function that returns threads in waiting flag */
  Int_t *GetWaiting() const { return fWaiting; }

  /** @brief Function returning the number of monitored features */
  Int_t GetMonFeatures() const;

  /** @brief Check if a monitoring feature is enabled */
  bool IsMonitored(EGeantMonitoringType feature) const;

  /** @brief Enable monitoring a feature */
  void SetMonitored(EGeantMonitoringType feature, bool flag = true);

  /** @brief Function that returns main basket scheduler */
  GeantScheduler *GetScheduler() const { return fScheduler; }

  /** @brief Function that returns scheduler locker */
  condition_locker &GetSchLocker() { return fSchLocker; }

  /** @brief Function that returns garbage collector locker */
  condition_locker &GetGbcLocker() { return fGbcLocker; }

  /**
   * @brief Function that create workload manager instance
   * 
   * @param nthreads Number of threads (by default 0)
   */
  static WorkloadManager *Instance(Int_t nthreads = 0);

  /** @brief Function that check if buffer is flushed */
  bool IsFlushed() const { return fFlushed; }

  /** @brief Function that check if worker queue is filling */
  bool IsFilling() const { return fFilling; }

  /** @brief Function that check stop flag */
  bool IsStopped() const { return fStopped; }
  
  /**
   * @brief Function that provide stop process by setting Stop flag = True
   */
  void Stop() { fStopped = kTRUE; }

  /**
   * @brief Function that provides buffer flushing
   * 
   * @param flag Flag for buffer flushing
   */
  void SetFlushed(bool flag) { fFlushed = flag; }

  /** @brief Function that returns basket generation */
  Int_t GetBasketGeneration() const { return fBasketGeneration; }

  /** @brief Print function */
  void Print(Option_t *option = "") const;

  /**
   * @brief  Function that set task broker
   * 
   * @param broker Broker to be set 
   */
  void SetTaskBroker(TaskBroker *broker);

  /** @brief Function that return minimum number of tracks in a basket to trigger transport */
  Int_t GetNminThreshold() const { return fNminThreshold; }

  /**
   * @brief Function that set minimum number of tracks in the basket to trigger transport
   * 
   * @param thr Thread for minimum value of threshold
   */
  void SetNminThreshold(Int_t thr) { fNminThreshold = thr; }

  /** @brief Function that provides start process of threads  */
  void StartThreads();

  /** @brief Function that provides join process of threads */
  void JoinThreads();

  /**
   * @brief Function that provides main scheduler
   *  
   * @param arg Arguments to be passed in the function
   */
  static void *MainScheduler(void *arg);

  /**
   * @brief Function that provides garbage collector thread
   * 
   * @param arg Arguments to be passed in the function
   */
  static void *GarbageCollectorThread(void *arg);

  /**
   * @brief Function that provides monitoring thread
   * 
   * @param arg Arguments to be passed in the function
   */
  static void *MonitoringThread(void *arg);

  /**
   * @brief Function that provides transporting tracks
   * 
   * @param arg Arguments to be passed in the function
   */
  static void *TransportTracks(void *arg);

  /**
   * @brief Function that provides transport tracks in coprocessor
   * 
   * @param arg Arguments to be passed in the function
   */
  static void *TransportTracksCoprocessor(void *arg);

  /** @brief Function that provides waiting of workers */
  void WaitWorkers();

private:
  /**
   * @brief Copy constructor for WorkloadManager
   * @details Still not implemented
   */
  WorkloadManager(const WorkloadManager &);

  /**
   * @brief Operator =
   * @todo Still not implemented
   */
  WorkloadManager &operator=(const WorkloadManager &);

  ClassDef(WorkloadManager, 0) // The work manager class.
};
#endif
