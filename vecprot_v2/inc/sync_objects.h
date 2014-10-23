#ifndef GEANT_SYNCOBJECTS
#define GEANT_SYNCOBJECTS
#include <deque>
#include <cassert>
#if __cplusplus >= 201103L
#include <atomic>
#endif
#include "TCondition.h"
#include "TMutex.h"

using namespace std;

class TStopwatch;
class TObject;

//______________________________________________________________________________
struct TimeCounter {
   Int_t       nthreads;
   TStopwatch *timer;
   Double_t    stamp;
   Double_t    realtime[100];
   
   TimeCounter(bool measure_time);
   TimeCounter(const TimeCounter &); // not implemented
   ~TimeCounter();
   TimeCounter &operator =(const TimeCounter &); // not implemented
   TimeCounter &operator++();
   TimeCounter &operator--();
   void Print();
};         

//______________________________________________________________________________
class concurrent_queue
{
private:
   deque<TObject*>   the_queue;
   mutable TMutex    the_mutex;
   TCondition        the_condition_variable;
   TimeCounter      *the_counter;
   Int_t             nobjects;
   Int_t             npriority;
   
   concurrent_queue(const concurrent_queue&); // not implemented
   concurrent_queue &operator=(const concurrent_queue&); // not implemented
public:
   concurrent_queue(bool counter=false);
   ~concurrent_queue();
   
   Int_t              assigned_workers() const {return the_counter->nthreads;}
   void               push(TObject *data, Bool_t priority=false);
   Int_t              size() const;
   Int_t              size_async() const {return nobjects;}
   Bool_t             empty() const;
   Bool_t             empty_async() const {return nobjects == 0;}
   TObject*           wait_and_pop();
   TObject*           wait_and_pop_max(UInt_t nmax, UInt_t &n, TObject **array);
   void               pop_many(UInt_t n, TObject **array);
   Int_t              size_priority() const {return npriority;}
   Int_t              size_objects() const {return nobjects;}
   void               Print();
};

// Reference counted atomic pointer
//______________________________________________________________________________
template <class T>
class ref_ptr {
public:
   T*                 fObjPtr;     // Object pointer
#if __cplusplus >= 201103L
   std::atomic_flag   fAcqLock;    // memory barrier to lock acquiring
   std::atomic_flag   fXcgLock;    // memory barrier to lock exchanging
   std::atomic_int    fRef;        // ref counter
#endif   
public:
   ref_ptr() : fObjPtr(0),fAcqLock(false),fXcgLock(false),fRef(0) {}
   ref_ptr(T* objptr) : fObjPtr(objptr),fAcqLock(false),fXcgLock(false),fRef(0) {}
   
   void               set(T* obj) {fObjPtr = obj;}
   T*                 acquire();
   void               clear_acq() {fAcqLock.clear();}
   void               clear_xcg() {fXcgLock.clear();}
   void               release();
   void               release_and_wait();
   bool               replace_and_wait(T* expected, T* newptr);
};

template <class T>
inline T* ref_ptr<T>::acquire()
{
// Acquire the pointer while issuing a memory barrier for other possible clients
  while (fAcqLock.test_and_set()) {};    // barrier here for other threads
//__mutexed code start
  T *cpy = fObjPtr; // private snapshot
  fRef++;
//__mutexed code end  
  fAcqLock.clear();
  return cpy;
}  

template <class T>
inline void ref_ptr<T>::release()
{
// Release the copy of the pointer
   fRef--;
   assert(fRef.load()>=0);
}

template <class T>
inline void ref_ptr<T>::release_and_wait()
{
// Release the copy of the pointer, but wait for everybody else to do it
// This acts as a spinlock
   fRef--;
   while (fRef.load() > 0) {};
   assert(fRef.load()>=0);
}

template <class T>
inline bool ref_ptr<T>::replace_and_wait(T* expected, T* newptr)
{
// Atomically replace the pointer only if the old equals expected
// and only when it is not used anymore. Returns true if the replacement
// succeeded, in which case acquiring is blocked. The lock must be released
// by the calling thread using clear() afterwards
  if (fXcgLock.test_and_set()) {
     // someone already passed for the same object
     release();
     while (fXcgLock.test_and_set()) {};
     fXcgLock.clear();
     assert(fRef.load()>=0);
     return false;
  }
  // If someone has stolen the pointer, exit   
  if (fObjPtr != expected) {
     // you thief!
     // Won't rely on expected pointer again
     release();
     assert(fRef.load()>=0);
     fXcgLock.clear();
     return false;
  }   
  // Thou shalt not pass
  fObjPtr = newptr;
  // go on, there's nothing to steal anymore
  // Lock acquiring new references
  while (fAcqLock.test_and_set()) {};
  release_and_wait();
  // now you're all mine...
  fAcqLock.clear();
  fXcgLock.clear();
  return true;
}

//______________________________________________________________________________
// Work pipe. It has a dequeue<T> plus a function Process(T*) provided by the 
//            user
//______________________________________________________________________________

class basepipe {
protected:
   int                 nobjects;   // Number of objects in the queue
   int                 npriority;  // Number of prioritized objects
   int                 priority;   // Assigned pipe priority
   
public:
   basepipe() : nobjects(0), npriority(0), priority(0) {}
   virtual ~basepipe() {}

   bool                empty() const {return (nobjects==0);}
   int                 size_priority() const {return npriority;}
   int                 size() const {return nobjects;}   
   virtual void        process() = 0;
};   

template <class T>
class workpipe : public basepipe {
   typedef void* (*ProcessFunc_t)(T *data);
   ProcessFunc_t       the_function; // Fuction to process this data
   deque<T*>           the_queue;  // Double-ended data queue
   mutable TMutex      the_mutex;  // General mutex for the queue
   TCondition          the_condition_variable; // Condition

   bool                process_if_any();
public:
   workpipe(ProcessFunc_t func): basepipe(), the_function(func), the_queue(), the_mutex(), the_condition_variable(&the_mutex) {}
   ~workpipe() {}
   void                push(T *data, bool priority=false);
   void                process();
};

template <class T>
void workpipe<T>::push(T *data, bool pr)
{
// Push an pointer of type T* in the queue. If pushed with priority, the pointer
// is put at the bask of the queue, otherwise to the front.
   the_mutex.Lock();
   nobjects++;
   if (pr) {the_queue.push_back(data); npriority++;}
   else          the_queue.push_front(data);
   the_condition_variable.Signal();
   the_mutex.UnLock();
}

template <class T>
void workpipe<T>::process()
{
// Gets the back object from the queue. Wait if none is available. Call the
// user Process function for the first available object.
   the_mutex.Lock();
   while(the_queue.empty()) the_condition_variable.Wait();        
   T *popped_value = the_queue.back();
   the_queue.pop_back();
   nobjects--;
   if (npriority>0) npriority--;
   the_mutex.UnLock();
   // Execute the user function
   the_function(popped_value);
}

template <class T>
bool workpipe<T>::process_if_any()
{
// If any object is available in the queue, pop it and call the
// user Process function for it. Returns true if the processing was called
   the_mutex.Lock();
   if (!nobjects) {
      the_mutex.UnLock();
      return false;
   }   
   T *popped_value = the_queue.back();
   the_queue.pop_back();
   nobjects--;
   if (npriority>0) npriority--;
   the_mutex.UnLock();
   // Execute the user function
   the_function(popped_value);
   return true;
}

//class workflow {
   
#endif
