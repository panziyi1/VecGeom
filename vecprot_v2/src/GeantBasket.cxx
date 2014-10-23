#include "TThread.h"
#include "GeantBasket.h"
#include "globals.h"
#include "GeantTrack.h"
#include "GeantPropagator.h"
#include "GeantScheduler.h"
#include "PhysicsProcess.h"
#include "WorkloadManager.h"

#include "TThread.h"
#include "TArrayI.h"
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNavigator.h"

ClassImp(GeantBasket)

//______________________________________________________________________________
GeantBasket::GeantBasket()
            :TObject(),
             fManager(0),
             fTracksIn(),
             fTracksOut(),
             fAddingOp(0)
{
// dummy ctor.
}

//______________________________________________________________________________
GeantBasket::GeantBasket(Int_t size, GeantBasketMgr *mgr)
            :TObject(),
             fManager(mgr),
             fTracksIn(size),
             fTracksOut(size),
             fAddingOp(0)
{
// ctor.
}

//______________________________________________________________________________
GeantBasket::~GeantBasket()
{
// dtor.
}
   
//______________________________________________________________________________
void GeantBasket::AddTrack(GeantTrack &track)
{
// Add a new track to this basket;
//   fTracksIn.AddTrack(track);
//   assert(fAddingOp.load()>0);
   fTracksIn.AddTrackSync(track);
}

//______________________________________________________________________________
void GeantBasket::AddTrack(GeantTrack_v &tracks, Int_t itr)
{
// Add track from a track_v array.
// Has to work concurrently
//   fTracksIn.AddTrack(tracks, itr, kTRUE);
//   assert(fAddingOp.load()>0);
   fTracksIn.AddTrackSync(tracks, itr);
}

//______________________________________________________________________________
void GeantBasket::AddTracks(GeantTrack_v &tracks, Int_t istart, Int_t iend)
{
// Add multiple tracks from a track_v array. Not used
//   assert(fAddingOp.load()>0);
   fTracksIn.AddTracks(tracks, istart, iend, kTRUE);
}
   
//______________________________________________________________________________
void GeantBasket::Clear(Option_t *option)
{
// Clear basket;
   SetMixed(kFALSE);
   fTracksIn.Clear(option);
   fTracksOut.Clear(option);
}   

//______________________________________________________________________________
Bool_t GeantBasket::Contains(Int_t evstart, Int_t nevents) const
{
// Checks if any of the array of tracks belongs to the given event.
   return fTracksIn.Contains(evstart, nevents);
}      

//______________________________________________________________________________
TGeoVolume *GeantBasket::GetVolume() const
{
// Returns volume for this basket
   return fManager->GetVolume();
}   

//______________________________________________________________________________
void GeantBasket::Print(Option_t *) const
{
// Print basket content.
   Printf("*** basket %s: ninput=%3d   noutput=%3d", GetName(), GetNinput(), GetNoutput());
}

//______________________________________________________________________________
void GeantBasket::PrintTrack(Int_t /*itr*/, Bool_t /*input*/) const
{
// Print a given track.
}

//______________________________________________________________________________
void GeantBasket::Recycle()
{
// Recycle the basket to the volume scheduler.
   fManager->RecycleBasket(this);   
}

//______________________________________________________________________________
// Basket manager for a given volume. Holds a list of free baskets stored in a
// concurrent queue
//______________________________________________________________________________

ClassImp(GeantBasketMgr)

//______________________________________________________________________________
GeantBasketMgr::GeantBasketMgr(GeantScheduler *sch, TGeoVolume *vol, Int_t number)
                  :TGeoExtension(),
                   fScheduler(sch),
                   fVolume(vol),
                   fNumber(number),
                   fBcap(0),
                   fThreshold(0),
                   fNbaskets(0),
                   fNused(0),
                   fCBasket(0),
                   fPBasket(0),
                   fLock(),
                   fBaskets(2<<15),
                   fFeeder(0),
                   fMutex()
{
// Constructor
   fBcap = GeantPropagator::Instance()->fNperBasket + 1;
   SetCBasket(GetNextBasket());
   SetPBasket(GetNextBasket());
}

//______________________________________________________________________________
GeantBasketMgr::~GeantBasketMgr()
{
// Clean up
//   delete fCBasket;
//   delete fPBasket;
}   

//______________________________________________________________________________
GeantBasket *GeantBasketMgr::StealAndReplace(atomic_basket &current)
{
// Steal the current pointer content and replace with a new basket from the 
// pool. If the operation succeeds, returns the released basket which is now
// thread local
   GeantBasket *newb, *oldb;
   // Backup cbasket. If a steal happened it does not matter
   newb = GetNextBasket();  // backup new basket to avoid leak
   oldb = current.load(); // both local vars, so identical content
   // Check if fPBasket pointer was stolen by other thread
   while (fLock.test_and_set(std::memory_order_acquire));
   Bool_t stolen = current.compare_exchange_strong(oldb, newb);
   fLock.clear(std::memory_order_release);   
   if (stolen) {
      // No steal: current points now to newbasket which takes all new 
      // tracks. We can push the backed-up old basket to the queue, but
      // only AFTER any ongoing track addition finishes.
      while (oldb->IsAddingOp()) {}; // new AddTrack go all to newbasket
      // Someone may still have a copy of the old basket  and haven't
      // started yet adding tracks to it
      // We are finally the owner of oldb
//      assert (!oldb->IsAddingOp());
      return oldb;
   } else {
      // Current basket stolen by other thread: we need to recicle the new one 
      // and ignore the backed-up old basket which was injected by the other thread
//      assert(!newb->IsAddingOp());
      newb->Recycle();
   }
   return 0;
}   

//______________________________________________________________________________
GeantBasket *GeantBasketMgr::StealAndPin(atomic_basket &current)
{
// The method pins non-atomically the basket to the caller while adding the
// fAddingOp flag. It makes sure that no StealAndReplace operation happened
// on the pinned basket.
   GeantBasket *oldb = 0;
   while (1) {
      // the 2 lines below should be atomically coupled
      while (fLock.test_and_set(std::memory_order_acquire));
      oldb = current.load();
      oldb->LockAddingOp();
      fLock.clear(std::memory_order_release);
      // If no steal happened in between return basket pointer
      if (oldb == current.load()) return oldb;
      oldb->UnLockAddingOp();
   }
}

//______________________________________________________________________________
Bool_t GeantBasketMgr::StealMatching(atomic_basket &global, GeantBasket *content)
{
// Steal the global basket if it has the right matching content
   // Prepare replacement
   GeantBasket *newb = GetNextBasket();
   while (fLock.test_and_set(std::memory_order_acquire));
   Bool_t stolen = global.compare_exchange_strong(content, newb);
   fLock.clear(std::memory_order_release);
   if (stolen) {
      // Basket stolen
      while (content->IsAddingOp()) {};
      return kTRUE;
   } else {
      newb->Recycle();
   }
   return kFALSE;
}         

//______________________________________________________________________________
Int_t GeantBasketMgr::AddTrack(GeantTrack_v &trackv, Int_t itr, Bool_t priority)
{
// Copy directly from a track_v a track to the basket manager.
// Has to work concurrently
   // Atomically pin the basket for the adding operation
   GeantBasket *oldb = StealAndPin(priority?fPBasket:fCBasket);
   // Now basket matches fP(C)Basket content and has the adding flag set  
   oldb->AddTrack(trackv, itr);
   if (oldb->GetNinput() >= fThreshold) {
      oldb->UnLockAddingOp();
      if (StealMatching(priority?fPBasket:fCBasket, oldb)) { 
         // we fully own now oldb
//         assert(!oldb->IsAddingOp());
         fFeeder->push(oldb, priority);         
         return 1;
      } else return 0;   
   }
   oldb->UnLockAddingOp();
   return 0;      
}

//______________________________________________________________________________
Int_t GeantBasketMgr::AddTrack(GeantTrack &track, Bool_t priority)
{
// Add a track to the volume basket manager. If the track number reaches the
// threshold, the basket is added to the feeder queue and replaced by an empty 
// one. The feeder must be defined beforehand. Returns the number of dispatched
// baskets
// Has to work concurrently
   // Atomically pin the basket for the adding operation
   GeantBasket *oldb = StealAndPin(priority?fPBasket:fCBasket);
   // Now basket matches fP(C)Basket content and has the adding flag set  
   oldb->AddTrack(track);
   if (oldb->GetNinput() >= fThreshold) {
      oldb->UnLockAddingOp();
      if (StealMatching(priority?fPBasket:fCBasket, oldb)) { 
         // we fully own now oldb
//         assert(!oldb->IsAddingOp());
         fFeeder->push(oldb, priority);         
         return 1;
      } else return 0;   
   }
   oldb->UnLockAddingOp();
   return 0;      
}

//______________________________________________________________________________
Int_t GeantBasketMgr::CollectPrioritizedTracks(Int_t evmin, Int_t evmax)
{
// Move current basket tracks to priority one. 
// *** NONE *** This should be done for all basket managers only once when 
// starting prioritizing an event range.
   // Lock and swap containers
   if (!GetPBasket()->GetNinput() && !GetCBasket()->GetNinput()) return 0;
//   Printf("=== CollectPrioritized");
//   Print();
   GeantBasket *pbasket=0, *cbasket=0, *basket;
   Int_t npush = 0;
   // We want to steal fPBasket and fCBasket
   while (pbasket==0) pbasket = StealAndReplace(fPBasket);
   while (cbasket==0) cbasket = StealAndReplace(fCBasket);
   // pbasket and cbasket are now thread local
   GeantTrack_v &tracks = cbasket->GetInputTracks();
   Int_t ntracks = tracks.GetNtracks();
   // Dump cbasket into pbasket if it contains tracks with priority and
   // inject the latter
   for (Int_t itr=0; itr<ntracks; itr++) {
      if (tracks.fEventV[itr]>=evmin && tracks.fEventV[itr]<=evmax) {
         pbasket->GetInputTracks().AddTracks(tracks, 0, ntracks-1);
         tracks.Clear();
         cbasket->Recycle();
//         assert(!pbasket->IsAddingOp());
         fFeeder->push(pbasket, kTRUE);
//         Printf("   -> pushed %d tracks", pbasket->GetNinput());
//         Print();
         return 1;
      }
   }
   // pbasket may contain priority tracks -> inject
   if (pbasket->GetNinput()) {
//      Printf("   -> pushed %d tracks", pbasket->GetNinput());
//      assert(!pbasket->IsAddingOp());
      fFeeder->push(pbasket, kTRUE);
      npush++;
   } else {
      pbasket->Recycle();
   }   
   // if cbasket empty -> just recycle
   if (cbasket->GetNinput() == 0) {
      cbasket->Recycle();
//         Print();
      return npush;
   }   
   // cbasket has to be merged back into fCBasket. Most profitable is to
   // exchange the pointer held by fCBasket (which is most likely empty if
   // no other thread was adding on top) with cbasket.
   cbasket = fCBasket.exchange(cbasket);
   while (cbasket->IsAddingOp()) {}; // wait possible adding to finish
   // Now add cbasket content on top of fCBasket
   GeantTrack_v &ctracks = cbasket->GetInputTracks();
   ntracks = ctracks.GetNtracks();
   for (Int_t itr=0; itr<ntracks; itr++) {
      basket = StealAndPin(fCBasket);
      basket->AddTrack(ctracks, itr);
      basket->UnLockAddingOp();
   }   
   ctracks.Clear();
   cbasket->Recycle();
//   Print();
   return npush;
}   

//______________________________________________________________________________
Int_t GeantBasketMgr::FlushPriorityBasket()
{
// Flush the baskets containing tracks. Returns the number of dispatched baskets.
   if (!GetPBasket()->GetNinput()) return 0;
// Just steal the priority basket
   GeantBasket *pbasket=0;
   while (pbasket==0) pbasket = StealAndReplace(fPBasket);
   if (pbasket->GetNinput()) {
//      assert(!pbasket->IsAddingOp());
      fFeeder->push(pbasket, kTRUE);
      return 1;
   }   
   pbasket->Recycle();
   return 0;
}   

//______________________________________________________________________________
Int_t GeantBasketMgr::GarbageCollect()
{
// Copy all priority tracks to the current basket and flush to queue
   GeantBasket *pbasket=0, *cbasket=0;
   if (GetPBasket()->GetNinput()) {
      // We want to steal fPBasket and fCBasket
      while (pbasket==0) pbasket = StealAndReplace(fPBasket);
      while (cbasket==0) cbasket = StealAndReplace(fCBasket);
      // pbasket and cbasket are now thread local -> we can copy all tracks
      GeantTrack_v &tracks = pbasket->GetInputTracks();
      Int_t ntracks = pbasket->GetNinput();
      if (ntracks) {
         // we own cbasket
         cbasket->GetInputTracks().AddTracks(tracks, 0, ntracks-1);
         tracks.Clear();
      }
      pbasket->Recycle();
//      assert(!cbasket->IsAddingOp());
      fFeeder->push(cbasket, kFALSE);
      return 1;
   } else {
      // We want to steal fCBasket
      if (GetCBasket()->GetNinput()) {
         while (cbasket==0) cbasket = StealAndReplace(fCBasket);         
         Int_t ntracks = cbasket->GetNinput();
         if (ntracks) {
//            assert(!cbasket->IsAddingOp());
            fFeeder->push(cbasket, kFALSE);
            return 1;
         }
         cbasket->Recycle();
      }
   }
   return 0;
}      

//______________________________________________________________________________
GeantBasket *GeantBasketMgr::GetNextBasket()
{
// Returns next empy basket if any available, else create a new basket.
//   GeantBasket *next = fBaskets.try_pop();
   GeantBasket *next = 0;
   Bool_t pulled = fBaskets.dequeue(next);
   if (!pulled) {
      next = new GeantBasket(fBcap, this);
      // === critical section if atomics not supported ===
      fNbaskets++;
      fNused++;
      // === end critical section ===
   } else {
      // === critical section if atomics not supported ===
      fNused++;
      // === end critical section ===
   }
   return next;
}

//______________________________________________________________________________
void GeantBasketMgr::RecycleBasket(GeantBasket *b)
{
// Recycle a basket.
//   assert(!b->GetNinput());
//   assert(!b->IsAddingOp()); 
   b->Clear();
   if (b->GetInputTracks().Capacity() < fThreshold) {
      b->GetInputTracks().Resize(fThreshold);
      // Resize also the output array if needed
      if (b->GetOutputTracks().Capacity() < fThreshold)
        b->GetOutputTracks().Resize(fThreshold); 
   }
   if (!fBaskets.enqueue(b)) {
      Printf("Fatal error: exceeded the size of the bounded queue for basket mgr: %s",
      b->GetName());
      exit(1);
   }
   fNused--;
}   

//______________________________________________________________________________
void GeantBasketMgr::Print(Option_t *) const
{
// Print info about the basket content.
   Printf("Bsk_mgr %s: priority: in=%d out=%d     current: in=%d out=%d", GetName(),
          GetPBasket()->GetNinput(), GetPBasket()->GetNoutput(),
          GetCBasket()->GetNinput(), GetCBasket()->GetNoutput());
}

