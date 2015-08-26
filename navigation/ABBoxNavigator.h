/*
 * ABBoxNavigator.h
 *
 *  Created on: 24.04.2015
 *      Author: swenzel
 */

#pragma once

#include "base/Global.h"

#include "volumes/PlacedVolume.h"
#include "base/Vector3D.h"
#include "management/GeoManager.h"
#include "navigation/NavigationState.h"
#include "base/Transformation3D.h"
#include "volumes/kernel/BoxImplementation.h"
#include "navigation/SimpleNavigator.h"

#ifdef VECGEOM_VC
#include "backend/vc/Backend.h"
#endif

#include <vector>
#include <cassert>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {


// A navigator using aligned bounding box = ABBox (hierarchies) to quickly find
// potential hit targets.
// This navigator goes into the direction of "voxel" navigators used in Geant4
// and ROOT. Checking single-rays against a set of aligned bounding boxes can be done
// in a vectorized fashion.
class ABBoxNavigator
{

public:
  VECGEOM_CUDA_HEADER_BOTH
  VECGEOM_INLINE
  ABBoxNavigator(){}


  int GetHitCandidates( LogicalVolume const * lvol,
          Vector3D<Precision> const & point,
          Vector3D<Precision> const & dir,
          ABBoxManager::ABBoxContainer_t const & corners, int size,
          ABBoxManager::HitContainer_t & hitlist
  ) const;

  int GetHitCandidates_v( LogicalVolume const * lvol,
            Vector3D<Precision> const & point,
            Vector3D<Precision> const & dir,
            ABBoxManager::ABBoxContainer_v const & corners, int size,
            ABBoxManager::HitContainer_t & hitlist
  ) const;


  void GetSafetyCandidates_v( Vector3D<Precision> const & /* point */,
                            ABBoxManager::ABBoxContainer_v const & /* corners */,
                            int size,
                            ABBoxManager::HitContainer_t & /* boxsafetypairs */,
                            Precision upperlimit ) const;


  // convert index to physical daugher
  VPlacedVolume const * LookupDaughter( LogicalVolume const *lvol, int id ) const {
      assert( id >= 0 && "access with negative index");
      assert( id < lvol->GetDaughtersp()->size() && "access beyond size of daughterlist ");
      return lvol->GetDaughtersp()->operator []( id );
  }

   /**
   * A function to navigate ( find next boundary and/or the step to do )
   */
   VECGEOM_CUDA_HEADER_BOTH
   VECGEOM_INLINE
   void FindNextBoundaryAndStep( Vector3D<Precision> const & /* global point */,
                          Vector3D<Precision> const & /* global dir */,
                          NavigationState const & /* currentstate */,
                          NavigationState & /* newstate */,
                          Precision const & /* proposed physical step */,
                          Precision & /*step*/
                         ) const;

   /**
    * A function to get back the safe distance; given a NavigationState object and a current global point
    * point
    */
   VECGEOM_CUDA_HEADER_BOTH
   VECGEOM_INLINE
   Precision GetSafety( Vector3D<Precision> const & /*global_point*/,
               NavigationState const & /* currentstate */
   ) const;

   // NOTE: there is no vector interface here yet --> this is part of SimpleNavigator


}; // end of class declaration

//#define VERBOSE
void
ABBoxNavigator::FindNextBoundaryAndStep( Vector3D<Precision> const & globalpoint,
                                          Vector3D<Precision> const & globaldir,
                                          NavigationState     const & currentstate,
                                          NavigationState           & newstate,
                                          Precision           const & pstep,
                                          Precision                 & step
                                        ) const
{
   // this information might have been cached in previous navigators??
#ifdef VERBOSE
    static int counter = 0;
    if( counter % 1 == 0 )
    std::cerr << counter << " " << globalpoint << " \n";

    counter++;
#endif

   Transformation3D m;
   currentstate.TopMatrix(m);
   Vector3D<Precision> localpoint=m.Transform(globalpoint);
   Vector3D<Precision> localdir=m.TransformDirection(globaldir);

   VPlacedVolume const * currentvolume = currentstate.Top();
   int nexthitvolume = -1; // means mother

   // StepType st = kPhysicsStep; // physics or geometry step
   step = currentvolume->DistanceToOut( localpoint, localdir, pstep );

   // NOTE: IF STEP IS NEGATIVE HERE, SOMETHING IS TERRIBLY WRONG. WE CAN TRY TO HANDLE THE SITUATION
   // IN TRYING TO PROPOSE THE RIGHT LOCATION IN NEWSTATE AND RETURN
   // I WOULD MUCH FAVOUR IF THIS WAS DONE OUTSIDE OF THIS FUNCTION BY THE USER
    if( step < 0. )
    {
       // TODO: instead of directly exiting we could see whether we hit a daughter
       // which is usally a logic thing to do
      // std::cerr << "negative DO\n";
     //  step = 0.;
     //  currentstate.CopyTo(&newstate);
     //  newstate.Pop();
     //  SimpleNavigator nav;
     //  nav.RelocatePointFromPath( localpoint, newstate );
      // return;
        step = kInfinity;
    }

   // if( step > 1E20 )
   //     std::cerr << "infinite DO\n";
   // TODO: compare steptoout and physics step and take minimum



   // do a quick and vectorized search using aligned bounding boxes
   // obtains a sorted container ( vector or list ) of hitboxstructs
   LogicalVolume const * currentlvol = currentstate.Top()->GetLogicalVolume();

#ifdef VERBOSE
   std::cerr << " I am in " << currentlvol->GetLabel() << "\n";
#endif
   if( currentlvol->GetDaughtersp()->size() > 0 ){
#ifdef VERBOSE
       std::cerr << " searching through " << currentlvol->GetDaughtersp()->size() << " daughters\n";
#endif
     ABBoxManager::HitContainer_t & hitlist = ABBoxManager::Instance().GetAllocatedHitContainer();
//       hitlist.clear();
       int size;
//       ABBoxManager::ABBoxContainer_t bboxes1 =  ABBoxManager::Instance().GetABBoxes( currentlvol , size );
//       GetHitCandidates( currentlvol,
//                         localpoint,
//                         localdir,
//                         bboxes1,
//                        size, hitlist );
#ifdef VERBOSE
       int c1 = hitlist.size();
      std::cerr << hitlist << "\n";
#endif
       hitlist.clear();
       ABBoxManager::ABBoxContainer_v bboxes =  ABBoxManager::Instance().GetABBoxes_v( currentlvol , size );
            GetHitCandidates_v( currentlvol,
                          localpoint,
                          localdir,
                          bboxes,
                          size, hitlist );
#ifdef VERBOSE
            int c2 = hitlist.size();
        std::cerr << hitlist << "\n";
        std::cerr << " hitting scalar " << c1 << " vs vector " << c2 << "\n";
 if( c1 != c2 )
     std::cerr << "HUHU " << c1 << " " << c2;
        #endif

        // sorting the histlist
        ABBoxManager::sort( hitlist, ABBoxManager::HitBoxComparatorFunctor() );

        // assumption: here hitlist is sorted in ascending distance order
#ifdef VERBOSE
        std::cerr << " hitting " << hitlist.size() << " boundary boxes\n";
#endif
        for( auto hitbox : hitlist )
        {
             VPlacedVolume const * candidate = LookupDaughter( currentlvol, hitbox.first );

            // only consider those hitboxes which are within potential reach of this step
            if( ! ( step < hitbox.second )) {
            //      std::cerr << "checking id " << hitbox.first << " at box distance " << hitbox.second << "\n";
             if( hitbox.second < 0 ){
                bool checkindaughter = candidate->Contains( localpoint );
                if( checkindaughter == true ){
                    // need to relocate
                    step = 0;
                    nexthitvolume = hitbox.first;
                    // THE ALTERNATIVE WOULD BE TO PUSH THE CURRENT STATE AND RETURN DIRECTLY
                    break;
                }
            }
            Precision ddistance = candidate->DistanceToIn( localpoint, localdir, step );
#ifdef VERBOSE
            std::cerr << "distance to " << candidate->GetLabel() << " is " << ddistance << "\n";
#endif
            nexthitvolume = (ddistance < step) ? hitbox.first : nexthitvolume;
            step      = (ddistance < step) ? ddistance  : step;
        }
      else
      {
          break;
      }
   }
   }

   // now we have the candidates
   // try
   currentstate.CopyTo(&newstate);

   // is geometry further away than physics step?
   // not necessarily true
   if(step > pstep)
   {
       assert( true && "impossible state");
       // don't need to do anything
       step = pstep;
       newstate.SetBoundaryState( false );
       return;
   }
   newstate.SetBoundaryState( true );

   assert( step >= 0 && "step negative");

   if( step > 1E30 )
     {
      //std::cout << "WARNING: STEP INFINITY; should never happen unless outside\n";
           //InspectEnvironmentForPointAndDirection( globalpoint, globaldir, currentstate );

           // set step to zero and retry one level higher
           step = 0;
           newstate.Pop();
           return;
      }

      if( step < 0. )
      {
        //std::cout << "WARNING: STEP NEGATIVE\n";
        //InspectEnvironmentForPointAndDirection( globalpoint, globaldir, currentstate );
         step = 0.;
      }

   // TODO: this is tedious, please provide operators in Vector3D!!
   // WE SHOULD HAVE A FUNCTION "TRANSPORT" FOR AN OPERATION LIKE THIS
   Vector3D<Precision> newpointafterboundary = localdir;
   newpointafterboundary*=(step + 1e-6);
   newpointafterboundary+=localpoint;

   if( nexthitvolume != -1 ) // not hitting mother
   {
      // continue directly further down
      VPlacedVolume const * nextvol = LookupDaughter( currentlvol, nexthitvolume );
      Transformation3D const * trans = nextvol->GetTransformation();

      SimpleNavigator nav;
      nav.LocatePoint( nextvol, trans->Transform(newpointafterboundary), newstate, false );
      assert( newstate.Top() != currentstate.Top() && " error relocating when entering ");
      return;
   }
   else // hitting mother
   {
      SimpleNavigator nav;
      nav.RelocatePointFromPath( newpointafterboundary, newstate );

      // can I push particle ?
      // int correctstep = 0;
      while( newstate.Top() == currentstate.Top() )
      {
     //     newstate.Print();
     //     step+=1E-6;
     //     SimpleNavigator nav;
     //     newstate.Clear();
     //     nav.LocatePoint( GeoManager::Instance().GetWorld(), globalpoint + (step)*globaldir, newstate, true );
     //     std::cerr << "correcting " << correctstep << " remaining dist to out "
      //              << currentvolume->DistanceToOut( localpoint + step*localdir, localdir, pstep )
      //              << " " << currentvolume->Contains( localpoint + step*localdir )
      //    << " " << currentvolume->SafetyToIn( localpoint + step*localdir )
      //    << " " << currentvolume->SafetyToOut( localpoint + step*localdir ) << "\n";
      //    currentvolume->PrintType();

      //    correctstep++;
       //   std::cerr << "Matrix error " << const_cast<NavigationState &> ( currentstate ).CalcTransformError( globalpoint, globaldir );
        newstate.Pop();
      }
//      if( newstate.Top() == currentstate.Top() )
//      {
//         std::cerr << "relocate failed; trying to locate from top for step " << step << "\n";
//         newstate.Clear();
//         SimpleNavigator nav;
//         nav.LocatePoint( GeoManager::Instance().GetWorld(), globalpoint + (step+1E-6)*globaldir, newstate, true );
//         //  std::cerr << "newstate top " << newstate.Top()->GetLabel() << "\n";
//      }
//      if( newstate.Top() == currentstate.Top() )
//      {
//         SimpleNavigator nav;
//         nav.InspectEnvironmentForPointAndDirection( globalpoint, globaldir, currentstate );
//      }
      assert( newstate.Top() != currentstate.Top() && " error relocating when leaving ");
   }
}

Precision ABBoxNavigator::GetSafety(Vector3D<Precision> const & globalpoint,
                            NavigationState const & currentstate) const
{
   // this information might have been cached already ??
   Transformation3D m;
   currentstate.TopMatrix(m);
   Vector3D<Precision> localpoint=m.Transform(globalpoint);

   VPlacedVolume const *currentvol = currentstate.Top();
   double safety = currentvol->SafetyToOut( localpoint );
   double safetysqr = safety*safety;

   // safety to bounding boxes
   LogicalVolume const *lvol = currentvol->GetLogicalVolume();
   if( safety > 0. && lvol->GetDaughtersp()->size() > 0 ){
       ABBoxManager &instance = ABBoxManager::Instance();
       ABBoxManager::HitContainer_t &boxsafetylist = instance.GetAllocatedHitContainer();

       int size;
       boxsafetylist.clear();

       ABBoxManager::ABBoxContainer_v bboxes =  instance.GetABBoxes_v( lvol , size );
       // calculate squared bounding box safeties in vectorized way
       GetSafetyCandidates_v(localpoint, bboxes, size, boxsafetylist, safetysqr);

       // sorting the list
       ABBoxManager::sort( boxsafetylist, ABBoxManager::HitBoxComparatorFunctor() );

       // at this moment boxsafetylist only contains
       // elements whose "bounding box" safetysqr is smaller than safetytooutsqr and which hence have to be checked
#ifdef VERBOSE
       std::cerr << "boxsafetylist has " << boxsafetylist.size() << " candidates \n";
#endif
       for( auto boxsafetypair : boxsafetylist ){
       if( boxsafetypair.second < safetysqr ){
	 //	 std::cerr << " id " << boxsafetypair.first << " safetysqr " << boxsafetypair.second << "\n"; 
         VPlacedVolume const *candidate = LookupDaughter( lvol, boxsafetypair.first );
         auto candidatesafety = candidate->SafetyToIn( localpoint );
#ifdef VERBOSE
         if( candidatesafety*candidatesafety > boxsafetypair.second && boxsafetypair.second > 0 )
             std::cerr << "real safety smaller than boxsafety \n";
#endif
         if( candidatesafety < safety ){
           safety = candidatesafety;
           safetysqr = safety*safety;
         }
         else { // this box has a safety which is larger than the best known safety so we can stop here
           break;
         }
       }
       }}
       return safety;
 }


} } // End global namespace
