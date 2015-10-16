/*
 * SC15Navigators.h
 *
 *  File implementating specialized navigators for the SimpleTracker.cpp example
 *  in this directory
 *
 *  Created on: Oct 14, 2015
 *      Author: swenzel
 */

#include "VNavigator.h"

#include "volumes/Box.h"
#include "volumes/Tube.h"
#include "navigation/NavigationState.h"
#include "volumes/PlacedVolume.h"
#include "backend/Backend.h"


#pragma once

#define INSTANCE(name)                \
    static VNavigator *Instance() {   \
      static name instance;           \
      return &instance;}

#define CLASSNAME(name)\
        static constexpr char *gClassNameString = #name;

#define BOILERPLATE(name)     \
    CLASSNAME(name)
//    INSTANCE(name)





// a navigator for the world volume
// this navigator knows that it is in a BoxWorld and that it has 1 tube daughter
template <bool fastcopy=false>
class WorldNavigator : public VNavigatorHelper<WorldNavigator<fastcopy>> {

public:
    BOILERPLATE(WorldNavigator)

//    static VNavigator *Instance() {
//      static WorldNavigator<fastcopy> instance;
//      return &instance;}
    INSTANCE(WorldNavigator)

    // a templated kernel to be used for scalar and vector interface
    template <typename Backend>
    static void WorldNavKernel(UnplacedBox const &mother, UnplacedTube const &daughter,
                               Vector3D<typename Backend::precision_v> const &pos,
                               Vector3D<typename Backend::precision_v> const &dir,
                               typename Backend::precision_v &distance, typename Backend::precision_v &distance1) {
          BoxImplementation<translation::kIdentity, rotation::kIdentity>::DistanceToOut<Backend>(
              mother, pos, dir, kInfinity, distance);

          TubeImplementation<translation::kIdentity, rotation::kIdentity, TubeTypes::UniversalTube>::DistanceToIn<Backend>(
              daughter, Transformation3D::kIdentity, pos, dir, kInfinity, distance1);
    }


    virtual double ComputeStepAndPropagatedState(Vector3D<double> const & globalpoint,
                                                 Vector3D<double> const & globaldir,
                                                 double pstep,
                                                 NavigationState const & in_state,
                                                 NavigationState & out_state) const override {
        // NO NEED TO TRANSFORM HERE

        // DISTANCETOOUT
        // actually this could be lazily done after distance to daughter because we are convex
        auto unplaced = in_state.Top()->GetUnplacedVolume();
        double distance;

        // DISTANCETODAUGHTER
        // daughter is a very specific tube:
        auto daughters= in_state.Top()->GetLogicalVolume()->GetDaughtersp();
        auto daughter = (*daughters)[0];
        double distance1;
        WorldNavKernel<kScalar>(*static_cast<UnplacedBox const*>(unplaced),
                *static_cast<UnplacedTube const*>(daughter->GetUnplacedVolume()), globalpoint,globaldir,distance,distance1);

        // we could speed this up if we knew the depth of this volume
        if(!fastcopy)
        {
          // we could speed this up if we knew the depth of this volume
          in_state.CopyTo(&out_state);
        }
        else {
          // do a very fast copy using the precomputed size of a NavigationState at this depth
          in_state.CopyToFixedSize<NavigationState::SizeOf(1)>(&out_state);
        }

        // do a very efficient relocation depending on distance
        if(distance1 < distance){
            out_state.Push(daughter);
            return distance1;
        }
        else // leaving mother ( have to be careful - we could actually leave the world volume )
        {
            out_state.Pop();
            return distance;
        }
    }

    // vector interface
  virtual void ComputeStepsAndPropagatedStates(SOA3D<double> const & globalpoints,
                                                 SOA3D<double> const & globaldirs,
                                                 double const *psteps,
                                                 NavStatePool const & in_states,
                                                 NavStatePool & out_states, double *out_steps) const override {
    using Backend = kVc;
    using Real_v = Backend::precision_v;
    int offset = globalpoints.size() - globalpoints.size() % Real_v::Size;
    auto unplaced = in_states[0]->Top()->GetUnplacedVolume();
    auto daughters= in_states[0]->Top()->GetLogicalVolume()->GetDaughtersp();
    auto daughter = (*daughters)[0];

    // vector part
    for (auto i = 0; i < offset; i += Real_v::Size) {
      Real_v distance, distance1;
      Vector3D<Real_v> p(Real_v(globalpoints.x() + i), Real_v(globalpoints.y() + i), Real_v(globalpoints.z() + i));
      Vector3D<Real_v> d(Real_v(globaldirs.x() + i), Real_v(globaldirs.y() + i), Real_v(globaldirs.z() + i));

      WorldNavKernel<Backend>(*static_cast<UnplacedBox const *>(unplaced),
                              *static_cast<UnplacedTube const *>(daughter->GetUnplacedVolume()), p, d,
                              distance, distance1);

      Min(distance,distance1).store(out_steps + i);
      for(unsigned int j=0;j<Real_v::Size;++j){
        // relocation here ( or in separate loop ? )
        if (!fastcopy) {
          in_states[i+j]->CopyTo(out_states[i+j]);
        } else {
          in_states[i+j]->CopyToFixedSize<NavigationState::SizeOf(1)>(out_states[i+j]);
        }

        // take a decision where to go
        if(distance1[j]<distance[j]) out_states[i+j]->Push(daughter);
        else
            out_states[i+j]->Pop();

      }
    }
    // tail part
    for (unsigned int i = offset; i < globalpoints.size(); ++i) {
      double distance, distance1;
      WorldNavKernel<kScalar>(*static_cast<UnplacedBox const *>(unplaced),
                              *static_cast<UnplacedTube const *>(daughter->GetUnplacedVolume()), globalpoints[i], globaldirs[i],
                              distance, distance1);

      // relocation here ( or in separate loop ? )
      if (!fastcopy) {
        in_states[i]->CopyTo(out_states[i]);
      } else {
        in_states[i]->CopyToFixedSize<NavigationState::SizeOf(1)>(out_states[i]);
      }
      // take a decision where to go
      if(distance1<distance) out_states[i]->Push(daughter);
            else
                out_states[i]->Pop();

      out_steps[i] = Min(distance1,distance);
    }
  }

};


// a navigator for the tube layers
// this navigator knows that it is in a tube and that it has 1 tube daughter
// we can use knowledge of the depth of the volume to speedup memcopy of NavigationStates
template <int Depth = -1>
class LayerNavigator : public VNavigatorHelper<LayerNavigator<Depth>> {
public:
    static VNavigator *Instance() {
          static LayerNavigator<Depth> instance;
          return &instance;}


    BOILERPLATE(LayerNavigator)

    // a templated kernel to be used for scalar and vector interface
    template <typename Backend>
    static void LayerNavKernel(UnplacedTube const &mother, UnplacedTube const &daughter,
                               Vector3D<typename Backend::precision_v> const &pos,
                               Vector3D<typename Backend::precision_v> const &dir,
                               typename Backend::precision_v &distance, typename Backend::precision_v &distance1) {
      TubeImplementation<translation::kIdentity, rotation::kIdentity, TubeTypes::UniversalTube>::DistanceToOut<Backend>(
          mother, pos, dir, kInfinity, distance);

      TubeImplementation<translation::kIdentity, rotation::kIdentity, TubeTypes::UniversalTube>::DistanceToIn<Backend>(
          daughter, Transformation3D::kIdentity, pos, dir, kInfinity, distance1);
    }

    virtual double ComputeStepAndPropagatedState(Vector3D<double> const &globalpoint, Vector3D<double> const &globaldir,
                                                 double pstep, NavigationState const &in_state,
                                                 NavigationState &out_state) const override {
      // NO NEED TO TRANSFORM HERE

      // DISTANCETOOUT
      // actually this could be done after distancetodaughter because we are convex
      auto unplaced = in_state.Top()->GetUnplacedVolume();
      double distance;
      double distance1;
      // daughter is a very specific tube:
      auto daughters= in_state.Top()->GetLogicalVolume()->GetDaughtersp();
      auto daughter = (*daughters)[0];
     
      LayerNavKernel<kScalar>(*static_cast<UnplacedTube const*>(unplaced), *static_cast<UnplacedTube const*>(daughter->GetUnplacedVolume()), globalpoint, globaldir, distance, distance1);

      if(Depth<0) /* depth is unknown */
      {
        in_state.CopyTo(&out_state);
      } else {
        // do a very fast copy using the precomputed size of a NavigationState at this depth
        in_state.CopyToFixedSize<NavigationState::SizeOf(Depth)>(&out_state);
      }

      // do a very efficient relocation depending on distance
      if(distance1 < distance){
          out_state.Push(daughter);
          return distance1;
      }
      else // leaving mother ( have to be careful - we could actually leave the world volume )
      {
          out_state.Pop();
          return distance;
      }
  }

    // vector interface
  virtual void ComputeStepsAndPropagatedStates(SOA3D<double> const & globalpoints,
                                               SOA3D<double> const & globaldirs,
                                               double const *psteps,
                                               NavStatePool const & in_states,
                                               NavStatePool & out_states, double *out_steps) const override {
    using Backend = kVc;
    using Real_v = Backend::precision_v;
    int offset = globalpoints.size() - globalpoints.size() % Real_v::Size;
    auto unplaced = in_states[0]->Top()->GetUnplacedVolume();
    auto daughters= in_states[0]->Top()->GetLogicalVolume()->GetDaughtersp();
    auto daughter = (*daughters)[0];

    // vector part
    for (auto i = 0; i < offset; i += Real_v::Size) {
      Real_v distance, distance1;
      Vector3D<Real_v> p(Real_v(globalpoints.x() + i), Real_v(globalpoints.y() + i), Real_v(globalpoints.z() + i));
      Vector3D<Real_v> d(Real_v(globaldirs.x() + i), Real_v(globaldirs.y() + i), Real_v(globaldirs.z() + i));

      LayerNavKernel<Backend>(*static_cast<UnplacedTube const *>(unplaced),
                              *static_cast<UnplacedTube const *>(daughter->GetUnplacedVolume()), p, d,
                              distance, distance1);

      Min(distance,distance1).store(out_steps + i);
      for(auto j=0;j<Real_v::Size;++j){
        // relocation here ( or in separate loop ? )
        if (Depth < 0) {
          in_states[i+j]->CopyTo(out_states[i+j]);
        } else {
          in_states[i+j]->CopyToFixedSize<NavigationState::SizeOf(Depth)>(out_states[i+j]);
        }

        // take a decision where to go
        if(distance1[j]<distance[j]) out_states[i+j]->Push(daughter);
        else
            out_states[i+j]->Pop();

      }
    }
    // tail part
    for (unsigned int i = offset; i < globalpoints.size(); ++i) {
      double distance, distance1;
      LayerNavKernel<kScalar>(*static_cast<UnplacedTube const *>(unplaced),
                              *static_cast<UnplacedTube const *>(daughter->GetUnplacedVolume()), globalpoints[i],
                              globaldirs[i], distance, distance1);

      // relocation here ( or in separate loop ? )
      if (Depth < 0) {
        in_states[i]->CopyTo(out_states[i]);
      } else {
        in_states[i]->CopyToFixedSize<NavigationState::SizeOf(Depth)>(out_states[i]);
      }
      // take a decision where to go
      if(distance1<distance) out_states[i]->Push(daughter);
            else
                out_states[i]->Pop();

      out_steps[i] = Min(distance1,distance);
    }
  }
};

// a navigator for the final tube ( no more daughters )
// this navigator knows that it is in a tube and that it has no further daughters
template <int Depth=-1>
class InnerMostTubeNavigator : public VNavigatorHelper<InnerMostTubeNavigator<Depth>> {

public:
    BOILERPLATE(InnerMostTubeNavigator)
static VNavigator *Instance() {
  static InnerMostTubeNavigator instance;
  return &instance;}

    virtual double ComputeStepAndPropagatedState(Vector3D<double> const & globalpoint,
                                       Vector3D<double> const & globaldir,
                                       double pstep,
                                       NavigationState const & in_state,
                                       NavigationState & out_state) const override {
      // NO NEED TO TRANSFORM HERE

      // DISTANCETOOUT
      auto unplaced = in_state.Top()->GetUnplacedVolume();
      double distance;
      TubeImplementation<translation::kIdentity, rotation::kIdentity, TubeTypes::UniversalTube>::DistanceToOut<kScalar>(
          *(UnplacedTube *)(unplaced), globalpoint, globaldir, kInfinity, distance);

      if (Depth < 0) {
        in_state.CopyTo(&out_state);
      } else {
        in_state.CopyToFixedSize<NavigationState::SizeOf(Depth)>(&out_state);
      }

      // always leaving mother ( have to be careful - we could actually leave the world volume )
      out_state.Pop();

      return distance;
    }

   virtual void ComputeStepsAndPropagatedStates(SOA3D<double> const & globalpoints,
                                             SOA3D<double> const & globaldirs,
                                             double const *psteps,
                                             NavStatePool const & in_states,
                                             NavStatePool & out_states, double *out_steps) const override { 
     using Backend = kVc;
     using Real_v = Backend::precision_v;
     int offset = globalpoints.size() - globalpoints.size() % Real_v::Size;
     auto unplaced = in_states[0]->Top()->GetUnplacedVolume();
     // vector part
     for(auto i=0; i<offset; i+=Real_v::Size){
       Real_v distance;
       Vector3D<Real_v> p(Real_v( globalpoints.x() + i ), Real_v( globalpoints.y() + i ), Real_v( globalpoints.z() + i ) ); 
       Vector3D<Real_v> d(Real_v( globaldirs.x() + i ), Real_v( globaldirs.y() + i ), Real_v( globaldirs.z() + i ) );
       TubeImplementation<translation::kIdentity, rotation::kIdentity,
                          TubeTypes::UniversalTube>::DistanceToOut<Backend>(*(UnplacedTube *)(unplaced), p, d,
                                                                            kInfinity, distance);

       distance.store(out_steps + i);

       for (auto j = 0; j < Real_v::Size; ++j) {
         // relocation here ( or in separate loop ? )
         if (Depth < 0) {
           in_states[i+j]->CopyTo(out_states[i+j]);
         } else {
           in_states[i+j]->CopyToFixedSize<NavigationState::SizeOf(Depth)>(out_states[i+j]);
         }
         out_states[i + j]->Pop();
       }
     }
     // tail part
     for( unsigned int i=offset; i < globalpoints.size(); ++i ){
       double distance;
       TubeImplementation<translation::kIdentity, rotation::kIdentity,
                          TubeTypes::UniversalTube>::DistanceToOut<kScalar>(*(UnplacedTube *)(unplaced),
                                                                            globalpoints[i], globaldirs[i], kInfinity,
                                                                            distance);

       // relocation here ( or in separate loop ? )
       if (Depth < 0) {
         in_states[i]->CopyTo(out_states[i]);
       } else {
         in_states[i]->CopyToFixedSize<NavigationState::SizeOf(Depth)>(out_states[i]);
       }
       out_steps[i] = i;

        // compute outstates
        out_states[i]->Pop();
     }
   }
};


