// this file is part of a demonstrator for SC15
// it mimics the possibilities of future VecGeom navigation
// author: Sandro Wenzel (sandro.wenzel@cern); October 2015

// some VecGeom includes
#include "base/SOA3D.h"
#include "navigation/NavigationState.h"
#include "navigation/NavStatePool.h"
#include "base/Vector3D.h"

using namespace vecgeom;

//! base class defining basic interface for navigation ( hit-detection )
//! sub classes implement optimized algorithms for logical volumes
class VNavigator {

public:
  //! computes the step (distance) to the next object in the geometry hierarchy obtained
  //! by propagating with step along the ray
  //! the next object could be after a boundary and does not necessarily coincide with the object
  //! hit by the ray

    //! this methods transforms the global coordinates into local ones usually calls more specialized methods
    //! like the hit detection on local coordinates
  virtual double ComputeStepAndPropagatedState(Vector3D<double> const & /*globalpoint*/,
                                             Vector3D<double> const & /*globaldir*/,
                                             double /*(physics) step limit */,
                                             NavigationState const & /*in_state*/,
                                             NavigationState & /*out_state*/) const = 0;

  virtual double ComputeStepAndHittingBoundaryForLocalPoint(Vector3D<double> const & /*localpoint*/,
                                                            Vector3D<double> const & /*globaldir*/,
                                                            double /*(physics) step limit */,
                                                            NavigationState const & /*in_state*/,
                                                            NavigationState & /*out_state*/) const
  {
    return 0.0;
  }

  // for vector navigation
  //! this methods transforms the global coordinates into local ones usually calls more specialized methods
  //! like the hit detection on local coordinates
  virtual void ComputeStepsAndPropagatedStates(SOA3D<double> const & /*globalpoints*/,
                                             SOA3D<double> const & /*globaldirs*/,
                                             double const * /*(physics) step limits */,
                                             NavStatePool const & /*in_states*/,
                                             NavStatePool & /*out_states*/, double * /*out_steps*/) const {};

protected:
  // a common relocate method ( to calculate propagated states after the boundary )
  virtual void Relocate(Vector3D<double> /*localpoint*/, NavigationState const & /*in_state*/,
                        NavigationState & /*out_state*/) const = 0;

  // a common function to be used by all navigators to ensure consistency in transporting points
  // after a boundary
  VECGEOM_INLINE
  Vector3D<double> MovePointAfterBoundary( Vector3D<double> const &localpoint, Vector3D<double> const &dir, double step ) const {
      const double extra=1E-6;
      return localpoint + (step+extra)*dir;
  }


public:
  virtual ~VNavigator() {};

  // get name of implementing class
  virtual const char *GetName() const = 0;
};

//! template class providing a standard implementation for
//! some interfaces in VSafetyEstimator (using the CRT pattern)
template <typename Impl>
class VNavigatorHelper : public VNavigator {

public:
    // a generic implementation
    // may be completely specialized in child-navigators
    virtual
    double
    ComputeStepAndPropagatedState(Vector3D<double> const & globalpoint,
                                Vector3D<double> const & globaldir,
                                double step_limit,
                                NavigationState const & in_state,
                                NavigationState & out_state) const override {
      /*
        // calculate local point from global point
      Transformation3D m;
      in_state.TopMatrix(m);
      auto localpoint = m.Transform(globalpoint);
      auto localdir = m.TransformDirection(globaldir);
      // "suck in" algorithm from Impl and treat hit detection in local coordinates
      double step = ((Impl*)this)->Impl::ComputeStepAndHittingBoundaryForLocalPoint(
              localpoint, localdir,
              step_limit, in_state, out_state);

      // step was physics limited
      if( ! out_state.IsOnBoundary() )
        return step;

      // otherwise if necessary do a relocation

      // try relocation to refine out_state to correct location after the boundary
      ((Impl*)this)->Impl::Relocate( MovePointAfterBoundary( localpoint, localdir, step ), in_state, out_state);
      return step;
    */
      return 0.0;
    }

protected:
  // a common relocate method ( to calculate propagated states after the boundary )
  virtual void Relocate(Vector3D<double> pointafterboundary, NavigationState const & in_state,
                        NavigationState & out_state) const override {
/*
    not implemented since not needed here
*/
  }

public:
  static const char *GetClassName() { return Impl::gClassNameString; }

  virtual const char *GetName() const override { return GetClassName(); }
}; // end class VNavigatorHelper


