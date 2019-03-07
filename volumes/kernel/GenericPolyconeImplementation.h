/// @file GenericPolyconeImplementation.h
/// @author Raman Sehgal (raman.sehgal@cern.ch)

#ifndef VECGEOM_VOLUMES_KERNEL_GENERICPOLYCONEIMPLEMENTATION_H_
#define VECGEOM_VOLUMES_KERNEL_GENERICPOLYCONEIMPLEMENTATION_H_

#include "base/Vector3D.h"
#include "volumes/GenericPolyconeStruct.h"
#include "volumes/kernel/GenericKernels.h"
#include <VecCore/VecCore>
#include "volumes/kernel/CoaxialConesImplementation.h"
#include "volumes/kernel/shapetypes/ConeTypes.h"

#include <cstdio>

namespace vecgeom {

VECGEOM_DEVICE_FORWARD_DECLARE(struct GenericPolyconeImplementation;);
VECGEOM_DEVICE_DECLARE_CONV(struct, GenericPolyconeImplementation);

inline namespace VECGEOM_IMPL_NAMESPACE {

class PlacedGenericPolycone;
template <typename T>
struct GenericPolyconeStruct;
class UnplacedGenericPolycone;

struct GenericPolyconeImplementation {

  using PlacedShape_t    = PlacedGenericPolycone;
  using UnplacedStruct_t = GenericPolyconeStruct<double>;
  using UnplacedVolume_t = UnplacedGenericPolycone;

  VECCORE_ATT_HOST_DEVICE
  static void PrintType()
  {
    //  printf("SpecializedGenericPolycone<%i, %i>", transCodeT, rotCodeT);
  }

  template <typename Stream>
  static void PrintType(Stream &st, int transCodeT = translation::kGeneric, int rotCodeT = rotation::kGeneric)
  {
    st << "SpecializedGenericPolycone<" << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintImplementationType(Stream &st)
  {
    (void)st;
    // st << "GenericPolyconeImplementation<" << transCodeT << "," << rotCodeT << ">";
  }

  template <typename Stream>
  static void PrintUnplacedType(Stream &st)
  {
    (void)st;
    // TODO: this is wrong
    // st << "UnplacedGenericPolycone";
  }

  template <typename Real_v, bool ForInside>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void GenericKernelForASection(UnplacedStruct_t const &unplaced, int isect,
                                       Vector3D<Real_v> const &polyconePoint,
                                       typename vecCore::Mask_v<Real_v> &secFullyInside,
                                       typename vecCore::Mask_v<Real_v> &secFullyOutside)
  {

    // using namespace PolyconeTypes;
    using namespace ConeTypes;

    if (isect < 0) {
      secFullyInside  = false;
      secFullyOutside = true;
      return;
    }

    GenericPolyconeSection const &sec = unplaced.GetSection(isect);
    Vector3D<Precision> secLocalp     = polyconePoint - Vector3D<Precision>(0, 0, sec.fShift);
#ifdef POLYCONEDEBUG
    std::cout << " isect=" << isect << "/" << unplaced.GetNSections() << " secLocalP=" << secLocalp
              << ", secShift=" << sec.fShift << " sec.fSolid=" << sec.fSolid << std::endl;
    if (sec.fSolid) sec.fSolid->Print();
#endif

    CoaxialConesImplementation::template GenericKernelForContainsAndInside<Real_v, typename vecCore::Mask_v<Real_v>,
                                                                           ForInside>(*sec.fCoaxialCones, secLocalp,
                                                                                      secFullyInside, secFullyOutside);
  }

  template <typename Real_v, typename Bool_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void Contains(UnplacedStruct_t const &genericPolycone, Vector3D<Real_v> const &point, Bool_v &inside)
  {
    Bool_v unused(false), outside(false);
    GenericKernelForContainsAndInside<Real_v, Bool_v, false>(genericPolycone, point, unused, outside);
    inside = !outside;
  }

  // BIG QUESTION: DO WE WANT TO GIVE ALL 3 TEMPLATE PARAMETERS
  // -- OR -- DO WE WANT TO DEDUCE Bool_v, Index_t from Real_v???
  template <typename Real_v, typename Inside_t>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void Inside(UnplacedStruct_t const &genericPolycone, Vector3D<Real_v> const &point, Inside_t &inside)
  {

    using Bool_v       = vecCore::Mask_v<Real_v>;
    using InsideBool_v = vecCore::Mask_v<Inside_t>;
    Bool_v completelyinside(false), completelyoutside(false);
    GenericKernelForContainsAndInside<Real_v, Bool_v, true>(genericPolycone, point, completelyinside,
                                                            completelyoutside);
    inside = EInside::kSurface;
    vecCore::MaskedAssign(inside, (InsideBool_v)completelyoutside, Inside_t(EInside::kOutside));
    vecCore::MaskedAssign(inside, (InsideBool_v)completelyinside, Inside_t(EInside::kInside));
  }

  template <typename Real_v, typename Bool_v, bool ForInside>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void GenericKernelForContainsAndInside(UnplacedStruct_t const &unplaced, Vector3D<Real_v> const &localPoint,
                                                Bool_v &completelyInside, Bool_v &completelyOutside)
  {
    /* TODO : Logic to check where the point is inside or not.
    **
    ** if ForInside is false then it will only check if the point is outside,
    ** and is used by Contains function
    **
    ** if ForInside is true then it will check whether the point is inside or outside,
    ** and if neither inside nor outside then it is on the surface.
    ** and is used by Inside function
    */

    typedef Bool_v Bool_t;

    int indexLow  = unplaced.GetSectionIndex(localPoint.z() - kTolerance);
    int indexHigh = unplaced.GetSectionIndex(localPoint.z() + kTolerance);
    if (indexLow < 0 && indexHigh < 0) {
      completelyOutside = true;
      return;
    }
    if (indexLow < 0 && indexHigh == 0) {
      // Check location in section 0 and return
      GenericKernelForASection<Real_v, ForInside>(unplaced, 0, localPoint, completelyInside, completelyOutside);
      return;
    }
    if (indexHigh < 0 && indexLow == (unplaced.GetNSections() - 1)) {
      // Check location in section N-1 and return
      GenericKernelForASection<Real_v, ForInside>(unplaced, (unplaced.GetNSections() - 1), localPoint, completelyInside,
                                                  completelyOutside);

      return;
    }
    if (indexLow >= 0 && indexHigh >= 0) {
      if (indexLow == indexHigh) {
        // Check location in section indexLow and return
        GenericKernelForASection<Real_v, ForInside>(unplaced, indexLow, localPoint, completelyInside,
                                                    completelyOutside);

        return;
      } else {

        Bool_t secInLow = false, secOutLow = false;
        Bool_t secInHigh = false, secOutHigh = false;

        GenericPolyconeSection const &sectionLow  = unplaced.GetSection(indexLow);
        GenericPolyconeSection const &sectionHigh = unplaced.GetSection(indexHigh);

        Bool_t onRing(false);
        onRing |= (CoaxialConesImplementation::template IsOnRing<Real_v, false>(
                       *sectionLow.fCoaxialCones, localPoint - Vector3D<Precision>(0, 0, sectionLow.fShift)) ||
                   CoaxialConesImplementation::template IsOnRing<Real_v, true>(
                       *sectionHigh.fCoaxialCones, localPoint - Vector3D<Precision>(0, 0, sectionHigh.fShift)));

        GenericKernelForASection<Real_v, ForInside>(unplaced, indexLow, localPoint, secInLow, secOutLow);
        GenericKernelForASection<Real_v, ForInside>(unplaced, indexHigh, localPoint, secInHigh, secOutHigh);
        Bool_t surfLow  = !secInLow && !secOutLow;
        Bool_t surfHigh = !secInHigh && !secOutHigh;

        if (surfLow && surfHigh) {
          completelyInside = !onRing; // true;
          return;
        } else {
          // else if point is on surface of only one of the two sections then point is actually on surface , the default
          // case,
          // so no need to check

          // What needs to check is if it is outside both ie. Outside indexLow section and Outside indexHigh section
          // then it is certainly outside
          if (secOutLow && secOutHigh) {
            completelyOutside = true;
            return;
          }
        }
      }
    }
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void DistanceToIn(UnplacedStruct_t const &polycone, Vector3D<Real_v> const &point,
                           Vector3D<Real_v> const &direction, Real_v const &stepMax, Real_v &distance)
  {
    // using namespace PolyconeTypes;
    Vector3D<Real_v> p = point;
    Vector3D<Real_v> v = direction;

#ifdef POLYCONEDEBUG
    std::cout << "Polycone::DistToIn() (spot 1): point=" << point << ", dir=" << direction << ", localPoint=" << p
              << ", localDir=" << v << "\n";
#endif

    // TODO: add bounding box check maybe??

    distance      = kInfLength;
    int increment = (v.z() > 0) ? 1 : -1;
    if (std::fabs(v.z()) < kTolerance) increment = 0;
    int index = polycone.GetSectionIndex(p.z());
    if (index == -1) index = 0;
    if (index == -2) index = polycone.GetNSections() - 1;

    do {
      // now we have to find a section
      GenericPolyconeSection const &sec = polycone.GetSection(index);

#ifdef POLYCONEDEBUG
      std::cout << "Polycone::DistToIn() (spot 2):"
                << " index=" << index << " NSec=" << polycone.GetNSections() << " &sec=" << &sec << " - secPars:"
                << " secOffset=" << sec.fShift << " Dz=" << sec.fSolid->GetDz() << " Rmin1=" << sec.fSolid->GetRmin1()
                << " Rmin2=" << sec.fSolid->GetRmin2() << " Rmax1=" << sec.fSolid->GetRmax1()
                << " Rmax2=" << sec.fSolid->GetRmax2() << " -- calling Cone::DistToIn()...\n";
#endif

      CoaxialConesImplementation::template DistanceToIn<Real_v>(
          *sec.fCoaxialCones, p - Vector3D<Precision>(0, 0, sec.fShift), v, stepMax, distance);

#ifdef POLYCONEDEBUG
      std::cout << "Polycone::DistToIn() (spot 3):"
                << " distToIn() = " << distance << "\n";
#endif

      if (distance < kInfLength || !increment) break;
      index += increment;
    } while (index >= 0 && index < polycone.GetNSections());
    return;
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void DistanceToOut(UnplacedStruct_t const &polycone, Vector3D<Real_v> const &point,
                            Vector3D<Real_v> const &dir, Real_v const &stepMax, Real_v &distance)
  {
    distance            = kInfLength;
    Vector3D<Real_v> pn = point;
    Precision dist      = 0.;
    int increment       = (dir.z() > 0) ? 1 : -1;

    // specialization for N==1??? It should be a cone in the first place
    if (polycone.GetNSections() == 1) {
      const GenericPolyconeSection &section = polycone.GetSection(0);

      CoaxialConesImplementation::template DistanceToOut<Real_v>(
          *section.fCoaxialCones, point - Vector3D<Precision>(0, 0, section.fShift), dir, stepMax, distance);

      return;
    }

    int indexLow  = polycone.GetSectionIndex(point.z() - kTolerance);
    int indexHigh = polycone.GetSectionIndex(point.z() + kTolerance);
    // std::cout <<"IndexLow : " << indexLow <<" : IndexHigh : " << indexHigh << std::endl;
    int index = 0;
    if (indexLow < 0 && indexHigh < 0) {
      distance = -1;
      return;
    }

    if (indexLow < 0 && indexHigh >= 0 && dir.z() < 0.) {
      index                                 = indexHigh;
      const GenericPolyconeSection &section = polycone.GetSection(index);
      CoaxialConesImplementation::template DistanceToOut<Real_v>(*section.fCoaxialCones, pn, dir, stepMax, dist);
      distance = dist;
      return;
    }

    if (indexLow > 0 && indexHigh < 0 && dir.z() > 0.) {
      index                                 = indexLow;
      const GenericPolyconeSection &section = polycone.GetSection(index);
      CoaxialConesImplementation::template DistanceToOut<Real_v>(*section.fCoaxialCones, pn, dir, stepMax, dist);
      distance = dist;
      return;
    }

    Inside_t inside;
    Precision totalDistance = 0.;
    int count               = 0;
    do {

      if (indexLow >= 0 && indexHigh >= 0) {
        count++;
        index                                 = indexLow;
        const GenericPolyconeSection &section = polycone.GetSection(index);
        pn.z() -= section.fShift;
        // CoaxialConesImplementation::template Inside<Real_v>(*section.fCoaxialCones, point - Vector3D<Precision>(0, 0,
        // section.fShift), inside);
        CoaxialConesImplementation::template Inside<Real_v>(*section.fCoaxialCones, pn, inside);
        if (inside == EInside::kOutside) {
          if (count == 1) {
            distance = -1;
            return;
          } else {
            distance = totalDistance;
            return;
          }
        } else {
          // pn.z() -= section.fShift;
          CoaxialConesImplementation::template DistanceToOut<Real_v>(*section.fCoaxialCones, pn, dir, stepMax, dist);
          if (dist < 0.) break;
          totalDistance += dist;
          pn += dir * dist;
          pn.z() += section.fShift;
        }
        indexLow += increment;
        indexHigh += increment;
      }
      //} while(indexLow > -1 && indexHigh < polycone.GetNSections());//end of do-while
    } while (indexLow > -1 && indexLow < polycone.GetNSections()); // end of do-while
    distance = totalDistance;
    return;
  }

#if (0)
  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void DistanceToOut(UnplacedStruct_t const &polycone, Vector3D<Real_v> const &point,
                            Vector3D<Real_v> const &dir, Real_v const &stepMax, Real_v &distance)
  {
    // using namespace PolyconeTypes;
    distance            = kInfLength;
    Vector3D<Real_v> pn = point;

    // specialization for N==1??? It should be a cone in the first place
    if (polycone.GetNSections() == 1) {
      const GenericPolyconeSection &section = polycone.GetSection(0);

      CoaxialConesImplementation::template DistanceToOut<Real_v>(
          *section.fCoaxialCones, point - Vector3D<Precision>(0, 0, section.fShift), dir, stepMax, distance);

      return;
    }

    int indexLow  = polycone.GetSectionIndex(point.z() - kTolerance);
    int indexHigh = polycone.GetSectionIndex(point.z() + kTolerance);
    std::cout << "IndexLow : " << indexLow << " : " << indexHigh << std::endl;
    int index = 0;

    // section index is -1 when out of left-end
    // section index is -2 when beyond right-end

    if (indexLow < 0 && indexHigh < 0) {
      distance = -1;
      return;
    } else if (indexLow < 0 && indexHigh >= 0) {
      index                                 = indexHigh;
      const GenericPolyconeSection &section = polycone.GetSection(index);

      Inside_t inside;
      //      ConeImplementation<ConeTypes::UniversalCone>::Inside<Real_v>(
      //          *section.fSolid, point - Vector3D<Precision>(0, 0, section.fShift), inside);
      CoaxialConesImplementation::template Inside<Real_v>(*section.fCoaxialCones,
                                                          point - Vector3D<Precision>(0, 0, section.fShift), inside);
      if (inside == EInside::kOutside) {
        distance = -1;
        return;
      }
    } else if (indexLow != indexHigh && (indexLow >= 0)) {
      // we are close to an intermediate Surface, section has to be identified
      const GenericPolyconeSection &section = polycone.GetSection(indexLow);

      Inside_t inside;
      //      ConeImplementation<ConeTypes::UniversalCone>::Inside<Real_v>(
      //        *section.fSolid, point - Vector3D<Precision>(0, 0, section.fShift), inside);

      CoaxialConesImplementation::template Inside<Real_v>(*section.fCoaxialCones,
                                                          point - Vector3D<Precision>(0, 0, section.fShift), inside);

      if (inside == EInside::kOutside) {
        index = indexHigh;
      } else {
        index = indexLow;
      }
    } else {
      index = indexLow;
      if (index < 0) index = polycone.GetSectionIndex(point.z());
    }
    if (index < 0) {
      distance = 0.;
      return;
    }
    // Added
    else {
      const GenericPolyconeSection &section = polycone.GetSection(index);

      Inside_t inside;
      //      ConeImplementation<ConeTypes::UniversalCone>::Inside<Real_v>(
      //          *section.fSolid, point - Vector3D<Precision>(0, 0, section.fShift), inside);
      CoaxialConesImplementation::template Inside<Real_v>(*section.fCoaxialCones,
                                                          point - Vector3D<Precision>(0, 0, section.fShift), inside);
      if (inside == EInside::kOutside) {
        distance = -1;
        return;
      }
    }

    Precision totalDistance = 0.;
    Precision dist          = 0.;
    int increment           = (dir.z() > 0) ? 1 : -1;
    if (std::fabs(dir.z()) < kTolerance) increment = 0;

    // int counter = 0 ;
    Precision localCheckZ = point.z();
    std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
    // What is the relevance of istep?
    int istep             = 0;
    bool anotherBreakCond = true;
    do {

      const GenericPolyconeSection &section = polycone.GetSection(index);

      if ((totalDistance != 0) || (istep < 2)) {
        pn = point + totalDistance * dir; // point must be shifted, so it could eventually get into another solid
        pn.z() -= section.fShift;
        Inside_t inside;
        //        ConeImplementation<ConeTypes::UniversalCone>::Inside<Real_v>(*section.fSolid, pn, inside);
        CoaxialConesImplementation::template Inside<Real_v>(*section.fCoaxialCones, pn, inside);

        if (inside == EInside::kOutside) {
          break;
        }
      } else
        pn.z() -= section.fShift;

      istep++;

      // ConeImplementation<ConeTypes::UniversalCone>::DistanceToOut<Real_v>(*section.fSolid, pn, dir, stepMax, dist);
      CoaxialConesImplementation::template DistanceToOut<Real_v>(*section.fCoaxialCones, pn, dir, stepMax, dist);
      if (dist == -1) return;

      // Section Surface case
      if (std::fabs(dist) < 0.5 * kTolerance) {
        int index1 = index;
        if ((index > 0) && (index < polycone.GetNSections() - 1)) {
          index1 += increment;
        } else {
          if ((index == 0) && (increment > 0)) index1 += increment;
          if ((index == polycone.GetNSections() - 1) && (increment < 0)) index1 += increment;
        }

        Vector3D<Precision> pte                = point + (totalDistance + dist) * dir;
        const GenericPolyconeSection &section1 = polycone.GetSection(index1);
        pte.z() -= section1.fShift;
        Vector3D<Precision> localp;
        Inside_t inside22;
        // ConeImplementation<ConeTypes::UniversalCone>::Inside<Real_v>(*section1.fSolid, pte, inside22);
        CoaxialConesImplementation::template Inside<Real_v>(*section1.fCoaxialCones, pte, inside22);
        if (inside22 == 3 || (increment == 0)) {
          break;
        }
      } // end if surface case

      // counter++;
      std::cout << "localCheckZ before moving : " << localCheckZ << " : dist : " << dist << std::endl;
      localCheckZ += dist * dir.z();
      std::cout << "localCheckZ : " << localCheckZ << " : fZs(index+1) : " << polycone.fZs[index + 1] << std::endl;

      totalDistance += dist;
      std::cout << "TotalDist : " << totalDistance << std::endl;
      index += increment;

      // anotherBreakCond = (index < polycone.GetNSections()-1) && (localCheckZ == polycone.fZs[index+1]);
      // anotherBreakCond = (localCheckZ == polycone.fZs[index+1]);

    } while (increment != 0 && index >= 0 && index < polycone.GetNSections() &&
             anotherBreakCond /* && localCheckZ == polycone.fZs[indexHigh+1] && indexLow==indexHigh*/);
    std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;

    distance = totalDistance;

    return;
  }
#endif

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void SafetyToIn(UnplacedStruct_t const &polycone, Vector3D<Real_v> const &point, Real_v &safety)
  {

    Vector3D<Real_v> p = point;
    int index          = polycone.GetSectionIndex(p.z());

    bool needZ = false;
    if (index < 0) {
      needZ = true;
      if (index == -1) index = 0;
      if (index == -2) index = polycone.GetNSections() - 1;
    }
    Precision minSafety               = 0; //= SafetyFromOutsideSection(index, p);
    GenericPolyconeSection const &sec = polycone.GetSection(index);
    // safety to current segment
    if (needZ) {
      CoaxialConesImplementation::template SafetyToIn<Real_v>(*sec.fCoaxialCones,
                                                              p - Vector3D<Precision>(0, 0, sec.fShift), safety);
    } else

      CoaxialConesImplementation::template SafetyToIn<Real_v>(*sec.fCoaxialCones,
                                                              p - Vector3D<Precision>(0, 0, sec.fShift), safety);

    if (safety < kTolerance) return;
    minSafety       = safety;
    Precision zbase = polycone.fZs[index + 1];
    // going right
    for (int i = index + 1; i < polycone.GetNSections(); ++i) {
      Precision dz = polycone.fZs[i] - zbase;
      if (dz >= minSafety) break;

      GenericPolyconeSection const &sect = polycone.GetSection(i);
      CoaxialConesImplementation::template SafetyToIn<Real_v>(*sect.fCoaxialCones,
                                                              p - Vector3D<Precision>(0, 0, sect.fShift), safety);
      if (safety < minSafety) minSafety = safety;
    }

    // going left if this is possible
    if (index > 0) {
      zbase = polycone.fZs[index - 1];
      for (int i = index - 1; i >= 0; --i) {
        Precision dz = zbase - polycone.fZs[i];
        if (dz >= minSafety) break;
        GenericPolyconeSection const &sect = polycone.GetSection(i);

        CoaxialConesImplementation::template SafetyToIn<Real_v>(*sect.fCoaxialCones,
                                                                p - Vector3D<Precision>(0, 0, sect.fShift), safety);

        if (safety < minSafety) minSafety = safety;
      }
    }
    safety = minSafety;

    return;
  }

  template <typename Real_v>
  VECGEOM_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  static void SafetyToOut(UnplacedStruct_t const &polycone, Vector3D<Real_v> const &point, Real_v &safety)
  {
    typedef typename vecCore::Mask_v<Real_v> Bool_v;
    Bool_v compIn(false), compOut(false);
    GenericKernelForContainsAndInside<Real_v, Bool_v, true>(polycone, point, compIn, compOut);
    if (compOut) {
      safety = -1;
      return;
    }

    int index = polycone.GetSectionIndex(point.z());
    if (index < 0) {
      safety = -1;
      return;
    }

    GenericPolyconeSection const &sec = polycone.GetSection(index);

    Vector3D<Real_v> p = point - Vector3D<Precision>(0, 0, sec.fShift);
    // ConeImplementation<ConeTypes::UniversalCone>::SafetyToOut<Real_v>(*sec.fSolid, p, safety);
    CoaxialConesImplementation::template SafetyToOut<Real_v>(*sec.fCoaxialCones, p, safety);

    Precision minSafety = safety;
    if (minSafety == kInfLength) {
      safety = 0.;
      return;
    }
    if (minSafety < kTolerance) {
      safety = 0.;
      return;
    }

    Precision zbase = polycone.fZs[index + 1];
    for (int i = index + 1; i < polycone.GetNSections(); ++i) {
      Precision dz = polycone.fZs[i] - zbase;
      if (dz >= minSafety) break;
      GenericPolyconeSection const &sect = polycone.GetSection(i);
      p                                  = point - Vector3D<Precision>(0, 0, sect.fShift);

      // ConeImplementation<ConeTypes::UniversalCone>::SafetyToIn<Real_v>(*sect.fSolid, p, safety);
      CoaxialConesImplementation::template SafetyToIn<Real_v>(*sect.fCoaxialCones, p, safety);

      if (safety < minSafety) minSafety = safety;
    }

    if (index > 0) {
      zbase = polycone.fZs[index - 1];
      for (int i = index - 1; i >= 0; --i) {
        Precision dz = zbase - polycone.fZs[i];
        if (dz >= minSafety) break;
        GenericPolyconeSection const &sect = polycone.GetSection(i);
        p                                  = point - Vector3D<Precision>(0, 0, sect.fShift);

        // ConeImplementation<ConeTypes::UniversalCone>::SafetyToIn<Real_v>(*sect.fSolid, p, safety);
        CoaxialConesImplementation::template SafetyToIn<Real_v>(*sect.fCoaxialCones, p, safety);

        if (safety < minSafety) minSafety = safety;
      }
    }

    safety = minSafety;
    return;
  }
  /*
    template <typename Real_v>
    VECGEOM_FORCE_INLINE
    VECCORE_ATT_HOST_DEVICE
    static Vector3D<Real_v> NormalKernel(UnplacedStruct_t const &genericPolycone, Vector3D<Real_v> const &point,
                                         typename vecCore::Mask_v<Real_v> &valid)
    {
      // Computes the normal on a surface and returns it as a unit vector
      Vector3D<Real_v> normal(0., 0., 0.);
      return normal;
    }
  */
};
} // namespace VECGEOM_IMPL_NAMESPACE
} // namespace vecgeom

#endif // VECGEOM_VOLUMES_KERNEL_GENERICPOLYCONEIMPLEMENTATION_H_
