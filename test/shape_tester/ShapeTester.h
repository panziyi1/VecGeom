//
// Definition of the batch solid test
//

#ifndef ShapeTester_hh
#define ShapeTester_hh

#include <iostream>
#include <sstream>

#include "VUSolid.hh"
#include "UUtils.hh"

const double kApproxEqualTolerance = 1E-6;

typedef struct sShapeTesterErrorList {
    std::string sMessage;
    int sNUsed;
    struct sShapeTesterErrorList *sNext;
  } ShapeTesterErrorList;

class ShapeTester {

//ALL MEMBER FUNCTIONS FIRST
public:
  ShapeTester();
  ~ShapeTester();

  int Run(VUSolid *testVolume);
  int RunMethod(VUSolid *testVolume, std::string fMethod1);
  inline void SetFilename(const std::string &newFilename) { fFilename = newFilename; }
  inline void SetMaxPoints(const int newMaxPoints) { fMaxPoints = newMaxPoints; }
  inline void SetMethod(const std::string &newMethod) { fMethod = newMethod; }
  inline void SetInsidePercent(const double percent) { fInsidePercent = percent; }
  inline void SetOutsidePercent(const double percent) { fOutsidePercent = percent; }
  inline void SetEdgePercent(const double percent) { fEdgePercent = percent; }
  inline void SetOutsideMaxRadiusMultiple(const double percent) { fOutsideMaxRadiusMultiple = percent; }
  inline void SetOutsideRandomDirectionPercent(const double percent) { fOutsideRandomDirectionPercent = percent; }
  inline void SetNewSaveValue(const double tolerance) { fMinDifference = tolerance; }
  inline void SetSaveAllData(const bool safe) { fIfSaveAllData = safe; }
  inline void SetRunAllTests(const bool safe) { fIfMoreTests = safe; }
  void SetFolder(const std::string &newFolder);
  void SetVerbose(int verbose) { fVerbose = verbose; }
  inline int GetMaxPoints() const { return fMaxPoints; }
  inline UVector3 GetPoint(int index) { return fPoints[index]; }
  inline void SetNumberOfScans(int num) { fGNumberOfScans = num; }

  /* Keep this Function as public to allow, if somebody just want
   * to do the Convention Check
   */
  bool RunConventionChecker(VUSolid *testVolume);


private:
  void SetDefaults();
  int SaveVectorToExternalFile(const std::vector<double> &vector, const std::string &fFilename);
  int SaveVectorToExternalFile(const std::vector<UVector3> &vector, const std::string &fFilename);
  int SaveLegend(const std::string &fFilename);
  int SaveDifLegend(const std::string &fFilename);
  int SaveDoubleResults(const std::string &fFilename);
  int SaveDifDoubleResults(const std::string &fFilename);
  int SaveVectorResults(const std::string &fFilename);
  int SaveDifVectorResults(const std::string &fFilename);
  int SaveDifVectorResults1(const std::string &fFilename);
  std::string PrintCoordinates(const UVector3 &vec, const std::string &delimiter, int precision = 4);
  std::string PrintCoordinates(const UVector3 &vec, const char *delimiter, int precision = 4);
  void PrintCoordinates(std::stringstream &ss, const UVector3 &vec, const std::string &delimiter, int precision = 4);
  void PrintCoordinates(std::stringstream &ss, const UVector3 &vec, const char *delimiter, int precision = 4);

  template <class T>
  void VectorDifference(const std::vector<T> &first, const std::vector<T> &second, std::vector<T> &result);

  void VectorToDouble(const std::vector<UVector3> &vectorUVector, std::vector<double> &vectorDouble);
  void BoolToDouble(const std::vector<bool> &vectorBool, std::vector<double> &vectorDouble);
  int CountDoubleDifferences(const std::vector<double> &differences);
  int CountDoubleDifferences(const std::vector<double> &differences, const std::vector<double> &values1,
                             const std::vector<double> &values2);

  void FlushSS(std::stringstream &ss);
  void Flush(const std::string &s);

  UVector3 GetPointOnOrb(double r);
  UVector3 GetRandomDirection();

  int TestConsistencySolids();
  int TestInsidePoint();
  int TestOutsidePoint();
  int TestSurfacePoint();

  int TestNormalSolids();

  int TestSafetyFromInsideSolids();
  int TestSafetyFromOutsideSolids();
  int ShapeSafetyFromInside(int max);
  int ShapeSafetyFromOutside(int max);

  void PropagatedNormal(const UVector3 &point, const UVector3 &direction, double distance, UVector3 &normal);
  void PropagatedNormalU(const UVector3 &point, const UVector3 &direction, double distance, UVector3 &normal);
  int TestDistanceToInSolids();
  int TestAccuracyDistanceToIn(double dist);
  int ShapeDistances();
  int TestFarAwayPoint();
  int TestDistanceToOutSolids();
  int ShapeNormal();
  int TestXRayProfile();
  int XRayProfile(double theta = 45, int nphi = 15, int ngrid = 1000, bool useeps = true);
  int Integration(double theta = 45, double phi = 45, int ngrid = 1000, bool useeps = true, int npercell = 1,
                  bool graphics = true);
  double CrossedLength(const UVector3 &point, const UVector3 &dir, bool useeps);
  void CreatePointsAndDirections();
  void CreatePointsAndDirectionsSurface();
  void CreatePointsAndDirectionsEdge();
  void CreatePointsAndDirectionsInside();
  void CreatePointsAndDirectionsOutside();

  void CompareAndSaveResults(const std::string &fMethod, double resG, double resR, double resU);
  int SaveResultsToFile(const std::string &fMethod);
  void SavePolyhedra(const std::string &fMethod);
  double MeasureTest(int (ShapeTester::*funcPtr)(int), const std::string &fMethod);
  double NormalizeToNanoseconds(double time);

  int TestMethod(int (ShapeTester::*funcPtr)());
  int TestMethodAll();

  inline double RandomRange(double min, double max) {
    double rand = min + (max - min) * UUtils::Random();
    return rand;
  }
  inline void GetVectorUSolids(UVector3 &point, const std::vector<UVector3> &afPoints, int index) {
    const UVector3 &p = afPoints[index];
    point.Set(p.x(), p.y(), p.z());
  }

  inline void SetVectorUSolids(const UVector3 &point, std::vector<UVector3> &afPoints, int index) {
    UVector3 &p = afPoints[index];
    p.Set(point.x(), point.y(), point.z());
  }

  inline double RandomIncrease() {
    double tolerance = VUSolid::Tolerance();
    double rand = -1 + 2 * UUtils::Random();
    double dif = tolerance * 0.1 * rand;
    return dif;
  }

  /* Private functions for Convention Checker, These functions never need
   * to be called from Outside the class
   */
  void PrintConventionMessages();
  void GenerateConventionReport();
  void SetupConventionMessages();
  bool ShapeConventionChecker();
  bool ShapeConventionSurfacePoint();
  bool ShapeConventionInsidePoint();
  bool ShapeConventionOutsidePoint();
  void SetNumDisp(int);
  bool ApproxEqual(const double x, const double y);
  //Return true if the 3vector check is approximately equal to target
  template <class Vec_t>
  bool ApproxEqual(const Vec_t &check, const Vec_t &target);


protected:
  UVector3 GetRandomPoint() const;
  double GaussianRandom(const double cutoff) const;
  void ReportError(int *nError, UVector3 &p, UVector3 &v, double distance,
                   std::string comment); //, std::ostream &fLogger );
  void ClearErrors();
  int CountErrors() const;


//ALL DATA MEMBERS
protected:
  int fMaxPoints; //Maximum num. of points to be generated for ShapeTester Tests.
  int fVerbose; // Variable to set verbose
  double fInsidePercent; // Percentage of inside points
  double fOutsidePercent; // Percentage of outside points
  double fEdgePercent; // Percentage of edge points
  double fOutsideMaxRadiusMultiple; // Range of outside points
  double fOutsideRandomDirectionPercent; // Percentage of outside random direction

  // XRay profile statistics
  int fGNumberOfScans;  // data member to store the number of different scan angle used for XRay profile
  double fGCapacitySampled; // data member to store calculated capacity
  double fGCapacityAnalytical; // data member to store analytical capacity
  double fGCapacityError; // data member to store error between above two.

  std::string fMethod; // data member to store the name of method to be executed

  ShapeTesterErrorList *fErrorList; // data member to store the list of errors

private:
  std::vector<UVector3> fPoints; //STL vector to store the points generated for various tests of ShapeTester
  std::vector<UVector3> fDirections; //STL vector to store the directions generated for corresponding points.

  VUSolid *fVolumeUSolids; // Data member to store pointer to the shape under test
  std::string fVolumeString; // data member to store the name of volume;

  std::vector<UVector3> fResultVectorUSolids; // stl vector for storing the vector results
  std::vector<UVector3> fResultVectorDifference; // stl vector for storing the vector difference
  std::vector<double> fResultDoubleUSolids; // stl vector for storing the double results
  std::vector<double> fResultDoubleDifference; // stl vector for storing the double difference
  std::vector<bool> fResultBoolUSolids;  // stl vector for storing the bool results.
  std::vector<bool> fResultBoolDifference; // stl vector for storing the bool difference.

  int fOffsetSurface; // offset of surface points
  int fOffsetInside; // offset of inside points
  int fOffsetOutside; // offset of outside points
  int fOffsetEdge; // offset of edge points
  int fMaxPointsInside; // Number of inside points
  int fMaxPointsOutside; // Number of outside points
  int fMaxPointsSurface; // Number of surface points
  int fMaxPointsEdge; // Number of edge points

  std::ostream *fLog; //Pointer to the directory storing all the log file for different tests

  std::string fFolder; // Name of the log folder
  std::string fFilename; // name of the file name depending on the method under test

  // Save only differences
  bool fIfSaveAllData; // save alldata, big files
  bool fIfMoreTests; // do all additional tests,
  // take more time, but not affect performance measures
  bool fIfDifUSolids; // save differences of Geant4 with Usolids or with ROOT
  double fMinDifference; // save data, when difference is bigger that min
  bool fDefinedNormal; // bool variable to skip normal calculation if it does not exist in the shape
  bool fIfException; // data memeber to abort ShapeTester if any error found

  //Added data member required for convention checker
  std::vector<std::string> fConventionMessage; //STL vector for convention error messages.
  int fScore; // an error code generate if conventions not followed, 0 mean convenetion followed.
  int fNumDisp; // number of points to be displayed in case a shape is not following conventions.

};

#endif
