// classes for dictionary generation

#include "base/Transformation3D.h"
#include "base/Vector.h"
#include "management/GeoManager.h"

#include "navigation/VNavigator.h"
#include "navigation/VLevelLocator.h"
#include "navigation/VSafetyEstimator.h"
#include "navigation/SimpleLevelLocator.h"


#include "volumes/UnplacedBox.h"
#include "volumes/UnplacedParaboloid.h"
#include "volumes/UnplacedParallelepiped.h"
#include "volumes/UnplacedSphere.h"
#include "volumes/UnplacedTrd.h"
#include "volumes/UnplacedGenTrap.h"
#include "volumes/UnplacedTrapezoid.h"
#include "volumes/UnplacedHype.h"
#include "volumes/UnplacedOrb.h"
#include "volumes/UnplacedCutTube.h"

#include "volumes/SpecializedPlacedVolImplHelper.h"
#include "volumes/LogicalVolume.h"
#include "volumes/PlacedVolume.h"
#include "volumes/PlacedBox.h"
#include "volumes/kernel/BoxImplementation.h"
#include "volumes/Box.h"
#include "volumes/Paraboloid.h"