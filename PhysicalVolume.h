/*
 * PhysicalVolume.h
 *
 *  Created on: Nov 8, 2013
 *      Author: swenzel
 */

#ifndef PHYSICALVOLUME_H_
#define PHYSICALVOLUME_H_

#include <list>
#include "TransformationMatrix.h"


// pure abstract class
class PhysicalVolume
{
	protected:

		TransformationMatrix *matrix; // placement matrix with respect to containing volume
		std::list<PhysicalVolume> daughterVolumes; // list or vector?

		// something like a logical volume id

		// I am not sure that this is appropriate
		// maybe we should do PhysicalShape and PhysicalVolume as different classes
		LogicalVolume const *logicalvol;


	public:

		virtual double DistanceToIn(Vector3D const &, Vector3D const &) const = 0;
		virtual double DistanceToOut(Vector3D const &, Vector3D const &) const = 0;

		// delete all explicit constructors etc...
		
		
		//
		virtual ~PhysicalVolume();

		// add factory methods
		void AddDaughter( PhysicalVolume const * vol );
		LogicalVolume const * getLogicalVolume(){return logicalvol;}
};


#endif /* PHYSICALVOLUME_H_ */
