/*
 * TransformationMatrix.h
 *
 *  Created on: Nov 8, 2013
 *      Author: swenzel
 */

#ifndef TRANSFORMATIONMATRIX_H_
#define TRANSFORMATIONMATRIX_H_

#include "Vector3D.h"
#include <cmath>

typedef int TranslationIdType;
typedef int RotationIdType;

class TransformationMatrix
{
private:
    // this can be varied depending on template specialization
	double trans[3];
	double rot[9];

	template<RotationIdType rid=-1>
	inline
	void
	emitrotationcode(Vector3D const &, Vector3D &) const;


public:
	static
	inline
	RotationIdType
	getfootprint(double const *r){
		int footprint=0;
		// count zero entries and give back footprint which classifies them
		for(int i=0;i<9;i++)
		{
			if( r[i]==0. )
				footprint+=i*i*i; // cubic power identifies cases uniquely
		}
		return footprint;
	}


	// constructor
	TransformationMatrix(double const *t, double const *r)
	{
	    trans[0]=t[0];
	    trans[1]=t[1];
		trans[2]=t[2];
		for(int i=0;i<9;i++)
			rot[i]=r[i];
		// we need to check more stuff ( for instance that product along diagonal is +1)
	}

	inline
	static TranslationIdType GetTranslationIdType(double const *t)
	{
		if( t[0]==0 && t[1]==0 && t[2]==0 )
			return 0;
		else
			return 1;
	}



	inline
	static RotationIdType GetRotationIdType(double const *r)
	{
		// this function has to be generated by some script
		if( r[0]==1 && r[1]==0 && r[2]==0 && r[3]==0 && r[4]==1 && r[5]==0 && r[6]==0 && r[7]==0 && r[8]==1 )
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}


	template <TranslationIdType tid=-1, RotationIdType rid=-1>
	inline
	void
	MasterToLocal(Vector3D const &, Vector3D &) const;

	// to transform real vectors, we don't care about translation
	template <RotationIdType rid=-1>
	inline
	void
	MasterToLocalVec(Vector3D const &, Vector3D &) const;

	// why not ??
	// inline Vector3D MasterToLocal(Vector3D const &) const;

	// for a Vc vector: not so easy
	// inline MasterToLocal(Vc::double) const;


	friend class PhysicalVolume;
	friend class PhysicalBox;
};


template <TranslationIdType tid=-1, RotationIdType rid=-1>
void
TransformationMatrix::MasterToLocal(Vector3D const & master, Vector3D & local) const
{
	if( tid==0 && rid ==0 ) // this means identity
	{
		local.x = master.x;
		local.y = master.y;
		local.z = master.z;
	}
	else if( tid != 0 && rid == 0 ) // tid == 1 means we have
	{
		local.x = master.x + trans[0];
		local.y = master.y + trans[1];
		local.z = master.z + trans[2];
	}
	else if( tid == 0 && rid!=0 ) // pure rotation
	{
		emitrotationcode<rid>(master,local);
	}
	else( tid != 0 && rid!=0 ) // both rotation and translation
	{
		Vector3D mt;
		mt.x = master.x + trans[0];
		mt.y=  master.y + trans[1];
		mt.z = master.z + trans[2];
		emitrotationcode<rid>(mt, local);
	}
}

template <RotationIdType rid=-1>
void
TransformationMatrix::MasterToLocalVec(Vector3D const & master, Vector3D & local ) const
{
	MasterToLocal<0, rid>(master, local);
}

template <RotationIdType rid=-1>
inline
void
TransformationMatrix::emitrotationcode(Vector3D const & mt, Vector3D & local) const
{
	if(rid==252){
	     local.x = mt.x*rot[0];
	     local.y = mt.y*rot[4]+mt.y*rot[7];
	     local.z = mt.z*rot[5]+mt.z*rot[8];
	     return;
	}
	if(rid==405){
	     local.x = mt.x*rot[3];
	     local.y = mt.y*rot[1]+mt.y*rot[7];
	     local.z = mt.z*rot[2]+mt.z*rot[8];
	     return;
	}
	if(rid==882){
	     local.x = mt.x*rot[6];
	     local.y = mt.y*rot[1]+mt.y*rot[4];
	     local.z = mt.z*rot[2]+mt.z*rot[5];
	     return;
	}
	if(rid==415){
	     local.x = mt.x*rot[3]+mt.x*rot[6];
	     local.y = mt.y*rot[1];
	     local.z = mt.z*rot[5]+mt.z*rot[8];
	     return;
	}
	if(rid==496){
	     local.x = mt.x*rot[0]+mt.x*rot[6];
	     local.y = mt.y*rot[4];
	     local.z = mt.z*rot[2]+mt.z*rot[8];
	     return;
	}
	if(rid==793){
	     local.x = mt.x*rot[0]+mt.x*rot[3];
	     local.y = mt.y*rot[7];
	     local.z = mt.z*rot[2]+mt.z*rot[5];
	     return;
	}
	if(rid==638){
	     local.x = mt.x*rot[3]+mt.x*rot[6];
	     local.y = mt.y*rot[4]+mt.y*rot[7];
	     local.z = mt.z*rot[2];
	     return;
	}
	if(rid==611){
	     local.x = mt.x*rot[0]+mt.x*rot[6];
	     local.y = mt.y*rot[1]+mt.y*rot[7];
	     local.z = mt.z*rot[5];
	     return;
	}
	if(rid==692){
	     local.x = mt.x*rot[0]+mt.x*rot[3];
	     local.y = mt.y*rot[1]+mt.y*rot[4];
	     local.z = mt.z*rot[8];
	     return;
	}
	if(rid==720){
	     local.x = mt.x*rot[0];
	     local.y = mt.y*rot[4];
	     local.z = mt.z*rot[8];
	     return;
	}
	if(rid==828){
	     local.x = mt.x*rot[0];
	     local.y = mt.y*rot[7];
	     local.z = mt.z*rot[5];
	     return;
	}
	if(rid==756){
	     local.x = mt.x*rot[3];
	     local.y = mt.y*rot[1];
	     local.z = mt.z*rot[8];
	     return;
	}
	if(rid==918){
	     local.x = mt.x*rot[3];
	     local.y = mt.y*rot[7];
	     local.z = mt.z*rot[2];
	     return;
	}
	if(rid==954){
	     local.x = mt.x*rot[6];
	     local.y = mt.y*rot[1];
	     local.z = mt.z*rot[5];
	     return;
	}
	if(rid==1008){
	     local.x = mt.x*rot[6];
	     local.y = mt.y*rot[4];
	     local.z = mt.z*rot[2];
	     return;
	}
	local.x = mt.x*rot[0]+mt.x*rot[3]+mt.x*rot[6];
	local.y = mt.y*rot[1]+mt.y*rot[4]+mt.y*rot[7];
	local.z = mt.z*rot[2]+mt.z*rot[5]+mt.z*rot[8];
}


#endif /* TRANSFORMATIONMATRIX_H_ */
