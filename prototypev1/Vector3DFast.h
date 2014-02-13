/*
 * Vector3DFast.h
 *
 *  Created on: Feb 7, 2014
 *      Author: swenzel
 */

#ifndef VECTOR3DFAST_H_
#define VECTOR3DFAST_H_

#include <Vc/Vc>
#include "Vc/vector.h"
#include "Vc/Memory"
#include <iostream>

class Vector3DFast
{
private:
	typedef Vc::Vector<double> base_t;
	typedef typename base_t::Mask mask_t;
	Vc::Memory<base_t, 3 > internalVcmemory;

	inline
	void SetVector(int i, Vc::Vector<double > const & v)
	{
		internalVcmemory.vector(i)=v;
	}

public:
	inline
	Vector3DFast( ) : internalVcmemory() {

		}

	inline
	Vector3DFast( Vector3DFast const & rhs ) : internalVcmemory() {
		for( int i=0; i < 1 + 3/Vc::Vector<double>::Size; i++ )
			{
				this->internalVcmemory.vector(i)=rhs.internalVcmemory.vector(i);
			}
	}


	inline
	Vector3DFast( double x,  double y,  double z) : internalVcmemory() {
			SetX(x);SetY(y);SetZ(z);
		}


	inline
	__attribute__((always_inline))
	Vector3DFast & operator+=( Vector3DFast const & rhs )
	{
		this->internalVcmemory+=rhs.internalVcmemory;
		return *this;
	}


	inline
	__attribute__((always_inline))
	Vector3DFast & operator*=( double s )
	{
		this->internalVcmemory*=s;
		return *this;
	}


	inline
		__attribute__((always_inline))
		Vector3DFast & operator-=( Vector3DFast const & rhs )
		{
			this->internalVcmemory-=rhs.internalVcmemory;
			return *this;
		}



	// assignment operator
	inline
	__attribute__((always_inline))
	Vector3DFast & operator=( Vector3DFast const & rhs )
	{
		for( int i=0; i < 1 + 3/Vc::Vector<double>::Size; i++ )
		{
			this->internalVcmemory.vector(i)=rhs.internalVcmemory.vector(i);
			return *this;
		}
	}


	inline
	__attribute__((always_inline))
	Vector3DFast & operator*=( Vector3DFast const & rhs )
	{
		/*
			for( int i=0; i < 1 + 3/Vc::Vector<double>::Size; i++ )
			{
				//		std::cerr << "adding vector " << i <<  std::endl;
				this->internalVcmemory.vector(i) *= rhs.internalVcmemory.vector(i);
			}
				*/
			this->internalVcmemory*=rhs.internalVcmemory;
			return *this;
	}


	inline
	__attribute__((always_inline))
	Vector3DFast & operator/=( Vector3DFast const & rhs )
	{
			this->internalVcmemory/= rhs.internalVcmemory;
			return *this;
	}

	inline
	__attribute__((always_inline))
	double ScalarProduct( Vector3DFast const & rhs ) const
	{
		double sum=0.;
		Vc::Vector<double> s(Vc::Zero);
		for( int i=0; i < 1 + 3/Vc::Vector<double>::Size; i++ )
		{
			//	std::cerr << "adding vector " << i <<  std::endl;
			Vc::Vector<double> tmp1 = this->internalVcmemory.vector(i);
			Vc::Vector<double> tmp2 = rhs.internalVcmemory.vector(i);
			s+=tmp1*tmp2;
		}
		return s.sum();
	}


	inline double GetX() const {
		return internalVcmemory[ 0 ];
	}

	inline double GetY() const {
		return internalVcmemory[ 1 ];
	}

	inline double GetZ() const {
		return internalVcmemory[ 2 ];

		//return internalVcmemory.vector( 2 / Vc::Vector<double>::Size  )[ 2 % Vc::Vector<double>::Size ];
	}

	inline void SetX(double x)  {
		internalVcmemory[0]=x;}

	inline void SetY(double y)  {
		internalVcmemory[1]=y;}

	inline void SetZ(double z)  {
		internalVcmemory[2]=z;}


	inline
	__attribute__((always_inline))
	Vector3DFast Abs() const
	{
		Vector3DFast tmp;
		for( int i=0; i < 1 + 3/Vc::Vector<double>::Size; i++ )
		{
			base_t v = this->internalVcmemory.vector(i);
			tmp.internalVcmemory.vector(i) = Vc::abs( v );
		}
		return tmp;
	}

	inline
	__attribute__((always_inline))
	// a starting comparison operator
	bool IsAnySmallerZero( ) const
	{
		bool result = false;
		Vector3DFast tmp;
		for( int i=0; i < 1 + 3/Vc::Vector<double>::Size; i++ )
		{
			mask_t m = this->internalVcmemory.vector(i) < 0;
			result |= ! m.isEmpty();
		}
		return result;
	}

	inline
	__attribute__((always_inline))
	// a starting comparison operator
	bool IsAnySmallerThan( Vector3DFast const & rhs ) const
	{
		bool result = false;
		for( int i=0; i < 1 + 3/Vc::Vector<double>::Size; i++ )
		{
			base_t v1,v2;
			v1 = this->internalVcmemory.vector(i);
			v2 = rhs.internalVcmemory.vector(i);
			mask_t m =  v1 < v2;
			result |= ! m.isEmpty();
// could early return here ...

		}
		return result;
	}

	inline
	__attribute__((always_inline))
	// a starting comparison operator
	bool IsAnyLargerThan( Vector3DFast const & rhs ) const
	{
			bool result = false;
			for( int i=0; i < 1 + 3/Vc::Vector<double>::Size; i++ )
			{
				base_t v1,v2;
				v1 = this->internalVcmemory.vector(i);
				v2 = rhs.internalVcmemory.vector(i);
				mask_t m =  v1 > v2;
				result |= ! m.isEmpty();
// could early return here
			}
			return result;
	}


};


inline
Vector3DFast const operator-(Vector3DFast const & lhs, Vector3DFast const & rhs)
{
	Vector3DFast tmp(lhs);
	tmp-=rhs;
	return tmp;
}


inline
Vector3DFast const operator+(Vector3DFast const & lhs, Vector3DFast const & rhs)
{
	Vector3DFast tmp(lhs);
	tmp+=rhs;
	return tmp;
}


inline
Vector3DFast const operator*(double lhs, Vector3DFast const & rhs)
{
	Vector3DFast tmp(rhs);
	tmp*=lhs;
	return tmp;
}

inline
Vector3DFast const operator*(Vector3DFast const & lhs, double rhs)
{
	Vector3DFast tmp(lhs);
	tmp*=rhs;
	return tmp;
}


inline
static
std::ostream& operator<<( std::ostream& stream, Vector3DFast const & vec )
{
	stream << "{ " << vec.GetX() << " , " << vec.GetY() << " , " << vec.GetZ() << " } ";
	return stream;
}

// to try something with Vc memory
class FastTransformationMatrix
{
// storage
private:
    // the translation
	Vector3DFast trans;

	// the columns of the rotation matrix ( as vectors )
	Vector3DFast rotcol1;
	Vector3DFast rotcol2;
	Vector3DFast rotcol3;

	// the rows of the rotation matrix ( as vectors )
	Vector3DFast rotrow1;
	Vector3DFast rotrow2;
	Vector3DFast rotrow3;

public:
	FastTransformationMatrix() : trans(), rotcol1(), rotcol2(), rotcol3(),
		rotrow1(), rotrow2(), rotrow3()
	{};

		//template <TranslationIdType tid, RotationIdType rid>
	inline
	void
	MasterToLocal(Vector3DFast const & master, Vector3DFast & local) const
	{
		Vector3DFast tmp = master-trans;
		local = tmp.GetX()*rotrow1 + tmp.GetY()*rotrow2 + tmp.GetZ()*rotrow3;
	}

	inline
	void
	LocalToMaster(Vector3DFast const & local, Vector3DFast & master) const
	{
		Vector3DFast tmp = local+trans;
		master = tmp.GetX()*rotcol1 + tmp.GetY()*rotcol2 + tmp.GetZ()*rotcol3;
	}

	void setAngles(double phi, double theta, double psi)
	{
		double degrad = M_PI/180.;
		double sinphi = sin(degrad*phi);
		double cosphi = cos(degrad*phi);
		double sinthe = sin(degrad*theta);
		double costhe = cos(degrad*theta);
		double sinpsi = sin(degrad*psi);
		double cospsi = cos(degrad*psi);

		// this is setting the rowvectors
		rotrow1.SetX( cospsi*cosphi - costhe*sinphi*sinpsi );
		rotrow1.SetY( -sinpsi*cosphi - costhe*sinphi*cospsi );
		rotrow1.SetZ( sinthe*sinphi );

		rotrow2.SetX( cospsi*sinphi + costhe*cosphi*sinpsi );
		rotrow2.SetY( -sinpsi*sinphi + costhe*cosphi*cospsi );
		rotrow2.SetZ( -sinthe*cosphi );

		rotrow3.SetX( sinpsi*sinthe );
		rotrow3.SetY( cospsi*sinthe );
		rotrow3.SetZ( costhe );

		rotcol1.SetX( rotrow1.GetX() ); rotcol1.SetY( rotrow2.GetX() ); rotcol1.SetZ( rotrow3.GetX() );
		rotcol2.SetX( rotrow1.GetY() ); rotcol2.SetY( rotrow2.GetY() ); rotcol2.SetZ( rotrow3.GetY() );
		rotcol3.SetX( rotrow1.GetZ() ); rotcol3.SetY( rotrow2.GetZ() ); rotcol3.SetZ( rotrow3.GetZ() );
	}


	FastTransformationMatrix(double tx, double ty, double tz, double phi, double theta, double psi)
	{
		trans.SetX(tx);
		trans.SetY(ty);
		trans.SetZ(tz);

		setAngles(phi, theta, psi);
	}

	void print() const
	{
		// print the matrix in 4x4 format
		printf("%12.11f\t%12.11f\t%12.11f    Tx = %10.6f\n", rotcol1.GetX(), rotcol2.GetX(), rotcol3.GetX(), trans.GetX());
		printf("%12.11f\t%12.11f\t%12.11f    Ty = %10.6f\n", rotcol1.GetY(), rotcol2.GetY(), rotcol3.GetY(), trans.GetY());
		printf("%12.11f\t%12.11f\t%12.11f    Tz = %10.6f\n", rotcol1.GetZ(), rotcol2.GetZ(), rotcol3.GetZ(), trans.GetZ());
	}
};

#endif /* VECTOR3DFAST_H_ */
