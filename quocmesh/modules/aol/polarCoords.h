#ifndef __POLARCOORDS_H
#define __POLARCOORDS_H

#include <aol.h>
#include <quoc.h>

// ******************************************************************
// polarCoords.h still only contains a class for calculating with
// 3d polar-coords, i.e. spherical coordinates.
// ******************************************************************

// neither Nemitz nor Schwen are convinced that this class works correctly ...

namespace aol {

//! 3d-vector in spherical coordinates
template <typename RealType> class SphericalVec {
protected:
  RealType _radius;           // length of the vector
  RealType _theta;            // angle in x,y-plane, in [0,2Pi)
  RealType _phi;              // angle in z - direction, between (-Pi/2, Pi/2)


public:

  // constructors, the first fills with 0's, the second gets the spherical coordinates
  SphericalVec() {
    _theta = _phi = _radius = 0.;
  }

  SphericalVec ( RealType theta, RealType phi, RealType radius ) {
    _radius = radius;
    _theta = theta;
    _phi = phi;
  }

  //! Copy-constructor
  SphericalVec ( const SphericalVec<RealType> &rhs ) {
    _radius = rhs.getRadius();
    _phi    = rhs.getPhi();
    _theta  = rhs.getTheta();
  }

  //! operator=
  SphericalVec<RealType>& operator= ( const SphericalVec<RealType> &rhs ) {
    _radius = rhs.getRadius();
    _phi    = rhs.getPhi();
    _theta  = rhs.getTheta();
    return *this;
  }



  // --------------------------- methods ----------------------------------------


  // get cartesic coordinates and calculate the spherical ones
  void initByCartesianCoords ( const aol::Vec3<RealType> arg ) {
    initByCartesianCoords ( arg[0], arg[1], arg[2] );
  }

  // get cartesic coordinates and calculate the spherical ones
  void initByCartesianCoords ( const RealType x, const RealType y, const RealType z ) {
    _radius = sqrt ( x * x + y * y + z * z );
    if ( _radius != 0 ) {
      _phi      = asin ( z / _radius );
      if ( sqrt ( x * x + y * y ) != 0.0 )
        _theta  = asin ( y / sqrt ( x * x + y * y ) );
      else
        _theta  = 0.0;

      // angle theta has to be in [0,2pi) => maybe add some offset
      if ( x <= 0 && y > 0 )  _theta += 0.5 * aol::NumberTrait<RealType>::pi;
      if ( x <= 0 && y <= 0 ) _theta += 1.5 * aol::NumberTrait<RealType>::pi;
      if ( x > 0 && y <= 0 )  _theta += 2.0 * aol::NumberTrait<RealType>::pi;
    } else {
      _phi = _theta = 0.;
    }
  }

  aol::Vec3<RealType> getCartesianCoords() const {
    aol::Vec3<RealType> result;
    result[0] = _radius * cos ( _phi ) * cos ( _theta );
    result[1] = _radius * cos ( _phi ) * sin ( _theta );
    result[2] = _radius * sin ( _phi );
    return result;
  }


  const RealType getRadius() const {
    return _radius;
  }
  const RealType getTheta()  const {
    return _theta;
  }
  const RealType getPhi()    const {
    return _phi;
  }

  void setRadius ( RealType radius ) {
    _radius = radius;
  }
  void setTheta ( RealType theta ) {
    _theta  = theta;
  }
  void setPhi ( RealType phi ) {
    _phi    = phi;
  }

  //! Equality operator.
  bool operator== ( const SphericalVec &Co ) const {
    if ( _radius != Co.getRadius() || _theta != Co.getTheta() || _phi != Co.getPhi() )
      return false;
    else
      return true;
  }

  //! Inequality operator.
  bool operator!= ( const SphericalVec &Co ) const {
    return !operator== ( Co );
  }

  //! normalize the vector
  void normalize() {
    _radius = 1.;
  }

  //! Clear the Vec.
  void setZero() {
    _radius = _theta = _phi = 0.;
  }


  void print () const {
    cerr << "(" << _radius << "," << _theta << "," << _phi << ")";
  }

};

// derivative computation methods formerly in namespace
// aol::polarTransform and only implemented for double moved to
// deprecated/

} // end namespace aol

#endif



