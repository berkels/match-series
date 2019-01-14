#ifndef __MULTIDOBJECT_H
#define __MULTIDOBJECT_H

#include <quoc.h>
#include <gridBase.h>

namespace qc {

// forward declaration
class GridStructure;

//! basis class for template specialization
template < qc::Dimension Dim >
class MultiDObject;


//! basis class for three-dimensional objects that have a size in x, y and z direction
template<>
class MultiDObject<qc::QC_3D> {
protected:
  int _numX, _numY, _numZ;

public:
  MultiDObject ( ) : _numX ( 0 ), _numY ( 0 ), _numZ ( 0 ) {
  }

  MultiDObject ( const int NumX, const int NumY, const int NumZ ) : _numX ( NumX ), _numY ( NumY ), _numZ ( NumZ ) {
  }

  template< class Struct3d >
  explicit MultiDObject ( const Struct3d &Strc ) : _numX ( Strc.getNumX() ), _numY ( Strc.getNumY() ), _numZ ( Strc.getNumZ() ) {
  }

  // copy constructor, assignment operator should do the correct thing

  virtual ~MultiDObject ( ) {
  }

  //! size comparison
  bool operator== ( const MultiDObject<qc::QC_3D> & other ) {
    return ( ( _numX == other._numX ) && ( _numY == other._numY ) && ( _numZ == other._numZ ) );
  }

  bool operator!= ( const MultiDObject<qc::QC_3D> & other ) {
    return ( ! ( operator== ( other ) ) );
  }

  inline int getNumX() const {
    return( _numX );
  }

  inline int getNumY() const {
    return( _numY );
  }

  inline int getNumZ() const {
    return( _numZ );
  }

  aol::Vec3<int> getSize() const {
    return ( aol::Vec3<int> ( _numX, _numY, _numZ ) );
  }

  inline int size() const {
    return ( _numX * _numY * _numZ );
  }

  //! return size in all directions if equal, throw exception otherwise
  int getNumXYZ ( ) const {
    if ( ( this->getNumX() != this->getNumY() ) || ( this->getNumY() != this->getNumZ() ) )
      throw aol::Exception ( "qc::MultiDObject<qc::QC_3D>::getNumXYZ() may not be called for non-cubic arrays", __FILE__, __LINE__ );
    return ( this->getNumX() );
  }

protected:
  //! compute index for 1D representation
  inline int oneDIndex ( const int ix, const int iy, const int iz ) const {
    return ( qc::ILexCombine3 ( ix, iy, iz, _numX, _numY ) );
  }

#ifdef BOUNDS_CHECK
  void boundsCheck ( const int ix, const int iy, const int iz, const char* methodName, const char* fileName, const int lineNumber ) const {
    if ( ! ( ( ix < _numX ) && ( iy < _numY ) && ( iz < _numZ ) && ( ix >= 0 ) && ( iy >= 0 ) && ( iz >= 0 ) ) ) {
      char errmsg[1024];
      sprintf( errmsg, "%s: %d %d %d are out of bounds (%d %d %d)!", methodName, ix, iy, iz, _numX, _numY, _numZ );
      throw aol::OutOfBoundsException( errmsg, fileName, lineNumber );
    }
  }
#endif

};


//! basis class for template specialization
template < qc::Dimension Dim >
class MultiDStorageObject;

//! basis class for three-dimensional objects that have a size in x, y and z direction and that permit changing the size via reallocate and/or resize
template<>
class MultiDStorageObject<qc::QC_3D> : public MultiDObject < qc::QC_3D > {
public:
  MultiDStorageObject ( const int NumX, const int NumY, const int NumZ ) : MultiDObject<qc::QC_3D> ( NumX, NumY, NumZ ) {
  }

  template< class Struct3d >
  explicit MultiDStorageObject ( const Struct3d &Strc ) : MultiDObject<qc::QC_3D> ( Strc ) {
  }

  virtual void reallocate ( const int, const int, const int ) {
    throw aol::UnimplementedCodeException ( "qc::MultiDStorageObject<qc::QC_3D>::reallocate ( const int, const int, const int ) must be overloaded on derived classes.", __FILE__, __LINE__ );
  }

  virtual void resize ( const int, const int, const int ) {
    throw aol::UnimplementedCodeException ( "qc::MultiDStorageObject<qc::QC_3D>::resize ( const int, const int, const int ) must be overloaded on derived classes.", __FILE__, __LINE__ );
  }

  void reallocate ( const GridStructure &Grid );

  template< typename Structure >
  void reallocate ( const Structure &Struc ) {
    this->reallocate ( Struc.getNumX(), Struc.getNumY(), Struc.getNumZ() );
  }

  void resize ( const GridStructure &Grid );

  template< typename Structure >
  void resize ( const Structure &Struc ) {
    resize ( Struc.getNumX(), Struc.getNumY(), Struc.getNumZ() );
  }

protected:
  void changeSizeTo ( const int Nx, const int Ny, const int Nz ) {
    _numX = Nx;
    _numY = Ny;
    _numZ = Nz;
  }
};

}

#endif
