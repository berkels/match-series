#ifndef __ARRAYEXTENSIONS_H
#define __ARRAYEXTENSIONS_H

#include <array.h>
#include <indexMapper.h>
#include <multiDObject.h>
#include <vectorExtensions.h>

namespace qc {

/** Basis class for template specialization */
template< typename DataType, qc::Dimension dim >
class AArray {};


/**
 * \author Berkels
 */
template< typename DataType >
class AArray < DataType, qc::QC_1D > : public aol::RandomAccessContainer < DataType > {
public:
  template< class Struct1d >
  explicit AArray ( const Struct1d &Strc ) : aol::RandomAccessContainer < DataType > ( Strc.getNumX() ) {
    cerr << "AArray < DataType, qc::QC_1D > has not been tested yet! Please don't remove this warning before testing it properly!\n";
  }

  //! Return reference to entry at position Pos
  DataType& getRef ( const CoordType &Pos ) {
    return (*this)[ Pos[0] ];
  }

  aol::Vec<1, int> getSize ( ) const { return aol::Vec<1, int> ( this->size() ); }
};

/** Abstract 2D array class that can store anything that can be put in a std::vector.
 *  AArray < DataType, qc::QC_2D > encapsulates the index arithmetics and provides getRef returning references to the instances stored here.
 *  It is intended for storing "big" objects, so there is not get/set which would require copying or assigning the contained objects.
 *  \author Berkels, Schwen
 */
template< typename DataType >
class AArray < DataType, qc::QC_2D > {
protected:
  aol::Vec2<int>           _size;
  aol::RandomAccessContainer < DataType > _data;

public:
  //! Standard constructor
  AArray ( ) : _data ( 0 ) {
    //  cerr << "AArray < DataType, qc::QC_2D > has not been tested yet! Please don't remove this warning before testing it properly!\n";
  }

  AArray ( const int NumX, const int NumY ) : _size ( NumX, NumY ), _data ( _size[0] * _size[1] ) {
    cerr << "AArray < DataType, qc::QC_2D > has not been tested yet! Please don't remove this warning before testing it properly!\n";
  }

  explicit AArray ( const aol::Vec<2, int> &Size ) : _size( Size ), _data ( Size.prod() ) {}

  template< class Struct2d >
  explicit AArray ( const Struct2d &Strc ) : _size ( Strc.getNumX(), Strc.getNumY() ), _data ( _size[0] * _size[1] ) {
    cerr << "AArray < DataType, qc::QC_2D > has not been tested yet! Please don't remove this warning before testing it properly!\n";
  }

  //! Copy constructor
  explicit AArray ( const AArray<DataType, qc::QC_2D> &Other ) : _size ( Other._size ), _data ( Other._data ) {
    cerr << "AArray < DataType, qc::QC_2D > has not been tested yet! Please don't remove this warning before testing it properly!\n";
  }

  ~AArray ( ) {}

  //! Assignment operator
  AArray<DataType, qc::QC_2D>& operator= ( const AArray<DataType, qc::QC_2D> &Other ) {
    _size = Other._size;
    _data = Other._data;

    return ( *this );
  }

  bool operator== ( const AArray<DataType, qc::QC_2D> &Other ) {
#ifdef BOUNDS_CHECK
    if ( _size != Other._size )
      throw aol::Exception ( "qc::AArray<DataType, qc::QC_2D>::operator==: should not compare AArrays of different size", __FILE__, __LINE__ );
#endif
    return ( this->_data == Other._data );
  }

  bool operator!= ( const AArray<DataType, qc::QC_2D> &Other ) {
    return ! ( operator== ( Other ) );
  }

public:
  //! Change size and clear data. Note that a resize preserving data is nontrivial.
  void reallocate ( const int NumX, const int  NumY ) {
    _size[0] = NumX;
    _size[1] = NumY;
    resetEntries();
  }

  template< typename Structure >
  void reallocate ( const Structure &struc ) {
    reallocate ( struc.getNumX(), struc.getNumY() );
  }

  //! reallocate to old size, i.e. deletes old contents and creates new contained objects by their standard constructor
  void resetEntries ( ) {
    _data.reallocate( _size[0] * _size[1] );
  }

  //! Return reference to entry at position ( I, J )
  DataType& getRef ( const int I, const int J ) {
#ifdef BOUNDS_CHECK
    boundsCheck ( I, J, "qc::AArray::get Index out of bounds", __FILE__, __LINE__ );
#endif
    return ( _data[ qc::ILexCombine2 ( I, J, _size[0] ) ] );
  }

  //! Return const reference to entry at position ( I, J )
  const DataType& getRef ( const int I, const int J ) const {
#ifdef BOUNDS_CHECK
    boundsCheck ( I, J, "qc::AArray::get Index out of bounds", __FILE__, __LINE__ );
#endif
    return ( _data[ qc::ILexCombine2 ( I, J, _size[0] ) ] );
  }

  //! Return reference to entry at position pos
  DataType& getRef ( const aol::Vec<2, int> &Pos ) {
    return getRef ( Pos[0], Pos[1] );
  }

  //! Return reference to entry at position Pos
  DataType& getRef ( const CoordType &Pos ) {
    return getRef ( Pos[0], Pos[1] );
  }

  //! Return const reference to entry at position Pos
  const DataType& getRef ( const aol::Vec<2, int> &Pos ) const {
    return getRef ( Pos[0], Pos[1] );
  }

  //! Return const reference to entry at position Pos
  const DataType& getRef ( const CoordType &Pos ) const {
    return getRef ( Pos[0], Pos[1] );
  }

  //! Access operator like for a vector
  DataType& operator[] ( const int I ) {
#ifdef BOUNDS_CHECK
    boundsCheck ( I, "qc::AArray::operator[] Index out of bounds", __FILE__, __LINE__ );
#endif
    return ( _data[I] );
  }

  //! Const access operator like for a vector
  const DataType& operator[] ( const int I ) const {
#ifdef BOUNDS_CHECK
    boundsCheck ( I, "qc::AArray::operator[] Index out of bounds", __FILE__, __LINE__ );
#endif
    return ( _data[I] );
  }

  int size ( ) const { return static_cast<int> ( _data.size() ); }
  int getNumX ( ) const { return ( _size[0] ); }
  int getNumY ( ) const { return ( _size[1] ); }
  aol::Vec2<int> getSize ( ) const { return _size; }

  //! assignment necessary for each entry, this may be slow for big DataType
  void setAll ( const DataType &value ) {
    for ( unsigned int I = 0; I < _data.size(); ++I ) {
      _data[I] = value;
    }
  }

#ifdef BOUNDS_CHECK
protected:
  inline void boundsCheck ( const int I, const int J, const char* Msg, const char* Fi, const int Li ) const {
    const bool isIn = ( I >= 0 && I < _size[0] && J >= 0 && J < _size[1] );
    if ( !isIn ) {
      char errmsg[1024];
      sprintf ( errmsg, "%s %d %d", Msg, I, J );
      throw aol::OutOfBoundsException ( errmsg, Fi, Li );
    }
  }

  inline void boundsCheck ( const int I, const char* Msg, const char* Fi, const int Li ) const {
    const bool isIn = I >= 0 && I < static_cast<int> ( _data.size() );
    if ( !isIn ) {
      char errmsg[1024];
      sprintf ( errmsg, "%s %d", Msg, I );
      throw aol::OutOfBoundsException ( errmsg, Fi, Li );
    }
  }
#endif
};


/** Abstract 3D array class that can store anything that can be put in an aol::RandomAccessContainer.
 *  AArray < DataType, qc::QC_3D > encapsulates the index arithmetics and provides getRef returning references to the instances stored here.
 *  It is intended for storing "big" objects, so there is not get/set which would require copying or assigning the contained objects.
 *  \author Schwen
 */
template< typename DataType >
class AArray < DataType, qc::QC_3D > : public qc::MultiDStorageObject < qc::QC_3D > {
protected:
  aol::RandomAccessContainer < DataType > _data;

public:
  //! Standard constructor
  AArray ( ) : MultiDStorageObject<qc::QC_3D> ( 0, 0, 0 ), _data ( 0 ) {}

  AArray ( const int numX, const int numY, const int numZ ) : MultiDStorageObject<qc::QC_3D> ( numX, numY, numZ ), _data ( numX * numY * numZ ) {}

  explicit AArray ( const aol::Vec3<int> &Size ) : MultiDStorageObject<qc::QC_3D> ( Size[0], Size[1], Size[2] ), _data ( Size.prod() ) {}

  template< class Struct3d >
  explicit AArray ( const Struct3d &strc ) : MultiDStorageObject<qc::QC_3D> ( strc ), _data ( strc.getNumX() * strc.getNumY() * strc.getNumZ() ) {}

  //! Copy constructor
  explicit AArray ( const AArray<DataType, qc::QC_3D> &other ) : MultiDStorageObject<qc::QC_3D> ( other ), _data ( other._data ) {}

  ~AArray ( ) {}

  //! Assignment operator
  AArray<DataType, qc::QC_3D>& operator= ( const AArray<DataType, qc::QC_3D> &other ) {
    MultiDObject<qc::QC_3D>::operator= ( other );
    _data = other._data;

    return ( *this );
  }

  bool operator== ( const AArray<DataType, qc::QC_3D> &other ) {
#ifdef BOUNDS_CHECK
    if ( MultiDObject<qc::QC_3D>::operator!= ( other ) )
      throw aol::Exception ( "qc::AArray<DataType, qc::QC_3D>::operator==: should not compare AArrays of different size", __FILE__, __LINE__ );
#endif
    return ( this->_data == other._data );
  }

  bool operator!= ( const AArray<DataType, qc::QC_3D> &other ) {
    return ! ( operator== ( other ) );
  }

public:
  //! Change size and clear data. Note that a resize preserving data is nontrivial.
  void reallocate ( const int NumX, const int  NumY, const int NumZ ) {
    this->changeSizeTo ( NumX, NumY, NumZ );
    resetEntries();
  }

  using qc::MultiDStorageObject<qc::QC_3D>::reallocate;

  //! reallocate to old size, i.e. deletes old contents and creates new contained objects by their standard constructor
  void resetEntries ( ) {
    _data.reallocate ( this->size() );
  }

  //! Return reference to entry at position ( i, j, k )
  DataType& getRef ( const int i, const int j, const int k ) {
#ifdef BOUNDS_CHECK
    this->boundsCheck ( i, j, k, "qc::AArray<qc::QC_3D>::get Index out of bounds", __FILE__, __LINE__ );
#endif
    return ( _data[ this->oneDIndex ( i, j, k ) ] );
  }

  //! Return const reference to entry at position ( i, j, k )
  const DataType& getRef ( const int i, const int j, const int k ) const {
#ifdef BOUNDS_CHECK
    this->boundsCheck ( i, j, k, "qc::AArray<qc::QC_3D>::get Index out of bounds", __FILE__, __LINE__ );
#endif
    return ( _data[ this->oneDIndex ( i, j, k ) ] );
  }

  //! Return reference to entry at position pos
  DataType& getRef ( const aol::Vec3<int> &pos ) {
#ifdef BOUNDS_CHECK
    this->boundsCheck ( pos[0], pos[1], pos[2], "qc::AArray<qc::QC_3D>::get Index out of bounds", __FILE__, __LINE__ );
#endif
    return ( _data[ this->oneDIndex ( pos[0], pos[1], pos[2] ) ] );
  }

  //! Return reference to entry at position pos
  DataType& getRef ( const CoordType &pos ) {
#ifdef BOUNDS_CHECK
    this->boundsCheck ( pos[0], pos[1], pos[2], "qc::AArray<qc::QC_3D>::get Index out of bounds", __FILE__, __LINE__ );
#endif
    return ( _data[ this->oneDIndex ( pos[0], pos[1], pos[2] ) ] );
  }

  //! Return const reference to entry at position pos
  const DataType& getRef ( const aol::Vec3<int> &pos ) const {
#ifdef BOUNDS_CHECK
    this->boundsCheck ( pos[0], pos[1], pos[2], "qc::AArray<qc::QC_3D>::get Index out of bounds", __FILE__, __LINE__ );
#endif
    return ( _data[ this->oneDIndex ( pos[0], pos[1], pos[2] ) ] );
  }

  //! Return const reference to entry at position pos
  const DataType& getRef ( const CoordType &pos ) const {
#ifdef BOUNDS_CHECK
    this->boundsCheck ( pos[0], pos[1], pos[2], "qc::AArray<qc::QC_3D>::get Index out of bounds", __FILE__, __LINE__ );
#endif
    return ( _data[ this->oneDIndex ( pos[0], pos[1], pos[2] ) ] );
  }

  //! Access operator like for a vector
  DataType& operator[] ( const int i ) {
#ifdef BOUNDS_CHECK
    this->boundsCheck ( i, "qc::AArray<qc::QC_3D>::operator[] Index out of bounds", __FILE__, __LINE__ );
#endif
    return ( _data[i] );
  }

  //! Access operator like for a vector
  const DataType& operator[] ( const int i ) const {
#ifdef BOUNDS_CHECK
    this->boundsCheck ( i, "qc::AArray<qc::QC_3D>::operator[] Index out of bounds", __FILE__, __LINE__ );
#endif
    return ( _data[i] );
  }

  // do not implement methods
  // DataType get ( const int i, const int j, const int k ) const;
  // DataType get ( const aol::Vec3<int> &pos ) const;

  //! assignment necessary for each entry, this may be slow for big DataType
  void setAll ( const DataType &value ) {
    for ( int i = 0; i < _data.size(); ++i ) {
      _data[i] = value;
    }
  }

  /** Reallocates current AAray and pastes Ohter into it.
   *  If Other is smaller in any space direction, it is centered, if Other is bigger in any space direction, only its central values are copied.
   */
  void padFrom ( const AArray<DataType, qc::QC_3D> &Other ) {
    qc::CoordType
      offsetb ( ( this->getSize() - Other.getSize() ) / 2 ),
      offsett ( ( this->getSize() - Other.getSize() ) - ( this->getSize() - Other.getSize() ) / 2 ) ; // integer division and conversion on purpose

    this->reallocate ( this->getNumX(), this->getNumY(), this->getNumZ() );

    for ( int z = aol::Max ( aol::NumberTrait<short>::zero, offsetb[2] ); z < aol::Min ( this->getNumX(), this->getNumX() - offsett[0] ); ++z ) {
      for ( int y = aol::Max ( aol::NumberTrait<short>::zero, offsetb[1] ); y < aol::Min ( this->getNumY(), this->getNumY() - offsett[1] ); ++y ) {
        for ( int x = aol::Max ( aol::NumberTrait<short>::zero, offsetb[0] ); x < aol::Min ( this->getNumZ(), this->getNumZ() - offsett[2] ); ++x ) {
          this->getRef ( x, y, z ) = Other.getRef ( x - offsetb[0], y - offsetb[1], z - offsetb[2] );
        }
      }
    }
  }

#ifdef BOUNDS_CHECK
protected:
  inline void boundsCheck ( const int i, const char* msg, const char* fi, const int li ) const {
    const bool isIn = i >= 0 && i < static_cast<int> ( _data.size() );
    if ( !isIn ) {
      char errmsg[1024];
      sprintf ( errmsg, "%s %d", msg, i );
      throw aol::OutOfBoundsException ( errmsg, fi, li );
    }
  }

  using MultiDObject<qc::QC_3D>::boundsCheck;
#endif
};


/** Basis class for template specialization */
template< typename DataType, typename ArrayType, qc::Dimension dim >
class RectangularContainer {};


// 2D version not implemented yet.


/** A RectangularContainer is a container class for using brick-shaped subsets of nD data structures (e. g. ScalarArray, BitArray).
 *  RectangularContainer basically encapsulates the index arithmetics (index shifts) and prevents waste of memory for the unused out-of-brick entries.
 *  \author Schwen
 */

template< typename DataType, typename ArrayType >
class RectangularContainer < DataType, ArrayType, qc::QC_3D > {
protected:
  aol::Vec3<int> _lower, _upper;
  ArrayType _data;

public:
  //! Standard Constructor creating empty brick
  RectangularContainer () : _lower ( 0, 0, 0 ), _upper ( 0, 0, 0 ), _data ( _upper - _lower ) {}

  //! Copy constructor relying on correct implementation of copy constructor of ArrayType
  RectangularContainer ( const RectangularContainer<DataType, ArrayType, qc::QC_3D> &other, aol::CopyFlag copyFlag = aol::DEEP_COPY ) : _lower ( other._lower ), _upper ( other._upper ), _data ( other._data, copyFlag ) {}

  //! Assignment Operator. Sizes must match.
  RectangularContainer<DataType, ArrayType, qc::QC_3D>& operator= ( const RectangularContainer<DataType, ArrayType, qc::QC_3D> &other ) {
    if ( _lower == other._lower && _upper == other._upper )
      _data = other._data;
    else
      throw aol::Exception ( "RectangularContainer<DataType, ArrayType, qc::QC_3D>::operator=: incompatible RectangularContainers", __FILE__, __LINE__ );
    return ( *this );
  }

  //! This is probably the most useful generic constructor
  RectangularContainer ( const aol::Vec3<int> &min, const aol::Vec3<int> &max ) :
    _lower ( min ), _upper ( max ), _data ( _upper[0] - _lower[0], _upper[1] - _lower[1], _upper[2] - _lower[2] ) {}

  //! Extract given brick from given array.
  RectangularContainer ( const aol::Vec3<int> &min, const aol::Vec3<int> &max, const ArrayType &array ) :
      _lower ( min ), _upper ( max ), _data ( _upper[0] - _lower[0], _upper[1] - _lower[1], _upper[2] - _lower[2] ) {
    for ( int k = _lower[2]; k < _upper[2]; ++k )
      for ( int j = _lower[1]; j < _upper[1]; ++j )
        for ( int i = _lower[0]; i < _upper[0]; ++i )
          _data.set ( i - _lower[0], j - _lower[1], k - _lower[2], array.get ( i, j, k ) ); // Relies on bounds checking in the array.
  }

public:
  //! Resize RectangularContainer, preserving contents.
  // \todo think about how to speed this up
  void resize (  const aol::Vec3<int> &min, const aol::Vec3<int> &max, const DataType default_for_unset = aol::NumberTrait<DataType>::zero ) {
    if ( ! ( min == _lower && max == _upper ) ) { // if nontrivial change
      const aol::Vec3<int> old_lower = _lower, old_upper = _upper;
      _lower = min; _upper = max;
      ArrayType old_data ( _data, aol::DEEP_COPY );
      _data.reallocate ( _upper[0] - _lower[0] , _upper[1] - _lower[1], _upper[2] - _lower[2] );
      _data.setAll ( default_for_unset ); // and afterwards overwrite known values

      for ( int i = aol::Max ( _lower[0], old_lower[0] ); i < aol::Min ( _upper[0], old_upper[0] ); ++i ) {
        for ( int j = aol::Max ( _lower[1], old_lower[1] ); j < aol::Min ( _upper[1], old_upper[1] ); ++j ) {
          for ( int k = aol::Max ( _lower[2], old_lower[2] ); k < aol::Min ( _upper[2], old_upper[2] ); ++k ) {
            const DataType value = old_data.get ( i - old_lower[0], j - old_lower[1], k - old_lower[2] );
            _data.set ( i - _lower[0], j - _lower[1], k - _lower[2], value );
          }
        }
      }
    }
  }

  //! Resize to match size of other RectangularContainer
  void resize ( const RectangularContainer<DataType, ArrayType, qc::QC_3D> &other, const DataType default_for_unset = aol::NumberTrait<DataType>::zero ) {
    resize ( other._lower, other._upper, default_for_unset );
  }

  //! change size of RectangularContainer, deleting contents
  void reallocate (  const aol::Vec3<int> &min, const aol::Vec3<int> &max ) {
    _lower = min; _upper = max;
    _data.reallocate (  _upper[0] - _lower[0] , _upper[1] - _lower[1], _upper[2] - _lower[2] );
  }

  //! Resize to match size of other RectangularContainer
  void reallocate ( const RectangularContainer<DataType, ArrayType, qc::QC_3D> &other ) {
    reallocate ( other._lower, other._upper );
  }

  //! Copy data from corresponding positions from src
  void getBrickDataFrom ( const ArrayType &src ) {
    for ( int i = _lower[0]; i < _upper[0] ; ++i )
      for ( int j = _lower[1]; j < _upper[1] ; ++j )
        for ( int k = _lower[2]; k < _upper[2] ; ++k )
          _data.set ( i - _lower[0], j - _lower[1], k - _lower[2], src.get ( i, j, k ) );
  }

  //! Copy data from corresponding positions from src
  void writeBrickDataTo ( ArrayType &dst ) const {
    for ( int i = _lower[0]; i < _upper[0] ; ++i )
      for ( int j = _lower[1]; j < _upper[1] ; ++j )
        for ( int k = _lower[2]; k < _upper[2] ; ++k )
          dst.set ( i, j, k, _data.get ( i - _lower[0], j - _lower[1], k - _lower[2] ) );
  }



  DataType get ( const int i, const int j, const int k ) const {
#ifdef BOUNDS_CHECK
    boundsCheck ( i, j, k, "qc::RectangularContainer::get: Index out of bounds. ", __FILE__, __LINE__ );
#endif
    return ( _data.get ( i - _lower[0], j - _lower[1], k - _lower[2] ) );
  }

  DataType get ( const aol::Vec3<int> &pos ) const {
    return ( get ( pos[0], pos[1], pos[2] ) );
  }

  DataType get ( const CoordType &pos ) const {
    return ( get ( pos[0], pos[1], pos[2] ) );
  }

  const DataType& getRef ( const int i, const int j, const int k ) const {
#ifdef BOUNDS_CHECK
    boundsCheck ( i, j, k, "qc::RectangularContainer::get: Index out of bounds. ", __FILE__, __LINE__ );
#endif
    return ( _data.getRef ( i - _lower[0], j - _lower[1], k - _lower[2] ) );
  }

  const DataType& getRef ( const aol::Vec3<int> &pos ) const {
    return ( getRef ( pos[0], pos[1], pos[2] ) );
  }

  const DataType& getRef ( const CoordType &pos ) const {
    return ( getRef ( pos[0], pos[1], pos[2] ) );
  }


  DataType& getRef ( const int i, const int j, const int k ) {
#ifdef BOUNDS_CHECK
    boundsCheck ( i, j, k, "qc::RectangularContainer::get: Index out of bounds. ", __FILE__, __LINE__ );
#endif
    return ( _data.getRef ( i - _lower[0], j - _lower[1], k - _lower[2] ) );
  }

  DataType& getRef ( const aol::Vec3<int> &pos ) {
    return ( getRef ( pos[0], pos[1], pos[2] ) );
  }

  DataType& getRef ( const CoordType &pos ) {
    return ( getRef ( pos[0], pos[1], pos[2] ) );
  }


  DataType operator[] ( const int i ) const {
    return ( _data[i] );
  }

  void set ( const int i, const int j, const int k, const DataType value ) {
#ifdef BOUNDS_CHECK
    boundsCheck ( i, j, k, "qc::RectangularContainer::set: Index out of bounds. ", __FILE__, __LINE__ );
#endif
    _data.set ( i - _lower[0], j - _lower[1], k - _lower[2], value );
  }

  void set ( const aol::Vec3<int> &pos, const DataType value ) {
    set ( pos[0], pos[1], pos[2], value );
  }

  void set ( const CoordType &pos, const DataType value ) {
    set ( pos[0], pos[1], pos[2], value );
  }

  void setAll ( const DataType value ) {
    _data.setAll ( value );
  }

  void setZero ( ) {
    _data.setZero();
  }

  aol::Vec3<int> getLower ( ) const {
    return ( _lower );
  }

  aol::Vec3<int> getUpper ( ) const {
    return ( _upper );
  }

  void dump ( ) const {
    cerr << "[ " << _lower[0] << ", " << _upper[0] << " ] x [ " << _lower[1] << ", " << _upper[1] << " ] x [ " << _lower[2] << ", " << _upper[2] << " ]" << endl;
    _data.dump();
  }

  void printSlices ( ) const {
    cerr << "[ " << _lower[0] << ", " << _upper[0] << " ] x [ " << _lower[1] << ", " << _upper[1] << " ] x [ " << _lower[2] << ", " << _upper[2] << " ]" << endl;
    _data.printSlices();
  }

  RectangularContainer<DataType, ArrayType, qc::QC_3D>& operator+= ( const RectangularContainer<DataType, ArrayType, qc::QC_3D> &Arg ) {
    // This obviously does not make sense for bools ...
    _data += Arg._data;
    return ( *this );
  }

  RectangularContainer<DataType, ArrayType, qc::QC_3D>& operator-= ( const RectangularContainer<DataType, ArrayType, qc::QC_3D> &Arg ) {
    // This obviously does not make sense for bools ...
    _data -= Arg._data;
    return ( *this );
  }

  int numOccurence ( const DataType value ) const {
    return ( _data.numOccurence ( value ) );
  }

  DataType norm() const {
    // This obviously does not make sense for bools ...
    return ( static_cast<DataType> ( _data.norm() ) );
  }

  int size ( ) const {
    return ( _data.size() );
  }

  ArrayType& getContainedRef ( ) {
    return ( _data );
  }

  const ArrayType& getContainedRef ( ) const {
    return ( _data );
  }

  inline bool containsPoint ( const int i, const int j, const int k ) const {
    return ( ( i >= _lower[0] ) && ( i < _upper[0] ) && ( j >= _lower[1] ) && ( j < _upper[1] ) && ( k >= _lower[2] ) && ( k < _upper[2] ) );
  }

  inline bool containsPoint ( const aol::Vec3<int> &pos ) {
    return ( containsPoint ( pos[0], pos[1], pos[2] ) );
  }

  inline bool containsPoint ( const CoordType &pos ) {
    return ( containsPoint ( pos[0], pos[1], pos[2] ) );
  }

#ifdef BOUNDS_CHECK
protected:
  inline bool boundsCheck ( const int i, const int j, const int k, const char* msg, const char* fi, const int li ) const {
    const bool isIn = containsPoint ( i, j, k );
    if ( !isIn ) {
      char errmsg[1024];
      sprintf ( errmsg, "%s %d %d %d", msg, i, j, k );
      throw aol::OutOfBoundsException ( errmsg, fi, li );
    }
    return ( isIn );
  }
#endif
};



} // end namespace

#endif
