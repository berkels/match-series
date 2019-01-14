#ifndef __BITVECTOR_H
#define __BITVECTOR_H

#include <aol.h>

namespace aol {

/** Class for compressed storage of bools in a vector-type data structure.
 *
 *  The size need not be a multiple of 8 so there may be some excess bits at the end.
 *  These need to be kept zero.
 */
class BitVector {
private:
  int   _size;
  unsigned char *_pFieldData;

public:
  BitVector ( ) : _size ( 0 ), _pFieldData ( NULL ) {}

  explicit BitVector ( int Size ) {
    _size = Size;
    alloc();
    setZero();
  }

  explicit BitVector ( int Size, const bool Value ) {
    _size = Size;
    alloc ( );
    setAll ( Value );
  }

  explicit BitVector ( const BitVector &other ) : _size ( other._size ) {
    alloc();
    memcpy ( _pFieldData, other._pFieldData, convertBitVectorLengthToCharLength ( _size ) ); // sic: memcpy( dest, source, n_bytes)
  }
  
  explicit BitVector ( const char *filename ) : _size ( 0 ), _pFieldData ( NULL )  {
    loadFromFile( filename );
  }

  virtual ~BitVector() {
    dealloc();
  }

  inline bool operator[] ( const int Index ) const {
#ifdef BOUNDS_CHECK
    if ( Index >= _size )  {
      throw aol::Exception ( "aol:BitVector::operator[]: Index out of range", __FILE__, __LINE__ );
    }
#endif
    const int element = Index / ( sizeof ( unsigned char ) * 8 );
    return ( _pFieldData[element] >> Index % ( sizeof ( unsigned char ) * 8 ) ) & 1;
  }

  inline bool get ( const int Index ) const {
      return ( *this ) [Index];
    }

  inline void set ( const int Index, const bool Value ) {
#ifdef BOUNDS_CHECK
    if ( Index >= _size )  {
      throw aol::Exception ( "aol:BitVector::set: Index out of range", __FILE__, __LINE__ );
    }
#endif

    int element = Index / ( sizeof ( unsigned char ) * 8 );

    if ( Value == true ) {
      _pFieldData[element] |= 1 << Index % ( sizeof ( unsigned char ) * 8 );
    } else if ( Value == false ) {
      _pFieldData[element] &= ~ ( 1 << Index % ( sizeof ( unsigned char ) * 8 ) );
    }
  }

  //! Set all entries to zero
  void setZero() {
    memset ( _pFieldData, 0, convertBitVectorLengthToCharLength ( _size ) );
  }


  inline int size() const {
    return static_cast<int> ( _size );
  }

  //! sets all entries to true
  void setAll() {
    setAll ( true );
  }

  void setAll ( const bool val ) {
    int ll = _size / ( sizeof ( unsigned char ) * 8 );
    memset ( _pFieldData, ( val == true ? 0xFF : 0 ), ll );

    // remaining entries, possibly not full last byte may be set to true (else resizing won't work) correctly
    for ( int i = 8 * ll; i < _size ; ++i ) {
      set ( i, val );
    }
  }

  ostream& dump ( ostream & out = cout ) const;

  int numOccurence ( const bool val ) const {
    //! \todo implement this more efficiently
    int ret = 0;
    for ( int i = 0; i < _size; ++i ) {
      if ( this->get ( i ) == val )
        ++ret;
    }

    return ( ret );
  }

  int numTrue() const {
    return ( numOccurence ( true ) );
  }

  int numFalse() const {
    return ( numOccurence ( false ) );
  }

  //! Sizes must match (if necessary, resize before assigning)
  BitVector& operator= ( const BitVector &other ) {
    if ( this == &other ) return ( *this ); // beware of self-assignment

    if ( this->_size != other._size ) {
      throw aol::Exception ( "aol::BitVector::operator= dimensions do not match!", __FILE__, __LINE__ );
    }

    // Beware of self-assignment
    if ( _pFieldData != other._pFieldData ) {
      memcpy ( _pFieldData, other._pFieldData, convertBitVectorLengthToCharLength ( _size ) ); // sic: memcpy( dest, source, n_bytes)
    }

    return ( *this );
  }

  //! comparison for equality
  bool operator== ( const BitVector &other ) const {
    if ( this->_size != other._size ) {
      throw aol::Exception ( "aol::BitVector::operator== dimensions do not match!", __FILE__, __LINE__ );
    }

    // We can only check the fully used bytes with memcmp, the remaining ones have to be checked manually.
    // Or we need to require that the unused part of the last byte always has to have a certain value.
    const int numBytesToCheck = this->_size / ( sizeof ( unsigned char ) * 8 );
    if ( memcmp ( _pFieldData, other._pFieldData, numBytesToCheck ) != 0 )
      return ( false );

    for ( int i = numBytesToCheck * sizeof ( unsigned char ) * 8; i < this->_size; ++i ) {
      if ( this->get ( i ) != other.get( i ) ) {
        return ( false );
      }
    }

    return ( true );
  }

  //! comparison for inequality
  bool operator!= ( const BitVector &other ) const {
    return ( ! this->operator== ( other ) );
  }

  //! component-wise AND
  BitVector & operator &= ( const BitVector & other ) {
    if ( this == &other ) return ( *this ); // nothing to do

    if ( this->_size != other._size ) {
      throw aol::Exception ( "aol::BitVector::operator &= dimensions do not match!", __FILE__, __LINE__ );
    }

    for ( int i = 0; i < this->_size; ++i )
      this->set ( i, this->get(i) && other[i] );

    return ( *this );
  }

  //! component-wise OR
  BitVector & operator |= ( const BitVector & other ) {
    if ( this == &other ) return ( *this ); // nothing to do

    if ( this->_size != other._size ) {
      throw aol::Exception ( "aol::BitVector::operator |= dimensions do not match!", __FILE__, __LINE__ );
    }

    for ( int i = 0; i < this->_size; ++i )
      this->set ( i, this->get(i) || other[i] );

    return ( *this );
  }


  /** resize the BitVector and do not delete contents.
   */
  virtual void resize ( const int Length ) {
    const int
      lCurrent = convertBitVectorLengthToCharLength ( _size ),
      lNew = convertBitVectorLengthToCharLength ( Length );

    unsigned char *tempData = new unsigned char [ lNew ];
    memset ( tempData, 0, lNew ); // make sure all memory (after end of copied data) is set to 0

    if ( _pFieldData != NULL ) {
      memcpy ( tempData, _pFieldData, sizeof ( unsigned char ) * aol::Min ( lCurrent, lNew ) );
      dealloc ( );
    }

    _size = Length;
    _pFieldData = tempData;

    // set values outside range in last partial byte to zero
    const int bitsToZero = sizeof ( unsigned char ) * 8 * lNew - Length;

    unsigned char bitMask = 0xFF;
    for ( int i = 7; i >= 8 - bitsToZero; --i ) {
      bitMask ^= ( 1 << i );
    }
    _pFieldData[ lNew - 1 ] &= bitMask;
  }

  /** Change size of the BitVector and set it to zero.
   */
  virtual void reallocate ( const int Length ) {
    dealloc();
    _size = Length;
    alloc();
    this->setZero();
  }

  void invert ( ) {
    for ( int i = 0; i < size(); ++i ) {
      bool temp = get ( i );
      set ( i, !temp );
    }
  }

  template< typename RealType, typename VectorType >
  void setNonzeroPatternFrom ( const VectorType &vec ) {
#ifdef BOUNDS_CHECK
    if ( _size != vec.size() ) {
      cerr << _size << " " << vec.size() << endl;
      throw aol::Exception( "aol::BitVector::setNonzeroPatternFrom: incorrect size of vector", __FILE__, __LINE__ );
    }
#endif
    for ( int i = 0; i < _size; ++i ) {
      set ( i, vec[i] != aol::NumberTrait<RealType>::zero );
    }
  }

  template <typename VectorType>
  void thresholdFrom ( const VectorType &Vec, const typename VectorType::DataType ThresholdValue ) {
    for ( int i = 0; i < _size; ++i ) {
      set ( i, Vec[i] > ThresholdValue );
    }
  }

  /** Compute 32 bit checksum */
  unsigned int crc32 ( ) const {
    return ( aol::crc32 ( _pFieldData, convertBitVectorLengthToCharLength ( _size ) ) );
  }

  //! lossless save in binary quocmesh format
  void saveToFile ( const char *filename ) const;

  //! load from file in quocmesh format, no conversion of data types
  void loadFromFile ( const char *filename );

protected:
  inline int convertBitVectorLengthToCharLength ( const int BitVectorLength ) const {
    return ( ( BitVectorLength % ( sizeof ( unsigned char ) * 8 ) == 0 ) ? 0 : 1 ) + BitVectorLength / ( sizeof ( unsigned char ) * 8 );
  }

  //! alloc does not clear
  void alloc ( ) {
    _pFieldData = new unsigned char [ convertBitVectorLengthToCharLength ( _size ) ];
  }

  void dealloc ( ) {
    if ( _pFieldData != NULL ) {
      delete[] _pFieldData;
      _pFieldData = NULL;
    }
  }

  //! do not allow setting any other type than bool (implemented explicitly above)
  template< typename AnyType >
  void set ( const int,  const AnyType & );

  //! do not allow setting any other type than bool
  template< typename AnyType >
  void setAll ( const AnyType & );

}; // end class BitVector


//! Shift BitVector into stream.
ostream &operator<< ( ostream &os, const BitVector &BVec );


} // namespace end


#endif
