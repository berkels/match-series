#ifndef __ARRAY_H
#define __ARRAY_H

#include <bitVector.h>
#include <quoc.h>
#include <gridBase.h>
#include <vec.h>
#include <qmException.h>
#include <bitArray.h>

namespace qc {

template <typename DataType>
class Array : public aol::Vector<DataType> {
public:
  //! Standard constructor creating array of dimensions 0 x 0 x 1
  Array ( ) :
      aol::Vector<DataType> ( 0 * 0 * 1 ),
      numX ( 0 ), numY ( 0 ), numZ ( 1 ) {
    init();
  }

  Array ( int NumX, int NumY, int NumZ, DataType *Data, aol::CopyFlag copyFlag = aol::FLAT_COPY ) : // default: flat copy!
      aol::Vector<DataType> ( Data, NumX * NumY * NumZ, copyFlag ),
      numX ( NumX ), numY ( NumY ), numZ ( NumZ ) {
    init();
  }

  /*! Generating an array of dimension \f$NumX\cdot NumY\cdot NumZ\f$.
     */
  Array ( int NumX, int NumY, int NumZ = 1 ) :
      aol::Vector<DataType> ( NumX * NumY * NumZ ),
      numX ( NumX ), numY ( NumY ), numZ ( NumZ ) {
    init();
  }

  template <typename GridType>
  Array ( const aol::Vector<DataType> &Vec,
          const GridType &Grid, aol::CopyFlag copyFlag = aol::FLAT_COPY ) : // default: flat copy!
      aol::Vector<DataType> ( Vec, copyFlag ),
    numX ( Grid.getNumX() ), numY ( Grid.getNumY() ), numZ ( Grid.getNumZ() ) {
    if ( static_cast<int> ( Vec.size() ) != numX*numY*numZ ) {
      char error[1024];
      sprintf ( error, "Array::Array( const aol::Vector<DataType> &, const GridType & ): Vectorlength = %d should be equal to size of array, which is %d", Vec.size(), numX*numY*numZ );
      throw ( aol::Exception ( error, __FILE__, __LINE__ ) );
    }
    init ();
  }

  Array ( const aol::Vector<DataType> &Vec, int NumX, int NumY, int NumZ = 1, aol::CopyFlag copyFlag = aol::FLAT_COPY ) : // default: flat copy!
  aol::Vector<DataType> ( Vec, copyFlag ),
  numX ( NumX ), numY ( NumY ), numZ ( NumZ ) {
    if ( static_cast<int> ( Vec.size() ) != numX*numY*numZ ) {
      char error[1024];
      sprintf ( error, "Array::Array( const aol::Vector<DataType> &, int, int, int ): Vectorlength = %d should be equal to size of array, which is %d", Vec.size(), numX*numY*numZ );
      throw ( aol::Exception ( error, __FILE__, __LINE__ ) );
    }
    init ();
  }

  Array ( int Width, DataType *Data, aol::CopyFlag copyFlag = aol::FLAT_COPY ) : // default: flat copy
      aol::Vector<DataType> ( Data, Width * Width * Width, copyFlag ),
      numX ( Width ), numY ( Width ), numZ ( Width ) {
    init();
  }

  /** This is the copy constructor.
   * @param Ref determines whether a flat copy (reference only) or a deep copy (entries copied, default) shall be made
   */
  explicit Array ( const Array<DataType> &Arr, aol::CopyFlag copyFlag = aol::DEEP_COPY ) :
      aol::Vector<DataType> ( Arr, copyFlag ),
      numX ( Arr.numX ),     numY ( Arr.numY ),     numZ ( Arr.numZ ),     _offset ( Arr._offset ) {
    // do not call init!
  }

  explicit Array ( int Width ) :
    aol::Vector<DataType> ( Width * Width * Width ) ,
    numX ( Width ), numY ( Width ), numZ ( Width ) {
    init();
  }

  explicit Array ( const GridStructure &Grid ) {
    numZ = Grid.getNumZ();
    numY = Grid.getNumY();
    numX = Grid.getNumX();
    reallocate ( numX, numY, numZ );
    init();
  }

  explicit Array ( const CoordType &Dim ) :
    aol::Vector<DataType> ( Dim.x() * Dim.y() * Dim.z() ) ,
    numX ( Dim.x() ), numY ( Dim.y() ), numZ ( Dim.z() ) {
    init();
  }

  virtual ~Array() {}

  //! Change size of the array, initializing full array with zero.
  void reallocate ( const int NumX, const int NumY, const int NumZ = 1 ) {
    if ( ( static_cast<int64_t> ( NumX ) * NumY * NumZ ) > std::numeric_limits<int>::max() )
      throw aol::Exception ( "Data too big for indexing with int!", __FILE__, __LINE__ );

    aol::Vector<DataType>::reallocate ( NumX * NumY * NumZ );
    numX = NumX;
    numY = NumY;
    numZ = NumZ;
    init ();
    this->setZero();
  }

  //! Inline-function for index-calculation from coordinates
  inline int index ( int X, int Y, int Z ) const;

  inline int index ( int X, int Y ) const;

  inline int index ( const CoordType &Coords ) const {
    return index ( Coords.x(), Coords.y(), Coords.z() );
  }

  inline int index ( const aol::Vec<2, short> &Coords ) const {
    return index ( Coords[0], Coords[1] );
  }

  const aol::Vec3<int>& getOffset() const {
    return _offset;
  }

  int getIndexOffset ( int dx, int dy, int dz ) const {
    return dx + dy * numX + dz * numX * numY;
  }

  Array<DataType> &operator= ( const Array<DataType> &Other ) {
    if ( this == &Other ) { // correct self-assignment
      return ( *this );
    }

    if ( this->numX != Other.numX ||  this->numY != Other.numY || this->numZ != Other.numZ ) {
      throw aol::Exception ( "trying to copy incompatible arrays.\n", __FILE__, __LINE__ );
    }
    aol::Vector<DataType>::operator= ( Other );
    return *this;
  }


  Array<DataType> &operator= ( const aol::Vector<DataType> &Data ) {
    aol::Vector<DataType>::operator= ( Data );
    return *this;
  }

  /** Returns the array entry at Coords */
  DataType get ( const CoordType &Coords ) const {
    return get ( Coords.x(), Coords.y(), Coords.z() );
  }

  DataType get ( const aol::Vec<2, short> &Coords ) const {
    return get ( Coords[0], Coords[1] );
  }
  
  DataType get ( const aol::Vec<2, int> &Coords ) const {
    return get ( Coords[0], Coords[1] );
  }

  /** Returns the array entry at the position (X,Y,Z) */
  DataType get ( int X, int Y, int Z ) const {
    return this->_pData[ Array<DataType>::index ( X, Y, Z ) ];
  }

  DataType get ( int X, int Y ) const {
    //    cerr << X << ", " << Y << " :::";
    return this->_pData[ Array<DataType>::index ( X, Y ) ];
  }

  using aol::Vector<DataType>::get;

  DataType& getReference ( int X, int Y, int Z ) {
    return this->_pData[ index ( X, Y, Z ) ];
  }

  DataType& getReference ( int X, int Y ) {
    return this->_pData[ index ( X, Y ) ];
  }

  DataType& getReference ( const CoordType &Coords ) {
    return getReference ( Coords.x(), Coords.y(), Coords.z() );
  }

  DataType& getReference ( const aol::Vec<2, short> &Coords ) {
    return getReference ( Coords[0], Coords[1] );
  }


  //! Returns array entry with periodic boundary clipping.
  DataType getPeriodic ( int X, int Y, int Z ) const {
    if ( X < 0 ) X += numX;
    else if ( X > numX - 1 ) X = X % numX;
    if ( Y < 0 ) Y += numY;
    else if ( Y > numY - 1 ) Y = Y % numY;
    if ( Z < 0 ) Z += numZ;
    else if ( Z > numZ - 1 ) Z = Z % numZ;
    return get ( X, Y, Z );
  }

  //! Returns array entry with forcing to the boundary.
  DataType getClip ( int X, int Y, int Z ) const {
    if ( X < 0 ) X = 0;
    else if ( X > numX - 1 ) X = numX - 1;
    if ( Y < 0 ) Y = 0;
    else if ( Y > numY - 1 ) Y = numY - 1;
    if ( Z < 0 ) Z = 0;
    else if ( Z > numZ - 1 ) Z = numZ - 1;
    return get ( X, Y, Z );
  }

  void set ( const CoordType &Coords, DataType value ) {
    set ( Coords.x(), Coords.y(), Coords.z(), value );
  }

  //! Sets a value of the array.
  void set ( int X, int Y, int Z, DataType value ) {
    this->_pData[ index ( X, Y, Z ) ] = value;
  }

  void set ( const aol::Vec<2, short> &Coords, DataType value ) {
    set ( Coords[0], Coords[1], value );
  }
  
  void set ( const aol::Vec<2, int> &Coords, DataType value ) {
    set ( Coords[0], Coords[1], value );
  }

  //! Sets a value of the array.
  void set ( int X, int Y, DataType value ) {
    this->_pData[ index ( X, Y ) ] = value;
  }

  using aol::Vector<DataType>::set;

  void add ( const CoordType &Coords, DataType value ) {
    add ( Coords.x(), Coords.y(), Coords.z(), value );
  }

  void add ( int X, int Y, int Z, DataType value ) {
    this->_pData[ index ( X, Y, Z ) ] += value;
  }

  void add ( const aol::Vec<2, short> &Coords, DataType value ) {
    add ( Coords[0], Coords[1], value );
  }

  void add ( int X, int Y, DataType value ) {
    this->_pData[ index ( X, Y ) ] += value;
  }

  using aol::Vector<DataType>::add;

  void add ( const Array<DataType> &AddedArray, DataType Factor ) {
    aol::Vector<DataType>::addMultiple ( AddedArray, Factor );
  }

  //! Returns width.
  int getNumX() const {
    return numX;
  }

  //! Returns height.
  int getNumY() const {
    return numY;
  }

  //! Returns depth.
  int getNumZ() const {
    return numZ;
  }

  aol::Vec3<int> getSize() const {
    return ( aol::Vec3<int> ( numX, numY, numZ ) );
  }

  /** Set the extent of this array in the given Vec<int> and furthermore return
   *  the dimension 2.
   */
  Dimension getDimensions ( aol::Vec3<int> &d ) const {
    d.set ( numX, numY, numZ );
    if ( numZ == 1 ) return QC_2D;
    return QC_3D;
  }

  //! checks whether both classes have equal dimensions
  int equalDimensions ( const Array<DataType> &Other ) const {
    return ( Other.getNumX() == numX &&
             Other.getNumY() == numY &&
             Other.getNumZ() == numZ );
  }

  //! copy a slice from *this (3d) into Dest (2d)
  void getSlice ( const Comp Comp, int Index, Array<DataType> &Dest ) const;

  //! copy a slice from Scr (2d) into *this (3d)
  void putSlice ( const Comp Comp, int Index, const Array<DataType> &Src );

  void printSlices ( ostream& out = cerr ) const {
    aol::MixedFormat printer ( 4, 4 );
    for ( int z = 0; z < numZ; ++z ) {
      for ( int y = 0; y < numY; ++y ) {
        for ( int x = 0; x < numX; ++x ) {
          out << printer.operator() ( get( x, y, z ) );
        }
        out << endl;
      }
      out << endl;
    }
    out << endl;
  }


private:
  //! Initialize some variables, do not touch data!
  void init() {
    // test if index range is large enough
    if ( !aol::productWillFit ( numX, numY ) ||
         !aol::productWillFit ( numX * numY, numZ ) ) {
      stringstream err;
      err << "Warning: array of size " << numX << " * "
          << numY << " * " << numZ << " cannot be indexed!";
      throw aol::OutOfBoundsException ( err.str(), __FILE__, __LINE__ );
    }
    _offset.set ( 1, numX, numY*numX );
  }


protected:
  int numX, numY, numZ;
  aol::Vec3<int> _offset;

  static const unsigned int DIFF_STENCIL = 3;
  static const float DIFF_STD_SIGMA; // float, independent of data type
  static const unsigned int MED_FILTER_WIDTH = 3;
  static const unsigned int MED_FILTER_SIZE = 3 * 3 * 3;

#ifdef BOUNDS_CHECK
  inline bool boundsCheck ( const int i, const int j, const int k, const char* msg, const char* fi, const int li ) const {
    const bool isIn = ( i >= 0 && i < numX &&
                        j >= 0 && j < numY &&
                        k >= 0 && k < numZ );
    if ( !isIn ) {
      char errmsg[1024];
      sprintf( errmsg, "%s %d %d %d (upper bounds: %d %d %d)", msg, i, j, k, numX, numY, numZ );
      throw aol::OutOfBoundsException ( errmsg, fi, li );
    }
    return ( isIn );
  }
#endif

private:
  void resize ( const int ) {
    throw aol::Exception("qc::Array::resize( int ) preserving contents does not make sense.", __FILE__, __LINE__ );
  }

  void reallocate ( const int ) {
    throw aol::Exception("qc::Array::reallocate( int ) preserving contents does not make sense.", __FILE__, __LINE__ );
  }

  void pushBack ( const DataType ) {
    throw aol::Exception("qc::Array::pushBack( DataType ) does not make sense.", __FILE__, __LINE__ );
  }

  void pushBackValues ( const aol::Vector<DataType> & ) {
    throw aol::Exception("qc::Array::pushBackValues( aol::Vector<DataType> ) does not make sense.", __FILE__, __LINE__ );
  }
};

template<typename DataType>
int Array<DataType>::index ( int X, int Y, int Z ) const {
#ifdef BOUNDS_CHECK
  if ( !boundsCheck ( X, Y, Z, "qc::Array::index: Index out of bounds", __FILE__, __LINE__ ) ) return 0;
#endif

  return ( qc::ILexCombine3 ( X, Y, Z, numX, numY ) );
}

template<typename DataType>
int Array<DataType>::index ( int X, int Y ) const {
#ifdef BOUNDS_CHECK
  if ( !boundsCheck ( X, Y, 0, "qc::Array::index out of bounds", __FILE__, __LINE__ ) ) return 0;
#endif

  return ( qc::ILexCombine2 ( X, Y, numX ) );
}


}

#endif
