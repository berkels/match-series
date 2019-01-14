#ifndef __BITARRAY_H
#define __BITARRAY_H

#include <bitVector.h>
#include <quoc.h>
#include <gridBase.h>
#include <vec.h>
#include <qmException.h>
#include <multiDObject.h>

namespace qc {

/** Basis class for template specialization
 */
template < qc::Dimension Dim >
class BitArray;

/**
 *  \author Berkels
 */
template<>
class BitArray<qc::QC_1D> : public aol::BitVector {
public:
  //! Construct from GridSize
  explicit BitArray ( const GridSize<QC_1D> &GridSize ) : BitVector ( GridSize.getNumX() ) {}
};

/** An extension of aol::BitVector that provides access functions similar to 2D qc::Arrays allowing access via (x,y).
 *  \author Horn
 */
template<>
class BitArray<qc::QC_2D> : public aol::BitVector {
protected:
  int _numX, _numY;


public:
  BitArray ( ) : BitVector( 0 ), _numX( 0 ), _numY ( 0 ) {
  }

  //! Construct from given dimensions
  BitArray ( const int NumX, const int NumY ) :
      BitVector ( NumX * NumY ),
      _numX ( NumX ),
      _numY ( NumY ) {}

  //! Construct from GridSize
  explicit BitArray ( const GridSize<QC_2D> &gridSize ) :
      BitVector ( gridSize.getNumX() * gridSize.getNumY() ),
      _numX ( gridSize.getNumX() ),
      _numY ( gridSize.getNumY() ) {}

  //! Construct from GridStructure
  explicit BitArray ( const GridStructure &gridStructure ) :
      BitVector ( gridStructure.getNumX() * gridStructure.getNumY() ),
      _numX ( gridStructure.getNumX() ),
      _numY ( gridStructure.getNumY() ) {}

  //! Construct from Vec2
  explicit BitArray ( const aol::Vec2<int> & Size ) :
      BitVector ( Size[0] * Size[1] ),
      _numX ( Size[0] ),
      _numY ( Size[1] ) {}

  //! Construct from Vec3, ignoring last component
  explicit BitArray ( const aol::Vec3<int> & Size ) :
      BitVector ( Size[0] * Size[1] ),
      _numX ( Size[0] ),
      _numY ( Size[1] ) {}

  //! Copy constructor
  explicit BitArray ( const BitArray<qc::QC_2D> &other, aol::CopyFlag copyFlag = aol::DEEP_COPY ) :
      BitVector ( other ),
      _numX ( other._numX ),
      _numY ( other._numY ) {
    if ( copyFlag != aol::DEEP_COPY )
      throw aol::Exception("qc::BitArray<qc::QC_2D> can only be DEEP_COPied via copy constructor so far", __FILE__, __LINE__ );
  }

  using BitVector::get;
  using BitVector::set;

  inline void set ( const int ix, const int iy, const bool Value ) {
#ifdef BOUNDS_CHECK
    boundsCheck ( ix, iy, "set" );
#endif
    BitVector::set ( ar_index ( ix, iy ), Value );
  }

  inline void set ( const CoordType &Coords, const bool Value ) {
#ifdef BOUNDS_CHECK
    boundsCheck ( Coords.x(), Coords.y(), "set" );
#endif
    BitVector::set ( ar_index ( Coords.x(), Coords.y() ), Value );
  }

  inline void set ( const aol::Vec<2, int> &Coords, const bool Value ) {
#ifdef BOUNDS_CHECK
    boundsCheck ( Coords[0], Coords[1], "set" );
#endif
    BitVector::set ( ar_index ( Coords[0], Coords[1] ), Value );
  }

  inline bool get ( const int ix, const int iy ) const {
#ifdef BOUNDS_CHECK
    boundsCheck ( ix, iy, "get" );
#endif
    return ( ( *this ) [ ar_index ( ix, iy ) ] );
  }

  inline bool get ( const CoordType &pos ) const {
#ifdef BOUNDS_CHECK
    boundsCheck ( pos[0], pos[1], "get" );
#endif
    return ( ( *this ) [ ar_index ( pos[0], pos[1] ) ] );
  }

  inline bool get ( const aol::Vec<2, int> &pos ) const {
#ifdef BOUNDS_CHECK
    boundsCheck ( pos[0], pos[1], "get" );
#endif
    return ( ( *this ) [ ar_index ( pos[0], pos[1] ) ] );
  }

  inline bool get ( const aol::Vec<2, short> &pos ) const {
    return get ( pos[0], pos[1] );
  }

  int getNumX() const {
    return( _numX );
  }

  int getNumY() const {
    return( _numY );
  }

  aol::Vec2<int> getSize() const {
    return ( aol::Vec2<int> ( _numX, _numY ) );
  }

  void reallocate ( const int numX, const int numY ) {
    aol::BitVector::reallocate ( numX * numY );
    _numX = numX;
    _numY = numY;
  }

  template< typename Structure >
  void reallocate ( const Structure &struc ) {
    reallocate ( struc.getNumX(), struc.getNumY() );
  }

  /**
   * Returns "true", if all vertices of the element have the value "true".
   *
   * \author Wirth
   */
  inline bool elementTrue ( const qc::Element &El ) const {
    return ( get ( El[0], El[1] ) && get ( El[0] + 1, El[1] ) && get ( El[0], El[1] + 1 ) && get ( El[0] + 1, El[1] + 1 ) );
  }

  void print ( ostream &out = cerr ) const {
    for (int x = 0; x < getNumX(); ++x) {
      for (int y = 0; y < getNumY(); ++y)
        out << (get ( x, y ) ? "*" : " ");
      out << endl;
    }
  }

  /**
   * Saves the BitArray<qc::QC_2D> as ASCII PGM (0 <> false, 255 <> true).
   *
   * \author Berkels
   */
  void save ( const char* FileName ) const;

  /**
   * Dilates the mask given by the "true" pixels by one pixel in
   * horizontal and vertical direction.
   *
   * \author Berkels
   */
  void dilateByOne ( );

  /**
   * Dilates the mask given by the "true" pixels by one pixel in direction Dir.
   *
   * \author Wirth
   */
  void dilateByOne( short Dir ) {
    for ( int x = 0; x < _numX; ++x )
      for ( int y = 0; y < _numY; ++y ) {
        CoordType coord( x, y, 0 );
        coord[Dir]--;
        coord[Dir] = aol::Max( coord[Dir], static_cast<short>( 0 ) );
        if ( get( x, y ) )
          set( coord, true );
      }
    CoordType maxCoord( _numX - 1, _numY - 1, 0 );
    for ( int x = _numX - 1; x >= 0; --x )
      for ( int y = _numY - 1; y >= 0; --y ) {
        CoordType coord( x, y, 0 );
        coord[Dir]++;
        coord[Dir] = aol::Min( coord[Dir], maxCoord[Dir] );
        if ( get( x, y ) )
          set( coord, true );
      }
  }
  
  /**
   * Erodes the mask given by the "true" pixels by one pixel in
   * horizontal and vertical direction.
   *
   * \author Mevenkamp
   */
  void erodeByOne ( );
  
  //! Erodes the mask given by the "true" pixels by one pixel in horizontal and vertical direction Size times
  void erodeBy ( short Size ) {
    for ( short it = 0; it < Size ; ++it )
      erodeByOne ( );
  }
  
  //! Dilates the mask given by the "true" pixels by one pixel in horizontal and vertical direction Size times
  void dilateBy ( short Size ) {
    for ( short it = 0; it < Size ; ++it )
      dilateByOne ( );
  }

  //! lossless save in binary quocmesh format
  void saveToFile ( const char *filename ) const;

  //! load from file in quocmesh format, no conversion of data types
  void loadFromFile ( const char *filename );
  
  void crop ( const aol::Vec<2, int> &CropStart, const aol::Vec<2, int> &CropSize )  {
    qc::BitArray<qc::QC_2D> copy ( *this, aol::DEEP_COPY );
    const int
      cpyStartX = aol::Clamp<int> ( CropStart[0], 0, _numX-1 ),
      cpyStartY = aol::Clamp<int> ( CropStart[1], 0, _numY-1 ),
      cpyNumX = aol::Min ( _numX - cpyStartX, CropSize[0] ),
      cpyNumY = aol::Min ( _numY - cpyStartY, CropSize[1] );
    
    this->reallocate ( CropSize[0], CropSize[1] );

    for ( int y = 0; y < cpyNumY; ++y ) {
      for ( int x = 0; x < cpyNumX; ++x ) { {
        this->set ( x, y, copy.get ( cpyStartX + x, cpyStartY + y ) );
      }
      }
    }
  }

protected:
  inline int ar_index ( const int ix, const int iy ) const {
    return ( iy * _numX + ix );
  }

  void resize ( const int newNumX, const int newNumY ) {
    qc::BitArray<qc::QC_2D> copyOfMe ( *this, aol::DEEP_COPY );
    const int
      cpyNumX = aol::Min ( _numX, newNumX ),
      cpyNumY = aol::Min ( _numY, newNumY );

    this->reallocate ( newNumX, newNumY ); // sic: reallocate (deleting old contents), then copy appropriate values

    for ( int y = 0; y < cpyNumY; ++y ) {
      for ( int x = 0; x < cpyNumX; ++x ) { {
        this->set ( x, y,  copyOfMe.get ( x, y ) );
        }
      }
    }
  }

  void resize ( const GridStructure & Grid ) {
    QUOC_ASSERT ( Grid.getDimOfWorld() == QC_2D );
    resize ( Grid.getNumX(), Grid.getNumY() );
  }

#ifdef BOUNDS_CHECK
  inline void boundsCheck ( const int ix, const int iy, const char* where ) const {
    if ( ! ( ( ix < _numX ) && ( iy < _numY ) && ( ix >= 0 ) && ( iy >= 0 ) ) ) {
      char errmsg[1024];
      sprintf( errmsg, "qc::BitArray<qc::QC_2D>::%s %d %d are out of bounds (%d %d)!", where, ix, iy, _numX, _numY );
      throw aol::OutOfBoundsException( errmsg, __FILE__, __LINE__ );
    }
  }
#endif

private:
  void resize ( const int ) {
    throw aol::Exception("qc::BitArray<qc::QC_2D>::resize( int ) preserving contents does not make sense.", __FILE__, __LINE__ );
  }

  void reallocate ( const int ) {
    throw aol::Exception("qc::BitArray<qc::QC_2D>::reallocate( int ) preserving contents does not make sense.", __FILE__, __LINE__ );
  }

  template< typename AnyType >
  void set ( const int, const int, const AnyType& );

  template< typename AnyType >
  void set ( const CoordType &, const AnyType& );

  template< typename AnyType >
  void set ( const aol::Vec2<int>, const AnyType& );

  // class BitArray<qc::QC_2D> end
};


/** An extension of aol::BitVector that provides access functions
 *  similar to 3D qc::Arrays allowing access via (x,y,z).
 *  \author Schwen
 */
template<>
class BitArray<qc::QC_3D> : public aol::BitVector, public qc::MultiDStorageObject<qc::QC_3D> {

public:
  BitArray ( ) : BitVector( 0 ), qc::MultiDStorageObject<qc::QC_3D> ( 0, 0, 0 ) {
  }

  //! Construct from given dimensions
  BitArray ( const int NumX, const int NumY, const int NumZ ) :
    aol::BitVector ( NumX * NumY * NumZ ),
    qc::MultiDStorageObject<qc::QC_3D> ( NumX, NumY, NumZ ) {
  }

  explicit BitArray ( const aol::Vec3<int> &Size ) :
    aol::BitVector ( Size[0] * Size[1] * Size[2] ),
    qc::MultiDStorageObject<qc::QC_3D> ( Size[0], Size[1], Size[2] ) {
  }

  //! Construct from gridSize
  explicit BitArray ( const GridSize<QC_3D> &gridSize ) :
    aol::BitVector ( gridSize.getNumX() * gridSize.getNumY() * gridSize.getNumZ() ),
    qc::MultiDStorageObject<qc::QC_3D> ( gridSize.getNumX(), gridSize.getNumY(), gridSize.getNumZ() ) {
  }

  //! Construct from gridStructure
  explicit BitArray ( const GridStructure &gridStructure ) :
    aol::BitVector ( gridStructure.getNumX() * gridStructure.getNumY() * gridStructure.getNumZ() ),
    qc::MultiDStorageObject<qc::QC_3D> ( gridStructure.getNumX(), gridStructure.getNumY(), gridStructure.getNumZ() ) {
  }

  //! Copy constructor
  explicit BitArray ( const BitArray<qc::QC_3D> &other, aol::CopyFlag copyFlag = aol::DEEP_COPY ) :
    aol::BitVector ( other ),
    qc::MultiDStorageObject<qc::QC_3D> ( other ) {
    if ( copyFlag != aol::DEEP_COPY )
      throw aol::Exception("qc::BitArray<qc::QC_3D> can only be DEEP_COPied via copy constructor so far", __FILE__, __LINE__ );
    }

  using aol::BitVector::get;
  using aol::BitVector::set;

  inline void set ( const int ix, const int iy, const int iz, const bool value ) {
#ifdef BOUNDS_CHECK
    boundsCheck ( ix, iy, iz, "qc::BitArray<qc::QC_3D>::set", __FILE__, __LINE__ );
#endif
    aol::BitVector::set ( this->oneDIndex ( ix, iy, iz ), value );
  }

  inline void set ( const CoordType &pos, const bool value ) {
#ifdef BOUNDS_CHECK
    boundsCheck ( pos.x(), pos.y(), pos.z(), "qc::BitArray<qc::QC_3D>::set", __FILE__, __LINE__ );
#endif
    aol::BitVector::set ( this->oneDIndex ( pos.x(), pos.y(), pos.z() ), value );
  }

  inline void set ( const aol::Vec3<int> &pos, const bool value ) {
#ifdef BOUNDS_CHECK
    boundsCheck ( pos.x(), pos.y(), pos.z(), "qc::BitArray<qc::QC_3D>::set", __FILE__, __LINE__ );
#endif
    aol::BitVector::set ( this->oneDIndex ( pos.x(), pos.y(), pos.z() ), value );
  }

  inline bool get ( const int ix, const int iy, const int iz ) const {
#ifdef BOUNDS_CHECK
    boundsCheck ( ix, iy, iz, "qc::BitArray<qc::QC_3D>::get", __FILE__, __LINE__ );
#endif
    return ( ( *this ) [ this->oneDIndex ( ix, iy, iz ) ] );
  }

  inline bool get ( const CoordType &pos ) const {
#ifdef BOUNDS_CHECK
    boundsCheck ( pos.x(), pos.y(), pos.z(), "qc::BitArray<qc::QC_3D>::get", __FILE__, __LINE__ );
#endif
    return ( ( *this ) [ this->oneDIndex ( pos.x(), pos.y(), pos.z() ) ] );
  }

  inline bool get ( const aol::Vec3<int> &pos ) const {
#ifdef BOUNDS_CHECK
    boundsCheck ( pos.x(), pos.y(), pos.z(), "qc::BitArray<qc::QC_3D>::get", __FILE__, __LINE__ );
#endif
    return ( ( *this ) [ this->oneDIndex ( pos.x(), pos.y(), pos.z() ) ] );
  }

  //! copy slice from *this into Dest
  void putSlice ( const Comp comp, int NumComp, BitArray<qc::QC_2D> & Dest ) const {
    switch ( comp ) {
      case QC_X:
        QUOC_ASSERT ( getNumY() == Dest.getNumX() );
        QUOC_ASSERT ( getNumZ() == Dest.getNumY() );
        for (int y = 0; y < getNumY(); ++y)
          for (int z = 0; z < getNumZ(); ++z)
            Dest.set ( y, z, get ( NumComp, y, z ) );
        break;

      case QC_Y:
        QUOC_ASSERT ( getNumX() == Dest.getNumX() );
        QUOC_ASSERT ( getNumZ() == Dest.getNumY() );
        for (int x = 0; x < getNumX(); ++x)
          for (int z = 0; z < getNumZ(); ++z)
            Dest.set ( x, z, get ( x, NumComp, z ) );
        break;

      case QC_Z:
        QUOC_ASSERT ( getNumX() == Dest.getNumX() );
        QUOC_ASSERT ( getNumY() == Dest.getNumY() );
        for (int x = 0; x < getNumX(); ++x)
          for (int y = 0; y < getNumY(); ++y)
            Dest.set ( x, y, get ( x, y, NumComp ) );
        break;
        
      default:
        throw aol::UnimplementedCodeException ( "Unknown qc::Comp", __FILE__, __LINE__ );
    }
  }

  //! copy a slice from Src into *this
  void getSlice ( const Comp comp, int NumComp, const BitArray<qc::QC_2D> & Src ) {
    switch ( comp ) {
      case QC_X:
        QUOC_ASSERT ( getNumY() == Src.getNumX() );
        QUOC_ASSERT ( getNumZ() == Src.getNumY() );
        for (int y = 0; y < getNumY(); ++y)
          for (int z = 0; z < getNumZ(); ++z)
            set ( NumComp, y, z, Src.get ( y, z ) );
        break;

      case QC_Y:
        QUOC_ASSERT ( getNumX() == Src.getNumX() );
        QUOC_ASSERT ( getNumZ() == Src.getNumY() );
        for (int x = 0; x < getNumX(); ++x)
          for (int z = 0; z < getNumZ(); ++z)
            set ( x, NumComp, z, Src.get ( x, z ) );
        break;

      case QC_Z:
        QUOC_ASSERT ( getNumX() == Src.getNumX() );
        QUOC_ASSERT ( getNumY() == Src.getNumY() );
        for (int x = 0; x < getNumX(); ++x)
          for (int y = 0; y < getNumY(); ++y)
            set ( x, y, NumComp, Src.get ( x, y ) );
        break;
        
      default:
        throw aol::UnimplementedCodeException ( "Unknown qc::Comp", __FILE__, __LINE__ );
    }
  }

  /**
   * Returns "true", if all vertices of the element have the value "true".
   *
   * \author Wirth
   */
  inline bool elementTrue ( const qc::Element &El ) const {
    return ( get ( El[0], El[1], El[2] ) && get ( El[0] + 1, El[1], El[2] ) && get ( El[0], El[1] + 1, El[2] ) && get ( El[0] + 1, El[1] + 1, El[2] ) &&
             get ( El[0], El[1], El[2] + 1 ) && get ( El[0] + 1, El[1], El[2] + 1 ) && get ( El[0], El[1] + 1, El[2] + 1 ) && get ( El[0] + 1, El[1] + 1, El[2] + 1 ) );
  }

  /**
   * Saves the BitArray<qc::QC_3D> as (0 <> false, 255 <> true).
   *
   * \author Wirth
   */
  void save ( const char* FileName ) const;

  void printSlices ( ostream &out = cerr ) const {
    for ( int z = 0; z < this->_numZ; ++z ) {
      for ( int y = 0; y < this->_numY; ++y ) {
        for ( int x = 0; x < this->_numX; ++x ) {
          out << get( x, y, z ) ;
        }
        out << endl;
      }
      out << endl;
    }
    out << endl;
  }

  /**
   * Dilates the mask given by the "true" pixels by one pixel in direction Dir.
   *
   * \author Wirth
   */
  void dilateByOne( short Dir ) {
    for ( int x = 0; x < this->_numX; ++x )
      for ( int y = 0; y < this->_numY; ++y )
        for ( int z = 0; z < this->_numZ; ++z ) {
          CoordType coord( x, y, z );
          coord[Dir]--;
          coord[Dir] = aol::Max( coord[Dir], static_cast<short>( 0 ) );
          if ( get( x, y, z ) )
            set( coord, true );
        }
    CoordType maxCoord( this->_numX - 1, this->_numY - 1, this->_numZ - 1 );
    for ( int x = this->_numX - 1; x >= 0; --x )
      for ( int y = this->_numY - 1; y >= 0; --y )
        for ( int z = this->_numZ - 1; z >= 0; --z ) {
          CoordType coord( x, y, z );
          coord[Dir]++;
          coord[Dir] = aol::Min( coord[Dir], maxCoord[Dir] );
          if ( get( x, y, z ) )
            set( coord, true );
        }
  }

  /**
   * Dilates the mask given by the "true" pixels by one pixel in all directions.
   *
   * \author Wirth
   */
  void dilateByOne() {
    for ( int i = 0; i < 3; ++i )
      dilateByOne( i );
  }

  //! Change size of BrickArray<qc::QC_3D> deleting old contents
  void reallocate ( const int numX, const int numY, const int numZ ) {
    aol::BitVector::reallocate ( numX * numY * numZ );
    this->changeSizeTo ( numX, numY, numZ );
  }

  using qc::MultiDStorageObject<qc::QC_3D>::size;
  using qc::MultiDStorageObject<qc::QC_3D>::reallocate;

  //! Change size of BitArray<qc::QC_3D> preserving old contents at correct positions
  void resize ( const int newNumX, const int newNumY, const int newNumZ ) {
    qc::BitArray<qc::QC_3D> copyOfMe ( *this, aol::DEEP_COPY );
    const int
      cpyNumX = aol::Min ( this->_numX, newNumX ),
      cpyNumY = aol::Min ( this->_numY, newNumY ),
      cpyNumZ = aol::Min ( this->_numZ, newNumZ );

    this->reallocate ( newNumX, newNumY, newNumZ ); // sic: reallocate (deleting old contents), then copy appropriate values

    for ( int z = 0; z < cpyNumZ; ++z ) {
      for ( int y = 0; y < cpyNumY; ++y ) {
        for ( int x = 0; x < cpyNumX; ++x ) {
          this->set ( x, y, z,  copyOfMe.get ( x, y, z ) );
        }
      }
    }
  }

  using qc::MultiDStorageObject<qc::QC_3D>::resize;

  //! fill a finer-level BitArray which "interpolates" this BitArray
  //! A node is set to "true" only if it is contained in the coarser
  //! array and is marked "true" there or if all its neighbours (l_infty
  //! distance 1) are marked "true".
  //! \note That means that prolongating a boolean array is not symmetric.
  //! If a node should only be marked "false" if all its coarse-grid
  //! neighbours are marked "false", you have to invert the coarse array,
  //! prolongate, and re-invert both arrays afterwards.
  //!
  //! \warning untested
  //!
  //! \author von Deylen
  void prolongateTo ( qc::BitArray<qc::QC_3D> & finerArray ) const;

  //! copies only nodes that are present in both arrays
  void restrictTo ( qc::BitArray<qc::QC_3D> & coarserArray ) const;

  //! lossless save in binary quocmesh format
  void saveToFile ( const char *filename ) const;

  //! load from file in quocmesh format, no conversion of data types
  void loadFromFile ( const char *filename );

private:
  void resize ( const int ) {
    throw aol::Exception("qc::BitArray<qc::QC_3D>::resize( int ) preserving contents does not make sense.", __FILE__, __LINE__ );
  }

  void reallocate ( const int ) {
    throw aol::Exception("qc::BitArray<qc::QC_3D>::reallocate( int ) preserving contents does not make sense.", __FILE__, __LINE__ );
  }

  template< typename AnyType >
  void set ( const int, const int, const int, const AnyType& );

  template< typename AnyType >
  void set ( const CoordType &, const AnyType& );

  template< typename AnyType >
  void set ( const aol::Vec3<int>, const AnyType& );

  // class BitArray<qc::QC_3D> end
};

/**
 * Trait that allows to convert an int into a BitArray of the corresponding dimension.
 *
 * \author Berkels
 */
template <int dim>
class BitArrayTrait {};

template <>
class BitArrayTrait<1> {
public:
  typedef qc::BitArray<qc::QC_1D> ArrayType;
};

template <>
class BitArrayTrait<2> {
public:
  typedef qc::BitArray<qc::QC_2D> ArrayType;
};

template <>
class BitArrayTrait<3> {
public:
  typedef qc::BitArray<qc::QC_3D> ArrayType;
};

}

#endif
