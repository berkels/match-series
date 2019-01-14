#ifndef __INDEXMAPPER_H
#define __INDEXMAPPER_H

#include <array.h>
#include <gridSize.h>

namespace qc {

/** On-the-fly index mapper using the usual inverse-lexicographical index mapping used in the quocmeshes: (x, y, z) is mapped to Nx*Ny*z + Nx*y + x
 *  \author Schwen
 */

template <Dimension Dim>
class OTFILexMapper;

// OTFILexMapper<qc::QC_1D> is a special case and derived from qc::FastILexMapper<qc::QC_1D>.

template <>
class OTFILexMapper<qc::QC_2D> {
protected:
  int      _Nx, _Ny;

public:
  //! standard constructor, somewhat useless ...
  OTFILexMapper ( ) : _Nx ( 0 ), _Ny ( 0 ) {
  }

  explicit OTFILexMapper ( const int Nx ) : _Nx ( Nx ), _Ny ( Nx ) {
  }

  explicit OTFILexMapper ( const int Nx, const int Ny ) : _Nx ( Nx ), _Ny ( Ny ) {
  }


  //! Constructor taking grid for node indexing or ScalarArray<QC_3D>
  template< class Struc2d >
  explicit OTFILexMapper ( const Struc2d &Str ) : _Nx ( Str.getNumX() ), _Ny ( Str.getNumY() ) {
  }

  //! constructor from GridSize
  explicit OTFILexMapper ( const GridSize<QC_2D> & gridSize ) : _Nx ( gridSize.getNumX() ), _Ny ( gridSize.getNumY() ) {
  }

  //! copy constructor
  OTFILexMapper ( const OTFILexMapper<qc::QC_2D> &other ) : _Nx ( other._Nx ), _Ny ( other._Ny ) {
  }

  //! assignment operator
  OTFILexMapper<qc::QC_2D>& operator= ( const OTFILexMapper<qc::QC_2D> &other ) {
    _Nx = other._Nx;
    _Ny = other._Ny;
    return ( *this );
  }

public:
  //! Change the size of the corresponding grid
  void resize ( const int Nx, const int Ny ) {
    _Nx = Nx;
    _Ny = Ny;
  }

  //! Compute global index from components
  inline int getGlobalIndex ( const int x, const int y ) const {
#ifdef BOUNDS_CHECK
    boundsCheck ( x, y, "qc::OTFILexMapper<qc::QC_2D>::getGlobalIndex: coordinates out of bounds: %d %d %d %d", __FILE__, __LINE__ );
#endif
    return ( qc::ILexCombine2 ( x, y, _Nx ) );
  }

  inline int getGlobalIndex ( const CoordType &Coord ) const {
    return ( getGlobalIndex ( Coord[0], Coord[1] ) );
  }

  //! get global index
  inline int operator() ( const int x, const int y ) const {
    return ( getGlobalIndex ( x, y ) );
  }

  //! get global index
  inline int operator() ( const CoordType &Coord ) const {
    return ( getGlobalIndex ( Coord[0], Coord[1] ) );
  }

   //! Split global index into components
  inline void splitGlobalIndex ( const int global, int &x, int &y ) const {
    x = global % _Nx;
    y = global / _Nx;
#ifdef BOUNDS_CHECK
    boundsCheck ( x, y, "qc::OTFILexMapper<qc::QC_2D>::splitGlobalIndex: coordinates out of bounds: %d %d %d %d", __FILE__, __LINE__ );
#endif
  }

   //! Split global index into components
  inline void splitGlobalIndex ( const int global, CoordType &Coord ) const {
    Coord[0] = global % _Nx;
    Coord[1] = global / _Nx;
#ifdef BOUNDS_CHECK
    boundsCheck ( Coord[0], Coord[1], "qc::OTFILexMapper<qc::QC_2D>::splitGlobalIndex: coordinates out of bounds: %d %d %d %d", __FILE__, __LINE__ );
#endif
  }

  //! \warning returning object is slower than filling one passed by reference
  inline CoordType splitGlobalIndex ( const int global ) const {
    const qc::CoordType Coord ( global % _Nx, global / _Nx, 0 );
#ifdef BOUNDS_CHECK
    boundsCheck ( Coord[0], Coord[1], "qc::OTFILexMapper<qc::QC_2D>::splitGlobalIndex: coordinates out of bounds: %d %d %d %d", __FILE__, __LINE__ );
#endif
    return ( Coord );
  }

protected:
#ifdef BOUNDS_CHECK
  inline void boundsCheck ( const int x, const int y, const char* msg, const char* fi, const int li ) const {
    if ( x < 0 || x >= _Nx || y < 0 || y >= _Ny ) {
      char errmsg[1024];
      sprintf( errmsg, msg, x, y, _Nx, _Ny );
      throw aol::Exception ( errmsg, fi, li );
    }
  }
#endif

};


template <>
class OTFILexMapper<qc::QC_3D> {
protected:
  int _Nx, _Ny, _Nz;

public:
  //! standard constructor, somewhat useless ...
  OTFILexMapper ( ) : _Nx ( 0 ), _Ny ( 0 ), _Nz ( 0 ) {
  }

  explicit OTFILexMapper ( const int Nx ) : _Nx ( Nx ), _Ny ( Nx ), _Nz ( Nx ) {
  }

  explicit OTFILexMapper ( const int Nx, const int Ny, const int Nz ) : _Nx ( Nx ), _Ny ( Ny ), _Nz ( Nz ) {
  }

  //! Constructor taking grid for node indexing or ScalarArray<QC_3D>
  template< class Struct3d >
  explicit OTFILexMapper ( const Struct3d &Str ) : _Nx ( Str.getNumX() ), _Ny ( Str.getNumY() ), _Nz ( Str.getNumZ() ) {
  }

  //! constructor from GridSize
  explicit OTFILexMapper ( const GridSize<QC_3D> & gridSize ) : _Nx ( gridSize.getNumX() ), _Ny ( gridSize.getNumY() ), _Nz ( gridSize.getNumZ() ) {
  }

  //! copy constructor
  OTFILexMapper ( const OTFILexMapper<qc::QC_3D> &other ) : _Nx ( other._Nx ), _Ny ( other._Ny ), _Nz ( other._Nz ) {
  }

  //! assignment operator
  OTFILexMapper<qc::QC_3D>& operator= ( const OTFILexMapper<qc::QC_3D> &other ) {
    _Nx = other._Nx;
    _Ny = other._Ny;
    _Nz = other._Nz;

    return ( *this );
  }

public:
  //! Change the size of the corresponding array
  void resize ( const int Nx, const int Ny, const int Nz ) {
    _Nx = Nx;
    _Ny = Ny;
    _Nz = Nz;
  }

  //! Compute global index from components
  inline int getGlobalIndex ( const int x, const int y, const int z ) const {
#ifdef BOUNDS_CHECK
    boundsCheck ( x, y, z, "qc::OTFILexMapper<qc::QC_3D>::getGlobalIndex: coordinates out of bounds: %d %d %d %d %d %d", __FILE__, __LINE__ );
#endif
    return ( qc::ILexCombine3 ( x, y, z, _Nx, _Ny ) );
  }

  //! Compute global index from components
  inline int getGlobalIndex ( const CoordType &Coord ) const {
    return ( getGlobalIndex ( Coord[0], Coord[1], Coord[2] ) );
  }

  //! Compute global index from components
  inline int getGlobalIndex ( const aol::Vec3<int> &Coord ) const {
    return ( getGlobalIndex ( Coord[0], Coord[1], Coord[2] ) );
  }

  //! get global index
  inline int operator() ( const int x, const int y, const int z ) const {
    return ( getGlobalIndex ( x, y, z ) );
  }

  //! get global index
  inline int operator() ( const CoordType &Coord ) const {
    return ( getGlobalIndex ( Coord[0], Coord[1], Coord[2] ) );
  }

  //! get global index
  inline int operator() ( const aol::Vec3<int> &Coord ) const {
    return ( getGlobalIndex ( Coord[0], Coord[1], Coord[2] ) );
  }

   //! Split global index into components
  inline void splitGlobalIndex ( const int global, int &x, int &y, int &z ) const {
    x = global %_Nx;
    y = ( global / _Nx ) % _Ny;
    z = global / ( _Nx * _Ny );
#ifdef BOUNDS_CHECK
    boundsCheck ( x, y, z, "qc::OTFILexMapper<qc::QC_3D>::splitGlobalIndex: coordinates out of bounds: %d %d %d %d %d %d", __FILE__, __LINE__ );
#endif
  }

  //! Split global index into components
  inline void splitGlobalIndex ( const int global, CoordType &Coord ) const {
    Coord[0] = global %_Nx;
    Coord[1] = ( global / _Nx ) % _Ny;
    Coord[2] = global / ( _Nx * _Ny );
#ifdef BOUNDS_CHECK
    boundsCheck ( Coord[0], Coord[1], Coord[2], "qc::OTFILexMapper<qc::QC_3D>::splitGlobalIndex: coordinates out of bounds: %d %d %d %d %d %d", __FILE__, __LINE__ );
#endif
  }

  //! \warning returning object is slower than filling one passed by reference
  inline CoordType splitGlobalIndex ( const int global ) const {
    const qc::CoordType Coord ( global %_Nx, ( global / _Nx ) % _Ny, global / ( _Nx * _Ny ) );
#ifdef BOUNDS_CHECK
    boundsCheck ( Coord[0], Coord[1], Coord[2], "qc::OTFILexMapper<qc::QC_3D>::splitGlobalIndex: coordinates out of bounds: %d %d %d %d %d %d", __FILE__, __LINE__ );
#endif
    return ( Coord );
  }

protected:
#ifdef BOUNDS_CHECK
  inline void boundsCheck ( const int x, const int y, const int z, const char* msg, const char* fi, const int li ) const {
    if ( x < 0 || x >= _Nx || y < 0 || y >= _Ny || z < 0 || z >= _Nz ) {
      char errmsg[1024];
      sprintf( errmsg, msg, x, y, z, _Nx, _Ny, _Nz );
      throw aol::Exception ( errmsg, fi, li );
    }
  }
#endif

};



/** Fast index mapper using the usual inverse-lexicographical index mapping used in the quocmeshes: (x, y, z) is mapped to Nx*Ny*z + Nx*y + x
 *  Fast refers to the fact that the multiplications are stored in a lookup table.
 *  \author Schwen
 */

template <Dimension Dim>
class FastILexMapper;

template <>
class FastILexMapper<qc::QC_1D> {
protected:
  int      _Nx;
  int _nmap[2];

public:
  //! standard constructor, somewhat useless ...
  //! The class does not have an extra copy constructor or
  //! assignment operator because deep copy of all member
  //! variables ist exactly the desired behaviour.
  FastILexMapper ( ) {
    resize ( 0 );
  }

  explicit FastILexMapper ( const int Nx ) {
    resize ( Nx );
  }

  //! Constructor taking grid for node indexing or ScalarArray<QC_1D>
  template< class Struct1d >
  explicit FastILexMapper ( const Struct1d &Str ) {
    resize ( Str.getNumX() );
  }

  //! constructor from GridSize
  explicit FastILexMapper ( const GridSize<QC_1D> & gridSize ) {
    resize ( gridSize.getNumX() );
  }

  //! Change the size of the grid
  void resize ( const int Nx ) {
    _Nx = Nx;

    _nmap[0] = 0;
    _nmap[1] = 1;
  }

  void resize ( const qc::CoordType &Coord ) {
    resize ( Coord[0] );
  }

  //! Compute global index from components
  inline int getGlobalIndex ( const int x ) const {
#ifdef BOUNDS_CHECK
    boundsCheck ( x, "qc::FastILexMapper<qc::QC_1D>::getGlobalIndex: coordinates out of bounds: %d %d", __FILE__, __LINE__ );
#endif
    return ( x );
  }

  inline int getGlobalIndex ( const CoordType &Coord ) const {
    return ( getGlobalIndex ( Coord[0] ) );
  }

  //! get global index
  inline int operator() ( const int x ) const {
    return ( getGlobalIndex ( x ) );
  }

  //! get global index
  inline int operator() ( const CoordType &Coord ) const {
    return ( getGlobalIndex ( Coord[0] ) );
  }

   //! Split global index into components
  inline void splitGlobalIndex ( const int global, int &x ) const {
    x = global;
#ifdef BOUNDS_CHECK
    boundsCheck ( x, "qc::FastILexMapper<qc::QC_1D>::splitGlobalIndex: coordinates out of bounds: %d %d", __FILE__, __LINE__ );
#endif
  }

  //! Split global index into components
  inline void splitGlobalIndex ( const int global, CoordType &Coord ) const {
    Coord[0] = global;
#ifdef BOUNDS_CHECK
    boundsCheck ( Coord[0], "qc::FastILexMapper<qc::QC_1D>::splitGlobalIndex: coordinates out of bounds: %d %d", __FILE__, __LINE__ );
#endif
  }


  //! \warning returning object is slower than filling one passed by reference
  inline CoordType splitGlobalIndex ( const int global ) const {
    const qc::CoordType Coord ( global, 0, 0 );
#ifdef BOUNDS_CHECK
    boundsCheck ( Coord[0], "qc::FastILexMapper<qc::QC_1D>::splitGlobalIndex: coordinates out of bounds: %d %d", __FILE__, __LINE__ );
#endif
    return ( Coord );
  }

  int localToGlobal ( const CoordType & Coord, int LocalIndex ) const {
    return getGlobalIndex ( Coord ) + _nmap[LocalIndex];
  }

protected:
#ifdef BOUNDS_CHECK
  inline void boundsCheck ( const int x, const char* msg, const char* fi, const int li ) const {
    if ( x < 0 || x >= _Nx ) {
      char errmsg[1024];
      sprintf( errmsg, msg, x, _Nx);
      throw aol::Exception ( errmsg, fi, li );
    }
  }
#endif

};

// qc::FastILexMapper<qc::QC_1D> doesn't have a lookup table, so it is actually OTF and can be used as OTFILexMapper<qc::QC_1D>.
template <>
class OTFILexMapper<qc::QC_1D> : public FastILexMapper<qc::QC_1D> {
public:
  template< class Struc1d >
  explicit OTFILexMapper ( const Struc1d &Str ) : FastILexMapper<qc::QC_1D> ( Str ) { }
};

template <>
class FastILexMapper<qc::QC_2D> {
protected:
  int      _Nx, _Ny;
  vector<int> _YLookup;
  int _nmap[4];

public:
  //! standard constructor, somewhat useless ...
  //! The class does not have an extra copy constructor or
  //! assignment operator because deep copy of all member
  //! variables ist exactly the desired behaviour.
  FastILexMapper ( ) {
    resize ( 0, 0 );
  }

  explicit FastILexMapper ( const int Nx ) {
    resize ( Nx, Nx );
  }

  FastILexMapper ( const int Nx, const int Ny ) {
    resize ( Nx, Ny );
  }

  //! Constructor taking grid for node indexing or ScalarArray<QC_2D>
  template< class Struct2d >
  explicit FastILexMapper ( const Struct2d &Str ) {
    resize ( Str.getNumX(), Str.getNumY() );
  }

  //! constructor from GridSize
  explicit FastILexMapper ( const GridSize<QC_2D> & gridSize ) {
    resize ( gridSize.getNumX(), gridSize.getNumY() );
  }

  //! Change the size of the grid
  void resize ( const int Nx, const int Ny ) {
    _Nx = Nx;
    _Ny = Ny;
    _YLookup.resize ( Ny );
    for ( int i = 0; i < Ny; ++i )
      _YLookup[i] = i * Nx;

    _nmap[0] = 0;
    _nmap[1] = 1;
    _nmap[2] = 0 + Nx;
    _nmap[3] = 1 + Nx;
  }

  template< class Struct2d >
  void resize ( const Struct2d &Str ) {
    resize ( Str.getNumX(), Str.getNumY() );
  }

  void resize ( const qc::CoordType &Coo ) {
    resize ( Coo[0], Coo[1] );
  }

  //! Compute global index from components
  inline int getGlobalIndex ( const int x, const int y ) const {
#ifdef BOUNDS_CHECK
    boundsCheck ( x, y, "qc::FastILexMapper<qc::QC_2D>::getGlobalIndex: coordinates out of bounds: %d %d %d %d %d", __FILE__, __LINE__ );
#endif
    return ( x + _YLookup[y] );
  }

  inline int getGlobalIndex ( const CoordType &Coord ) const {
    return ( getGlobalIndex ( Coord[0], Coord[1] ) );
  }

  //! get global index
  inline int operator() ( const int x, const int y ) const {
    return ( getGlobalIndex ( x, y ) );
  }

  //! get global index
  inline int operator() ( const CoordType &Coord ) const {
    return ( getGlobalIndex ( Coord[0], Coord[1] ) );
  }

   //! Split global index into components
  inline void splitGlobalIndex ( const int global, int &x, int &y ) const {
    x = global % _Nx;
    y = global / _Nx;
#ifdef BOUNDS_CHECK
    boundsCheck ( x, y, "qc::FastILexMapper<qc::QC_2D>::splitGlobalIndex: coordinates out of bounds: %d %d %d %d %d", __FILE__, __LINE__ );
#endif
  }

  //! Split global index into components
  inline void splitGlobalIndex ( const int global, CoordType &Coord ) const {
    Coord[0] = global % _Nx;
    Coord[1] = global / _Nx;
#ifdef BOUNDS_CHECK
    boundsCheck ( Coord[0], Coord[1], "qc::FastILexMapper<qc::QC_2D>::splitGlobalIndex: coordinates out of bounds: %d %d %d %d %d", __FILE__, __LINE__ );
#endif
  }


  //! \warning returning object is slower than filling one passed by reference
  inline CoordType splitGlobalIndex ( const int global ) const {
    const qc::CoordType Coord ( global % _Nx, global / _Nx, 0 );
#ifdef BOUNDS_CHECK
    boundsCheck ( Coord[0], Coord[1], "qc::OTFILexMapper<qc::QC_2D>::splitGlobalIndex: coordinates out of bounds: %d %d %d %d", __FILE__, __LINE__ );
#endif
    return ( Coord );
  }

  int localToGlobal ( const CoordType & Coord, int LocalIndex ) const {
    return getGlobalIndex ( Coord ) + _nmap[LocalIndex];
  }

protected:
#ifdef BOUNDS_CHECK
  inline void boundsCheck ( const int x, const int y, const char* msg, const char* fi, const int li ) const {
    if ( x < 0 || x >= _Nx || y < 0 || y >= _Ny || y >= static_cast<int>( _YLookup.size() ) ) {
      char errmsg[1024];
      sprintf( errmsg, msg, x, y, _Nx, _Ny, _YLookup.size() );
      throw aol::Exception ( errmsg, fi, li );
    }
  }
#endif

};


template <>
class FastILexMapper<qc::QC_3D> {
protected:
  int      _Nx, _Ny, _Nz;
  vector<int> _YLookup, _ZLookup;
  int _nmap[8];

public:
  //! standard constructor, somewhat useless ...
  FastILexMapper ( ) {
    resize ( 0, 0, 0 );
  }

  explicit FastILexMapper ( const int Nx ) : _YLookup( Nx ), _ZLookup( Nx ) {
    resize ( Nx, Nx, Nx );
  }

  FastILexMapper ( const int Nx, const int Ny, const int Nz ) : _YLookup( Ny ), _ZLookup( Nz ) {
    resize ( Nx, Ny, Nz );
  }

  //! Constructor taking grid for node indexing or ScalarArray<QC_3D>
  template< class Struc3d >
  explicit FastILexMapper ( const Struc3d &Str ) : _YLookup( Str.getNumY() ), _ZLookup( Str.getNumZ() ) {
    resize ( Str.getNumX(), Str.getNumY(), Str.getNumZ() );
  }

  //! constructor from GridSize
  explicit FastILexMapper ( const GridSize<QC_3D> & gridSize ) {
    resize ( gridSize.getNumX(), gridSize.getNumY(), gridSize.getNumZ() );
  }

public:
  template< class Struc3d >
  void resize ( const Struc3d &Str ) {
    resize ( Str.getNumX(), Str.getNumY(), Str.getNumZ() );
  }

  void resize ( const qc::CoordType &Coo ) {
    resize ( Coo[0], Coo[1], Coo[2] );
  }

  //! Change the size of the grid
  void resize ( const int Nx, const int Ny, const int Nz ) {
    _Nx = Nx;
    _Ny = Ny;
    _Nz = Nz;

    const int Nxy = Nx * Ny;

    if ( static_cast<int>( _YLookup.size() ) != ( Ny ) )
      _YLookup.resize ( Ny );

    if ( static_cast<int>( _ZLookup.size() ) != ( Nz ) )
      _ZLookup.resize ( Nz );

    for ( int i = 0; i < Ny; ++i )
      _YLookup[i] = i * Nx;

    for ( int i = 0; i < Nz; ++i )
      _ZLookup[i] = i * Nxy;

    _nmap[0] = 0;
    _nmap[1] = 1;
    _nmap[2] = 0 + Nx;
    _nmap[3] = 1 + Nx;
    _nmap[4] = 0      + Nxy;
    _nmap[5] = 1      + Nxy;
    _nmap[6] = 0 + Nx + Nxy;
    _nmap[7] = 1 + Nx + Nxy;
 }

  //! Compute global index from components
  inline int getGlobalIndex ( const int x, const int y, const int z ) const {
#ifdef BOUNDS_CHECK
    boundsCheck ( x, y, z, "qc::FastILexMapper<qc::QC_3D>::getGlobalIndex: coordinates out of bounds: %d %d %d %d %d %d %d %d", __FILE__, __LINE__ );
#endif
    return ( x + _YLookup[ y ] + _ZLookup [ z ] );
  }

  //! Compute global index from components
  inline int getGlobalIndex ( const CoordType &Coords ) const {
    return ( getGlobalIndex ( Coords[0], Coords[1], Coords[2] ) );
  }

  //! get global index
  inline int operator() ( const int x, const int y, const int z ) const {
    return ( getGlobalIndex ( x, y, z ) );
  }

  //! get global index
  inline int operator() ( const CoordType &Coord ) const {
    return ( getGlobalIndex ( Coord[0], Coord[1], Coord[2] ) );
  }

   //! Split global index into components
  inline void splitGlobalIndex ( const int global, int &x, int &y, int &z ) const {
    x = global %_Nx;
    y = ( global / _Nx ) % _Ny;
    z = global / ( _Nx * _Ny );
#ifdef BOUNDS_CHECK
    boundsCheck ( x, y, z, "qc::FastILexMapper<qc::QC_3D>::splitGlobalIndex: coordinates out of bounds: %d %d %d %d %d %d %d %d", __FILE__, __LINE__ );
#endif
  }

  //! Split global index into components
  inline void splitGlobalIndex ( const int global, CoordType &Coords ) const {
    Coords[0] = global %_Nx;
    Coords[1] = ( global / _Nx ) % _Ny;
    Coords[2] = global / ( _Nx * _Ny );
#ifdef BOUNDS_CHECK
    boundsCheck ( Coords[0], Coords[1], Coords[2], "qc::FastILexMapper<qc::QC_3D>::splitGlobalIndex: coordinates out of bounds: %d %d %d %d %d %d %d %d", __FILE__, __LINE__ );
#endif
  }

  //! \warning returning object is slower than filling one passed by reference
  inline CoordType splitGlobalIndex ( const int global ) const {
    const qc::CoordType Coord ( global %_Nx, ( global / _Nx ) % _Ny, global / ( _Nx * _Ny ) );
#ifdef BOUNDS_CHECK
    boundsCheck ( Coord[0], Coord[1], Coord[2], "qc::OTFILexMapper<qc::QC_3D>::splitGlobalIndex: coordinates out of bounds: %d %d %d %d %d %d", __FILE__, __LINE__ );
#endif
    return ( Coord );
  }

  int localToGlobal ( const CoordType & Coord, int LocalIndex ) const {
    return getGlobalIndex ( Coord ) + _nmap[LocalIndex];
  }

protected:
#ifdef BOUNDS_CHECK
  inline void boundsCheck ( const int x, const int y, const int z, const char* msg, const char* fi, const int li ) const {
    if ( x < 0 || x >= _Nx || y < 0 || y >= _Ny || y >= static_cast<int>( _YLookup.size() ) || z < 0 || z >= _Nz || z >= static_cast<int>( _ZLookup.size() ) ) {
      char errmsg[1024];
      sprintf( errmsg, msg, x, y, z, _Nx, _Ny, _Nz, _YLookup.size(), _ZLookup.size() );
      throw aol::Exception ( errmsg, fi, li );
    }
  }
#endif

};

template <Dimension Dim>
class FastQuadILexMapper {
public:

};


/****
 * \Brief Optimized version of the inverse-lexikographical indexmapper for 2D and cubic elements.
 * \author Droske
 */
template <>
class FastQuadILexMapper<qc::QC_2D> {
protected:
  const int      _depth, _width;
  struct nmap_T {
    int v[9];
  }
  ; // int[9] doesn't work with std::vector
  vector<int>    _ymap;
  vector<nmap_T> _nmap;

public:
  template< class QStruc2d >
  explicit FastQuadILexMapper ( const QStruc2d &Grid )
      : _depth ( Grid.getGridDepth() ), _width ( Grid.getWidth() ) {
    init();
  }

  inline int getGlobalIndex ( const CoordType &Coords ) const {
    return ( Coords.x() << 1 )  + _ymap[Coords.y() ];
  }

  inline int localToGlobal ( const Element &El, const int localIndex ) const {
    return getGlobalIndex ( El ) + _nmap[El.level() ].v[localIndex];
  }

protected:
  void init() {

    const int dofsPerRow = ( ( _width - 1 ) << 1 ) + 1;

    _ymap.resize ( _width );
    _nmap.resize ( _depth + 1 );

    for ( int j = 0; j < _width; j++ ) {
      _ymap[j] = 2 * dofsPerRow * j;
    }

    for ( int i = 0; i <= _depth; i++ ) {
      _nmap[i].v[0] = 0;
      _nmap[i].v[1] = 1;
      _nmap[i].v[2] = 2;

      _nmap[i].v[3] = dofsPerRow;
      _nmap[i].v[4] = dofsPerRow + 1;
      _nmap[i].v[5] = dofsPerRow + 2;

      _nmap[i].v[6] = 2 * dofsPerRow;
      _nmap[i].v[7] = 2 * dofsPerRow + 1;
      _nmap[i].v[8] = 2 * dofsPerRow + 2;
    }
  }
};



template <Dimension Dim>
class FastQuartILexMapper {
public:

};


/*
 * \Brief inverse-lexikographical indexmapper for bi-quartic elements in 2D on cubic elements.
 * \author Geihe
 * based on FastQuadILexMapper
 */
template <>
class FastQuartILexMapper<qc::QC_2D> {
protected:
  const int      _depth, _width;
  struct nmap_T {
    int v[25];
  };
  vector<int>    _ymap;
  vector<nmap_T> _nmap;

public:
  template< class QStruc2d >
  explicit FastQuartILexMapper ( const QStruc2d &Grid )
  : _depth ( Grid.getGridDepth() ), _width ( Grid.getWidth() ) {
    init();
  }

  inline int getGlobalIndex ( const CoordType &Coords ) const {
    return ( Coords.x() << 2 )  + _ymap[Coords.y() ];
  }

  inline int localToGlobal ( const Element &El, const int localIndex ) const {
    return getGlobalIndex ( El ) + _nmap[El.level() ].v[localIndex];
  }

protected:
  void init() {

    const int dofsPerRow = ( ( _width - 1 ) << 2 ) + 1;

    _ymap.resize ( _width );
    _nmap.resize ( _depth + 1 );

    for ( int j = 0; j < _width; j++ ) {
      _ymap[j] = 4 * dofsPerRow * j;
    }

    for ( int i = 0; i <= _depth; i++ ) {
      _nmap[i].v[0] = 0;
      _nmap[i].v[1] = 1;
      _nmap[i].v[2] = 2;
      _nmap[i].v[3] = 3;
      _nmap[i].v[4] = 4;

      _nmap[i].v[5] = dofsPerRow;
      _nmap[i].v[6] = dofsPerRow + 1;
      _nmap[i].v[7] = dofsPerRow + 2;
      _nmap[i].v[8] = dofsPerRow + 3;
      _nmap[i].v[9] = dofsPerRow + 4;

      _nmap[i].v[10] = 2 * dofsPerRow;
      _nmap[i].v[11] = 2 * dofsPerRow + 1;
      _nmap[i].v[12] = 2 * dofsPerRow + 2;
      _nmap[i].v[13] = 2 * dofsPerRow + 3;
      _nmap[i].v[14] = 2 * dofsPerRow + 4;

      _nmap[i].v[15] = 3 * dofsPerRow;
      _nmap[i].v[16] = 3 * dofsPerRow + 1;
      _nmap[i].v[17] = 3 * dofsPerRow + 2;
      _nmap[i].v[18] = 3 * dofsPerRow + 3;
      _nmap[i].v[19] = 3 * dofsPerRow + 4;

      _nmap[i].v[20] = 4 * dofsPerRow;
      _nmap[i].v[21] = 4 * dofsPerRow + 1;
      _nmap[i].v[22] = 4 * dofsPerRow + 2;
      _nmap[i].v[23] = 4 * dofsPerRow + 3;
      _nmap[i].v[24] = 4 * dofsPerRow + 4;
    }
  }
};



template <Dimension Dim>
class PeriodicFastILexMapper {
};

/****
 * \Brief inversely-lexikographical indexmapper for 2D with periodic boundary condition.
 */
template <>
class PeriodicFastILexMapper<QC_2D> {
protected:
  const int      _depth, _width;
  struct nmap_T {
    int v[4];
  };
  vector<int>    _ymap;
  vector<nmap_T> _nmap;

public:
  template< class QStruc2d >
  explicit PeriodicFastILexMapper ( const QStruc2d &Grid )
      : _depth ( Grid.getGridDepth() ), _width ( Grid.getWidth() ) {
    init();
  }

  PeriodicFastILexMapper ( const int Depth )
      : _depth ( Depth ), _width ( ( 1 << Depth ) + 1 ) {
    init();
  }

  inline int getGlobalIndex ( const CoordType &Coords ) const {
    return Coords.x() + _ymap[Coords.y() ];
  }

  inline int localToGlobal ( const Element &El, const int localIndex ) const {
    int _localNode = _nmap[El.level() ].v[localIndex];
    int _globalIndex = getGlobalIndex ( El );
    if ( El.x() + 1 == _width ) {
      if ( localIndex == 1 || localIndex == 3 ) {
        _globalIndex -= _width;
      }
    }
    if ( El.y() + 1 == _width ) {
      if ( localIndex == 2 || localIndex == 3 ) {
        _globalIndex -= _width * _width;
      }
    }
    return _globalIndex + _localNode;
  }

protected:
  void init() {
    _ymap.resize ( _width );
    _nmap.resize ( _depth + 1 );

    for ( int j = 0; j < _width; j++ ) {
      _ymap[j] = _width * j;
    }

    for ( int i = 0; i <= _depth; i++ ) {
      _nmap[i].v[0] = 0;
      _nmap[i].v[1] = 1;
      _nmap[i].v[2] = _width;
      _nmap[i].v[3] = _width + 1;
    }
  }
};



/****
 * \Brief inversely-lexikographical indexmapper for 3D with periodic boundary condition.
 */
template <>
class PeriodicFastILexMapper<QC_3D> {
protected:
  const int      _depth, _width;
  struct nmap_T {
    int v[8];
  };
  vector<int>    _ymap, _zmap;
  vector<nmap_T> _nmap;

public:
  template< class QStruc3d >
  explicit PeriodicFastILexMapper ( const QStruc3d &Grid )
      : _depth ( Grid.getGridDepth() ), _width ( Grid.getWidth() ) {
    init();
  }

  explicit PeriodicFastILexMapper ( const int Depth )
      : _depth ( Depth ), _width ( ( 1 << Depth ) + 1 ) {
    init();
  }


  int getGlobalIndex ( const CoordType &Coords ) const {
    return Coords.x() + _ymap[Coords.y() ] + _zmap[Coords.z() ];
  }

  int localToGlobal ( const Element &El, const int localIndex ) const {
    int _localNode = _nmap[El.level() ].v[localIndex];
    int _globalIndex = getGlobalIndex ( El );
    if ( El.x() + 1 == _width ) {
      if ( localIndex == 1 || localIndex == 3 || localIndex == 5 || localIndex == 7 ) {
        _globalIndex -= _width;
      }
    }
    if ( El.y() + 1 == _width ) {
      if ( localIndex == 2 || localIndex == 3 || localIndex == 6 || localIndex == 7 ) {
        _globalIndex -= _width * _width;
      }
    }
    if ( El.z() + 1 == _width ) {
      if ( localIndex > 3 ) {
        _globalIndex -= _width * _width * _width;
      }
    }
    return _globalIndex + _localNode;
  }

protected:
  void init() {
    _ymap.resize ( _width );
    _zmap.resize ( _width );
    _nmap.resize ( _depth + 1 );

    const int os = aol::Sqr ( _width );

    for ( int j = 0; j < _width; j++ ) {
      _ymap[j] = _width * j;
      _zmap[j] = os * j;
    }

    for ( int i = 0; i <= _depth; i++ ) {
      _nmap[i].v[0] = 0;
      _nmap[i].v[1] = 1;
      _nmap[i].v[2] = _width;
      _nmap[i].v[3] = _width + 1;
      _nmap[i].v[4] = 0 + os;
      _nmap[i].v[5] = 1 + os;
      _nmap[i].v[6] = _width + os;
      _nmap[i].v[7] = _width + 1 + os;
    }
  }
};

} // end namespace qc

#endif
