#ifndef __TRIANGMESH_H
#define __TRIANGMESH_H

#include <pointerClasses.h>
#include <meshWithData.h>
#include <vectorExtensions.h>
#include <multiVector.h>
#include <bzipiostream.h>
#include <smallMat.h>
#include <geom.h>
#include <rgbColorMap.h>

namespace aol {

/** Class for storage of triangle meshes.
 *  \todo bounds checking
 *  \author Schwen, von Deylen
 */
template< typename DataType, typename _TriangleType = Triangle<DataType> >
class TriangMesh {
public:
  typedef DataType RealType;
  typedef _TriangleType TriangleType;
  typedef std::vector< Vec3<int> > CellArray;

  //----------- iterators -------------------------------------------------------------------------------

  // --- triangle iterator ---
  class TriangleIterator {

    typedef TriangMesh<DataType, TriangleType> MeshT;

  protected:
    const MeshT &  _mesh;
    unsigned int _iter;

    mutable TriangleType _cur;
    mutable bool _filled;

  public:
    typedef TriangleIterator Self;
    typedef MeshT            BeginType;
    typedef MeshT            EndType;
    typedef TriangleType IteratedType;

    TriangleIterator ( const MeshT & mesh )
        : _mesh ( mesh ),
        _iter ( 0 ),
        _filled ( 0 ) {}

    bool operator== ( const EndType & ) const {
      return atEnd();
    }

    bool operator!= ( const EndType & ) const {
      return notAtEnd();
    }

    bool atEnd () const {
      return _iter == _mesh.getTriangs().size();
    }

    bool notAtEnd () const {
      return ! atEnd ();
    }

    Self& operator++( ) {
      ++_iter;
      _filled = false;
      return *this;
    }

    const IteratedType& operator*() const {
      if (!_filled)
        fillTriangle();
      return _cur;
    }

    const IteratedType* operator->() const {
      if (!_filled)
        fillTriangle();
      return &_cur;
    }

    int getIndex () const {
      return static_cast<int> ( _iter );
    }
    
    Vec3<int> getNodeIndices() const {
      return  ( _mesh.getTriangs() [_iter] );
    }
    
    int getNodeIndex( int localIdx ) const {
      return  ( _mesh.getTriangs() [_iter] )[localIdx];
    }

  protected:
    void fillTriangle( ) const {
      _cur.fillElement(_mesh,_iter);
      _filled = true;
    }
  };

  // --- node iterator ---
  class VertexIterator {

    typedef TriangMesh<DataType, _TriangleType> MeshT;

  protected:
    const MeshT &  _mesh;
    int            _iter;

  public:
    typedef VertexIterator Self;
    typedef MeshT          BeginType;
    typedef MeshT          EndType;
    typedef Vec3<DataType> IteratedType;

    VertexIterator ( const MeshT &mesh )
        : _mesh ( mesh ),
        _iter ( 0 ) {}

    bool operator== ( const EndType & ) const {
      return atEnd();
    }

    bool operator!= ( const EndType & ) const {
      return notAtEnd();
    }

    bool atEnd () const {
      return _iter == _mesh.getNumVertices ();
    }

    bool notAtEnd () const {
      return ! atEnd ();
    }

    Self& operator++( ) {
      ++_iter;
      return *this;
    }

    IteratedType operator*() const {
      return _mesh.getVertex ( _iter );
    }

    int getIndex () const {
      return _iter;
    }

    aol::Vec3<DataType> getCoords () const {
      return operator* ();
    }
  };
  //----------- end of iterators ------------------------------------------------------------------------

  typedef TriangleType ElementType;
  typedef TriangleIterator  ElementIteratorType;
  typedef VertexIterator    NodeIteratorType;

protected:
  MultiVector< DataType >  _vertex;  //!< coordinates of the vertices, stored as MultiVector because it is used for numerics in SurfMesh
  DeleteFlagPointer<CellArray> _triang;  //!< vertex indices of the triangles

  aol::MultiVector<DataType> _vertexData, _triangData;
  std::vector<std::string>   _vertexDataDescr, _triangDataDescr;

protected:
  mutable vector< aol::Vec3<int>      > _neighbour_;   //!< Do NOT access directly!
public:
  //! Create empty TriangMesh
  TriangMesh ( ) : _vertex ( 3, 0 ), _triang ( new CellArray, true ) {
  }

  //! Create empty TriangMesh and reserve memory for reserveNumVertex vertices and reserveNumTriang triangles
  TriangMesh ( const int reserveNumVertex, const int reserveNumTriang ) :
      _vertex ( 3, 0 ),
      _triang ( new CellArray, true ) {
    reserve ( reserveNumVertex, reserveNumTriang );
  }

  //! performs ony flat copy of geometry and topology vectors
  TriangMesh ( MultiVector<DataType> & vertex, CellArray & triang ) :
      _vertex ( vertex, FLAT_COPY ),
      _triang ( &triang, false ) {
  }

  TriangMesh ( const TriangMesh<DataType, TriangleType> & other ) :
      _vertex ( other._vertex ),
      _triang ( new CellArray ( *other._triang ), true ),
      _vertexData ( other._vertexData ),
      _triangData ( other._triangData ),
      _vertexDataDescr ( other._vertexDataDescr ),
      _triangDataDescr ( other._triangDataDescr ) {}
      
  TriangMesh ( string filename ) : _vertex ( 3, 0 ), _triang ( new CellArray, true ) {
    this->loadBasedOnSuffix(filename);
  }

  //! destructor
  virtual ~TriangMesh ( ) {
    // nothing to destroy
  }

  TriangMesh<DataType, TriangleType> & operator= ( const TriangMesh<DataType, TriangleType> & other ) {
    _vertex.reallocate ( other._vertex );
    _vertex = other._vertex;
    _triang.reset ( new CellArray ( *other._triang ), true );
    _vertexData.reallocate ( other._vertexData );
    _vertexData = other._vertexData;
    _triangData.reallocate ( other._triangData );
    _triangData = other._triangData;
    _vertexDataDescr = other._vertexDataDescr;
    _triangDataDescr = other._triangDataDescr;
    _neighbour_ = other._neighbour_;

    return *this;
  }

  // implicitly generated copy constructor (non-explicit) and assignment operator
  // do the right job, no extra implementation neccessary.
  // TriangMesh ( const TriangMesh<DataType> &other );
  // TriangMesh<DataType>& operator= ( const TriangMesh<DataType> &other );

public:
  int getNumVertices ( ) const {
    return ( _vertex[0].size() );
  }

  int getMaxNumVertices ( ) const {
    return ( _vertex[0].capacity() );
  }

  int getNumTriangs ( ) const {
    return ( static_cast<int> (  _triang->size() ) );
  }

  int getNumFaces ( ) const {
    return getNumTriangs ();
  }

  int getMaxNumTriangs ( ) const {
    return ( static_cast<int> (  _triang->capacity() ) );
  }

  //! add given number of new vertex data vectors, preserving old data
  void addVertexData ( const int numNewVertexDataVectors = 1 ) {
    for ( int i = 0; i < numNewVertexDataVectors; ++i )
      _vertexData.appendReference ( * ( new Vector<DataType> ( getNumVertices() ) ), true );
    resizeVertexDataDescription ( getNumVertexDataVectors() );
  }

  //! add new vector for vertex data with given name, preserving old data
  void addVertexData ( string vertexDataDescr ) {
    _vertexData.appendReference ( * ( new Vector<DataType> ( getNumVertices() ) ), true );
    _vertexDataDescr.push_back ( vertexDataDescr );
  }

  //! create given number of vectors for vertex data, deleting old vertex data
  void createVertexData ( const int numVertexDataVectors = 1 ) {
    resizeVertexDataDescription ( numVertexDataVectors );
    reallocateVertexData ( numVertexDataVectors );
  }

  //! create vectors for new vertex data for the given names, deleting old data
  void createVertexData ( const vector<string> & vertexDataDescr ) {
    _vertexDataDescr = vertexDataDescr;
    reallocateVertexData ( static_cast<int> ( vertexDataDescr.size() ) );
  }

  //! add given number of new vertex data vectors, preserving old data
  void addTriangData ( const int numNewTriangDataVectors = 1 ) {
    for ( int i = 0; i < numNewTriangDataVectors; ++i )
      _triangData.appendReference ( * ( new Vector<DataType> ( getNumTriangs() ) ), true );
    resizeTriangDataDescription ( getNumTriangDataVectors() );
  }

  //! add new vector for triangle data with given name, preserving old data
  void addTriangData ( string triangDataDescr ) {
    _triangData.appendReference ( * ( new Vector<DataType> ( getNumTriangs() ) ), true );
    _triangDataDescr.push_back ( triangDataDescr );
  }

  //! create given number of vectors for new triangle data for the given names, deleting old data
  void createTriangData ( const int numTriangDataVectors = 1 ) {
    resizeTriangDataDescription ( numTriangDataVectors );
    reallocateTriangData ( numTriangDataVectors );
  }

  //! create vectors for new triangle data for the given names, deleting old data
  void createTriangData ( const vector<string> & triangDataDescr ) {
    _triangDataDescr = triangDataDescr;
    reallocateTriangData ( static_cast<int> ( triangDataDescr.size() ) );
  }

  int getNumVertexDataVectors ( ) const {
    return _vertexData.numComponents();
  }

  int getNumTriangDataVectors ( ) const {
    return _triangData.numComponents();
  }

  bool hasVertexData ( ) const {
    return getNumVertexDataVectors ( ) != 0 ;
  }

  bool hasTriangData ( ) const {
    return getNumTriangDataVectors ( ) != 0 ;
  }

  //! reserve memory for given number of vertices and triangles
  void reserve ( const int reserveNumVertex, const int reserveNumTriang ) {
    for ( short i = 0; i < 3; ++i )
      _vertex[i].reserve ( reserveNumVertex );
    _triang->reserve ( reserveNumTriang );

    _vertexData.reserve ( getNumVertexDataVectors(), reserveNumVertex );
    _triangData.reserve ( getNumTriangDataVectors(), reserveNumTriang );

    _neighbour_.reserve ( reserveNumTriang );
  }

  //! change number of vertices and triangles
  void resize ( const int newNumVertex, const int newNumTriang ) {
    for ( short i = 0; i < 3; ++i )
      _vertex[i].resize ( newNumVertex );
    _triang->resize ( newNumTriang );

    _vertexData.resize ( getNumVertexDataVectors(), newNumVertex );
    _triangData.resize ( getNumTriangDataVectors(), newNumTriang );

    _neighbour_.resize ( newNumTriang );
  }

  //! change number of vertices and triangles, deleting old data
  void reallocate ( const int newNumVertex, const int newNumTriang ) {
    for ( short i = 0; i < 3; ++i )
      _vertex[i].reallocate ( newNumVertex );
    _triang->clear();
    _triang->resize ( newNumTriang );

    _vertexData.reallocate ( getNumVertexDataVectors(), newNumVertex );
    _triangData.reallocate ( getNumTriangDataVectors(), newNumTriang );

    _neighbour_.clear();
    _neighbour_.resize ( newNumTriang );
  }

  //! insert new vertex and return global index
  int pushBackVertex ( const Vec3<DataType> newVertex ) {
    for ( short i = 0; i < 3; ++i ) {
      _vertex[i].pushBack ( newVertex[i] );
    }

    for ( int i = 0; i < getNumVertexDataVectors(); ++i ) {
      _vertexData[i].pushBack ( aol::ZTrait<DataType>::zero );
    }
    
    return getNumVertices() - 1;
  }

  //! insert new triangle and return global index
  int pushBackTriang ( const Vec3<int> newTriang ) {
    _triang->push_back ( newTriang );

    for ( int i = 0; i < getNumTriangDataVectors(); ++i ) {
      _triangData[i].pushBack ( aol::ZTrait<DataType>::zero );
    }
    
    return getNumFaces() - 1;
  }

  const MultiVector<DataType> & getVertexCoords () const {
    return _vertex;
  }

  const CellArray & getTriangs () const {
    return *_triang;
  }

  DataType getVertexCoord ( const int numVertex, const int numCoord ) const {
    return _vertex[numCoord][numVertex];
  }

  Vec3<DataType> getVertex ( const int num ) const {
    return ( Vec3<DataType> ( _vertex[0][num], _vertex[1][num], _vertex[2][num] ) );
  }

  void setVertex ( const int num, const Vec3<DataType> Arg ) {
    for ( short i = 0; i < 3; ++i )
      _vertex[i][num] = Arg[i];
  }

  int getNeighbour ( const int elementID, const int acrossLocalNode ) const {
    return _neighbour_[elementID][acrossLocalNode];
  }

  void setNeighbour ( const int elementID, const int acrossLocalNode, const int value ) const {
    _neighbour_[elementID][acrossLocalNode] = value;
  }

  const aol::Vector<DataType>& getVertexData ( const int i = 0 ) const {
    return ( _vertexData[i] );
  }

  aol::Vector<DataType>& getVertexData ( const int i = 0 ) {
    return ( _vertexData[i] );
  }

  const string & getVertexDataDescr ( const int i = 0 ) const {
    return _vertexDataDescr[i];
  }

  string & getVertexDataDescr ( const int i = 0 ) {
    return _vertexDataDescr[i];
  }

  const Vec3<int> & getTriang ( const int num ) const {
    return ( ( *_triang ) [num] );
  }

  Vec3<int> & getTriang ( const int num ) {
    return ( ( *_triang ) [num] );
  }

  void setTriang ( const int num, const Vec3<int> Arg ) {
    ( *_triang ) [num] = Arg;
  }

  int getTriangNodeIdx ( const int num, const int localNode ) const {
    return (*_triang)[num][localNode];
  }

  void setTriangNodeIdx ( const int num, const int localNode, const int value ) {
    (*_triang)[num][localNode] = value;
  }

  const aol::Vector<DataType>& getTriangData ( const int i = 0 ) const {
    return ( _triangData[i] );
  }

  aol::Vector<DataType>& getTriangData ( const int i = 0 ) {
    return ( _triangData[i] );
  }

  const string & getTriangDataDescr ( const int i = 0 ) const {
    return _triangDataDescr[i];
  }

  string & getTriangDataDescr ( const int i = 0 ) {
    return _triangDataDescr[i];
  }

  //! save in UD ply format
  void saveAsUDPLY ( string filename ) const {
    saveAsStanfordOrUDPLY ( filename, UD_PLY );
  }

  //! save in Stanford ply format
  void saveAsPLY ( string filename, const bool Binary = false ) const {
    saveAsStanfordOrUDPLY ( filename, STANFORD_PLY, -1, Binary );
  }

  //! load from file in UD ply format
  void loadFromUDPLY ( string filename ) {
    loadFromStanfordOrUDPLY ( filename, UD_PLY );
  }

  //! load from file in Stanford ply format
  void loadFromPLY ( string filename ) {
    loadFromStanfordOrUDPLY ( filename, STANFORD_PLY );
  }

  //! load from file while trying to figure out the format from the suffix of the file name
  void loadBasedOnSuffix ( string filename ) {
    if ( fileNameEndsWith ( filename.c_str(), ".ply" ) )
      loadFromPLY ( filename );
    else if ( fileNameEndsWith ( filename.c_str(), ".udply" ) )
      loadFromUDPLY ( filename );
    else if ( fileNameEndsWith ( filename.c_str(), ".obj" ) )
      loadFromOBJ ( filename );
    else if ( fileNameEndsWith ( filename.c_str(), ".srf" ) )
      loadFromSRF ( filename );
    else if ( fileNameEndsWith ( filename.c_str(), ".stl" ) )
      loadFromSTL ( filename );
    else
      throw aol::Exception ( aol::strprintf( "aol::TriangMesh<DataType,TriangleType>::loadBasedOnSuffix(): unknown suffix (filename = %s)", filename.c_str() ).c_str(), __FILE__, __LINE__ );
  }

  void setZero ( ) {
    throw aol::Exception ( "aol::TriangMesh<DataType,TriangleType>::setZero() must not be called and exists for compatibility only", __FILE__, __LINE__ );
  }

  void clear ( ) {
    resize ( 0, 0 );
    _vertexData.reallocate ( 0, 0 );
    _triangData.reallocate ( 0, 0 );
  }

  Vec3<RealType> centerOfBoundingBox ( ) const {
    aol::Vec3<RealType> minXYZ, maxXYZ;
    computeBoundingBox ( minXYZ, maxXYZ );
    maxXYZ += minXYZ;
    maxXYZ *= .5;
    return maxXYZ;
  }

  //! scales the mesh by multiplying each coordinate with the respective factor
  void scaleSizeByFactor ( const aol::Vec3<DataType> &ScaleFactorXYZDirection ) {
    for ( int j = 0; j < 3; j++ ) {
      _vertex[j] *= ScaleFactorXYZDirection[j];
    }
  }

  //! scales the mesh by multiplying each coordinate with the same factor
  void scaleSizeByFactor ( const DataType Factor ) {
    for ( int j = 0; j < 3; j++ ) {
      _vertex[j] *= Factor;
    }
  }

  //! shifts the mesh by adding the respective offset to each coordinate
  void shiftByOffset ( const aol::Vec<3, DataType> &Offset ) {
    for ( int j = 0; j < 3; j++ )
      _vertex[j].addToAll ( Offset[j] );
  }

  template <typename ParametricDeformationType>
  void transformParametric ( const ParametricDeformationType &ParDef ) {
    Vec3<DataType> oldVertex, newVertex;
    for ( int i = 0; i < getNumVertices(); ++i ) {
      for ( short j = 0; j < 3; ++j )
        oldVertex[j] = _vertex[j][i];

      ParDef.evaluateDeformationOn01 ( ParDef, oldVertex, newVertex );

        for ( short j = 0; j < 3; ++j )
          _vertex[j][i] = newVertex[j];
    }
  }

  void clearNeighbours ( ) {
    _neighbour_.clear();
  }

protected:
  void reallocateVertexData ( const int numVertexDataVectors ) {
    _vertexData.reallocate ( numVertexDataVectors, 0 );
    for ( int i = 0; i < numVertexDataVectors; ++i ) {
      _vertexData[i].reserve ( getMaxNumVertices() );
      _vertexData[i].resize ( getNumVertices() );
    }
  }

  void reallocateTriangData ( const int numTriangDataVectors ) {
    _triangData.reallocate ( numTriangDataVectors, 0 );
    for ( int i = 0; i < numTriangDataVectors; ++i ) {
      _triangData[i].reserve ( getMaxNumTriangs() );
      _triangData[i].resize ( getNumTriangs() );
    }
  }

  void resizeVertexDataDescription ( const int numVertexDataVectors ) {
    if ( numVertexDataVectors > getNumVertexDataVectors() )
      for ( int i = getNumVertexDataVectors(); i < numVertexDataVectors; ++i )
        _vertexDataDescr.push_back ( strprintf ( "generic_scalar_values%02i", i ) );
    else
      _vertexDataDescr.resize ( numVertexDataVectors );
  }

  void resizeTriangDataDescription ( const int numTriangDataVectors ) {
    if ( numTriangDataVectors > getNumTriangDataVectors() )
      for ( int i = getNumTriangDataVectors(); i < numTriangDataVectors; ++i )
        _triangDataDescr.push_back ( strprintf ( "generic_scalar_values%02i", i ) );
    else
      _triangDataDescr.resize ( numTriangDataVectors );
  }

public:
  //-------------------------------------------------------
  //				member needed by FEOps!!!
  //-------------------------------------------------------
  
  virtual bool isAdaptive() const {
    return false;
  }  

  virtual int checkForHangingNode ( const ElementType &, int ) const {
    return -1;
  }
  
  int getGridDepth() const {
     throw aol::Exception ( "TriangMesh does not support getGridDepth()!", __FILE__, __LINE__ );
  }

  inline int getNumX() const {
    return 0;
  }

  void saveAsStanfordOrUDPLY ( string filename, aol::PLY_FORMAT format, const int Precision = -1, const bool Binary = false ) const {

  MeshWithData<TriangMesh<DataType,_TriangleType> > meshSaver ( *this );
  for ( int i = 0; i < getNumVertexDataVectors(); ++i )
    meshSaver.addData ( getVertexData ( i ), getVertexDataDescr ( i ), VERTEX_DATA );
  for ( int i = 0; i < getNumTriangDataVectors(); ++i )
    meshSaver.addData ( getTriangData ( i ), getTriangDataDescr ( i ), FACE_DATA );

  if ( Precision != -1 )
    meshSaver.setPrecisionTo ( Precision );

  if ( format == UD_PLY )
    meshSaver.saveAsUDPLY ( filename );
  else
    meshSaver.saveAsPLY ( filename, Binary );
}

  //! save in legacy VTK format
void saveAsLegacyVTK ( string filename ) const {

  MeshWithData<TriangMesh<DataType,_TriangleType> > meshSaver ( *this );
  for ( int i = 0; i < getNumVertexDataVectors(); ++i )
    meshSaver.addData ( getVertexData ( i ), getVertexDataDescr ( i ), VERTEX_DATA );
  for ( int i = 0; i < getNumTriangDataVectors(); ++i )
    meshSaver.addData ( getTriangData ( i ), getTriangDataDescr ( i ), FACE_DATA );

  meshSaver.saveAsLegacyVTK ( filename );
}

  //! save in VRML 2.0 format
  //! Note: .wrl file format can only store either ONE vertex color or ONE face color (which can be chosen by 'numOfDataVector' ).
void saveAsWRL ( string filename, int numOfDataVector = 0, colorTrans Colormap = WHITE_TO_BLUE ) const {

  if ( this->hasVertexData() && this->hasTriangData() )
    throw ( aol::Exception ( "aol::TriangMesh<DataType,_TriangleType>::saveAsWRL: .wrl can only save one data type!", __FILE__, __LINE__ ) );
  
  MeshWithData<TriangMesh<DataType,_TriangleType> > meshSaver ( *this );

  aol::MultiVector<DataType> colors;

  if ( this->hasVertexData() ){
    // create a color map for your scalar data
    aol::RGBColorMap<DataType> colormap ( getVertexData( numOfDataVector ).getMinValue() , getVertexData( numOfDataVector ).getMaxValue() , Colormap );
    // compute colors from scalar value
    colors.resize( 3, this->getNumVertices() );
    for( int i = 0; i < getVertexData( numOfDataVector ).size();++i ){
      aol::Vec3<RealType> colorRGB;
      colormap.scalarToColor( getVertexData(numOfDataVector)[i], colorRGB );
      for( int j = 0; j < 3; ++j )
        colors[j][i] = colorRGB[j];
    }
    meshSaver.addData ( colors, getVertexDataDescr ( numOfDataVector ), VERTEX_DATA );    
  }

  if ( this->hasTriangData() ){
    // create a color map for your scalar data
    aol::RGBColorMap<DataType> colormap ( getTriangData( numOfDataVector ).getMinValue() , getTriangData( numOfDataVector ).getMaxValue() , Colormap );
    // compute colors from scalar value
    colors.resize( 3, this->getNumTriangs() );
    for( int i = 0; i < getTriangData( numOfDataVector ).size();++i ){
      aol::Vec3<RealType> colorRGB;
      colormap.scalarToColor( getTriangData(numOfDataVector)[i], colorRGB );
      for( int j = 0; j < 3; ++j )
        colors[j][i] = colorRGB[j];
    }
    meshSaver.addData ( colors, getTriangDataDescr ( numOfDataVector ), FACE_DATA );    
  }

  meshSaver.saveAsWRL ( filename );
}

private:
  enum PLY_DATA_FORMAT { ASCII, BINARY_LITTLE_ENDIAN };

  template <typename EntryType>
  EntryType readPLYDataEntry ( std::istream &In, const PLY_DATA_FORMAT DataFormat ) const {
    switch ( DataFormat ) {
      case ASCII:
        {
          EntryType tmp;
          In >> tmp;
          return tmp;
        }
        break;
      case BINARY_LITTLE_ENDIAN:
        return aol::readBinaryData<EntryType, EntryType> ( In );
      default:
        throw aol::Exception ( "Unsupported data format!", __FILE__, __LINE__ );
        break;
    }
  }

public:

// \todo do some checking whether header contains expected data
void loadFromStanfordOrUDPLY ( string filename, aol::PLY_FORMAT format ) {

  aol::Bzipifstream file ( filename.c_str() );

  if ( !file )
    throw ( aol::Exception ( "aol::TriangMesh<DataType,_TriangleType>::loadFromStanfordOrUDPLY: could not open file for reading", __FILE__, __LINE__ ) );

  int numVertex = 0, numTriang = 0;
  string line;

  getline ( file, line );

  if ( line.find ( "ply" ) == string::npos  )
    throw aol::Exception ( "aol::TriangMesh<DataType,_TriangleType>::loadFromStanfordOrUDPLY: file not in ply format", __FILE__, __LINE__ );

  // Older versions of the read routine always assumed ASCII data, so use this as default.
  PLY_DATA_FORMAT dataFormat = ASCII;
  if ( format == STANFORD_PLY ) {
    while ( line.find ( "format " ) == string::npos )
      getline ( file, line );

    line.erase ( line.find ( "format " ), strlen ( "format " ) );
    if ( line.compare ( "ascii 1.0" ) == 0 )
      dataFormat = ASCII;
    else if ( line.compare ( "binary_little_endian 1.0" ) == 0 )
      dataFormat = BINARY_LITTLE_ENDIAN;
  }

  while ( line.find ( "element vertex" ) == string::npos )
    getline ( file, line );

  line.erase ( line.find ( "element vertex" ), strlen ( "element vertex" ) );
  numVertex = atoi ( line.c_str() );

  while ( line.find ( "property float" ) == string::npos )
    getline ( file, line );
  getline ( file, line ); // we know there are at least three lines: property float x, y, z.
  getline ( file, line );

  vector<string> vertexDataDescr;
  getline ( file, line );
  while ( line.find ( "property float" ) != string::npos ) {
    vertexDataDescr.push_back ( line.substr ( line.find ( "property float" ) + 15 ) );
    getline ( file, line );
  }
  createVertexData ( vertexDataDescr );

  while ( line.find ( "element face" ) == string::npos )
    getline ( file, line );

  line.erase ( line.find ( "element face" ), strlen ( "element face" ) );
  numTriang = atoi ( line.c_str() );

  while ( line.find ( "property list" ) == string::npos )
    getline ( file, line );

  vector<string> triangDataDescr;
  getline ( file, line );
  while ( line.find ( "property float" ) != string::npos ) {
    triangDataDescr.push_back ( line.substr ( line.find ( "property float" ) + 15 ) );
    getline ( file, line );
  }
  createTriangData ( triangDataDescr );

  while ( line.find ( "end_header" ) == string::npos )
    getline ( file, line );

  resize ( numVertex, numTriang );

  for ( int i = 0; i < numVertex; i++ ) {
    _vertex[0][i] = readPLYDataEntry<float> ( file, dataFormat );
    _vertex[1][i] = readPLYDataEntry<float> ( file, dataFormat );
    _vertex[2][i] = readPLYDataEntry<float> ( file, dataFormat );
    for ( int j = 0; j < getNumVertexDataVectors(); ++j )
      getVertexData ( j ) [i] = readPLYDataEntry<float> ( file, dataFormat );
  }

  int vertexCount;
  for ( int i = 0; i < numTriang; i++ ) {
    if ( format == STANFORD_PLY ) {
      // The vertex count is usually stored as unsigned char and needs special treatment.
      if ( dataFormat == ASCII )
        file >> vertexCount;
      else
        vertexCount = readPLYDataEntry<unsigned char> ( file, dataFormat );
      if ( vertexCount != 3 ) {
        throw aol::FileFormatException ( "TriangMesh::loadFromStanfordOrUDPLY(): file contains polygons with more that 3 vertices (not a Stanford PLY file?).", __FILE__, __LINE__ );
      }
    }
    ( *_triang ) [i][0] = readPLYDataEntry<int> ( file, dataFormat );
    ( *_triang ) [i][1] = readPLYDataEntry<int> ( file, dataFormat );
    ( *_triang ) [i][2] = readPLYDataEntry<int> ( file, dataFormat );
    for ( int j = 0; j < getNumTriangDataVectors(); ++j )
      getTriangData ( j ) [i] = readPLYDataEntry<float> ( file, dataFormat );
  }

  if ( dataFormat == ASCII ) {
    while ( file.eof() == false ) {
      char ch;
      file.get ( ch );
      // 9 is tab, 10 is LF (line feed), 13 is CR (carriage return), 32 is space
      if ( ( ch != 0 ) && ( ch != 9 ) && ( ch != 10 ) && ( ch != 13 ) && ( ch != 32 ) ) {
        throw aol::FileFormatException ( aol::strprintf ( "TriangMesh::loadFromStanfordOrUDPLY(): file \"%s\" contains too much data (reading a Stanford PLY file as UDPLY?).", filename.c_str() ) , __FILE__, __LINE__ );
      }
    }
  }
}

//! load from file in the Wavefront .obj file format. Currently only supports 'v' and 'f' lines, nothing else.
void loadFromOBJ ( string filename, bool mode2D = false ) {
  ifstream in ( filename.c_str() );
  int numOfCoords = mode2D ? 2 : 3;
  do
  {
    string temp;
    in >> temp;
    // skip line of comments
    if ( temp.compare ( "#" ) == 0 ){
      char line[256];
      in.getline ( line, 256 );
      continue;
    }
    // otherwise expect data
    if ( temp.compare ( "v" ) == 0 ) {
      aol::Vec3<DataType> vertex;
      for ( int i = 0; i < numOfCoords; ++i )
        in >> vertex[i];
      pushBackVertex ( vertex );
    }
    else if ( temp.compare ( "f" ) == 0 ) {
      aol::Vec3<int> triang;
      for ( int i = 0; i < 3; ++i ) {
        in >> triang[i];
        --triang[i];
      }
      pushBackTriang ( triang );
    }
    else if ( temp.length() > 0 )
      throw aol::FileFormatException ( "TriangMesh::loadFromOBJ(): Error parsing the input file.", __FILE__, __LINE__ );
  }
  while ( in.eof() == false );
}

//! load from file in the .vtk file format. Currently only loads geometric information.
void loadFromLegacyVTK( string filename ) {
  
  ifstream vtkFile ( filename.c_str() );  
  bool readVertices = false;
  bool readFaces = false;
    
  while ( !vtkFile.eof() ) {
    char line[256];
    vtkFile.getline ( line, 256 );
	
    // first: expect vertices
    if( !strncmp ( line, "POINTS ", strlen ( "POINTS " ) ) ){
      readVertices = true; 
      continue;
    }
	
    // second: expect triangles (i.e. starting with "3 " )
    if( !strncmp ( line, "POLYGONS ", strlen ( "POLYGONS " ) ) ){
      readVertices = false;
      readFaces = true;
      continue;      
    }

    // geometric information ends with the first line that does not start with "3 " anymore
    if( readFaces && strncmp ( line, "3 ", strlen ( "3 " ) ) )
      break;

    // read in the vertex coordinates and add it to the mesh
    if( readVertices ){
      aol::Vec3<RealType> vertex;            
      char * substring = strtok ( line, " " );
      for ( int i = 0; i < 3; i++ ){
        vertex[i] = atof( substring );
        substring = strtok (NULL, " ");
      }
      pushBackVertex ( vertex );
    }
         
    // read in the face and add it to the mesh
    if( readFaces ){
      aol::Vec3<int> triangle;          
      char * substring = strtok ( line, " " );
      for ( int i = 0; i < 3; i++ ){
        substring = strtok (NULL, " ");
        triangle[i] = atof( substring );	    
      }
      pushBackTriang( triangle );
    }
  }
  vtkFile.close();
}

//! load from file in the BrainVoyager SRF file format. Currently only supports version 4.
//! \author Berkels
void loadFromSRF ( string Filename ) {
  std::ifstream in ( Filename.c_str(), std::ios::binary );

  if ( aol::readBinaryData<float, float> ( in ) != 4 )
    throw aol::FileFormatException ( "Only SRF version 4 is supported!", __FILE__, __LINE__ );

  // The next int is reserved (should be 0 according to the documentation,
  // but seems to be 1 for the files I tested).
  aol::readBinaryData<int32_t, int> ( in );

  const int numVert = aol::readBinaryData<int32_t, int> ( in );
  const int numTriangs = aol::readBinaryData<int32_t, int> ( in );
  cerr << "numVert = " << numVert << endl;
  cerr << "numTriangs = " << numTriangs << endl;
  cerr << "MeshCenter = [ ";
  for ( int i = 0; i < 3; ++i )
    cerr << aol::readBinaryData<float, float> ( in ) << " ";
  cerr << "]\n";

  aol::MultiVector<DataType> vertices ( 3, numVert );
  for ( int i = 0; i < 3; ++i )
    aol::readBinaryData<float, DataType> ( in, vertices[i].getData(), numVert );

  aol::MultiVector<DataType> normals ( 3, numVert );
  for ( int i = 0; i < 3; ++i )
    aol::readBinaryData<float, DataType> ( in, normals[i].getData(), numVert );

  cerr << "Convex curvature color = RGBA [ ";
  for ( int i = 0; i < 4; ++i )
    cerr << aol::readBinaryData<float, float> ( in ) << " ";
  cerr << "]\n";
  cerr << "Concave curvature color = RGBA [ ";
  for ( int i = 0; i < 4; ++i )
    cerr << aol::readBinaryData<float, float> ( in ) << " ";
  cerr << "]\n";

  aol::Vector<DataType> meshColor ( numVert );
  aol::readBinaryData<float, DataType> ( in, meshColor.getData(), numVert );

  for ( int i = 0; i < numVert; ++i ) {
    const int numNearestNeighbors = aol::readBinaryData<int32_t, int> ( in );
    aol::Vector<int> neighbors ( numNearestNeighbors );
    aol::readBinaryData<int32_t, int> ( in, neighbors.getData(), numNearestNeighbors );
  }

  aol::RandomAccessContainer<aol::Vec3<int> > triangs ( numTriangs );
  for ( int i = 0; i < numTriangs; ++i )
    aol::readBinaryData<int32_t, int> ( in, triangs[i].getData(), 3 );

  if ( aol::readBinaryData<int32_t, int> ( in ) != 0 )
    cerr << "Warning: Reading triangle strip elements is not supported.\n" << endl;

  if ( aol::readBinaryData<char, char> ( in ) != 0 )
    cerr << "Warning: MTC files are not supported.\n" << endl;

  // If we read till the end of the file, but not beyond it, the stream will
  // still be good, but peeking the next char will fail.
  if ( in.good() )
  {
    in.peek();
    if ( !in.eof() )
      cerr << "Warning: File not fully parsed.\n";
  }
  else
    cerr << "Warning: Error parsing file.\n";

  reserve ( numVert, numTriangs );

  for ( int i = 0; i < numVert; ++i )
    pushBackVertex ( aol::Vec3<DataType> ( vertices[0][i], vertices[1][i], vertices[2][i] ) );

  for ( int i = 0; i < numVert; ++i )
    pushBackTriang ( triangs[i] );
}

  //! load from file in the STL (STereoLithography) file format.
  //! \note Does not merge duplicated vertices.
  //! \author Berkels
  void loadFromSTL ( string Filename ) {
    std::ifstream in ( Filename.c_str(), std::ios::binary );

    if ( in.good() == false )
      throw aol::FileException ( aol::strprintf ( "Cannot open file %s for reading", Filename.c_str() ).c_str(), __FILE__, __LINE__ );

    // Skip the 80 byte header.
    in.seekg ( 80 );
    const int numTriangs = aol::readBinaryData<uint32_t, int> ( in );

    for ( int i = 0; i < numTriangs; ++i ) {
      for ( int j = 0; j < 3; ++j )
        aol::readBinaryData<float, DataType> ( in );

      Vec3<DataType> newVertex;
      Vec3<int> newTriang;
      for ( int j = 0; j < 3; ++j ) {
        for ( int k = 0; k < 3; ++k )
          newVertex[k] = aol::readBinaryData<float, RealType> ( in );
        newTriang[j] = pushBackVertex ( newVertex );
      }
      pushBackTriang ( newTriang );
      aol::readBinaryData<int16_t, int> ( in );
    }
  }

  //! write triangle mesh in format readable by povray (without vertex or triangle data)
void saveAsPov ( string filename ) const {
  aol::Bzipofstream out ( filename.c_str() );
  out << "mesh2 {" << endl
  << "vertex_vectors {" << endl
  << getNumVertices();

  for ( int i = 0; i < getNumVertices(); ++i ) {
    out << ( i == 0 ? "" : "," ) << endl
    << "  < " << _vertex[0][i] << ", " << _vertex[1][i] << ", " << _vertex[2][i] << " >";
  }
  out << endl << "}" << endl

  << "face_indices {" << endl
  << getNumTriangs();

  for ( int i = 0; i < getNumTriangs(); ++i ) {
    out << ( i == 0 ? "" : "," ) << endl
    << "  < " << ( *_triang ) [i][0] << ", " << ( *_triang ) [i][1] << ", " << ( *_triang ) [i][2] << " >";
  }
  out << endl
  << "}" << endl
  << "}" << endl;
}

// save as OBJ
void saveAsOBJ ( string filename, int precision = 8 ) const {

  MeshWithData<TriangMesh<DataType,_TriangleType> > meshSaver ( *this );
  for ( int i = 0; i < getNumVertexDataVectors(); ++i )
    meshSaver.addData ( getVertexData ( i ), getVertexDataDescr ( i ), VERTEX_DATA );
  if( getNumTriangDataVectors() > 0 )
    cerr << "WARNING in aol::TriangMesh::saveAsOBJ(): obj format does not support any face data!" << endl;

  meshSaver.saveAsOBJ ( filename, precision );
}

// comparison of two TriangMeshes: tolerance refers to coordinates and data; triangulations may not have different indexing
bool isApproxEqual ( const TriangMesh<DataType, _TriangleType> &other, const DataType tolerance = 1e-6 ) const {
  if ( ( _vertexData.compareDim (other._vertexData) ) &&
       ( _triangData.compareDim (other._triangData) ) &&
       ( *_triang == *other._triang ) ) {

    aol::MultiVector<DataType> vertexDiff ( this->_vertex );
    vertexDiff -= other._vertex;
    if ( vertexDiff.norm() > tolerance )
      return ( false );

    aol::MultiVector<DataType> vertexDataDiff ( this->_vertexData );
    vertexDataDiff -= other._vertexData;
    if ( vertexDataDiff.norm() > tolerance )
      return ( false );

    aol::MultiVector<DataType> triangDataDiff ( this->_triangData );
    triangDataDiff -= other._triangData;
    if ( triangDataDiff.norm() > tolerance )
      return ( false );

    return ( true );

  } else {
    return ( false );
  }
}


  //! paste other TriangMesh into this one (note: treatment of vertex and triangle data need to be refined)
void pasteFrom ( const TriangMesh<DataType,_TriangleType> &other ) {
  const int
    oldNumVertices  = this->getNumVertices(),
    addlNumVertices = other.getNumVertices(),
    oldNumTriangs   = this->getNumTriangs(),
    addlNumTriangs  = other.getNumTriangs();

  this->reserve ( oldNumVertices + addlNumVertices, oldNumTriangs + addlNumTriangs );

  // copy vertices
  for ( int i = 0; i < addlNumVertices; ++i ) {
    this->pushBackVertex ( other.getVertex ( i ) );
  }

  // copy vertex data
  for ( int j = 0; j < this->getNumVertexDataVectors(); ++j ) {
    for ( int i = 0; i < addlNumVertices; ++i ) {
      ( this->getVertexData ( j ) ) [ oldNumVertices + i ] = ( other.getVertexData ( j ) ) [i];
    }
  }

  // copy triangles
  aol::Vec3<int> offsetTriang ( oldNumVertices, oldNumVertices, oldNumVertices );
  for ( int i = 0; i < addlNumTriangs; ++i ) {
    this->pushBackTriang ( other.getTriang ( i ) + offsetTriang );
  }

  // copy triang data
  for ( int j = 0; j < this->getNumTriangDataVectors(); ++j ) {
    for ( int i = 0; i < addlNumTriangs; ++i ) {
      ( this->getTriangData ( j ) ) [ oldNumTriangs + i ] = ( other.getTriangData ( j ) ) [i];
    }
  }

}

  //! returns the maximum width of the bounding box and minimum x,y,z-coordinates in MinXYZ and the analogue for MaxXYZ
DataType computeBoundingBox ( aol::Vec3<DataType> &MinXYZ, aol::Vec3<DataType> &MaxXYZ ) const {
  // compute minimum and maximum coordinates
  for ( int j = 0; j < 3; j++ ) {
    MinXYZ[j] = _vertex[j].getMinValue();
    MaxXYZ[j] = _vertex[j].getMaxValue();
  }
  // compute maximum width
  return ( ( MaxXYZ - MinXYZ ).getMaxValue() );
}

  //! transform vertices by a linear mapping
void transformLinear ( const Matrix33<DataType> & transfMatrix ) {
  Vec3<DataType> oldVertex, newVertex;
  for ( int i = 0; i < getNumVertices(); ++i ) {
    for ( short j = 0; j < 3; ++j ) {
      oldVertex[j] = _vertex[j][i];
    }
    newVertex = transfMatrix * oldVertex;
    for ( short j = 0; j < 3; ++j ) {
      _vertex[j][i] = newVertex[j];
    }
  }
}

  //! transform vertices by a mapping given in homogeneous coordinates
void transformHomogeneous ( const Matrix44<DataType> & transfMatrix ) {
  Vec3<DataType> oldVertex, newVertex;
  for ( int i = 0; i < getNumVertices(); ++i ) {
    for ( short j = 0; j < 3; ++j ) {
      oldVertex[j] = _vertex[j][i];
    }
    aol::TransformCartesianCoordinatesByHomogeneousMapping ( oldVertex, transfMatrix, newVertex );
    for ( short j = 0; j < 3; ++j ) {
      _vertex[j][i] = newVertex[j];
    }
  }
}


void makeNeighbour() const {
  
  int  *n_to_t, *num_of_n, *start_of_n;
  
  int i, j, k, l, n, v, node1, node2, node3, noe, nop;

  nop = getNumVertices();
  noe = getNumTriangs();
  
  if ( int(_neighbour_.size()) != noe )
    _neighbour_.resize ( noe );

  n_to_t     = new int[noe*3];
  num_of_n   = new int[nop];
  start_of_n = new int[nop+1];

  /* Wieviel Referenzen pro Knoten initialisieren */
  // iterate over all vertices and set num_of_n to zero, i.e. initialize
  for ( i = 0; i < nop; i++ ) num_of_n[i] = 0;

  /* Wieviel Referenzen pro Knoten zaehlen, nb initialisieren */
  // iterate over all triangles
  for ( i = 0; i < noe; i++ )
  // iterate over all neighbours, i.e. 3, since every triangle has 3 neighbours for closed meshes!!!
    for ( j = 0; j < 3; j++ )  /*  tr->vertex[i][j] ist der Knoten */
    {
    //
      num_of_n[getTriangNodeIdx ( i,j ) ]++;
      //nb[3*i+j] = -1;
      setNeighbour ( i, j, -1 );
    }
    
  /* Startindex fuer Knoten bestimmen */
  i = 0; start_of_n[i++] = 0;
  while ( i < nop ) {
    start_of_n[i] = start_of_n[i-1] + num_of_n[i-1];
    i++;
  }
  start_of_n[nop] = 3 * noe;

  /* Referenzen initialisieren */
  for ( i = 0; i < noe*3; i++ ) n_to_t[i] = -1;

  /* Referenzen erzeugen */
  for ( i = 0; i < noe; i++ )
    for ( j = 0; j < 3; j++ ) {
      k = start_of_n[getTriangNodeIdx ( i,j ) ];
      while ( n_to_t[k] > -1 ) k++;
      n_to_t[k] = i;
    }
  /* Nachbarn suchen */
  for ( i = 0; i < noe; i++ )
    for ( j = 0; j < 3; j++ ) {
      node1 = getTriangNodeIdx ( i, j        ); //vx[i][j];
      node2 = getTriangNodeIdx ( i, ( j + 1 ) % 3 ); //vx[i][(j+1) % 3];
      node3 = getTriangNodeIdx ( i, ( j + 2 ) % 3 ); //vx[i][(j+2) % 3];

      for ( k = start_of_n[node1];k < start_of_n[node1+1]; k++ ) {
        n = n_to_t[k]; /* Tetraeder */
        if ( i < n ) /* Nachbarschaft nur einmal setzen */
          for ( l = 0; l < 3; l++ ) {
            if ( node3 == getTriangNodeIdx ( n, l ) ) {
              //nb[3*i+((j+1) % 3)] = n;
              setNeighbour ( i, ( j + 1 ) % 3, n );
              for ( v = 0; v < 3; v++ )
                if ( v != l && getTriangNodeIdx ( n, v ) != node1 )
                  setNeighbour ( n, v, i );
              //nb[3*n+v] = i;
            }
            else {
              if ( node2 == getTriangNodeIdx ( n, l ) ) {
                //nb[3*i+((j+2) % 3)] = n;
                setNeighbour ( i, ( j + 2 ) % 3, n );
                for ( v = 0; v < 3; v++ )
                  if ( v != l && getTriangNodeIdx ( n, v ) != node1 )
                    setNeighbour ( n, v, i );
                //nb[3*n+v] = i;
              }
            }
          }
      }
   }

  delete[] n_to_t;
  delete[] num_of_n;
  delete[] start_of_n;
}



void makeOrientationConsistent() {
  // need neighbors
  if ( int(_neighbour_.size()) != this->getNumTriangs() )
    makeNeighbour();
  
  // true for all triangles T, whose neighboring triangles have already been oriented consistent with T
  aol::BitVector alreadyHandled( this->getNumTriangs() );
  // true for all triangles who have already been dealt with or who are already waiting in queue
  aol::BitVector inQueue( this->getNumTriangs() );
  // contains all triangles which are already oriented and whose neighbors will be dealt with next (may contain triangles twice for simplicity)
  std::queue<int> activeTriangles;
  activeTriangles.push( 0 );
  inQueue.set( 0, true );
  // the triangle whose neighbors are currently handled
  int currentTriangle;
  // while there are triangles left whose neighbors are not consistently oriented...
  while( !activeTriangles.empty() ){
    currentTriangle = activeTriangles.front();
    activeTriangles.pop();
    // deal with all three neighbors of currentTriangle, i.e. orient them and add them to the list to deal with their neighbors
    for ( int i = 0; i < 3; i++ ){
      int neighbor = getNeighbour( currentTriangle, i );
      if ( neighbor >= 0 && neighbor < this->getNumTriangs() && !alreadyHandled[neighbor] ){
        // compute the nodes "currentTriangle" and "neighbor" have in common
        int node1 = getTriangNodeIdx ( currentTriangle, ( i + 1 ) % 3 );
        int node2 = getTriangNodeIdx ( currentTriangle, ( i + 2 ) % 3 );
        // check whether common nodes occur in reversed order in "neighbor", if not, change order
        int j = 0;
        while ( getTriangNodeIdx ( neighbor, j ) != node2 )
          j++;
        if ( getTriangNodeIdx ( neighbor, ( j + 1 ) % 3 ) != node1 ){
          // change order of nodes
          int exchangeCache = getTriangNodeIdx ( neighbor, ( j + 1 ) % 3 );
          setTriangNodeIdx( neighbor, ( j + 1 ) % 3, getTriangNodeIdx ( neighbor, ( j + 2 ) % 3 ) );
          setTriangNodeIdx( neighbor, ( j + 2 ) % 3, exchangeCache );
          // change order of corresponding neighbours
          exchangeCache = getNeighbour( neighbor, ( j + 1 ) % 3 );
          setNeighbour( neighbor, ( j + 1 ) % 3, getNeighbour( neighbor, ( j + 2 ) % 3 ) );
          setNeighbour( neighbor, ( j + 2 ) % 3, exchangeCache );
        }
        if ( !inQueue[neighbor] ){
          activeTriangles.push( neighbor );
          inQueue.set( neighbor, true );
        }
      }
    }
    alreadyHandled.set( currentTriangle, true );
  }
}

  // end of class TriangMesh<DataType>
};


/** \brief Class for storage of triangle meshes that can be used for FE simulations.
 *  \author Heeren, Perl
 * 
 *  \note Equivalent to aol::TriangMesh<> but provides a function to fill boundary mask with boundary vertices
 */
template< typename DataType , typename TriangleType = TriangBaseElement< DataType > >
class FETriangMesh : public TriangMesh< DataType , TriangleType  > {
  
public:
  FETriangMesh ( ) : TriangMesh< DataType, TriangleType > () {}
  
  FETriangMesh (string filename){
    this->loadBasedOnSuffix(filename);
  }     
  
  typedef typename FETriangMesh<DataType, TriangleType >::TriangleIterator  ElementIteratorType;
  typedef typename FETriangMesh<DataType, TriangleType >::VertexIterator    NodeIteratorType;  
  
  void fillBoundaryMask( aol::BitVector& mask ) const {
    int num =  this->getNumVertices();
    mask.resize(num);
    mask.setAll( false );
    this->makeNeighbour();
    for( ElementIteratorType elementIte (*this); elementIte.notAtEnd(); ++elementIte )
      for( int i = 0; i < 3; i++ ) 
	if ( this->getNeighbour( elementIte.getIndex() , i ) == -1 ){
	  mask.set( elementIte->globNodeIdx( (i+2) % 3 ) , true);
	  mask.set( elementIte->globNodeIdx( (i+1) % 3 ) , true);
        }
  }  
};


}

#endif
