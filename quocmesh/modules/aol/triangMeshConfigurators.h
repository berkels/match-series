#ifndef __TRIANGMESHCONFIGURATORS_H
#define __TRIANGMESHCONFIGURATORS_H

#include <aol.h>
#include <configurators.h>
#include <FEOpInterface.h>
#include <solver.h>
#include <triangMesh.h>
#include <triangGaussQuadrature.h>


namespace aol {



//! BasefunctionSet for linear finite element functions on triangular meshes
//! NOTE: in order to use the gradient of the basis functions, 
//!       TriangleType has to provide a function that returns the inverted metric (first fund. form)
/*!
 * \author Droske, vDeylen, Heeren
 */
template <typename RealType, typename QuadRuleType, typename TriangleType>
class TriangMeshBaseFunctionSet : public aol::BaseFunctionSetInterface<RealType, aol::Vec3<RealType>,
  aol::Vec2<RealType>, 3, QuadRuleType, TriangMeshBaseFunctionSet<RealType, QuadRuleType, TriangleType> > {


  static RealType _b1   ( const aol::Vec2<RealType> &c ) { return 1. - c[0] - c[1]; }
  static RealType _b2   ( const aol::Vec2<RealType> &c ) { return c[0]; }
  static RealType _b3   ( const aol::Vec2<RealType> &c ) { return c[1]; }

  static RealType _d1_b1   ( const aol::Vec2<RealType> & ) { return - 1.; }
  static RealType _d1_b2   ( const aol::Vec2<RealType> & ) { return 1.; }
  static RealType _d1_b3   ( const aol::Vec2<RealType> & ) { return 0.; }

  static RealType _d2_b1   ( const aol::Vec2<RealType> & ) { return - 1.; }
  static RealType _d2_b2   ( const aol::Vec2<RealType> & ) { return 0.; }
  static RealType _d2_b3   ( const aol::Vec2<RealType> & ) { return 1.; }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const aol::Vec2<RealType> &RefCoord );
  BASIS_FUNC_TYPE _dbasis[2][3];
  BASIS_FUNC_TYPE _basis[3];

  const TriangleType *_triangle;

public:
  TriangMeshBaseFunctionSet(  ) : _triangle ( NULL ) {
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;

    _dbasis[0][0] = _d1_b1;
    _dbasis[0][1] = _d1_b2;
    _dbasis[0][2] = _d1_b3;

    _dbasis[1][0] = _d2_b1;
    _dbasis[1][1] = _d2_b2;
    _dbasis[1][2] = _d2_b3;
  }

  enum { numBaseFuncs = 3 };

  void setTriangle ( const TriangleType &T ) {
    _triangle = &T;
  }

  void evaluateGradient ( int BaseFuncNum, const aol::Vec2<RealType> &RefCoord, aol::Vec3<RealType> &Gradient ) const {
    // initialize vectors
    aol::Vec2<RealType> tmp, tmp2;
    // gradient at quad point in barycentric coords
    tmp[0] = _dbasis[0][BaseFuncNum] ( RefCoord );
    tmp[1] = _dbasis[1][BaseFuncNum] ( RefCoord );
    // change coordinate system
    _triangle->ginv(  ).mult ( tmp, tmp2 );
    // compute tangent vectors
    const aol::Vec3<RealType> dir0 = _triangle->edge(0,1);
    const aol::Vec3<RealType> dir1 = _triangle->edge(0,2);
    for ( int i = 0; i < 3; i++ ) {
      Gradient[i] = tmp2[0] * ( dir0[i] ) + tmp2[1] * ( dir1[i] );
    }
  }

  inline aol::Vec3<RealType> evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    aol::Vec3<RealType> g;
    evaluateGradient ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), g );
    return g;
  }

  RealType evaluate ( int BaseFuncNum, const aol::Vec2<RealType> &RefCoord ) const {
    return _basis[BaseFuncNum] ( RefCoord );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return evaluate ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );
  }
protected:

};


//! Mesh configurator for linear finite element functions on general triangular meshes, e.g. aol::TriangMesh<> or om::Trimesh<>.
/*!
 * \author Droske, vDeylen, Heeren
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename MeshType, typename _QuadType>
class TriangMeshConfigurator {
protected:
  const MeshType &_mesh;
public:
  typedef _RealType RealType;
  typedef _QuadType QuadType;

  static const qc::Dimension Dim = qc::QC_3D;
  enum { DomDim = 2 }; // ???

  static const qc::Dimension DimOfWorld = qc::QC_3D; // ??? currently not used, but should be correct anyway ...

  typedef RealType Real;
  typedef MeshType                                 InitType;               //!< that's the type that is needed by the constructor of the configurator
  typedef aol::Vec<2, RealType>                        DomVecType;
  typedef aol::Vec3<RealType>                          VecType;
  typedef aol::Mat<3, 3, RealType>                     MatType;
  typedef aol::Vector<RealType>                        VectorType;
  typedef aol::Vector<RealType>                        ArrayType;
  typedef aol::SparseMatrix<RealType>                  MatrixType;
  typedef aol::BitVector                               MaskType;
  typedef typename MeshType::ElementType           ElementType;
  typedef TriangMeshBaseFunctionSet<RealType, QuadType,
            ElementType>                               BaseFuncSetType;

  typedef typename MeshType::ElementIteratorType   ElementIteratorType;
  typedef typename MeshType::NodeIteratorType      NodeIteratorType;

  class DOFIterator : public NodeIteratorType {
    typedef TriangMeshConfigurator<_RealType, MeshType, _QuadType> ConfType;
  public:
    DOFIterator (const ConfType& conf ) : NodeIteratorType( conf.getInitializer() ){}
  };

  typedef DOFIterator DOFIteratorType;


  TriangMeshConfigurator ( const InitType &Mesh ) :
      _mesh ( Mesh ) {}

  //! returns the begin iterator of the grid
  inline const MeshType & begin( ) const {
    return _mesh;
  }

  //! returns the end iterator of the grid
  inline const MeshType & end( ) const {
    return _mesh;
  }

  const InitType& getInitializer( ) const { return this->_mesh; }



  mutable BaseFuncSetType _baseFuncSet;

  static const int maxNumLocalDofs = 3;

  inline int getNumLocalDofs ( const ElementType & ) const {
    return 3;
  }

  int getNumGlobalDofs( ) const {
    return this->_mesh.getNumVertices();
  }

  int maxNumQuadPoints( ) const {
    return QuadType::numQuadPoints;
  }


  const BaseFuncSetType& getBaseFunctionSet ( const ElementType &T ) const {
    _baseFuncSet.setTriangle ( T );
    return _baseFuncSet;
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const ElementType &T, int localIndex ) const {
    return T.globNodeIdx( localIndex );
  }

  //! \warning copied and pasted from QuocConfiguratorTraitMultiLin.
  //! I have absolutely no idea if this works or not.
  inline void localToGlobal ( const ElementType &T, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal ( T, localIndex0 );
    glob[1] = localToGlobal ( T, localIndex1 );
  }

  RealType vol ( const ElementType &T ) const {
    return T.area();
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    int num = getNumGlobalDofs();
    MatrixType *mat = new MatrixType ( num, num );
    return mat;
  }

  void fillBoundaryMask( MaskType& mask ) const {
    _mesh.fillBoundaryMask( mask );
  }

};


//! This class realizes the localToGlobal mapping when working with quadratic base functions,
//! i.e. we have 6 (local) degrees of freedom (dofs) per element.
/*!
 * \author Heeren, Perl, Simon
 */
template< typename MeshType >
class P2Mapper {
  typedef aol::Vec2<int> KeyType;
  typedef std::pair<KeyType, int> PairType;
  typedef std::pair<int, KeyType> InvPairType;
      
  const MeshType& _grid;
  int _globalDofs;
  
  std::map< KeyType, int > _map;        // maps triangle index and local node index to corresponding global node index
  std::map< KeyType, int > _edgeMap;    // maps global node indices of two nodes representing one edge to global index of dof on that edge 
  std::map< int, KeyType > _invEdgeMap; // maps one global node index to (i) itsef, if vertex; (ii) two adjacent global node indices, if on edge (neglecting order here!)
  
public:
  P2Mapper( const MeshType& Grid ) : _grid( Grid ) {
    initMapper();
  }
  
  int localToGlobal( const int& T, const int& loc ) const {
    return typename std::map< KeyType, int >::const_iterator( _map.find( KeyType( T, loc ) ) )->second;
  }
  
  int getGlobalDofs () const{
    return _globalDofs;
  }
  
  int getEdgeDOF( const int& v1, const int& v2 ) const {
    return typename std::map< KeyType, int >::const_iterator( _edgeMap.find( KeyType( v1, v2 ) ) )->second;
  }
  
  //
  KeyType getGlobalPairing( const int& globIdx ) const {
    KeyType res( globIdx, globIdx );
    typename std::map< int, KeyType >::const_iterator it = _invEdgeMap.find( globIdx );
    if( it != _invEdgeMap.end() )
      return it->second;
    return res;
  }
  
protected:
  void initMapper() {
    // first fill edge map to define global node indices for dofs on edges
    fillEdgeMaps();
    
    // add 6 local dofs for each element to the map
    for( typename MeshType::ElementIteratorType iter = _grid; iter.notAtEnd(); ++iter ){
      
      aol::Vec3<int> nodes ( iter.getNodeIndices() );
      int triangIdx = iter.getIndex();

      // ordering is (0,1,2,3,4,5) when going anti-clockwise
      _map.insert( PairType( KeyType(triangIdx,0), nodes[0] ) );
      _map.insert( PairType( KeyType(triangIdx,1), _edgeMap[ KeyType( min(nodes[0],nodes[1]), max(nodes[0],nodes[1]) ) ] ) );
      _map.insert( PairType( KeyType(triangIdx,2), nodes[1] ) );
      _map.insert( PairType( KeyType(triangIdx,3), _edgeMap[ KeyType( min(nodes[1],nodes[2]), max(nodes[1],nodes[2]) ) ] ) );
      _map.insert( PairType( KeyType(triangIdx,4), nodes[2] ) );
      _map.insert( PairType( KeyType(triangIdx,5), _edgeMap[ KeyType( min(nodes[2],nodes[0]), max(nodes[2],nodes[0]) ) ] ) );
      
    }
  }
  
  void fillEdgeMaps () {     
    _globalDofs = _grid.getNumVertices();
    // construct edge map and increase number of dofs
    for( typename MeshType::ElementIteratorType iter = _grid; iter.notAtEnd(); ++iter ){
      aol::Vec3<int> nodes ( iter.getNodeIndices() );
      for( int i = 0; i < 3; i++ ){
        typename std::map< KeyType, int >::const_iterator it = _edgeMap.find( KeyType( min(nodes[i],nodes[(i+1)%3]), max(nodes[i],nodes[(i+1)%3]) ) );
        if( it == _edgeMap.end() ){
          _edgeMap.insert( PairType( KeyType( min(nodes[i],nodes[(i+1)%3]), max(nodes[i],nodes[(i+1)%3]) ), _globalDofs ) );
                _invEdgeMap.insert( InvPairType( _globalDofs++, KeyType( nodes[i], nodes[(i+1)%3] ) ) ); // order in Vec2<> does not matter as it is not the key!
             }
      }
    }
  }
};


//! This class realizes the localToGlobal mapping when working with reduced cubic base functions,
//! i.e. the dof in the center of the triangle has been removed
//! Hence, we have 9 (local) degrees of freedom (dofs) per element, which are all on the boundary of that element
//! and oredered in anti-clockwise direction
/*!
 * \author Heeren, Perl, Simon
 */
template< typename MeshType >
class P3redMapper {
  typedef aol::Vec2<int> KeyType;
  typedef std::pair<KeyType, int> PairType;
  typedef std::pair<int, KeyType> InvPairType;
      
  const MeshType& _grid;
  int _globalDofs;
  
  std::map< KeyType, int > _map;        // maps triangle index and local node index to corresponding global node index
  std::map< KeyType, int > _edgeMap;    // maps global node indices of two nodes representing one edge to global index of dof on that edge 
  std::map< int, KeyType > _invEdgeMap; // maps one global node index to (i) itsef, if vertex; (ii) two adjacent global node indices, if on edge
  
public:
  P3redMapper( const MeshType& Grid ) : _grid( Grid ) {
    initMapper();
  }
  
  int localToGlobal( const int& T, const int& loc ) const {
    return typename std::map< KeyType, int >::const_iterator( _map.find( KeyType( T, loc ) ) )->second;
  }
  
  int getGlobalDofs () const{
    return _globalDofs;
  }
  
  int getEdgeDOF( const int& v1, const int& v2 ) const {
    return typename std::map< KeyType, int >::const_iterator( _edgeMap.find( KeyType( v1, v2 ) ) )->second;
  }
  
  //
  KeyType getGlobalPairing( const int& globIdx ) const {
    KeyType res( globIdx, globIdx );
    typename std::map< int, KeyType >::const_iterator it = _invEdgeMap.find( globIdx );
    if( it != _invEdgeMap.end() )
      return it->second;
    return res;
  }
  
protected:
  void initMapper() {
    // first fill edge map to define global node indices for dofs on edges
    fillEdgeMaps();
    
    // add 9 local dofs for each element to the map
    for( typename MeshType::ElementIteratorType iter = _grid; iter.notAtEnd(); ++iter ){
      
      aol::Vec3<int> nodes ( iter.getNodeIndices() );
      int triangIdx = iter.getIndex();

      // ordering is (0,1,2,3,4,5,6,7,8) when going anti-clockwise
      _map.insert( PairType( KeyType(triangIdx,0), nodes[0] ) );
      _map.insert( PairType( KeyType(triangIdx,1), _edgeMap[ KeyType( nodes[0], nodes[1] ) ] ) );
      _map.insert( PairType( KeyType(triangIdx,2), _edgeMap[ KeyType( nodes[1], nodes[0] ) ] ) );
      _map.insert( PairType( KeyType(triangIdx,3), nodes[1] ) );
      _map.insert( PairType( KeyType(triangIdx,4), _edgeMap[ KeyType( nodes[1], nodes[2] ) ] ) );
      _map.insert( PairType( KeyType(triangIdx,5), _edgeMap[ KeyType( nodes[2], nodes[1] ) ] ) );
      _map.insert( PairType( KeyType(triangIdx,6), nodes[2] ) );
      _map.insert( PairType( KeyType(triangIdx,7), _edgeMap[ KeyType( nodes[2], nodes[0] ) ] ) );
      _map.insert( PairType( KeyType(triangIdx,8), _edgeMap[ KeyType( nodes[0], nodes[2] ) ] ) );
      
    }
  }
  
  void fillEdgeMaps () {     
    _globalDofs = _grid.getNumVertices();
    // construct edge map and increase number of dofs
    for( typename MeshType::ElementIteratorType iter = _grid; iter.notAtEnd(); ++iter ){
      aol::Vec3<int> nodes ( iter.getNodeIndices() );
      for( int i = 0; i < 3; i++ ){
        typename std::map< KeyType, int >::const_iterator it = _edgeMap.find( KeyType( nodes[i], nodes[(i+1)%3] ) );
        if( it == _edgeMap.end() ){
          // first dof on edge (i,j) is closer to i, hence node i is first entry in key (j=i+1)
          _edgeMap.insert( PairType( KeyType( nodes[i], nodes[(i+1)%3] ), _globalDofs ) );
          //cerr << nodes[i] << ", " <<  nodes[(i+1)%3]  << " : " << _globalDofs << "\n";
          _invEdgeMap.insert( InvPairType( _globalDofs++, KeyType( nodes[i], nodes[(i+1)%3] ) ) );
          // second dof on edge (i,j) is closer to j, hence node j is first entry in key (j=i+1)
          _edgeMap.insert( PairType( KeyType( nodes[(i+1)%3], nodes[i] ), _globalDofs ) );
          //cerr << nodes[(i+1)%3] << ", " <<  nodes[i]  << " : " << _globalDofs << "\n";
          _invEdgeMap.insert( InvPairType( _globalDofs++, KeyType( nodes[(i+1)%3], nodes[i] ) ) );
        }
      }
    }
    //abort();
  }
};


//! BasefunctionSet for quadratic finite element functions on triangular meshes
//! NOTE: in order to use the gradient of the basis functions, 
//!       TriangleType has to provide a function that returns the inverted metric (first fund. form)
/*!
 * \author Heeren, Perl, Simon
 */
template <typename RealType, typename QuadRuleType, typename TriangleType>
class QuadTriangMeshBaseFunctionSet : public aol::BaseFunctionSetInterface<RealType, aol::Vec3<RealType>,
  aol::Vec2<RealType>, 6, QuadRuleType, TriangMeshBaseFunctionSet<RealType, QuadRuleType, TriangleType> > {

public:
  // ordering anti-clockwise beginning at node zero, i.e. _triangle-> globNodeIdx( 0 )
  static RealType _b1   ( const aol::Vec2<RealType> &c ) { return 2. * ( 1. - c[0] - c[1] ) * ( 1./2. - c[0] - c[1] ); }
  static RealType _b2   ( const aol::Vec2<RealType> &c ) { return 4. * c[0] * ( 1. - c[0] - c[1] ); }
  static RealType _b3   ( const aol::Vec2<RealType> &c ) { return c[0] * ( 2. * c[0] - 1. ); }
  static RealType _b4   ( const aol::Vec2<RealType> &c ) { return 4. * c[0] * c[1]; }
  static RealType _b5   ( const aol::Vec2<RealType> &c ) { return c[1] * ( 2. * c[1] - 1. ); }
  static RealType _b6   ( const aol::Vec2<RealType> &c ) { return 4. * c[1] * ( 1. - c[0] - c[1] ); }

  static RealType _d1_b1   ( const aol::Vec2<RealType> &c ) { return - 4. * ( 3./4. - c[0] - c[1] ); }
  static RealType _d1_b2   ( const aol::Vec2<RealType> &c ) { return 4. - 8.* c[0] - 4. * c[1]; }
  static RealType _d1_b3   ( const aol::Vec2<RealType> &c ) { return 4. * c[0] - 1.; }
  static RealType _d1_b4   ( const aol::Vec2<RealType> &c ) { return 4. * c[1]; }
  static RealType _d1_b5   ( const aol::Vec2<RealType> &  ) { return 0.; }
  static RealType _d1_b6   ( const aol::Vec2<RealType> &c ) { return -4. * c[1]; }

  static RealType _d2_b1   ( const aol::Vec2<RealType> &c ) { return - 4. * ( 3./4. - c[0] - c[1] ); }
  static RealType _d2_b2   ( const aol::Vec2<RealType> &c ) { return - 4. * c[0]; }
  static RealType _d2_b3   ( const aol::Vec2<RealType> &  ) { return 0.; }
  static RealType _d2_b4   ( const aol::Vec2<RealType> &c ) { return 4. * c[0]; }
  static RealType _d2_b5   ( const aol::Vec2<RealType> &c ) { return 4. * c[1] - 1.; }
  static RealType _d2_b6   ( const aol::Vec2<RealType> &c ) { return 4. - 4. * c[0] - 8. * c[1]; }
  
protected:
  typedef RealType ( *BASIS_FUNC_TYPE ) ( const aol::Vec2<RealType> &RefCoord );
  BASIS_FUNC_TYPE _dbasis[2][6];
  BASIS_FUNC_TYPE _basis[6];

  const TriangleType *_triangle;

public:
  QuadTriangMeshBaseFunctionSet(  ) : _triangle ( NULL ) {
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;
    _basis[4] = _b5;
    _basis[5] = _b6;

    _dbasis[0][0] = _d1_b1;
    _dbasis[0][1] = _d1_b2;
    _dbasis[0][2] = _d1_b3;
    _dbasis[0][3] = _d1_b4;
    _dbasis[0][4] = _d1_b5;
    _dbasis[0][5] = _d1_b6;

    _dbasis[1][0] = _d2_b1;
    _dbasis[1][1] = _d2_b2;
    _dbasis[1][2] = _d2_b3;
    _dbasis[1][3] = _d2_b4;
    _dbasis[1][4] = _d2_b5;
    _dbasis[1][5] = _d2_b6;
  }

  enum { numBaseFuncs = 6 };

  void setTriangle ( const TriangleType &T ) {
    _triangle = &T;
  }
  
  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return evaluate ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );
  }

  inline aol::Vec3<RealType> evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    aol::Vec3<RealType> g;
    evaluateGradient ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), g );
    return g;
  }

protected:
  
  RealType evaluate ( int BaseFuncNum, const aol::Vec2<RealType> &RefCoord ) const {
    return _basis[BaseFuncNum] ( RefCoord );
  }

  void evaluateGradient ( int BaseFuncNum, const aol::Vec2<RealType> &RefCoord, aol::Vec3<RealType> &Gradient ) const {
    //
    aol::Vec2<RealType> tmp, tmp2;
    aol::Vec3<RealType> dir0, dir1;
    //
    tmp[0] = _dbasis[0][BaseFuncNum] ( RefCoord );
    tmp[1] = _dbasis[1][BaseFuncNum] ( RefCoord );
    //
    _triangle->ginv().mult ( tmp, tmp2 );
    //
    dir0 = _triangle->edge(0,1);
    dir1 = _triangle->edge(0,2);
    for ( int i = 0; i < 3; i++ ) {
      Gradient[i] = tmp2[0] * ( dir0[i] ) + tmp2[1] * ( dir1[i] );
    }
  }
};


//! BasefunctionSet for quadratic finite element functions on triangular meshes
//! NOTE: in order to use the gradient of the basis functions, 
//!       TriangleType has to provide a function that returns the inverted metric (first fund. form)
/*!
 * \author Heeren, Perl, Simon
 */
template <typename RealType, typename QuadRuleType, typename TriangleType>
class ReducedCubicTriangMeshBaseFunctionSet : public aol::BaseFunctionSetInterface<RealType, aol::Vec3<RealType>,
  aol::Vec2<RealType>, 9, QuadRuleType, TriangMeshBaseFunctionSet<RealType, QuadRuleType, TriangleType> > {

public:
  // ordering anti-clockwise beginning at node zero, i.e. _triangle-> globNodeIdx( 0 )
  static RealType _b1   ( const aol::Vec2<RealType> &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 1./2. * c2 * ( 3. * c2 - 1. ) * ( 3. * c2 - 2. ) - 9./2. * c[0] * c[1] * c2;
  }
  static RealType _b2   ( const aol::Vec2<RealType> &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 9./2. * c2 * c[0] * ( 3. * c2 - 1. ) + 27./4. * c[0] * c[1] * c2;
  }
  static RealType _b3   ( const aol::Vec2<RealType> &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 9./2. * c2 * c[0] * ( 3. * c[0] - 1. ) + 27./4. * c[0] * c[1] * c2;
  }
  static RealType _b4   ( const aol::Vec2<RealType> &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 1./2. * c[0] * ( 3. * c[0] - 1. ) * ( 3. * c[0] - 2. ) - 9./2. * c[0] * c[1] * c2;
  }
  static RealType _b5   ( const aol::Vec2<RealType> &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 9./2. * c[0] * c[1] * ( 3. * c[0] - 1. ) + 27./4. * c[0] * c[1] * c2;
  }
  static RealType _b6   ( const aol::Vec2<RealType> &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 9./2. * c[1] * c[0] * ( 3. * c[1] - 1. ) + 27./4. * c[0] * c[1] * c2;
  }
  static RealType _b7   ( const aol::Vec2<RealType> &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 1./2. * c[1] * ( 3. * c[1] - 1. ) * ( 3. * c[1] - 2. ) - 9./2. * c[0] * c[1] * c2; 
  }
  static RealType _b8   ( const aol::Vec2<RealType> &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 9./2. * c[1] * c2 * ( 3. * c[1] - 1. ) + 27./4. * c[0] * c[1] * c2; 
  }
  static RealType _b9   ( const aol::Vec2<RealType> &c ) { 
    RealType c2 = 1. - c[0] - c[1];
    return 9./2. * c2 * c[1] * ( 3. * c2 - 1. ) + 27./4. * c[0] * c[1] * c2;
  }

  //
  static RealType _d1_b1   ( const aol::Vec2<RealType> &c ) { 
    return -11./2. + 18. * c[0] - 27./2. * Sqr(c[0]) + 27./2. * c[1] - 18. * c[0] * c[1] - 9. * Sqr(c[1]);
  }
  static RealType _d1_b2   ( const aol::Vec2<RealType> &c ) { 
    return 9. - 45. * c[0] + 81./2. * Sqr(c[0]) - 63./4. * c[1] + 81./2. * c[0] * c[1] + 27./4. * Sqr(c[1]);
  }
  static RealType _d1_b3   ( const aol::Vec2<RealType> &c ) { 
    return - 9./2. + 36. * c[0] - 81./2. * Sqr(c[0]) + 45./4. * c[1] - 81./2. * c[0] * c[1] - 27./4. * Sqr(c[1]);
  }
  static RealType _d1_b4   ( const aol::Vec2<RealType> &c ) { 
    return 1. - 9. * c[0] + 27./2. * Sqr(c[0]) - 9./2. * c[1] + 9. * c[0] * c[1] + 9./2. * Sqr(c[1]);
  }
  static RealType _d1_b5   ( const aol::Vec2<RealType> &c ) { 
    return 9./4. * c[1] + 27./2. * c[0] * c[1] - 27./4. * Sqr(c[1]);
  }
  static RealType _d1_b6   ( const aol::Vec2<RealType> &c ) { 
    return 9./4. * c[1] - 27./2. * c[0] * c[1] + 27./4. * Sqr(c[1]);
  }
  static RealType _d1_b7   ( const aol::Vec2<RealType> &c ) { 
    return - 9./2. * c[1] + 9. * c[0] * c[1] + 9./2. * Sqr(c[1]);
  }
  static RealType _d1_b8   ( const aol::Vec2<RealType> &c ) { 
    return 45./4. * c[1] - 27./2. * c[0] * c[1] - 81./4. * Sqr(c[1]);
  }
  static RealType _d1_b9   ( const aol::Vec2<RealType> &c ) { 
    return - 63./4. * c[1] + 27./2. * c[0] * c[1] + 81./4. * Sqr(c[1]);
  }

  //
  static RealType _d2_b1   ( const aol::Vec2<RealType> &c ) { 
    return -11./2. + 27./2. * c[0] - 9. * Sqr(c[0]) + 18. * c[1] - 18. * c[0] * c[1] - 27./2. * Sqr(c[1]);
  }
  static RealType _d2_b2   ( const aol::Vec2<RealType> &c ) { 
    return - 63./4. * c[0] + 81./4. * Sqr(c[0]) + 27./2. * c[0] * c[1] ;
  }
  static RealType _d2_b3   ( const aol::Vec2<RealType> &c ) { 
    return 45./4. * c[0] - 81./4. * Sqr(c[0]) - 27./2. * c[0] * c[1] ;
  }
  static RealType _d2_b4   ( const aol::Vec2<RealType> &c ) { 
    return -9./2. * c[0] + 9./2. * Sqr(c[0]) + 9. * c[0] * c[1] ;
  }
  static RealType _d2_b5   ( const aol::Vec2<RealType> &c ) { 
    return 9./4. * c[0] + 27./4. * Sqr(c[0]) - 27./2. * c[0] * c[1] ;
  }
  static RealType _d2_b6   ( const aol::Vec2<RealType> &c ) { 
    return 9./4. * c[0] - 27./4. * Sqr(c[0]) + 27./2. * c[0] * c[1] ;
  }
  static RealType _d2_b7   ( const aol::Vec2<RealType> &c ) { 
    return 1. - 9./2. * c[0] + 9./2. * Sqr(c[0]) - 9. * c[1] + 9. * c[0] * c[1] + 27./2. * Sqr(c[1]);
  }
  static RealType _d2_b8   ( const aol::Vec2<RealType> &c ) { 
    return - 9./2. + 45./4. * c[0] - 27./4. * Sqr(c[0]) + 36. * c[1] - 81./2. * c[0] * c[1] - 81./2. * Sqr(c[1]);
  }
  static RealType _d2_b9   ( const aol::Vec2<RealType> &c ) { 
    return 9. - 63./4. * c[0] + 27./4. * Sqr(c[0]) - 45. * c[1] + 81./2. * c[0] * c[1] + 81./2. * Sqr(c[1]);
  }
  
protected:
  typedef RealType ( *BASIS_FUNC_TYPE ) ( const aol::Vec2<RealType> &RefCoord );
  BASIS_FUNC_TYPE _dbasis[2][9];
  BASIS_FUNC_TYPE _basis[9];

  const TriangleType *_triangle;

public:
  ReducedCubicTriangMeshBaseFunctionSet(  ) : _triangle ( NULL ) {
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;
    _basis[4] = _b5;
    _basis[5] = _b6;
    _basis[6] = _b7;
    _basis[7] = _b8;
    _basis[8] = _b9;

    _dbasis[0][0] = _d1_b1;
    _dbasis[0][1] = _d1_b2;
    _dbasis[0][2] = _d1_b3;
    _dbasis[0][3] = _d1_b4;
    _dbasis[0][4] = _d1_b5;
    _dbasis[0][5] = _d1_b6;
    _dbasis[0][6] = _d1_b7;
    _dbasis[0][7] = _d1_b8;
    _dbasis[0][8] = _d1_b9;

    _dbasis[1][0] = _d2_b1;
    _dbasis[1][1] = _d2_b2;
    _dbasis[1][2] = _d2_b3;
    _dbasis[1][3] = _d2_b4;
    _dbasis[1][4] = _d2_b5;
    _dbasis[1][5] = _d2_b6;
    _dbasis[1][6] = _d2_b7;
    _dbasis[1][7] = _d2_b8;
    _dbasis[1][8] = _d2_b9;
  }

  enum { numBaseFuncs = 9 };

  void setTriangle ( const TriangleType &T ) {
    _triangle = &T;
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return evaluate ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );
  }
  
  inline aol::Vec3<RealType> evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    aol::Vec3<RealType> g;
    evaluateGradient ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), g );
    return g;
  }
  
protected:
  RealType evaluate ( int BaseFuncNum, const aol::Vec2<RealType> &RefCoord ) const {
    return _basis[BaseFuncNum] ( RefCoord );
  }
  
  void evaluateGradient ( int BaseFuncNum, const aol::Vec2<RealType> &RefCoord, aol::Vec3<RealType> &Gradient ) const {
    //
    aol::Vec2<RealType> tmp, tmp2;
    aol::Vec3<RealType> dir0, dir1;
    //
    tmp[0] = _dbasis[0][BaseFuncNum] ( RefCoord );
    tmp[1] = _dbasis[1][BaseFuncNum] ( RefCoord );
    //
    _triangle->ginv().mult ( tmp, tmp2 );
    //
    dir0 = _triangle->edge(0,1);
    dir1 = _triangle->edge(0,2);
    for ( int i = 0; i < 3; i++ ) {
      Gradient[i] = tmp2[0] * ( dir0[i] ) + tmp2[1] * ( dir1[i] );
    }
  }

};


/*! \brief BasefunctionSet for gradient (i.e. \f$ \theta = \nabla w \f$) in Discrete Kirchhoff Triangle (DKT)
 *  In particular, the element has to be aol::DKTPlateElement<RealType>.
 * 
 *  The 9 DOFs per element T are denoted by \f$ W = [ w_0 w_{x0} w_{y0} w_1 w_{x1} w_{y1} w_2 w_{x2} w_{y2} ] \f$,
 *  where "w_i" refers to the function value of "w" at the vertex in "T" with local index "i"
 *  and "w_{zi}" to the corresponding derivative with respect to "z".
 *  
 *  The base functions are given by vector valued functions \f$ \theta \f$,
 *  i.e. \f$ \theta = (\theta^1, \theta^2): T \rightarrow \R^2, (\xi, \eta) \mapsto \theta(\xi, \eta) \f$.
 *  We shall use the notation \f$ w(\xi,\eta) = ( H_x W, H_y W )^T \f$.
 * 
 *  In applications with DKT the most relevant object is the stiffnes matrix \f$ L_ij = \int C \epsilon[\theta_i] : \epsilon[\theta_j] dx \f$,
 *  where \f$ \epsilon[\theta] \in \R^{2,2} \f$ denotes the symmetrized gradient and \f$ C \f$ is the elastic tensor, 
 *  i.e. \f$ \epsilon[\theta]_ij = 1/2 ( \partial_i\theta^j + \partial_j\theta^i) \f$.
 *  However, the evaluation of the gradient of \f$ \theta \f$ returns \f$ ( \partial_1\theta^1, \partial_2\theta^2, \partial_1\theta^2 + \partial_2\theta^1 ) \in \R^3 \f$.
 * 
 *  As there is no direct need for evaluations of \f$ \theta \f$, we did not implement them so far.
 *
 *  \author Heeren, Perl, Simon
 */
template <typename RealType, typename QuadRuleType>
class DKTPlateGradientBaseFunctionSet : public aol::BaseFunctionSetInterface<RealType, aol::Vec3<RealType>, 
  aol::Vec2<RealType>, 9, QuadRuleType, TriangMeshBaseFunctionSet<RealType, QuadRuleType, aol::DKTPlateElement<RealType> > > {
    
  typedef aol::QuadTriangMeshBaseFunctionSet<RealType, QuadRuleType, aol::DKTPlateElement<RealType> > P2ShapeFnc;
  typedef aol::DKTPlateElement<RealType> TriangleType;
  
 public:
  // Numbering of local DOFs corresponds to order in W, i.e. 0 = w_0, 1 = w_{x0}, ..., 8 = w_{y2} 
  // Ordering of vertices anti-clockwise beginning at node zero, i.e. _triangle-> globNodeIdx( 0 ) 
  
  //!NOTE numbering of shape function of P2 differs from the one in the tex file!
  
  // Hx
  static RealType _tet11 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return 3.*( t->edgeXNormSqr(1) * P2ShapeFnc::_b6(c) - t->edgeXNormSqr(2) * P2ShapeFnc::_b2(c) ) / 2.;
  }
  
  static RealType _tet12 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return ( 4. * P2ShapeFnc::_b1(c) - t->edgeXXMin2YYNormSqr(2) * P2ShapeFnc::_b2(c) - t->edgeXXMin2YYNormSqr(1) * P2ShapeFnc::_b6(c) ) / 4.;
  }
  
  static RealType _tet13 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return -3.*( t->edgeXYNormSqr(2) * P2ShapeFnc::_b2(c) - t->edgeXYNormSqr(1) * P2ShapeFnc::_b6(c) ) / 4.;
  }
  
  static RealType _tet14 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return 3.*( t->edgeXNormSqr(2) * P2ShapeFnc::_b2(c) - t->edgeXNormSqr(0) * P2ShapeFnc::_b4(c) ) / 2.;
  } 
    
  static RealType _tet15 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return ( 4. * P2ShapeFnc::_b3(c) - t->edgeXXMin2YYNormSqr(2) * P2ShapeFnc::_b2(c) - t->edgeXXMin2YYNormSqr(0) * P2ShapeFnc::_b4(c) ) / 4.;
  }
  
  static RealType _tet16 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return -3.*( t->edgeXYNormSqr(2) * P2ShapeFnc::_b2(c) - t->edgeXYNormSqr(0) * P2ShapeFnc::_b4(c) ) / 4.;
  }
    
  static RealType _tet17 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return 3.*( t->edgeXNormSqr(0) * P2ShapeFnc::_b4(c) - t->edgeXNormSqr(1) * P2ShapeFnc::_b6(c) ) / 2.;
  } 
      
  static RealType _tet18 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return ( 4. * P2ShapeFnc::_b5(c) - t->edgeXXMin2YYNormSqr(0) * P2ShapeFnc::_b4(c) - t->edgeXXMin2YYNormSqr(1) * P2ShapeFnc::_b6(c) ) / 4.;
  }
    
  static RealType _tet19 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return -3.*( t->edgeXYNormSqr(0) * P2ShapeFnc::_b4(c) - t->edgeXYNormSqr(1) * P2ShapeFnc::_b6(c) ) / 4.;
  }

  //Hy
  static RealType _tet21 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return 3.*( t->edgeYNormSqr(1) * P2ShapeFnc::_b6(c) - t->edgeYNormSqr(2) * P2ShapeFnc::_b2(c) ) / 2.;
  }
  
  static RealType _tet22 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return -3.*( t->edgeXYNormSqr(2) * P2ShapeFnc::_b2(c) - t->edgeXYNormSqr(1) * P2ShapeFnc::_b6(c) ) / 4.;
  }
    
  static RealType _tet23 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return ( 4. * P2ShapeFnc::_b1(c) - t->edgeYYMin2XXNormSqr(2) * P2ShapeFnc::_b2(c) - t->edgeYYMin2XXNormSqr(1) * P2ShapeFnc::_b6(c) ) / 4.;
  }
  
  static RealType _tet24 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return 3.*( t->edgeYNormSqr(2) * P2ShapeFnc::_b2(c) - t->edgeYNormSqr(0) * P2ShapeFnc::_b4(c) ) / 2.;
  }
  
  static RealType _tet25 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return -3.*( t->edgeXYNormSqr(2) * P2ShapeFnc::_b2(c) - t->edgeXYNormSqr(0) * P2ShapeFnc::_b4(c) ) / 4.;
  }
    
  static RealType _tet26 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return ( 4. * P2ShapeFnc::_b3(c) - t->edgeYYMin2XXNormSqr(2) * P2ShapeFnc::_b2(c) - t->edgeYYMin2XXNormSqr(0) * P2ShapeFnc::_b4(c) ) / 4.;
  }
  
  static RealType _tet27 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return 3.*( t->edgeYNormSqr(0) * P2ShapeFnc::_b4(c) - t->edgeYNormSqr(1) * P2ShapeFnc::_b6(c) ) / 2.;
  }
  
  static RealType _tet28 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return -3.*( t->edgeXYNormSqr(0) * P2ShapeFnc::_b4(c) - t->edgeXYNormSqr(1) * P2ShapeFnc::_b6(c) ) / 4.;
  }
    
  static RealType _tet29 ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return ( 4. * P2ShapeFnc::_b5(c) - t->edgeYYMin2XXNormSqr(0) * P2ShapeFnc::_b4(c) - t->edgeYYMin2XXNormSqr(1) * P2ShapeFnc::_b6(c) ) / 4.;
  }
  
  
  // Hx,\xi
  static RealType _d1_tet11   ( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return - 6. * t->edgeXNormSqr(2) * ( 1. - 2. * c[0] ) + 6. * ( t->edgeXNormSqr(2) - t->edgeXNormSqr(1) ) * c[1];
  }
  
  static RealType _d1_tet12   ( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return -1. * t->edgeXXMin2YYNormSqr(2) * ( 1. - 2. * c[0] ) + 4. * ( c[0] + c[1] ) + ( t->edgeXXMin2YYNormSqr(1) + t->edgeXXMin2YYNormSqr(2) ) * c[1] - 3.;
  }
 
  static RealType _d1_tet13   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {   
    return -3. * t->edgeXYNormSqr(2) * ( 1. - 2. * c[0] ) + 3. * ( t->edgeXYNormSqr(1) + t->edgeXYNormSqr(2) ) * c[1];
  }

  static RealType _d1_tet14   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {   
    return 6. * t->edgeXNormSqr(2) * ( 1. - 2. * c[0] ) - 6. * ( t->edgeXNormSqr(0) + t->edgeXNormSqr(2) ) * c[1]  ;
  }  
 
  static RealType _d1_tet15   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return -1. * t->edgeXXMin2YYNormSqr(2) * ( 1. - 2. * c[0] ) + ( t->edgeXXMin2YYNormSqr(2) - t->edgeXXMin2YYNormSqr(0) ) * c[1] + 4.*c[0] - 1.;
  }
  
  static RealType _d1_tet16   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return -3. * t->edgeXYNormSqr(2) * (1. - 2. * c[0] ) + 3. * ( t->edgeXYNormSqr(2) - t->edgeXYNormSqr(0) ) * c[1];
  }
 
  static RealType _d1_tet17   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {   
    return 6. * ( t->edgeXNormSqr(0) + t->edgeXNormSqr(1) ) * c[1];
  } 

  static RealType _d1_tet18   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {   
    return ( t->edgeXXMin2YYNormSqr(1) - t->edgeXXMin2YYNormSqr(0) ) * c[1];
  }
  
  static RealType _d1_tet19   ( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return 3. * ( t->edgeXYNormSqr(1) - t->edgeXYNormSqr(0) ) * c[1];
  }

  // Hx,\eta
  static RealType _d2_tet11   ( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return 6. * t->edgeXNormSqr(1) * ( 1. - 2. * c[1] ) + 6. * ( t->edgeXNormSqr(2) - t->edgeXNormSqr(1) ) * c[0];
  }
  
  static RealType _d2_tet12   ( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return -1. * t->edgeXXMin2YYNormSqr(1) * ( 1. - 2. * c[1] ) + 4. * ( c[0] + c[1] ) + ( t->edgeXXMin2YYNormSqr(2) + t->edgeXXMin2YYNormSqr(1) ) * c[0]  - 3.;
  }    
   
  static RealType _d2_tet13   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return -3. * t->edgeXYNormSqr(1) * ( 1. - 2 * c[1] ) + 3.*( t->edgeXYNormSqr(2) + t->edgeXYNormSqr(1) ) * c[0];
  }  
 
  static RealType _d2_tet14   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return -6. * ( t->edgeXNormSqr(2) + t->edgeXNormSqr(0) ) * c[0];
  }  
  
  static RealType _d2_tet15   ( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( t->edgeXXMin2YYNormSqr(2) - t->edgeXXMin2YYNormSqr(0) ) * c[0];
  }
  
  static RealType _d2_tet16   ( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return 3. * ( t->edgeXYNormSqr(2) - t->edgeXYNormSqr(0) ) * c[0];
  }
  
  static RealType _d2_tet17   ( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return -6. * t->edgeXNormSqr(1) * ( 1. - 2. * c[1] ) + 6. * ( t->edgeXNormSqr(0) + t->edgeXNormSqr(1) ) * c[0];
  }
  
  static RealType _d2_tet18   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {   
    return -1. * t->edgeXXMin2YYNormSqr(1) * ( 1. - 2. * c[1] ) + 4. * c[1]   + ( t->edgeXXMin2YYNormSqr(1) - t->edgeXXMin2YYNormSqr(0) ) * c[0] - 1.;
  }
  
  static RealType _d2_tet19   ( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return -3. * t->edgeXYNormSqr(1) * ( 1. - 2. * c[1] ) + 3.*( t->edgeXYNormSqr(1) - t->edgeXYNormSqr(0) ) * c[0];
  }
  
  
  // Hy,\xi
  static RealType _d1_tet21   ( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return  -6. * t->edgeYNormSqr(2) * (1. - 2. * c[0]) + 6. * (t->edgeYNormSqr(2) - t->edgeYNormSqr(1)) * c[1];
  }
  
  static RealType _d1_tet22   ( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return -3. * t->edgeXYNormSqr(2) * ( 1. - 2. * c[0] ) + 3. * ( t->edgeXYNormSqr(2) + t->edgeXYNormSqr(1) ) * c[1];
  }
 
  static RealType _d1_tet23   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {   
    return  -1. * t->edgeYYMin2XXNormSqr(2) * ( 1. - 2. * c[0] ) + 4. * ( c[0] + c[1] ) + ( t->edgeYYMin2XXNormSqr(2) + t->edgeYYMin2XXNormSqr(1) ) * c[1] - 3.;
  }

  static RealType _d1_tet24   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {   
    return 6. * t->edgeYNormSqr(2) * ( 1. - 2. * c[0] ) - 6. * ( t->edgeYNormSqr(2) + t->edgeYNormSqr(0) ) * c[1];
  }  
 
  static RealType _d1_tet25   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return -3. * t->edgeXYNormSqr(2) * ( 1. - 2. * c[0] ) + 3. * ( t->edgeXYNormSqr(2) - t->edgeXYNormSqr(0) ) * c[1];
  }
  
  static RealType _d1_tet26   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return  -1. * t->edgeYYMin2XXNormSqr(2) * ( 1. - 2. * c[0] ) + ( t->edgeYYMin2XXNormSqr(2) - t->edgeYYMin2XXNormSqr(0) ) *  c[1] + 4. * c[0] - 1.;
  }
 
  static RealType _d1_tet27   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {   
    return 6. * ( t->edgeYNormSqr(0) + t->edgeYNormSqr(1) ) *  c[1];
  } 

  static RealType _d1_tet28   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {   
    return 3. * ( t->edgeXYNormSqr(1) - t->edgeXYNormSqr(0) ) * c[1];
  }
  
  static RealType _d1_tet29   ( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( t->edgeYYMin2XXNormSqr(1) - t->edgeYYMin2XXNormSqr(0) ) * c[1];
  }
  
  // Hy,\eta
  static RealType _d2_tet21   ( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return  6. * t->edgeYNormSqr(1) * ( 1. - 2. * c[1] ) + 6. * ( t->edgeYNormSqr(2) - t->edgeYNormSqr(1) ) * c[0];
  }
  
  static RealType _d2_tet22   ( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return  -3. * t->edgeXYNormSqr(1) *  ( 1. - 2. * c[1] ) + 3. * ( t->edgeXYNormSqr(2) + t->edgeXYNormSqr(1) ) * c[0];
  }
 
  static RealType _d2_tet23   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {   
    return -1. * t->edgeYYMin2XXNormSqr(1) * ( 1. - 2. * c[1] ) + 4.*( c[0] + c[1] ) + ( t->edgeYYMin2XXNormSqr(2) + t->edgeYYMin2XXNormSqr(1) ) * c[0] - 3.;
  }

  static RealType _d2_tet24   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {   
    return -6. * ( t->edgeYNormSqr(2) + t->edgeYNormSqr(0) ) * c[0];
  }  
 
  static RealType _d2_tet25   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return 3. * ( t->edgeXYNormSqr(2) - t->edgeXYNormSqr(0) ) * c[0];
  }
  
  static RealType _d2_tet26   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {  
    return ( t->edgeYYMin2XXNormSqr(2) - t->edgeYYMin2XXNormSqr(0) ) * c[0];
  }
 
  static RealType _d2_tet27   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {   
    return  - 6. * t->edgeYNormSqr(1) * ( 1. - 2. * c[1] ) + 6. * ( t->edgeYNormSqr(0) + t->edgeYNormSqr(1) ) * c[0];
  } 

  static RealType _d2_tet28   ( const TriangleType* t, const aol::Vec2<RealType> &c ) {   
    return  -3. * t->edgeXYNormSqr(1) * ( 1. - 2. * c[1] ) + 3.*( t->edgeXYNormSqr(1) - t->edgeXYNormSqr(0) ) * c[0];
  }
  
  static RealType _d2_tet29   ( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return -1. * t->edgeYYMin2XXNormSqr(1) * ( 1. - 2. * c[1] ) + ( t->edgeYYMin2XXNormSqr(1) - t->edgeYYMin2XXNormSqr(0) ) * c[0]  + 4. * c[1] - 1.;
  }

protected:
  //
  typedef RealType ( *BASIS_FUNC_TYPE ) ( const TriangleType* Triangle, const aol::Vec2<RealType> &RefCoord );
  BASIS_FUNC_TYPE _basis1[9];
  BASIS_FUNC_TYPE _basis2[9];
  BASIS_FUNC_TYPE _dbasis1[2][9];
  BASIS_FUNC_TYPE _dbasis2[2][9];

  const TriangleType *_triangle;

public:
  DKTPlateGradientBaseFunctionSet(  ) : _triangle ( NULL ) {
    // Hx
    _basis1[0] = _tet11;
    _basis1[1] = _tet12;
    _basis1[2] = _tet13;
    _basis1[3] = _tet14;
    _basis1[4] = _tet15;
    _basis1[5] = _tet16;
    _basis1[6] = _tet17;
    _basis1[7] = _tet18;
    _basis1[8] = _tet19;
    
    // Hy
    _basis2[0] = _tet21;
    _basis2[1] = _tet22;
    _basis2[2] = _tet23;
    _basis2[3] = _tet24;
    _basis2[4] = _tet25;
    _basis2[5] = _tet26;
    _basis2[6] = _tet27;
    _basis2[7] = _tet28;
    _basis2[8] = _tet29;
    
    // Hx,\xi
    _dbasis1[0][0] = _d1_tet11;
    _dbasis1[0][1] = _d1_tet12;
    _dbasis1[0][2] = _d1_tet13;
    _dbasis1[0][3] = _d1_tet14;
    _dbasis1[0][4] = _d1_tet15;
    _dbasis1[0][5] = _d1_tet16;
    _dbasis1[0][6] = _d1_tet17;
    _dbasis1[0][7] = _d1_tet18;
    _dbasis1[0][8] = _d1_tet19;
    
    // Hx,\eta
    _dbasis1[1][0] = _d2_tet11;
    _dbasis1[1][1] = _d2_tet12;
    _dbasis1[1][2] = _d2_tet13;
    _dbasis1[1][3] = _d2_tet14;
    _dbasis1[1][4] = _d2_tet15;
    _dbasis1[1][5] = _d2_tet16;
    _dbasis1[1][6] = _d2_tet17;
    _dbasis1[1][7] = _d2_tet18;
    _dbasis1[1][8] = _d2_tet19;
    
    // Hy,\xi
    _dbasis2[0][0] = _d1_tet21;
    _dbasis2[0][1] = _d1_tet22;
    _dbasis2[0][2] = _d1_tet23;
    _dbasis2[0][3] = _d1_tet24;
    _dbasis2[0][4] = _d1_tet25;
    _dbasis2[0][5] = _d1_tet26;
    _dbasis2[0][6] = _d1_tet27;
    _dbasis2[0][7] = _d1_tet28;
    _dbasis2[0][8] = _d1_tet29;
    
    // Hy, \eta
    _dbasis2[1][0] = _d2_tet21;
    _dbasis2[1][1] = _d2_tet22;
    _dbasis2[1][2] = _d2_tet23;
    _dbasis2[1][3] = _d2_tet24;
    _dbasis2[1][4] = _d2_tet25;
    _dbasis2[1][5] = _d2_tet26;
    _dbasis2[1][6] = _d2_tet27;
    _dbasis2[1][7] = _d2_tet28;
    _dbasis2[1][8] = _d2_tet29;
  }

  enum { numBaseFuncs = 9 };

  void setTriangle ( const TriangleType &T ) {
    _triangle = &T;
  }

  inline const RealType evaluate ( int /*BaseFuncNum*/, int /*QuadPoint*/ ) const {
    throw aol::Exception ( "DKTPlateGradientBaseFunctionSet::evaluate() does no make sense! Use DKTPlateGradientBaseFunctionSet::getTheta() instead!", __FILE__, __LINE__ );
    return 0.;
  }
  
  inline aol::Vec2<RealType> getTheta ( int BaseFuncNum, int QuadPoint ) const {
    return aol::Vec2<RealType>( _basis1[BaseFuncNum](_triangle, this->_quadRule.getRefCoord ( QuadPoint ) ), _basis2[BaseFuncNum](_triangle, this->_quadRule.getRefCoord ( QuadPoint ) ) );
  }
  
  void getTheta ( int BaseFuncNum, const aol::Vec2<RealType> &RefCoord, aol::Vec2<RealType>& theta ) const {
    theta[0] = _basis1[BaseFuncNum](RefCoord);
    theta[1] = _basis2[BaseFuncNum](RefCoord);
  }
  
  inline aol::Vec3<RealType> evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {    
    aol::Vec3<RealType> g;    
//     aol::Vec2<RealType> theta = getTheta( BaseFuncNum, QuadPoint );
//     aol::Vec2<RealType> tmp;
//     _triangle->ginv().mult ( theta, tmp );
//     for ( int i = 0; i < 2; i++ )
//       g[i] = tmp[0] * _triangle->edge(2)[i] - tmp[1] * _triangle->edge(1)[i];
//     g[2] = 0.;
    evaluateGradient ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), g );
    return g;
  }
  
protected:  
  void evaluateGradient ( int BaseFuncNum, const aol::Vec2<RealType> &RefCoord, aol::Vec3<RealType> &Gradient ) const {    
    // chain rule: D_x F(x(z)) = D_z F(x(z)) D_z x(z)^{-1}   [ as D_z (F \circ x)(z) = D_x F(x(z)) D_z x(z) ]
    // Here D_z F corresponds to H, however, we want to return D_x F = H * g^{-1}    [as Dx = g]
    // Note: gradient of \theta: T \rightarrow \R^2 is symmetric 2x2-matrix D\theta
    // return here  \grad\theta := (D\theta_11, D\theta_22, D\theta_12 + D\theta_21) \in \R^3
    
//     const aol::Matrix22<RealType>& ginv = _triangle->ginv();
    const aol::Matrix22<RealType>& metric = _triangle->metric();
    
    RealType Hx1( _dbasis1[0][BaseFuncNum](_triangle, RefCoord) ), Hx2( _dbasis1[1][BaseFuncNum](_triangle, RefCoord) );
    RealType Hy1( _dbasis2[0][BaseFuncNum](_triangle, RefCoord) ), Hy2( _dbasis2[1][BaseFuncNum](_triangle, RefCoord) );
    
//     Gradient[0] = ginv.get(0,0) * Hx1 + ginv.get(1,0) * Hx2;
//     Gradient[1] = ginv.get(0,1) * Hy1 + ginv.get(1,1) * Hy2;
//     Gradient[2] = ginv.get(0,1) * Hx1 + ginv.get(1,1) * Hx2 + ginv.get(0,0) * Hy1 + ginv.get(1,0) * Hy2;
    Gradient[0] = metric.get(0,0) * Hx1 + metric.get(1,0) * Hx2;
    Gradient[1] = metric.get(0,1) * Hy1 + metric.get(1,1) * Hy2;
    Gradient[2] = metric.get(0,1) * Hx1 + metric.get(1,1) * Hx2 + metric.get(0,0) * Hy1 + metric.get(1,0) * Hy2;
  }
};


//! Mesh configurator for general quadratic triangular meshes, e.g. aol::TriangMesh<> or om::Trimesh<>.
/*!
 * \author Heeren, Perl, Simon
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename MeshType, typename _QuadType>
class QuadTriangMeshConfigurator {
protected:
  const MeshType &_mesh;
  const P2Mapper<MeshType> _dofMapper;
 
public:
  // --- dof iterator ---
  class DOFIterator {

  protected:
    typedef QuadTriangMeshConfigurator<_RealType, MeshType, _QuadType> ConfType;
    const ConfType &  _conf;
    int               _iter;

  public:
    typedef DOFIterator Self;

    DOFIterator ( const ConfType &conf ) : _conf( conf ), _iter ( 0 ) {}

    template<typename EndType>
    bool operator== ( const EndType & ) const {
      return atEnd();
    }

    template<typename EndType>
    bool operator!= ( const EndType & ) const {
      return notAtEnd();
    }

    bool atEnd () const {
      return _iter == _conf.getNumGlobalDofs();
    }

    bool notAtEnd () const {
      return ! atEnd ();
    }

    Self& operator++( ) {
      ++_iter;
      return *this;
    }

    // expensive! Use only if necessary!
    Vec3<_RealType> operator*() const {
      Vec2<int> globIdx( _conf.getMapper().getGlobalPairing( _iter ) );
      Vec3<_RealType> coords( _conf.getInitializer().getVertex( globIdx[0] ) + _conf.getInitializer().getVertex( globIdx[1] ) );
      coords /= 2.;
      return coords;
    }

    int getIndex () const {
      return _iter;
    }

    aol::Vec3<_RealType> getCoords () const {
      return operator* ();
    }
  };
  
public:
  typedef _RealType RealType;
  typedef _QuadType QuadType;

  static const qc::Dimension Dim = qc::QC_3D;
  enum { DomDim = 2 }; // ???

  static const qc::Dimension DimOfWorld = qc::QC_3D; // ??? currently not used, but should be correct anyway ...

  typedef RealType Real;
  typedef MeshType                                 InitType;               //!< that's the type that is needed by the constructor of the configurator
  typedef aol::Vec<2, RealType>                        DomVecType;
  typedef aol::Vec3<RealType>                          VecType;
  typedef aol::Mat<3, 3, RealType>                     MatType;
  typedef aol::Vector<RealType>                        VectorType;
  typedef aol::Vector<RealType>                        ArrayType;
  typedef aol::SparseMatrix<RealType>                  MatrixType;
  typedef aol::BitVector                               MaskType;
  typedef typename MeshType::ElementType           ElementType;
  // different basus functions set
  typedef QuadTriangMeshBaseFunctionSet<RealType, QuadType,
            ElementType>                               BaseFuncSetType;

  typedef typename MeshType::ElementIteratorType   ElementIteratorType;
  typedef typename MeshType::NodeIteratorType      NodeIteratorType;
  typedef DOFIterator                              DOFIteratorType;

  QuadTriangMeshConfigurator ( const InitType &Mesh ) :
      _mesh ( Mesh ), _dofMapper( Mesh ) {}

  //! returns the begin iterator of the grid
  inline const MeshType & begin( ) const {
    return _mesh;
  }

  //! returns the end iterator of the grid
  inline const MeshType & end( ) const {
    return _mesh;
  }

  const InitType& getInitializer( ) const { return this->_mesh; }

  const P2Mapper<MeshType>& getMapper( ) const { return this->_dofMapper; }
  
  mutable BaseFuncSetType _baseFuncSet;

  static const int maxNumLocalDofs = 6;

  inline int getNumLocalDofs ( const ElementType & ) const {
    return 6;
  }

  int getNumGlobalDofs( ) const {
    return _dofMapper.getGlobalDofs();
  }

  int maxNumQuadPoints( ) const {
    return QuadType::numQuadPoints;
  }


  const BaseFuncSetType& getBaseFunctionSet ( const ElementType &T ) const {
    _baseFuncSet.setTriangle ( T );
    return _baseFuncSet;
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const ElementType &T, int localIndex ) const {
    return _dofMapper.localToGlobal( T.getIndex(), localIndex );
    //return T.globNodeIdx( localIndex );
  }

  //! \warning copied and pasted from QuocConfiguratorTraitMultiLin.
  //! I have absolutely no idea if this works or not.
  inline void localToGlobal ( const ElementType &T, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal ( T, localIndex0 );
    glob[1] = localToGlobal ( T, localIndex1 );
  }

  RealType vol ( const ElementType &T ) const {
    return T.area();
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    int num = getNumGlobalDofs();
    MatrixType *mat = new MatrixType ( num, num );
    return mat;
  }
  
  //! fill boundary mask
  void fillBoundaryMask( MaskType& mask ) const {
    mask.resize( _dofMapper.getGlobalDofs() );
    mask.setAll( false );
    _mesh.makeNeighbour();
    for( ElementIteratorType iter = _mesh; iter.notAtEnd(); ++iter ){
      for( int loc = 0; loc < 3; loc++ ){
        if ( _mesh.getNeighbour( iter.getIndex() , loc ) == -1 ){
                int v1( iter->globNodeIdx( (loc+1) % 3 ) ), v2( iter->globNodeIdx( (loc+2) % 3 ) );
          mask.set( v1, true);
          mask.set( v2, true);
                mask.set( _dofMapper.getEdgeDOF( min(v1,v2) , max(v1,v2) ), true );
        }
      }
    }
  }

};


//! \brief Mesh configurator for general triangular meshes with P3 reduced elements.
//!        This means, there is no degree of freedom in the center of the element
//!        Reasonable choices for MeshType are aol::FETriangMesh<> or om::Trimesh<>.
/*!
 * \author Heeren, Perl, Simon
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename MeshType, typename _QuadType>
class ReducedCubicTriangMeshConfigurator {
protected:
  const MeshType &_mesh;
  const P3redMapper<MeshType> _dofMapper;
 
public:
  // --- dof iterator ---
  class DOFIterator {

  protected:
    typedef ReducedCubicTriangMeshConfigurator<_RealType, MeshType, _QuadType> ConfType;
    const ConfType &  _conf;
    int               _iter;

  public:
    typedef DOFIterator Self;

    DOFIterator ( const ConfType &conf ) : _conf( conf ), _iter ( 0 ) {}

    template<typename EndType>
    bool operator== ( const EndType & ) const {
      return atEnd();
    }

    template<typename EndType>
    bool operator!= ( const EndType & ) const {
      return notAtEnd();
    }

    bool atEnd () const {
      return _iter == _conf.getNumGlobalDofs();
    }

    bool notAtEnd () const {
      return ! atEnd ();
    }

    Self& operator++( ) {
      ++_iter;
      return *this;
    }

    // expensive! Use only if necessary!
    Vec3<_RealType> operator*() const {
      Vec2<int> globIdx( _conf.getMapper().getGlobalPairing( _iter ) );
      Vec3<_RealType> coords( _conf.getInitializer().getVertex( globIdx[0] ) );
      if( globIdx[0] != globIdx[1] ){
        coords *= 2./3;
        coords.addMultiple( _conf.getInitializer().getVertex( globIdx[1] ), 1./3 );
      }
      return coords;
    }

    int getIndex () const {
      return _iter;
    }

    aol::Vec3<_RealType> getCoords () const {
      return operator* ();
    }
  };
  
public:
  typedef _RealType RealType;
  typedef _QuadType QuadType;

  static const qc::Dimension Dim = qc::QC_3D;
  enum { DomDim = 2 }; // ???

  static const qc::Dimension DimOfWorld = qc::QC_3D; // ??? currently not used, but should be correct anyway ...

  typedef RealType Real;
  typedef MeshType                                 InitType;               //!< that's the type that is needed by the constructor of the configurator
  typedef aol::Vec<2, RealType>                        DomVecType;
  typedef aol::Vec3<RealType>                          VecType;
  typedef aol::Mat<3, 3, RealType>                     MatType;
  typedef aol::Vector<RealType>                        VectorType;
  typedef aol::Vector<RealType>                        ArrayType;
  typedef aol::SparseMatrix<RealType>                  MatrixType;
  typedef aol::BitVector                               MaskType;
  typedef typename MeshType::ElementType           ElementType;
  // different basus functions set
  typedef ReducedCubicTriangMeshBaseFunctionSet<RealType, QuadType,
            ElementType>                               BaseFuncSetType;

  typedef typename MeshType::ElementIteratorType   ElementIteratorType;
  typedef typename MeshType::NodeIteratorType      NodeIteratorType;
  typedef DOFIterator                              DOFIteratorType;

  ReducedCubicTriangMeshConfigurator ( const InitType &Mesh ) :
      _mesh ( Mesh ), _dofMapper( Mesh ) {}

  //! returns the begin iterator of the grid
  inline const MeshType & begin( ) const {
    return _mesh;
  }

  //! returns the end iterator of the grid
  inline const MeshType & end( ) const {
    return _mesh;
  }

  const InitType& getInitializer( ) const { return this->_mesh; }

  const P3redMapper<MeshType>& getMapper( ) const { return this->_dofMapper; }
  
  mutable BaseFuncSetType _baseFuncSet;

  static const int maxNumLocalDofs = 9;

  inline int getNumLocalDofs ( const ElementType & ) const {
    return 9;
  }

  int getNumGlobalDofs( ) const {
    return _dofMapper.getGlobalDofs();
  }

  int maxNumQuadPoints( ) const {
    return QuadType::numQuadPoints;
  }


  const BaseFuncSetType& getBaseFunctionSet ( const ElementType &T ) const {
    _baseFuncSet.setTriangle ( T );
    return _baseFuncSet;
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const ElementType &T, int localIndex ) const {
    return _dofMapper.localToGlobal( T.getIndex(), localIndex );
  }

  //! \warning copied and pasted from QuocConfiguratorTraitMultiLin.
  //! I have absolutely no idea if this works or not.
  inline void localToGlobal ( const ElementType &T, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal ( T, localIndex0 );
    glob[1] = localToGlobal ( T, localIndex1 );
  }

  RealType vol ( const ElementType &T ) const {
    return T.area();
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    int num = getNumGlobalDofs();
    MatrixType *mat = new MatrixType ( num, num );
    return mat;
  }
  
  //! fill boundary mask
  void fillBoundaryMask( MaskType& mask ) const {
    mask.resize( _dofMapper.getGlobalDofs() );
    mask.setAll( false );
    _mesh.makeNeighbour();
    for( ElementIteratorType iter = _mesh; iter.notAtEnd(); ++iter ){
      for( int loc = 0; loc < 3; loc++ ){
        if ( _mesh.getNeighbour( iter.getIndex() , loc ) == -1 ){
    int v1( iter->globNodeIdx( (loc+2) % 3 ) ), v2( iter->globNodeIdx( (loc+1) % 3 ) );
          mask.set( v1, true);
          mask.set( v2, true);
          mask.set( _dofMapper.getEdgeDOF( v1, v2 ), true );
          mask.set( _dofMapper.getEdgeDOF( v2, v1 ), true );
        }
      }
    }
  }

};






class DKTPlateTriangMeshConfiguratorBase {};



//! Mesh configurator for Discrete Kirchhoff Triangle (DKT) for plates
/*!
 * \author Heeren, Perl, Simon
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename MeshType, typename _QuadType>
class DKTPlateTriangMeshConfigurator : DKTPlateTriangMeshConfiguratorBase {
protected:
  const MeshType &_mesh;

public:
  typedef _RealType RealType;
  typedef _QuadType QuadType;

  static const qc::Dimension Dim = qc::QC_3D;
  enum { DomDim = 2 }; // ???

  static const qc::Dimension DimOfWorld = qc::QC_3D; // ??? currently not used, but should be correct anyway ...

  typedef RealType Real;
  typedef MeshType                        InitType;  //!< that's the type that is needed by the constructor of the configurator
  typedef aol::Vec<2, RealType>           DomVecType;
  typedef aol::Vec3<RealType>             VecType;
  typedef aol::Mat<3, 3, RealType>        MatType;
  typedef aol::Vector<RealType>           VectorType;
  typedef aol::Vector<RealType>           ArrayType;
  typedef aol::SparseMatrix<RealType>     MatrixType;
  typedef aol::BitVector                  MaskType;
  typedef typename MeshType::ElementType  ElementType;
  //typedef aol::DKTPlateElement<RealType> ElementType;
  
  // different basis functions set
  typedef DKTPlateGradientBaseFunctionSet<RealType, QuadType> BaseFuncSetType;

  typedef typename MeshType::ElementIteratorType   ElementIteratorType;
  typedef typename MeshType::NodeIteratorType      NodeIteratorType;
  
  class DOFIterator : public NodeIteratorType {
    typedef DKTPlateTriangMeshConfigurator<_RealType, MeshType, _QuadType> ConfType;
  public:
    DOFIterator (const ConfType& conf ) : NodeIteratorType( conf.getInitializer() ){
      throw aol::Exception ( "DOFIterator in DKTPlateTriangMeshConfigurator is not implemented and returns only NodeIterator", __FILE__, __LINE__ );
    }
  };

  typedef DOFIterator DOFIteratorType;
  

  DKTPlateTriangMeshConfigurator ( const InitType &Mesh ) :
      _mesh ( Mesh ) {}

  //! returns the begin iterator of the grid
  inline const MeshType & begin( ) const {
    return _mesh;
  }

  //! returns the end iterator of the grid
  inline const MeshType & end( ) const {
    return _mesh;
  }

  const InitType& getInitializer( ) const { return this->_mesh; }
  
  mutable BaseFuncSetType _baseFuncSet;

  static const int maxNumLocalDofs = 9;

  inline int getNumLocalDofs ( const ElementType & ) const {
    return 9;
  }

  int getNumGlobalDofs( ) const {
    return 3 * _mesh.getNumVertices();
  }

  int maxNumQuadPoints( ) const {
    return QuadType::numQuadPoints;
  }

  const BaseFuncSetType& getBaseFunctionSet ( const ElementType &T ) const {
    _baseFuncSet.setTriangle ( T );
    return _baseFuncSet;
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const ElementType &T, int localIndex ) const {
    int aux = ( localIndex%3 == 0) ? 0 : (  ( (localIndex+2)%3 == 0) ? 1 : 2 );
    return aux * _mesh.getNumVertices() + T.globNodeIdx( std::floor(localIndex/3) );
  }

//   //! \warning copied and pasted from QuocConfiguratorTraitMultiLin.
//   //! I have absolutely no idea if this works or not.
//   inline void localToGlobal ( const ElementType &T, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
//     glob[0] = localToGlobal ( T, localIndex0 );
//     glob[1] = localToGlobal ( T, localIndex1 );
//   }

  RealType vol ( const ElementType &T ) const {
    return T.area();
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    int num = getNumGlobalDofs();
    MatrixType *mat = new MatrixType ( num, num );
    return mat;
  }
  
  //! fill boundary mask
  void fillBoundaryMaskClamped( MaskType& mask ) const {
    mask.resize( 3 * _mesh.getNumVertices() );
    mask.setAll( false );
    _mesh.makeNeighbour();
    for( ElementIteratorType iter = _mesh; iter.notAtEnd(); ++iter ){
      for( int loc = 0; loc < 3; loc++ ){
        if ( _mesh.getNeighbour( iter.getIndex() , loc ) == -1 ){
          int boundVertexIdx1( iter->globNodeIdx( (loc+1) % 3 ) ), boundVertexIdx2( iter->globNodeIdx( (loc+2) % 3 ) );
          for ( int i = 0; i < 3; i++ ){
            mask.set( i * _mesh.getNumVertices() + boundVertexIdx1, true);
            mask.set( i * _mesh.getNumVertices() + boundVertexIdx2, true);
          }
        }
      }
    }
  }
  
  //! fill boundary mask
  void fillBoundaryMaskSimplySupported( MaskType& mask ) const {
    mask.resize( 3 * _mesh.getNumVertices() );
    mask.setAll( false );
    _mesh.makeNeighbour();
    for( ElementIteratorType iter = _mesh; iter.notAtEnd(); ++iter ){
      for( int loc = 0; loc < 3; loc++ ){
        if ( _mesh.getNeighbour( iter.getIndex() , loc ) == -1 ){
          int boundVertexIdx1( iter->globNodeIdx( (loc+1) % 3 ) ), boundVertexIdx2( iter->globNodeIdx( (loc+2) % 3 ) );
          // w at the two boundary vertices boundVertexIdx1 and boundVertexIdx2
          mask.set( boundVertexIdx1, true);
          mask.set( boundVertexIdx2, true);
        }
      }
    }
  }

};























//! Basefunction for DKT-Funktions = \{ w \in P3Red , Dw continuous at nodes \}
/*!
 * \author Heeren, Perl, Simon
 */
template <typename RealType, typename QuadRuleType>
class DKTPlateBaseFunctionSet : public aol::BaseFunctionSetInterface<RealType, aol::Vec3<RealType>, 
  aol::Vec2<RealType>, 9, QuadRuleType, TriangMeshBaseFunctionSet<RealType, QuadRuleType, aol::DKTPlateElement<RealType> > > {
  
protected:
  typedef aol::DKTPlateElement<RealType> TriangleType;
  typedef aol::ReducedCubicTriangMeshBaseFunctionSet<RealType, QuadRuleType, TriangleType> P3RedShapeFnc;  
 
  // ordering anti-clockwise beginning at node zero, i.e. _triangle-> globNodeIdx( 0 )
  static RealType _w1( const TriangleType* /*t*/, const aol::Vec2<RealType> &c ) { 
    return ( 7.*P3RedShapeFnc::_b8(c) + 20.*P3RedShapeFnc::_b9(c) + 27.*P3RedShapeFnc::_b1(c) + 20.*P3RedShapeFnc::_b2(c) + 7.*P3RedShapeFnc::_b3(c) ) / 27.;     
  }
  
  static RealType _w2( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(2)[0]*P3RedShapeFnc::_b2(c) + 2.*t->edge(2)[0]*P3RedShapeFnc::_b3(c) - 2.*t->edge(1)[0]*P3RedShapeFnc::_b8(c) - 4.*t->edge(1)[0]*P3RedShapeFnc::_b9(c) ) / 27.;     
  }
  
  static RealType _w3( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(2)[1]*P3RedShapeFnc::_b2(c) + 2.*t->edge(2)[1]*P3RedShapeFnc::_b3(c) - 2.*t->edge(1)[1]*P3RedShapeFnc::_b8(c) - 4.*t->edge(1)[1]*P3RedShapeFnc::_b9(c) ) / 27.;     
  }
  
  static RealType _w4( const TriangleType* /*t*/, const aol::Vec2<RealType> &c ) { 
    return ( 7.*P3RedShapeFnc::_b2(c) + 20.*P3RedShapeFnc::_b3(c) + 27.*P3RedShapeFnc::_b4(c) + 20.*P3RedShapeFnc::_b5(c) + 7.*P3RedShapeFnc::_b6(c) ) / 27.;     
  }  
    
  static RealType _w5( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(0)[0]*P3RedShapeFnc::_b5(c) + 2.*t->edge(0)[0]*P3RedShapeFnc::_b6(c) - 2.*t->edge(2)[0]*P3RedShapeFnc::_b2(c) - 4.*t->edge(2)[0]*P3RedShapeFnc::_b3(c) ) / 27.;     
  }
  
  static RealType _w6( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(0)[1]*P3RedShapeFnc::_b5(c) + 2.*t->edge(0)[1]*P3RedShapeFnc::_b6(c) - 2.*t->edge(2)[1]*P3RedShapeFnc::_b2(c) - 4.*t->edge(2)[1]*P3RedShapeFnc::_b3(c) ) / 27.;     
  }  
    
  static RealType _w7( const TriangleType* /*t*/, const aol::Vec2<RealType> &c ) { 
    return ( 7.*P3RedShapeFnc::_b5(c) + 20.*P3RedShapeFnc::_b6(c) + 27.*P3RedShapeFnc::_b7(c) + 20.*P3RedShapeFnc::_b8(c) + 7.*P3RedShapeFnc::_b9(c) ) / 27.;     
  }  
    
  static RealType _w8( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(1)[0]*P3RedShapeFnc::_b8(c) + 2.*t->edge(1)[0]*P3RedShapeFnc::_b9(c) - 2.*t->edge(0)[0]*P3RedShapeFnc::_b5(c) - 4.*t->edge(0)[0]*P3RedShapeFnc::_b6(c) ) / 27.;     
  }
  
  static RealType _w9( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(1)[1]*P3RedShapeFnc::_b8(c) + 2.*t->edge(1)[1]*P3RedShapeFnc::_b9(c) - 2.*t->edge(0)[1]*P3RedShapeFnc::_b5(c) - 4.*t->edge(0)[1]*P3RedShapeFnc::_b6(c) ) / 27.;     
  }
  
  // derivative with respect to \xi
  static RealType _d1_w1( const TriangleType* /*t*/, const aol::Vec2<RealType> &c ) { 
    return ( 7.*P3RedShapeFnc::_d1_b8(c) + 20.*P3RedShapeFnc::_d1_b9(c) + 27.*P3RedShapeFnc::_d1_b1(c) + 20.*P3RedShapeFnc::_d1_b2(c) + 7.*P3RedShapeFnc::_d1_b3(c) ) / 27.;     
  }
  
  static RealType _d1_w2( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(2)[0]*P3RedShapeFnc::_d1_b2(c) + 2.*t->edge(2)[0]*P3RedShapeFnc::_d1_b3(c) - 2.*t->edge(1)[0]*P3RedShapeFnc::_d1_b8(c) - 4.*t->edge(1)[0]*P3RedShapeFnc::_d1_b9(c) ) / 27.;     
  }
  
  static RealType _d1_w3( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(2)[1]*P3RedShapeFnc::_d1_b2(c) + 2.*t->edge(2)[1]*P3RedShapeFnc::_d1_b3(c) - 2.*t->edge(1)[1]*P3RedShapeFnc::_d1_b8(c) - 4.*t->edge(1)[1]*P3RedShapeFnc::_d1_b9(c) ) / 27.;     
  }
  
  static RealType _d1_w4( const TriangleType* /*t*/, const aol::Vec2<RealType> &c ) { 
    return ( 7.*P3RedShapeFnc::_d1_b2(c) + 20.*P3RedShapeFnc::_d1_b3(c) + 27.*P3RedShapeFnc::_d1_b4(c) + 20.*P3RedShapeFnc::_d1_b5(c) + 7.*P3RedShapeFnc::_d1_b6(c) ) / 27.;     
  }  
    
  static RealType _d1_w5( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(0)[0]*P3RedShapeFnc::_d1_b5(c) + 2.*t->edge(0)[0]*P3RedShapeFnc::_d1_b6(c) - 2.*t->edge(2)[0]*P3RedShapeFnc::_d1_b2(c) - 4.*t->edge(2)[0]*P3RedShapeFnc::_d1_b3(c) ) / 27.;     
  }
  
  static RealType _d1_w6( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(0)[1]*P3RedShapeFnc::_d1_b5(c) + 2.*t->edge(0)[1]*P3RedShapeFnc::_d1_b6(c) - 2.*t->edge(2)[1]*P3RedShapeFnc::_d1_b2(c) - 4.*t->edge(2)[1]*P3RedShapeFnc::_d1_b3(c) ) / 27.;     
  }  
    
  static RealType _d1_w7( const TriangleType* /*t*/, const aol::Vec2<RealType> &c ) { 
    return ( 7.*P3RedShapeFnc::_d1_b5(c) + 20.*P3RedShapeFnc::_d1_b6(c) + 27.*P3RedShapeFnc::_d1_b7(c) + 20.*P3RedShapeFnc::_d1_b8(c) + 7.*P3RedShapeFnc::_d1_b9(c) ) / 27.;     
  }  
    
  static RealType _d1_w8( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(1)[0]*P3RedShapeFnc::_d1_b8(c) + 2.*t->edge(1)[0]*P3RedShapeFnc::_d1_b9(c) - 2.*t->edge(0)[0]*P3RedShapeFnc::_d1_b5(c) - 4.*t->edge(0)[0]*P3RedShapeFnc::_d1_b6(c) ) / 27.;     
  }
  
  static RealType _d1_w9( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(1)[1]*P3RedShapeFnc::_d1_b8(c) + 2.*t->edge(1)[1]*P3RedShapeFnc::_d1_b9(c) - 2.*t->edge(0)[1]*P3RedShapeFnc::_d1_b5(c) - 4.*t->edge(0)[1]*P3RedShapeFnc::_d1_b6(c) ) / 27.;     
  }
  
 // derivative with respect to \eta
  static RealType _d2_w1( const TriangleType* /*t*/, const aol::Vec2<RealType> &c ) { 
    return ( 7.*P3RedShapeFnc::_d2_b8(c) + 20.*P3RedShapeFnc::_d2_b9(c) + 27.*P3RedShapeFnc::_d2_b1(c) + 20.*P3RedShapeFnc::_d2_b2(c) + 7.*P3RedShapeFnc::_d2_b3(c) ) / 27.;     
  }
  
  static RealType _d2_w2( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(2)[0]*P3RedShapeFnc::_d2_b2(c) + 2.*t->edge(2)[0]*P3RedShapeFnc::_d2_b3(c) - 2.*t->edge(1)[0]*P3RedShapeFnc::_d2_b8(c) - 4.*t->edge(1)[0]*P3RedShapeFnc::_d2_b9(c) ) / 27.;     
  }
  
  static RealType _d2_w3( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(2)[1]*P3RedShapeFnc::_d2_b2(c) + 2.*t->edge(2)[1]*P3RedShapeFnc::_d2_b3(c) - 2.*t->edge(1)[1]*P3RedShapeFnc::_d2_b8(c) - 4.*t->edge(1)[1]*P3RedShapeFnc::_d2_b9(c) ) / 27.;     
  }
  
  static RealType _d2_w4( const TriangleType* /*t*/, const aol::Vec2<RealType> &c ) { 
    return ( 7.*P3RedShapeFnc::_d2_b2(c) + 20.*P3RedShapeFnc::_d2_b3(c) + 27.*P3RedShapeFnc::_d2_b4(c) + 20.*P3RedShapeFnc::_d2_b5(c) + 7.*P3RedShapeFnc::_d2_b6(c) ) / 27.;     
  }  
    
  static RealType _d2_w5( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(0)[0]*P3RedShapeFnc::_d2_b5(c) + 2.*t->edge(0)[0]*P3RedShapeFnc::_d2_b6(c) - 2.*t->edge(2)[0]*P3RedShapeFnc::_d2_b2(c) - 4.*t->edge(2)[0]*P3RedShapeFnc::_d2_b3(c) ) / 27.;     
  }
  
  static RealType _d2_w6( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(0)[1]*P3RedShapeFnc::_d2_b5(c) + 2.*t->edge(0)[1]*P3RedShapeFnc::_d2_b6(c) - 2.*t->edge(2)[1]*P3RedShapeFnc::_d2_b2(c) - 4.*t->edge(2)[1]*P3RedShapeFnc::_d2_b3(c) ) / 27.;     
  }  
    
  static RealType _d2_w7( const TriangleType* /*t*/, const aol::Vec2<RealType> &c ) { 
    return ( 7.*P3RedShapeFnc::_d2_b5(c) + 20.*P3RedShapeFnc::_d2_b6(c) + 27.*P3RedShapeFnc::_d2_b7(c) + 20.*P3RedShapeFnc::_d2_b8(c) + 7.*P3RedShapeFnc::_d2_b9(c) ) / 27.;     
  }  
    
  static RealType _d2_w8( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(1)[0]*P3RedShapeFnc::_d2_b8(c) + 2.*t->edge(1)[0]*P3RedShapeFnc::_d2_b9(c) - 2.*t->edge(0)[0]*P3RedShapeFnc::_d2_b5(c) - 4.*t->edge(0)[0]*P3RedShapeFnc::_d2_b6(c) ) / 27.;     
  }
  
  static RealType _d2_w9( const TriangleType* t, const aol::Vec2<RealType> &c ) { 
    return ( 4.*t->edge(1)[1]*P3RedShapeFnc::_d2_b8(c) + 2.*t->edge(1)[1]*P3RedShapeFnc::_d2_b9(c) - 2.*t->edge(0)[1]*P3RedShapeFnc::_d2_b5(c) - 4.*t->edge(0)[1]*P3RedShapeFnc::_d2_b6(c) ) / 27.;     
  }  

  //
  typedef RealType ( *BASIS_FUNC_TYPE ) ( const TriangleType* Triangle, const aol::Vec2<RealType> &RefCoord );
  BASIS_FUNC_TYPE _dbasis[2][9];
  BASIS_FUNC_TYPE _basis[9];

  const TriangleType *_triangle;

public:
  DKTPlateBaseFunctionSet(  ) : _triangle ( NULL ) {
    _basis[0] = _w1;
    _basis[1] = _w2;
    _basis[2] = _w3;
    _basis[3] = _w4;
    _basis[4] = _w5;
    _basis[5] = _w6;
    _basis[6] = _w7;
    _basis[7] = _w8;
    _basis[8] = _w9;

    _dbasis[0][0] = _d1_w1;
    _dbasis[0][1] = _d1_w2;
    _dbasis[0][2] = _d1_w3;
    _dbasis[0][3] = _d1_w4;
    _dbasis[0][4] = _d1_w5;
    _dbasis[0][5] = _d1_w6;
    _dbasis[0][6] = _d1_w7;
    _dbasis[0][7] = _d1_w8;
    _dbasis[0][8] = _d1_w9;

    _dbasis[1][0] = _d2_w1;
    _dbasis[1][1] = _d2_w2;
    _dbasis[1][2] = _d2_w3;
    _dbasis[1][3] = _d2_w4;
    _dbasis[1][4] = _d2_w5;
    _dbasis[1][5] = _d2_w6;
    _dbasis[1][6] = _d2_w7;
    _dbasis[1][7] = _d2_w8;
    _dbasis[1][8] = _d2_w9;
  }

  enum { numBaseFuncs = 9 };

  void setTriangle ( const TriangleType &T ) {
    _triangle = &T;
  }
  
  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return evaluate ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );
  }

  inline aol::Vec3<RealType> evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    aol::Vec3<RealType> g;
    evaluateGradient ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), g );
    return g;
  }
  
  RealType evaluate ( int BaseFuncNum, const aol::Vec2<RealType> &RefCoord ) const {
    return _basis[BaseFuncNum] ( _triangle, RefCoord );
  }

  void evaluateGradient ( int BaseFuncNum, const aol::Vec2<RealType> &RefCoord, aol::Vec3<RealType> &Gradient ) const {
    //
    aol::Vec2<RealType> tmp, tmp2;
    tmp[0] = _dbasis[0][BaseFuncNum] ( _triangle, RefCoord );
    tmp[1] = _dbasis[1][BaseFuncNum] ( _triangle, RefCoord );
    _triangle->ginv().mult ( tmp, tmp2 );

    for ( int i = 0; i < 2; i++ )
      Gradient[i] = tmp2[0] * _triangle->edge(2)[i] - tmp2[1] * _triangle->edge(1)[i];
    Gradient[2] = 0.;
  }
};





//! Mesh configurator for Discrete Kirchhoff Triangle (DKT) for plates
/*!
 * \author Heeren, Perl, Simon
 */
template <typename _RealType, typename MeshType, typename _QuadType>
class FullDKTPlateTriangMeshConfigurator : DKTPlateTriangMeshConfiguratorBase {
protected:
  const MeshType &_mesh;

public:
  typedef _RealType RealType;
  typedef _QuadType QuadType;

  static const qc::Dimension Dim = qc::QC_3D;
  enum { DomDim = 2 }; // ???

  static const qc::Dimension DimOfWorld = qc::QC_3D; // ??? currently not used, but should be correct anyway ...

  typedef RealType Real;
  typedef MeshType                        InitType;  //!< that's the type that is needed by the constructor of the configurator
  typedef aol::Vec<2, RealType>           DomVecType;
  typedef aol::Vec3<RealType>             VecType;
  typedef aol::Mat<3, 3, RealType>        MatType;
  typedef aol::Vector<RealType>           VectorType;
  typedef aol::Vector<RealType>           ArrayType;
  typedef aol::SparseMatrix<RealType>     MatrixType;
  typedef aol::BitVector                  MaskType;
  typedef typename MeshType::ElementType  ElementType;
  
  typedef DKTPlateBaseFunctionSet<RealType, QuadType> BaseFuncSetType;
  typedef DKTPlateGradientBaseFunctionSet<RealType, QuadType> ApproxGradientBaseFuncSetType;

  typedef typename MeshType::ElementIteratorType   ElementIteratorType;
  typedef typename MeshType::NodeIteratorType      NodeIteratorType;

  class DOFIterator : public NodeIteratorType {
    typedef FullDKTPlateTriangMeshConfigurator<_RealType, MeshType, _QuadType> ConfType;
  public:
    DOFIterator (const ConfType& conf ) : NodeIteratorType( conf.getInitializer() ){
      throw aol::Exception ( "DOFIterator in FullDKTPlateTriangMeshConfigurator is not implemented and returns only NodeIterator", __FILE__, __LINE__ );
    }
  };

  typedef DOFIterator DOFIteratorType;
  
  
  FullDKTPlateTriangMeshConfigurator ( const InitType &Mesh ) :
      _mesh ( Mesh ) {}

  //! returns the begin iterator of the grid
  inline const MeshType & begin( ) const {
    return _mesh;
  }

  //! returns the end iterator of the grid
  inline const MeshType & end( ) const {
    return _mesh;
  }

  const InitType& getInitializer( ) const { return this->_mesh; }
  
  mutable BaseFuncSetType _baseFuncSet;
  mutable ApproxGradientBaseFuncSetType _approxGradBaseFuncSet;

  static const int maxNumLocalDofs = 9;

  inline int getNumLocalDofs ( const ElementType & ) const {
    return 9;
  }

  int getNumGlobalDofs( ) const {
    return 3 * _mesh.getNumVertices();
  }

  int maxNumQuadPoints( ) const {
    return QuadType::numQuadPoints;
  }

  const BaseFuncSetType& getBaseFunctionSet ( const ElementType &T ) const {
    _baseFuncSet.setTriangle ( T );
    return _baseFuncSet;
  }
  
  const ApproxGradientBaseFuncSetType& getApproxGradientBaseFunctionSet ( const ElementType &T ) const {
    _approxGradBaseFuncSet.setTriangle ( T );
    return _approxGradBaseFuncSet;
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const ElementType &T, int localIndex ) const {
    int aux = ( localIndex%3 == 0) ? 0 : (  ( (localIndex+2)%3 == 0) ? 1 : 2 );
    return aux * _mesh.getNumVertices() + T.globNodeIdx( std::floor(localIndex/3) );
  }
  
  RealType vol ( const ElementType &T ) const {
    return T.area();
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    int num = getNumGlobalDofs();
    MatrixType *mat = new MatrixType ( num, num );
    return mat;
  }
  
  //! fill boundary mask
  void fillBoundaryMaskClamped( MaskType& mask ) const {
    mask.resize( 3 * _mesh.getNumVertices() );
    mask.setAll( false );
    _mesh.makeNeighbour();
    for( ElementIteratorType iter = _mesh; iter.notAtEnd(); ++iter ){
      for( int loc = 0; loc < 3; loc++ ){
        if ( _mesh.getNeighbour( iter.getIndex() , loc ) == -1 ){
          int boundVertexIdx1( iter->globNodeIdx( (loc+1) % 3 ) ), boundVertexIdx2( iter->globNodeIdx( (loc+2) % 3 ) );
          for ( int i = 0; i < 3; i++ ){
            mask.set( i * _mesh.getNumVertices() + boundVertexIdx1, true);
            mask.set( i * _mesh.getNumVertices() + boundVertexIdx2, true);
          }
        }
      }
    }
  }
  
  //! fill boundary mask
  void fillBoundaryMaskSimplySupported( MaskType& mask ) const {
    mask.resize( 3 * _mesh.getNumVertices() );
    mask.setAll( false );
    _mesh.makeNeighbour();
    for( ElementIteratorType iter = _mesh; iter.notAtEnd(); ++iter ){
      for( int loc = 0; loc < 3; loc++ ){
        if ( _mesh.getNeighbour( iter.getIndex() , loc ) == -1 ){
          int boundVertexIdx1( iter->globNodeIdx( (loc+1) % 3 ) ), boundVertexIdx2( iter->globNodeIdx( (loc+2) % 3 ) );
          // w at the two boundary vertices boundVertexIdx1 and boundVertexIdx2
          mask.set( boundVertexIdx1, true);
          mask.set( boundVertexIdx2, true);
        }
      }
    }
  }

};












} // namespace aol

#endif // __TRIANGMESHCONFIGURATORS_H
