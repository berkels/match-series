#ifndef __ADAPTIVETRIANGMESH_H
#define __ADAPTIVETRIANGMESH_H

#include <aol.h>
#include <vec.h>
#include <math.h>
#include <triangMeshConfigurators.h>
#include <triangMesh.h>


//#define DEBUGMODE

//! defines for local node and edge numbering
typedef short LocalIndex;
typedef int GlobalIndex;
const int IndexNotSet = -1;

//! different possible plot formats for AdaptiveFETriangMesh<>::plot2DProjection().
enum FormatType{ MATLAB, GNUPLOT, POSTSCRIPT, JPG, PLY };



/*! \brief Class for adaptive triangular meshes (refinement realized by bisection algorithm)
 *  \author Geihe, Heeren, vDeylen
 *  \todo atm makeOrientationConsistent is called after each refinement operation
 *  \todo gnuplotting should be done quocmesh style
 * 
 *  Works just as aol::FETriangMesh<>, ie. all classical FE operation are provided.
 * 
 *  Additional functions for refinement (via bisection): 
 *   - triangle (with global index t) to be refined has to be marked via AdaptiveFETriangMesh<>::mark( t ).
 *   - after marking all desired triangles, call AdaptiveFETriangMesh<>::refineMarkedTriangles().
 *   - if you want to refine all triangles simultaneously, call AdaptiveFETriangMesh<>::refineAll()
 *   - manual unmarking via AdaptiveFETriangMesh<>::unmark( t ) or AdaptiveFETriangMesh<>::unmarkAll( )
 * 
 *  Features: 
 *   - prolongation possible via linear interpolation (each node created by refinement is an edge midpoint, hence function value at new node is interpolation between the edge endpoints)
 *   - can plot a 2D-projection of the mesh with possibly showing neighbouring relations and marked elements in various FormatTypes (MATLAB, GNUPLOT, POSTSCRIPT, JPG, PLY)
 *   - provides a protected class DartIterator to simulate a half-edge data structure for efficient grid navigation.
 *  
 */
template< typename RealType, typename TriangleType = aol::TriangBaseElement<RealType> >
class AdaptiveFETriangMesh : public aol::FETriangMesh<RealType, TriangleType> {
  
public:    
  typedef typename aol::FETriangMesh<RealType, TriangleType>::NodeIteratorType NodeIteratorType;
  typedef typename aol::FETriangMesh<RealType, TriangleType>::ElementIteratorType ElementIteratorType;
  
protected:  
// half-edge iterator (only for internal refinement!)
class DartIterator{
  
   typedef AdaptiveFETriangMesh<RealType, TriangleType> MeshType;
    
  protected:
    const MeshType & _grid;
    GlobalIndex _triangle;
    LocalIndex  _node;
    LocalIndex  _edge;
	
  public:
    DartIterator( const MeshType & Grid, GlobalIndex triangleIndex, LocalIndex localEdgeIndex, LocalIndex localNodeIndex ) : _grid( Grid ), _triangle (triangleIndex), _node(localNodeIndex), _edge(localEdgeIndex){}    
    DartIterator( const MeshType & Grid, GlobalIndex triangleIndex, LocalIndex localEdgeIndex ) : _grid( Grid ), _triangle (triangleIndex), _node((localEdgeIndex + 1) % 3), _edge(localEdgeIndex){}

    void set ( const GlobalIndex triangle, const LocalIndex edge, const LocalIndex node )  {
      _triangle = triangle;
      _node     = node;
      _edge     = edge;
    }

    GlobalIndex getGlobalTriangleIndex() const{
      return _triangle;
    }

    GlobalIndex getGlobalNodeIndex( ) const{
      return _grid.getTriangNodeIdx( _triangle, _node );
    }

    LocalIndex getLocalNodeIndex() const {
      return _node;
    }

    LocalIndex getLocalEdgeIndex() const {
      return _edge;
    }

    LocalIndex  getNextNodeLocalIndex() const{
      return 3 - (_node + _edge);
    }

    LocalIndex  getNextEdgeLocalIndex() const{
      return 3 - (_node + _edge);
    }  

    GlobalIndex getNextNodeGlobalIndex() const {
      return _grid.getTriangNodeIdx( _triangle, getNextNodeLocalIndex() );
    }
    
    GlobalIndex getNextTriangleIndex() const {
      return _grid.getNeighbour( _triangle, _edge );
    }
 
    // returns local index of common node (first checking if "d" also refers to _triangle, but to a different node)
    template<typename DartIteratorType>
    LocalIndex getCommonNodeLocalIndex( const DartIteratorType & d ) const {
      if ( _triangle !=  d.getGlobalTriangleIndex() ){
        stringstream err_str;
        err_str << "DartIterator::getCommonNodeLocalIndex(): Passed DartIterator points "
                   "to triangle " << d.getGlobalTriangleIndex() << ", but this->triangle is " << _triangle << ".";
        throw aol::Exception ( err_str.str(), __FILE__, __LINE__ ); 
      }
      if ( _edge == d.getLocalEdgeIndex() ){
        stringstream err_str;
        err_str << "DartIterator::getCommonNodeLocalIndex(): Passed DartIterator and *this "
                   "point to the same edge with local GlobalIndex " << _edge << ".";
        throw aol::Exception ( err_str.str(), __FILE__, __LINE__ ); 
      }
      return 3 - ( getLocalEdgeIndex() + d.getLocalEdgeIndex() );
    }

    // returns global index of common node (also if "d" refers to a different triangle than _triangle )
    template<typename DartIteratorType>
    GlobalIndex getCommonNodeGlobalIndex( const DartIteratorType & d ) const {
      if ( getGlobalNodeIndex() == d.getGlobalNodeIndex() || getGlobalNodeIndex() == d.getNextNodeGlobalIndex() )
        return getGlobalNodeIndex();

      if ( getNextNodeGlobalIndex() == d.getGlobalNodeIndex() || getNextNodeGlobalIndex() == d.getNextNodeGlobalIndex() )
        return getNextNodeGlobalIndex();

      stringstream err_str;
      err_str << "DartIterator::getCommonNodeGlobalIndex(): Passed Dart d = " << d << " and *this = " << *this << " have no node in common.";
      throw aol::Exception ( err_str.str(), __FILE__, __LINE__ ); 
      return -1;
    }
    
    void print() const{
      cerr << "triangle = " << _triangle << ", node = " << _node << ", edge = " << _edge << endl;
    }
	

    // does _triangle have a neighbour across _edge?
    bool canFlipTriangle() const {
      return this->getNextTriangleIndex() != IndexNotSet ;
    }

    // moves to neighbouring triangle across _edge
    void flipTriangle();
    
    // moves to other node along _edge (inside current triangle)
    void flipNode() {
      _node = getNextNodeLocalIndex();
    }
    
    // moves to other edge with _node (inside current triangle)
    void flipEdge() {
      _edge = getNextEdgeLocalIndex();
    }

};

public:
  //! boundary iterator (using DartIterator to move along edges)
  //! @note starting from triangle 0 a boundary edge is searched and then followed; if there are several boundaries only one of them will be visited
  class BoundaryIterator
  {
  public:
    typedef AdaptiveFETriangMesh<RealType, TriangleType> MeshType;

    BoundaryIterator (const MeshType & grid);

    bool notAtEnd() const;
    BoundaryIterator & operator++();
    BoundaryIterator & operator--();

    GlobalIndex getStartIndex() const;
    GlobalIndex getEndIndex() const;

  protected:
    const MeshType & _grid;
    DartIterator     _dart;
    GlobalIndex      _startIndex;
    bool             _initialized;
  };

protected:  
  // store if faces are marked for refinement or not  
  aol::RandomAccessContainer<bool> _markedForRefinement;
  // each node n created by refinement is an edge midpoint, ie. the mean value of two nodes n_1 and n_2; _interpolationMap: n /mapsto (n_1, n_2)
  std::map< int, aol::Vec2<int> > _interpolationMap;
  
public:
  AdaptiveFETriangMesh() : aol::FETriangMesh<RealType, TriangleType >() {
    this->makeNeighbour();
  } 
  
  AdaptiveFETriangMesh ( string filename ) : aol::FETriangMesh<RealType, TriangleType >( filename ) , _markedForRefinement( this->getNumFaces() ) {
    this->makeNeighbour();
    unmarkAll();
  }   
  
  // mark for refinement
  void mark( int element ){
    if( !( element < this->getNumFaces() ) )
      throw aol::Exception ( "AdaptiveFETriangMesh<>::mark(): out of range!", __FILE__, __LINE__ );    
    _markedForRefinement[ element ] = true;
  }
  
  // unmark
  void unmark( int element ){  
    _markedForRefinement[ element ] = false;
  }
  
  // unmark all elements
  void unmarkAll(){
    for( int i = 0; i < _markedForRefinement.size(); i++ )
      _markedForRefinement[i] = false;
  }
  
  bool isMarkedForRefinement( int element ) const {
    return _markedForRefinement[element];
  }
  
  // when adding a new face we have to expand the array _markedForRefinement correspondingly
  int pushBackTriang ( const aol::Vec3<int> newTriang ) {
    _markedForRefinement.pushBack( false );
    return aol::FETriangMesh<RealType, TriangleType >::pushBackTriang( newTriang );
  }
  
  // refinement
  void subdivide();
  void refineAll();
  void refineMarkedTriangles();
  // prolongate by linear interpolation
  void prolongateLinearly( aol::Vector<RealType>& function ) const;
  void prolongateConst( aol::Vector<RealType>& function, const RealType constant = 0.0 ) const;
  
  // plot projection of *this to (x,y)-plane in desired format, possibly with encoding neighbouring relations and refinement markers
  void plot2DProjection( string filename = "output", FormatType format = POSTSCRIPT, bool neighbours = false, bool marked = false, double relWidthHeight = 1.0 ) const;  

protected:
  // internal refinement functions
  GlobalIndex refine( GlobalIndex triangleToRefine );
  GlobalIndex refine( const DartIterator& d);
  GlobalIndex refineOnlyThis( const DartIterator& d, GlobalIndex midpoint);
  // get longest edge index (starting search possibly with a preferred edge)
  LocalIndex getLongestEdgeIndex( GlobalIndex triangle ) const;
  // add edge midpoint on d.edge, update _interpolationMap and return global index of node
  GlobalIndex addEdgeMidpoint( const DartIterator& d );

};

#endif
