#ifndef __TRIANGLE_H
#define __TRIANGLE_H

#include <smallVec.h>
#include <smallMat.h>
#include <aol.h>

namespace aol {

// forward declaration
/*
template <typename DataType> class DefaultTriangleDataStorage;
template <typename DataType, typename StorageType = DefaultTriangleDataStorage<DataType> > class Triangle;
template <typename DataType, typename _TriangleType = Triangle<DataType> > class TriangMesh;
*/

/** Default storage class for data of triangles in 3D
 *  stores nodal positions only.
 *  \author Heeren
 */
template<typename DataType>
class DefaultTriangleDataStorage {

protected:
  Vec3<DataType> _nodes[3];

public:
  DefaultTriangleDataStorage(){}

  DefaultTriangleDataStorage ( const Vec3<DataType> &Node0, const Vec3<DataType> &Node1, const Vec3<DataType> &Node2 ) {
    _nodes[0] = Node0;  _nodes[1] = Node1;  _nodes[2] = Node2;
  }

  DefaultTriangleDataStorage ( const DefaultTriangleDataStorage<DataType> &rhs ) {
    _nodes[0] = rhs.getNode(0); _nodes[1] =  rhs.getNode(1); _nodes[2] =  rhs.getNode(2);
  }

  //! Copy the TriangleNum-th triangle from the TriangMesh argument.
  template <typename MeshType>
  DefaultTriangleDataStorage ( const MeshType &Mesh, const int TriangleNum ){
    for ( int i = 0; i < 3; ++i ) {
      const int idx = Mesh.getTriangNodeIdx ( TriangleNum, i );
      for ( int j = 0; j < 3; ++j ) {
        _nodes[i][j] = Mesh.getVertexCoord ( idx, j );
      }
    }
  }

  DefaultTriangleDataStorage& operator= ( const DefaultTriangleDataStorage<DataType> &rhs ) {
    _nodes[0] = rhs.getNode(0); _nodes[1] =  rhs.getNode(1); _nodes[2] =  rhs.getNode(2);
    return *this;
  }

  const Vec3<DataType>& getNode ( int i ) const {
    return _nodes[i];
  }

  Vec3<DataType>& getNode ( int i ) {
    return _nodes[i];
  }
  
  void setNode ( int i, const Vec3<DataType>& node ) {
     _nodes[i] = node;
   }

};

/** Pointer storage class for data of triangles in 3D
  *  stores pointers to nodal positions only.
  *  \author Heeren, Perl
  */
template<typename DataType>
class PointerTriangleDataStorage {

protected:
   const Vec3<DataType>* _nodes[3];

public:
   PointerTriangleDataStorage(){}

   PointerTriangleDataStorage ( const Vec3<DataType> &Node0, const Vec3<DataType> &Node1, const Vec3<DataType> &Node2 ) {
     _nodes[0] = &Node0;  _nodes[1] = &Node1;  _nodes[2] = &Node2;
   }

   PointerTriangleDataStorage ( const PointerTriangleDataStorage<DataType> &rhs ) {
     _nodes[0] = &rhs.getNode(0); _nodes[1] = &rhs.getNode(1); _nodes[2] =  &rhs.getNode(2);
   }

   //! Copy the TriangleNum-th triangle from the MeshType argument.
   template <typename MeshType>
   PointerTriangleDataStorage ( const MeshType &Mesh, const int TriangleNum );

   PointerTriangleDataStorage& operator= ( const PointerTriangleDataStorage<DataType> &rhs ) {
     _nodes[0] = &rhs.getNode(0); _nodes[1] = &rhs.getNode(1); _nodes[2] =  &rhs.getNode(2);
     return *this;
   }

   const Vec3<DataType>& getNode ( int i ) const {
     return *_nodes[i];
   }

   Vec3<DataType>& getNode ( int /*i*/ ) {
     throw aol::Exception ( "PointerTriangleDataStorage::getNode: ERROR!!", __FILE__, __LINE__ );
   }

   void setNode ( int i, const Vec3<DataType>& node ) {
     _nodes[i] = &node;
   }

};


/** Triangles in 3D
 *  \author Horn
 */
template<typename DataType, typename StorageType = DefaultTriangleDataStorage<DataType> >
class Triangle : public StorageType {

public:
  Triangle () {}

  Triangle ( const Vec3<DataType> &Node0, const Vec3<DataType> &Node1, const Vec3<DataType> &Node2 ) : StorageType( Node0, Node1, Node2 ){ }

  Triangle ( const Vec3<DataType> Node[3] ) : StorageType( Node[0], Node[1], Node[2] ){ }

  Triangle ( const Triangle<DataType, StorageType> &rhs ) : StorageType( rhs ) { }

  //! Copy the TriangleNum-th triangle from the TriangMesh argument.
  template <typename MeshType>
  Triangle ( const MeshType &Mesh, const int TriangleNum ) : StorageType( Mesh, TriangleNum ) { }

  // default destructor OK

  const Vec3<DataType>& operator[] ( int i ) const {
    return this->getNode(i);
  }

  Vec3<DataType>& operator[] ( int i ) {
    return this->getNode(i);
  }

 //  fill Element function
 template <typename MeshType>
  void fillElement( const MeshType& Mesh, int globalIdx ){
    /// fill Geometry corresponding to StorageType
    for ( int i = 0; i < 3; ++i )
	  this->setNode(i, Mesh.getVertex( Mesh.getTriangNodeIdx ( globalIdx , i ) ) );
  }

  void printNodes() const {
    cerr << "node[0] = " << this->getNode(0) << "  node[1] = " << this->getNode(1) << "  node[2] = " << this->getNode(2) << "  ";
  }

 
  /** calculates barycentric coordinate of the point chi with vectors dx1 and dx2 spanning the triangle
   *  \attention method name misleading
   */
  int barycenter ( DataType &lambda0,
                   DataType &lambda1,
                   DataType &lambda2,
                   const Vec3<DataType> &dx1,
                   const Vec3<DataType> &dx2,
                   const Vec3<DataType> &chi ) const {

    const DataType g11 = ( dx1 * dx1 );
    const DataType g12 = ( dx1 * dx2 );
    const DataType g22 = ( dx2 * dx2 );

    const DataType g = g11 * g22 - Sqr ( g12 );

#ifdef VERBOSE
    if ( g < 0.00001 ) {
      cerr << "g = " << g << endl;
    }
#endif

    const DataType g_11 = g22 / g;
    const DataType g_12 = -g12 / g;
    const DataType g_22 = g11 / g;

    const DataType alpha1 = chi * dx1;
    const DataType alpha2 = chi * dx2;
    lambda1 = ( g_11 * alpha1 ) + ( g_12 * alpha2 );
    lambda2 = ( g_12 * alpha1 ) + ( g_22 * alpha2 );
    lambda0 = aol::ZOTrait<DataType>::one - lambda1 - lambda2;

    const DataType EPSILON = 0.0;

    int type = 0;
    if ( lambda0 <= EPSILON ) type |= 4;
    if ( lambda1 <= EPSILON ) type |= 2;
    if ( lambda2 <= EPSILON ) type |= 1;
    return type;
  }


  // for chi given in the plane of the triangle,
  // computes the barycentric coordinates
  int coordToBary ( const aol::Vec3<DataType> &chi,
                    aol::Vec3<DataType> &bary ) const {

    aol::Vec3<DataType> p = chi;

    aol::Vec3<DataType> dx1 = this->getNode(1), dx2 = this->getNode(2);
    dx1 -= this->getNode(0);
    dx2 -= this->getNode(0);
    p -= this->getNode(0);

    return barycenter ( bary[0], bary[1], bary[2], dx1, dx2, p );
  }

  //! transform barycentric coordinates to cartesian coordinates
  void barToCartCoord ( const aol::Vec<2, DataType>& barCoord, aol::Vec<3, DataType>& CartCoord ) const {
    aol::Mat< 3, 2, DataType > transform;
    transform.setCol (0, this->getNode(0) - this->getNode(2));
    transform.setCol (1, this->getNode(1) - this->getNode(2));
    transform.mult ( barCoord, CartCoord);
    CartCoord += this->getNode(2);
  }

  Vec3<DataType> centerOfMass ( ) const {
    return  Vec3<DataType> ( ( this->getNode(0)[0] + this->getNode(1)[0] + this->getNode(2)[0] )/3.,
                             ( this->getNode(0)[1] + this->getNode(1)[1] + this->getNode(2)[1] )/3.,
                             ( this->getNode(0)[2] + this->getNode(1)[2] + this->getNode(2)[2] )/3. );
  }
  
  //! convention: local edge index corresponds to opposite local node index
  Vec3<DataType> getEdgeMidpoint ( int localEdgeIndex ) const {
    return  Vec3<DataType> ( ( this->getNode((localEdgeIndex+1)%3)[0] + this->getNode((localEdgeIndex+2)%3)[0] )/2.,
                             ( this->getNode((localEdgeIndex+1)%3)[1] + this->getNode((localEdgeIndex+2)%3)[1] )/2.,
                             ( this->getNode((localEdgeIndex+1)%3)[2] + this->getNode((localEdgeIndex+2)%3)[2] )/2. );
  }

  inline void barycenter ( Vec3<DataType> &center ) const {
    center = centerOfMass();
  }

  inline void centerOfMass ( Vec3<DataType> &center ) const {
    center = centerOfMass();
  }


  void weightedNormal ( aol::Vec3<DataType> &normal ) const {
    aol::Vec3<DataType> dx1 = this->getNode(1), dx2 = this->getNode(2);
    dx1 -= this->getNode(0);
    dx2 -= this->getNode(0);
    normal.makeCrossProduct ( dx1, dx2 );
  }

  aol::Vec3<DataType> weightedNormal ( ) const {
    aol::Vec3<DataType> ret;
    weightedNormal ( ret );
    return ( ret );
  }

  void normalizedNormal ( aol::Vec3<DataType> &normal ) const {
    weightedNormal ( normal );
    normal.normalize();
  }

  aol::Vec3<DataType> normalizedNormal ( ) const {
    aol::Vec3<DataType> ret;
    normalizedNormal ( ret );
    return ( ret );
  }

  inline DataType area ( ) const {
    Vec3<DataType> n;
    weightedNormal ( n );
    return n.norm() / 2.;
  }

  //! 2D volume of the triangle, i.e. the area of the triangle.
  //! \note necessary to be compatible with om::OpenMeshConfigurator
  //! \todo Shall we try to get rid of this?
  DataType vol( ) const {
    return area();
  }

  void getBoundingBox ( aol::Vec3<DataType> &min,
                        aol::Vec3<DataType> &max ) const {
    min = this->getNode(0);
    max = this->getNode(0);

    for ( int i = 1; i < 3; i++ ) {
      for ( int j = 0; j < 3; j++ ) {
        min[j] = aol::Min ( min[j], this->getNode(i)[j] );
        max[j] = aol::Max ( max[j], this->getNode(i)[j] );
      }
    }
  }

  //! \attention method name misleading
  void getBoundingBox ( Vec3<int> &min, Vec3<int> &max ) const {
    Vec3<DataType> rmin, rmax;
    getBoundingBox ( rmin, rmax );
    for ( int i = 0; i < 3; i++ ) {
      min[i] = static_cast<int> ( floor ( rmin[i] ) );
      max[i] = static_cast<int> ( ceil ( rmax[i] ) );
    }
  }

  Vec3<DataType> edge ( const int locNode1, const int locNode2 ) const {
    return Vec3<DataType> ( this->getNode(locNode2) - this->getNode(locNode1) );
  }

  DataType edgeLength ( const int id0, const int id1 ) const {
    return (  this->getNode(id0) - this->getNode(id1) ).norm();
  }

  DataType incircleRadius() const {
    return ( 2. * area() / ( edgeLength ( 0, 1 ) + edgeLength ( 0, 2 ) + edgeLength ( 1, 2 ) ) );
  }

  DataType circumcircleRadius() const {
    return ( ( edgeLength ( 0, 1 ) *edgeLength ( 0, 2 ) *edgeLength ( 1, 2 ) ) / ( 4.*area() ) );
  }

  DataType maxEdgeLength() const {
    return ( aol::Max ( edgeLength ( 0, 1 ), edgeLength ( 0, 2 ), edgeLength ( 1, 2 ) ) );
  }

  DataType minEdgeLength() const {
    return ( aol::Min ( edgeLength ( 0, 1 ), edgeLength ( 0, 2 ), edgeLength ( 1, 2 ) ) );
  }

  bool isDegenerate ( const DataType Tolerance = aol::NumberTrait<DataType>::zero ) const {
    aol::Vec3<DataType> wN;
    weightedNormal ( wN );
    return ( fabs ( wN.norm() ) <= Tolerance );
  }

  //! Calculates the unsigned distance form point to the triangle,
  //! if edge or point of triangle is more closely to the point than the distance to the plan, than the distance to edge / point is used.
  DataType calcDist ( const Vec3<DataType> &point ) const {
    Vec3<DataType> dx1 = this->getNode(1);
    dx1 -= this->getNode(0);
    Vec3<DataType> dx2 = this->getNode(2);
    dx2 -= this->getNode(0);

    Vec3<DataType> normal;
    normalizedNormal ( normal );
    return fabs ( calcDist ( point, normal ) );
  }

  DataType calcDist ( const Vec3<DataType> &point, const Vec3<DataType> &_normal ) const {
    DataType aux;
    return calcDist ( point, _normal, aux );
  }

  //! SignedDistToPlane will contain the signed (depending on the normal) distance of point to the plane through the triangle
  DataType calcDist ( const Vec3<DataType> &point, const Vec3<DataType> &_normal, DataType &SignedDistToPlane ) const {

    Vec3<DataType> dx1 = this->getNode(1);
    dx1 -= this->getNode(0);
    Vec3<DataType> dx2 = this->getNode(2);
    dx2 -= this->getNode(0);
    Vec3<DataType> chi;

    DataType l0, l1, l2;
    int flag;
    chi = calProject ( point, _normal );
    // chi is projected point on plane now

    chi -= this->getNode(0);
    flag = barycenter ( l0, l1, l2, dx1, dx2, chi );
    chi += this->getNode(0);

    // return euclidianDist( chi, point );

    Vec3<DataType> p;

    Vec3<DataType> aux( point );
    aux -= chi;
    SignedDistToPlane = aux * _normal;

    switch ( flag ) {
    case 0:
      return euclidianDist ( chi, point );
      break;
    case 1:
      projectToSegment ( chi, this->getNode(0), this->getNode(1), p );
      return euclidianDist ( p, point );
      break;
    case 2:
      projectToSegment ( chi, this->getNode(0), this->getNode(2), p );
      return euclidianDist ( p, point );
      break;
    case 4:
      projectToSegment ( chi, this->getNode(1), this->getNode(2), p );
      return euclidianDist ( p, point );
      break;
    case 3:
      return euclidianDist ( this->getNode(0), point );
      break;
    case 5:
      return euclidianDist ( this->getNode(1), point );
      break;
    case 6:
      return euclidianDist ( this->getNode(2), point );
      break;
    default:
      cerr << "no case" << endl;
      return -1;
    }
  }


  /** Project point to the plane of the triangle
   *  \attention method name misleading
   *  \attention duplicate method?
   */
  Vec3<DataType> projectPoint ( const aol::Vec3<DataType> &Point ) const {
    aol::Vec3<DataType> projection;
    projectToTriangPlane ( Point, projection );
    return ( projection );
  }

  /** normal has to be the normalized normal of the triangle, than calProject returns the projection of point into the plane in which the triangle lies
   *  \attention duplicate method?
   */
  Vec3<DataType> calProject ( const Vec3<DataType> &point,
                              const Vec3<DataType> &normal ) const {
    const Vec3<DataType> sp = point - this->getNode(0);
    const DataType fac = sp * normal;
    return ( point - fac * normal );
  }

  /** Project point to the plane of the triangle?
   *  \attention duplicate method?
   */
  void projectToTriangPlane ( const aol::Vec3<DataType> &p,
                              const aol::Vec3<DataType> &normal,
                              aol::Vec3<DataType> &proj ) const {
    proj = this->getNode(0) - p;
    const DataType tmp = proj * normal;
    proj = tmp * normal + p;
  }

  /** \attention duplicate method?
   */
  void projectToTriangPlane ( const aol::Vec3<DataType> &p,
                              aol::Vec3<DataType> &proj ) const {
    aol::Vec3<DataType> normal;
    weightedNormal ( normal );
    normal.normalize();
    projectToTriangPlane ( p, normal, proj );
  }

  //! does this method project to the edges of the triangle?
  inline void projectToSegment ( const Vec3<DataType> &point,
                                 const Vec3<DataType> &p1,
                                 const Vec3<DataType> &p2,
                                 Vec3<DataType> &projPt ) const {
    DataType l = ( point[0] - p2[0] ) * ( p1[0] - p2[0] ) + ( point[1] - p2[1] ) * ( p1[1] - p2[1] ) + ( point[2] - p2[2] ) * ( p1[2] - p2[2] );
    l /= Sqr ( p1[0] - p2[0] ) + Sqr ( p1[1] - p2[1] ) + Sqr ( p1[2] - p2[2] );
    if ( l >= 1. ) {
      projPt = p1;
    } else if ( l <= 0. ) {
      projPt = p2;
    } else {
      projPt = p1;
      projPt -= p2;
      projPt *= l;
      projPt += p2;
    }
  }

  //! \attention no size checking, what does this method do?
  inline void moment1 ( DataType* m ) const {
    DataType w1[3], w2[3], w3[3];

    w1[0] = ( this->getNode(0)[0] + this->getNode(1)[0] ) / 2;  w1[1] = ( this->getNode(0)[1] + this->getNode(1)[1] ) / 2; w1[2] = ( this->getNode(0)[2] + this->getNode(1)[2] ) / 2;
    w2[0] = ( this->getNode(1)[0] + this->getNode(2)[0] ) / 2;  w2[1] = ( this->getNode(1)[1] + this->getNode(2)[1] ) / 2; w2[2] = ( this->getNode(1)[2] + this->getNode(2)[2] ) / 2;
    w3[0] = ( this->getNode(0)[0] + this->getNode(2)[0] ) / 2;  w3[1] = ( this->getNode(0)[1] + this->getNode(2)[1] ) / 2; w3[2] = ( this->getNode(0)[2] + this->getNode(2)[2] ) / 2;

    // This does not make sense, if DataType is non floating point
    DataType a = static_cast<DataType>(1.0 / 3);

    m[0] = a * ( w1[0] * w1[0] + w2[0] * w2[0] + w3[0] * w3[0] );   //M[0,0]
    m[1] = a * ( w1[0] * w1[1] + w2[0] * w2[1] + w3[0] * w3[1] );   //M[0,1]
    m[2] = a * ( w1[0] * w1[2] + w2[0] * w2[2] + w3[0] * w3[2] );   //M[0,2]
    m[3] = a * ( w1[1] * w1[1] + w2[1] * w2[1] + w3[1] * w3[1] );   //M[1,1]
    m[4] = a * ( w1[1] * w1[2] + w2[1] * w2[2] + w3[1] * w3[2] );   //M[1,2]
    m[5] = a * ( w1[2] * w1[2] + w2[2] * w2[2] + w3[2] * w3[2] );   //M[2,2]
  }

  static bool isSelfTestOK() {
    cerr << "\nTesting distances from points to triangle ... ";
    bool isTestOK = true;

    Vec3<DataType> node0 ( 0., 0., 0. ), node1( 1., 0., 0. ), node2( 0., 1., 0. );
    Vec3<DataType> point[6];

    const aol::Triangle<DataType,StorageType> triangle ( node0, node1, node2 );
    point[0].set( 0.5, 0.5, 1. );
    point[1].set( 0.5, 0.5, -1. );
    point[2].set( 0.5, 0.5, 1. );
    point[3].set( 0., 0., 1. );
    point[4].set( 2., 0., 0. );
    point[5].set( 0., 2., 0. );
    for( int i = 0; i < 6; i++ )
      if ( triangle.calcDist( point[i] ) != 1. ) {
        cerr << "FAILED!" << endl;
        isTestOK = false;
      }
    if( isTestOK )
      cerr << "OK." <<endl;
    return isTestOK;
  }
  
  	// FEOps require this function!!!
	int level() const{
		return 0;
	}
	// FEOps require this function!!!
	int type() const{
		return 0;
	}

};



/** \brief Element for classical Finite Element simulations.
 *  \author Heeren, Perl
 * 
 *  As aol::Triangle<>, but has a relation to an underlying grid, i.e. 
 *   - it knows its global triangle index in the grid
 *   - it knows the global node indices of its nodes in the grid
 *   - it stores the inverse metric (which is needed several times when evaluating FE base function gradients)
 *   - it stores its area (which is needed several times when doing numerical quarature for assembling FE matrices)
 */
template< typename RealType, typename StorageType = DefaultTriangleDataStorage<RealType> >
class TriangBaseElement : public aol::Triangle<RealType, StorageType>{

  protected:   
    // global indices of element an nodes
    int _globIdx;
    aol::Vec<3, int> _globNodeIdx;
    // inverse metric
    aol::Matrix22<RealType> _ginv;
    // area
    RealType _area;
    
  public:
    TriangBaseElement() : aol::Triangle<RealType, StorageType>(), _globIdx(-1){}

    template<typename MeshType>
    TriangBaseElement( const MeshType& Mesh, int globalIdx ) : aol::Triangle<RealType, StorageType>(), _globIdx(globalIdx){
      fillElement( Mesh , globalIdx );
    }

    ~TriangBaseElement(){}

    //  fill Element function
    template<typename MeshType>
    void fillElement( const MeshType& Mesh, int globalIdx ){
      _globIdx = globalIdx;
      /// fill Geometry corresponding to StorageType
      for ( int i = 0; i < 3; ++i ){
	    _globNodeIdx[i] =  Mesh.getTriang(_globIdx)[i];
	    this->setNode( i , Mesh.getVertex( _globNodeIdx[i] ) );
	  }
	  /// compute Inverse Metric and save
	  setInverseMetric();
    }

    int globIdx(  ) const {
      return _globIdx;
    }
    // get and set functions
    int getIndex() const{
      return _globIdx;
    }

    int globNodeIdx(int localIndex) const{
      return _globNodeIdx[localIndex];
    }
    
    const aol::Matrix22<RealType>& ginv() const{
      return _ginv;
    }
  
    inline RealType area ( ) const { 
      return _area;
    }
  
    // For compatibility with qc::Element
    short level() const {
      throw aol::Exception ( "aol::TriangBaseElement does not have a level.", __FILE__, __LINE__ );
      return -1;
    }
    unsigned char type() const {
      throw aol::Exception ( "aol::TriangBaseElement does not have a type.", __FILE__, __LINE__ );
      return 23;
    }

  protected:
    void setInverseMetric() {
      aol::Vec3<RealType> dx1 = this->edge(0,1), dx2 = this->edge(0,2);
      RealType g11 = ( dx1 * dx1 );
      RealType g12 = ( dx1 * dx2 );
      RealType g22 = ( dx2 * dx2 );
      RealType detg = g11 * g22 - Sqr ( g12 );
      if ( detg < 1e-15 )
	    throw aol::Exception ( "aol::TriangBaseElement::setInverseMetric: determinant not positive.\n", __FILE__, __LINE__ );
      _ginv.set(0,0,g22/detg);
      _ginv.set(1,1,g11/detg);
      _ginv.set(1,0,-g12/detg);
      _ginv.set(0,1,-g12/detg);
      _area = 0.5 * sqrt( detg );
    }
 
};


/** Element for Discrete Kirchhoff Triangle (for plates) stores edges and squared edge lengths.
 * 
 *  Convention: local edge number corresponds to local node index of opposite node!
 *  Notation:   ei = e_i = x_{i-1} - x_{i+1} = x_{i+2} - x_{i+1}, i = 0,1,2
 *              ei = (ei_x, ei_y) \in \R^2
 * 
 *  \author Heeren, Perl
 */
template< typename RealType, typename StorageType = DefaultTriangleDataStorage<RealType> >
class DKTPlateElement : public aol::TriangBaseElement<RealType, StorageType> {
  
protected:
  Vec3<RealType> _edges[3];
  Vec3<RealType> _edgeLengthSqr;
  aol::Matrix22<RealType> _metric; //TODO this should be D (\phi) whereas ginv = (D phi^-1)^T (D \phi^-1)
  
public:
  DKTPlateElement() : aol::TriangBaseElement<RealType, StorageType>(){}

  template<typename MeshType>
  DKTPlateElement( const MeshType& Mesh, int globalIdx ) : aol::TriangBaseElement<RealType, StorageType>(){
    fillElement(Mesh, globalIdx);
  }
  
  ~DKTPlateElement(){}
  
  using aol::TriangBaseElement<RealType, StorageType>::edge;
  
  const Vec3<RealType>& edge ( const int locNode ) const {
    return _edges[locNode];
  }
  
  // get d = ei_x / |ei|^2, i = local node index
  inline RealType edgeXNormSqr ( const int locNode ) const {
    return _edges[locNode][0] / _edgeLengthSqr[locNode];
  }
  
  // get a = ei_y / |ei|^2, i = local node index
  inline RealType edgeYNormSqr ( const int locNode ) const {
    return _edges[locNode][1] / _edgeLengthSqr[locNode];
  }  
  
  // get b = ei_x*ei_y / |ei|^2, i = local node index
  inline RealType edgeXYNormSqr ( const int locNode ) const {
    return _edges[locNode][0]*_edges[locNode][1] / _edgeLengthSqr[locNode];
  }
  
  // get c = ( ei_x*ei_x - 2*ei_y*ei_y) / |ei|^2, i = local node index
  inline RealType edgeXXMin2YYNormSqr ( const int locNode ) const {
    return (_edges[locNode][0]*_edges[locNode][0] - 2.*_edges[locNode][1]*_edges[locNode][1]) / _edgeLengthSqr[locNode];
  }
  
  // get e = ( ei_y*ei_y - 2*ei_x*ei_x) / |ei|^2, i = local node index
  inline RealType edgeYYMin2XXNormSqr ( const int locNode ) const {
    return (_edges[locNode][1]*_edges[locNode][1] - 2.*_edges[locNode][0]*_edges[locNode][0]) / _edgeLengthSqr[locNode];
  }
  
  // get l^2 = |ei|^2, i = local node index
  inline RealType edgeLengthSqr( int locNode ) const {
    return _edgeLengthSqr[locNode];
  } 

  template <typename MeshType>
  void fillElement( const MeshType& Mesh, int globalIdx ){ 
    this->_globIdx = globalIdx;
    // fill Geometry corresponding to StorageType
    for ( int i = 0; i < 3; ++i ){
      this->_globNodeIdx[i] =  Mesh.getTriang(this->_globIdx)[i];
      this->setNode( i , Mesh.getVertex( this->_globNodeIdx[i] ) );
    }
    // fill edges and squared edge lengths
    for ( int i = 0; i < 3; ++i ){
      _edges[i] = this->getNode((i+2)%3) - this->getNode((i+1)%3);
      _edgeLengthSqr[i] = _edges[i].normSqr();
    }    
    // compute inverse Metric and save
    setMetric();
    this->setInverseMetric();
  }
  
  //!TODO check whether this corresponds to new edge ordering and remove warning!
  void print() const {
    cerr << "CAUTION: numbering might still correspond to old convention of local edge numbering!" << endl;
    cerr << "-----------------------------------------" << endl;
    for( int i = 0; i < 3; i++ )
      cerr << "x_" << i+1 << " = " << this->getNode(i) << endl;
    cerr << "-----------------------------------------" << endl;
    
    for( int i = 0; i < 3; i++ ){
      cerr << "x_" << i+1 << ((i+1)%3)+1 <<  " = " << _edges[i][0] << endl;
      cerr << "y_" << i+1 << ((i+1)%3)+1 <<  " = " << _edges[i][1] << endl;
      cerr << endl;
    }
    cerr << "-----------------------------------------" << endl;
    
    for( int i = 0; i < 3; i++ ){
      cerr << "a_" << i+1 << ((i+1)%3)+1 << " = " << edgeYNormSqr ( i ) << endl;
      cerr << "b_" << i+1 << ((i+1)%3)+1 << " = " << edgeXYNormSqr ( i )<< endl;
      cerr << "c_" << i+1 << ((i+1)%3)+1 << " = " << edgeXXMin2YYNormSqr ( i )<< endl;
      cerr << "d_" << i+1 << ((i+1)%3)+1 << " = " << edgeXNormSqr ( i )<< endl;
      cerr << "e_" << i+1 << ((i+1)%3)+1 << " = " << edgeYYMin2XXNormSqr ( i )<< endl << endl;
    }
    cerr << "-----------------------------------------" << endl;    
    cerr << "area = " << this->_area << endl;
    cerr << "-----------------------------------------" << endl;
  }
  
  
  const aol::Matrix22<RealType>& metric() const{
      return _metric;
  }
  
protected:  
//   void setInverseMetric() {
//     // g = [e2|-e1] = [ e2_0  (-e1_0); e2_1 (-e1_1)],  detg = e2_0*(-e1_1) - e2_1*(-e1_0)
//     RealType detg = _edges[1][0]*_edges[2][1] - _edges[2][0]*_edges[1][1];
//     if ( detg < 1e-22 ){
//       stringstream err_msg;
//         err_msg << "aol::DKTPlateElement::setInverseMetric: determinant of triangle " << this->_globIdx << " not positive.\n";
// 	err_msg << "Nodes = (" << this->getNode(0) << "), ("<< this->getNode(1) << "), ("<< this->getNode(2) << ")\n";
//       throw aol::Exception ( err_msg.str(), __FILE__, __LINE__ );
//     }
//     // g^{-1} = 1/detg [ (-e1_1) e1_0; (-e2_1) e2_0 ]
//     this->_ginv.set( 0, 0, -1.*_edges[1][1]/detg ); 
//     this->_ginv.set( 0, 1,     _edges[1][0]/detg );
//     this->_ginv.set( 1, 0, -1.*_edges[2][1]/detg );
//     this->_ginv.set( 1, 1,     _edges[2][0]/detg );
//     this->_area = 0.5 * detg;
//   }
  
  
  
  void setMetric() {
    // g = [e2|-e1] = [ e2_0  (-e1_0); e2_1 (-e1_1)],  detg = e2_0*(-e1_1) - e2_1*(-e1_0)
    RealType detg = _edges[1][0]*_edges[2][1] - _edges[2][0]*_edges[1][1];
    if ( detg < 1e-22 ){
      stringstream err_msg;
        err_msg << "aol::DKTPlateElement::setInverseMetric: determinant of triangle " << this->_globIdx << "is " << detg << " , hence not positive.\n";
        err_msg << "Nodes = (" << this->getNode(0) << "), ("<< this->getNode(1) << "), ("<< this->getNode(2) << ")\n";
      throw aol::Exception ( err_msg.str(), __FILE__, __LINE__ );
    }
    // g^{-1} = 1/detg [ (-e1_1) e1_0; (-e2_1) e2_0 ]
    this->_metric.set( 0, 0, -1.*_edges[1][1]/detg ); 
    this->_metric.set( 0, 1,     _edges[1][0]/detg );
    this->_metric.set( 1, 0, -1.*_edges[2][1]/detg );
    this->_metric.set( 1, 1,     _edges[2][0]/detg );
//     this->_area = 0.5 * detg;
  }
  
  
};

}

#endif
