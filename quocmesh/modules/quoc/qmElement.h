#ifndef __QMELEMENT_H
#define __QMELEMENT_H

#include <quoc.h>

namespace qc {

/** A class for 2D/3D grid elements.
 *  The local numbering of degrees of freedom is as follows:<p>
 *   In 2D:&nbsp; \f{picture}(85, 75)
 *  \put(8,4){\line(1,0){34}}
 *  \put(8,46){\line(1,0){34}}
 *  \put(4,9){\line(0,1){32}}
 *  \put(46,9){\line(0,1){32}}
 *  \put(0,0){\makebox(7,7)[c]{0}}
 *  \put(43,0){\makebox(7,7)[c]{1}}
 *  \put(0,43){\makebox(7,7)[c]{2}}
 *  \put(43,43){\makebox(7,7)[c]{3}}
 *  \put(55,10){\vector(1,0){20}}\put(75,8){x}
 *  \put(55,10){\vector(0,1){20}}\put(55,32){y}
 *  \f}
 *  &nbsp; In 3D: &nbsp; \f{picture}(110, 75)
 *  \put(8,4){\line(1,0){34}}
 *  \put(8,46){\line(1,0){34}}
 *  \put(4,9){\line(0,1){32}}
 *  \put(46,9){\line(0,1){32}}
 *  \put(0,0){\makebox(7,7)[c]{0}}
 *  \put(43,0){\makebox(7,7)[c]{1}}
 *  \put(0,43){\makebox(7,7)[c]{2}}
 *  \put(43,43){\makebox(7,7)[c]{3}}
 *  \put(28,24){\line(1,0){15}}\put(48,24){\line(1,0){14}}
 *  \put(28,66){\line(1,0){34}}
 *  \put(24,29){\line(0,1){15}}\put(24,48){\line(0,1){13}}
 *  \put(66,29){\line(0,1){32}}
 *  \put(20,20){\makebox(7,7)[c]{4}}
 *  \put(63,20){\makebox(7,7)[c]{5}}
 *  \put(20,63){\makebox(7,7)[c]{6}}
 *  \put(63,63){\makebox(7,7)[c]{7}}
 *  \put(7,7){\line(1,1){12}}
 *  \put(50,7){\line(1,1){12}}
 *  \put(7,50){\line(1,1){12}}
 *  \put(50,50){\line(1,1){12}}
 *  \put(80,10){\vector(1,0){20}}\put(100,8){x}
 *  \put(80,10){\vector(1,1){15}}\put(97,27){z}
 *  \put(80,10){\vector(0,1){20}}\put(80,32){y}
 *  \f}
 *  @ingroup Element
 */
class Element : public aol::Vec3<short> {
public:
  Element() : aol::Vec3<short> ( 0, 0, 0 ), leveL ( 0 ), typE ( 9 )  {
    //index = 0;
  }

  Element ( short X, short Y, short Z = 0, short Level = -1, unsigned char Type = 9, int /*Index*/ = 0 )
      : aol::Vec3<short> ( X, Y, Z ), leveL ( Level ), typE ( Type ) /*, index(Index) */{
  }


  Element ( const aol::Vec3<short> &Coord, short Level, unsigned char Type = 9 )
      : aol::Vec3<short> ( Coord ), leveL ( Level ), typE ( Type ) {}

  Element ( const aol::Vec2<short> &Coord, short Level, unsigned char Type = 9 )
      : aol::Vec3<short> ( Coord[0], Coord[1], 0 ), leveL ( Level ), typE ( Type ) {}

  const Element & getCubicElement () const {
    return *this;
  }

  Element & getCubicElement () {
    return *this;
  }

  void setType ( unsigned char Type ) {
    typE = Type;
  }

  void set ( short X, short Y, short Z, short Level, unsigned char Type ) {
    aol::Vec3<short>::set ( X, Y, Z );
    leveL = Level;
    typE = Type;
  }

  void set ( short X, short Y, short Z = 0, short Level = -1 ) {
    aol::Vec3<short>::set ( X, Y, Z );
    leveL = Level;
  }

  void set ( const aol::Vec3<short>& Other ) {
    aol::Vec3<short>::set ( Other );
  }

  void set ( const Element & Other ) {
    *this = Other;
  }

  void setLevel ( short Level ) {
    leveL = Level;
  }

  void writeCoordsTo ( short & X, short & Y, short & Z ) const {
    X = (*this)[0];
    Y = (*this)[1];
    Z = (*this)[2];
  }

  bool operator== ( const Element& El ) const {
    return ( static_cast<aol::Vec3<short> > ( El ) == static_cast<aol::Vec3<short> > ( *this ) &&
             leveL == El.leveL && typE == El.typE );
  }

  bool operator!= ( const Element& El ) const {
    return ( static_cast<aol::Vec3<short> > ( El ) != static_cast<aol::Vec3<short> > ( *this ) ||
             leveL != El.leveL || typE != El.typE );
  }

  short level() const {
    return leveL;
  }
  unsigned char type() const {
    return typE;
  }

  void print() const {
    fprintf ( stderr, "Element: x=%d, y=%d, z=%d, type=%d, level=%d\n", x(), y(), z(), type(), level() );
  }

  //int index;

protected:
  short leveL;
  unsigned char typE;
};


inline ostream &operator<< ( ostream &os, const Element &el ) {
  os << "Element [( " << el.x() << ", " << el.y() << ", " << el.z() << "), level: " << el.level() << ", type: " << static_cast<int> ( el.type() ) << "]";
  return os;
}


/**
 * \brief Base class of 2D/3D boundary element. It knows on which face of the boundary it lies and the corresponding normal vector.
 * \author Teusner
 **/
template <typename RealType>
class BoundaryFaceElementBase : public qc::Element {

public:
  //! \brief BoundaryFaceElementBase::X_LOWER_BOUNDARY denotes the y-z-plane with x = 0, BoundaryFaceElementBase::X_UPPER_BOUNDARY the y-z-plane with x = 1 and so on
  enum BoundaryFaceType { X_LOWER_BOUNDARY = 0, ///< y-z-plane with x = 0
                          X_UPPER_BOUNDARY = 1, ///< y-z-plane with x = 1
                          Y_LOWER_BOUNDARY = 2,
                          Y_UPPER_BOUNDARY = 3,
                          Z_LOWER_BOUNDARY = 4,
                          Z_UPPER_BOUNDARY = 5,
                          UNDEF_BOUNDARY = 42 };

  BoundaryFaceElementBase ( short X, short Y, short Z = 0, short Level = -1, unsigned char Type = 9, int Index = 0,
                            BoundaryFaceType BFType = BoundaryFaceElementBase::UNDEF_BOUNDARY ) :
    Element ( X, Y, Z, Level, Type, Index ),
    _bfType ( BFType ) {}

  BoundaryFaceElementBase ( Element & el, BoundaryFaceType BFType = BoundaryFaceElementBase::UNDEF_BOUNDARY ) :
    Element ( el, el.level(), el.type() ),
    _bfType ( BFType ) {}

  BoundaryFaceElementBase ( ) :
    Element ( ),
    _bfType ( BoundaryFaceElementBase::X_LOWER_BOUNDARY ) {}

  void set ( short X, short Y, short Z = 0, short Level = -1,
             BoundaryFaceType BFType = BoundaryFaceElementBase::UNDEF_BOUNDARY ) {
    qc::Element::set ( X, Y, Z, Level );
    _bfType = BFType;
  }

  BoundaryFaceType getBoundaryFaceType () const {
    return _bfType;
  }

  void setBoundaryFaceType ( BoundaryFaceType BFType ) {
    _bfType = BFType;
  }

protected:
  BoundaryFaceType _bfType;
};

template <typename RealType,qc::Dimension Dim>
class BoundaryFaceElement {};

template <typename RealType>
class BoundaryFaceElement<RealType,qc::QC_3D> :
  public BoundaryFaceElementBase<RealType> {

public:
  typedef typename BoundaryFaceElementBase<RealType>::BoundaryFaceType BoundaryFaceType;

  BoundaryFaceElement ( short X, short Y, short Z = 0, short Level = -1, unsigned char Type = 9, int Index = 0,
                        BoundaryFaceType BFType = BoundaryFaceElementBase<RealType>::UNDEF_BOUNDARY ) :
    BoundaryFaceElementBase<RealType> ( X, Y, Z, Level, Type, Index, BFType ) {}

  BoundaryFaceElement ( ) :
    BoundaryFaceElementBase<RealType> ( ) {}

  //! returns the outer normal corresponding to the face on which the element lies
  aol::Vec3<RealType> getNormal ( ) const {
    static const aol::Vec3<RealType> Lookup[6] = { aol::Vec3<RealType> ( -1., 0., 0. ),
                                                   aol::Vec3<RealType> ( 1., 0., 0. ),
                                                   aol::Vec3<RealType> ( 0., -1., 0. ),
                                                   aol::Vec3<RealType> ( 0., 1., 0. ),
                                                   aol::Vec3<RealType> ( 0., 0., -1. ),
                                                   aol::Vec3<RealType> ( 0., 0., 1. ) };
    return Lookup[ static_cast<int>( this->_bfType ) ];
  }

  //! The boundary face is a 2D object, but the element itself is a 3D object. This function gets a 2D reference coordinate and returns the corresponding 3D reference coordinate.
  aol::Vec3<RealType> getRefCoordOnElementFromRefCoordOnFace ( const aol::Vec2<RealType> RefCoord2D ) const {
    RealType rcEntries[4] = { 0., 1., RefCoord2D[0], RefCoord2D[1] };

    static const short Lookup[6][3] = { { 0, 2, 3 }, { 1, 2, 3 }, { 2, 0, 3 }, { 2, 1, 3 }, { 2, 3, 0 }, { 2, 3, 1 } };

    aol::Vec3<RealType> refCoord ( rcEntries [ Lookup [ static_cast<int>( this->_bfType ) ][ 0 ]],
                                   rcEntries [ Lookup [ static_cast<int>( this->_bfType ) ][ 1 ]],
                                   rcEntries [ Lookup [ static_cast<int>( this->_bfType ) ][ 2 ]] );
    return refCoord;
  }

  //! This function returns the coordinates of the boundary nodes, which lie on the given boundary face.
  aol::Vec3<RealType> getBoundaryNode (short i)  const {
    if ( i < 0 || i > 3 ) {
      throw aol::Exception ( "wrong choice of i, it must be 0 <= i <= 3", __FILE__, __LINE__ );
    }
    static const aol::Vec3<RealType> Lookup[6][4] = { { aol::Vec3<RealType> ( 0., 0., 0. ), aol::Vec3<RealType> ( 0., 1., 0. ), aol::Vec3<RealType> ( 0., 0., 1. ), aol::Vec3<RealType> ( 0., 1., 1. ) },
                                                      { aol::Vec3<RealType> ( 1., 0., 0. ), aol::Vec3<RealType> ( 1., 1., 0. ), aol::Vec3<RealType> ( 1., 0., 1. ), aol::Vec3<RealType> ( 1., 1., 1. ) },
                                                      { aol::Vec3<RealType> ( 0., 0., 0. ), aol::Vec3<RealType> ( 1., 0., 0. ), aol::Vec3<RealType> ( 0., 0., 1. ), aol::Vec3<RealType> ( 1., 0., 1. ) },
                                                      { aol::Vec3<RealType> ( 0., 1., 0. ), aol::Vec3<RealType> ( 1., 1., 0. ), aol::Vec3<RealType> ( 0., 1., 1. ), aol::Vec3<RealType> ( 1., 1., 1. ) },
                                                      { aol::Vec3<RealType> ( 0., 0., 0. ), aol::Vec3<RealType> ( 1., 0., 0. ), aol::Vec3<RealType> ( 0., 1., 0. ), aol::Vec3<RealType> ( 1., 1., 0. ) },
                                                      { aol::Vec3<RealType> ( 0., 0., 1. ), aol::Vec3<RealType> ( 1., 0., 1. ), aol::Vec3<RealType> ( 0., 1., 1. ), aol::Vec3<RealType> ( 1., 1., 1. ) } };
    aol::Vec3<RealType> node ( Lookup[ static_cast<int>( this->_bfType ) ][ i ] );

    return node;
  }

  //! This function converts the local nodes on the boundary face to local nodes on the element.
  short getLocalNodeOnBoundary ( short i ) const {
    if ( i < 0 || i > 3 ) {
      throw aol::Exception ( "wrong choice of i, it must be 0 <= i <= 3", __FILE__, __LINE__ );
    }
    static const short Lookup[6][4] = { {0, 2, 4, 6}, {1, 3, 5, 7}, {0, 1, 4, 5}, {2, 3, 6, 7}, {0, 1, 2, 3}, {4, 5, 6, 7} };

    short locNode = Lookup[ static_cast<int>( this->_bfType ) ][ i ];

    return locNode;
  }
};

template <typename RealType>
class BoundaryFaceElement<RealType,qc::QC_2D> :
  public BoundaryFaceElementBase<RealType> {

public:
  typedef typename BoundaryFaceElementBase<RealType>::BoundaryFaceType BoundaryFaceType;

  BoundaryFaceElement ( short X, short Y, short Z = 0, short Level = -1, unsigned char Type = 9, int Index = 0,
                        BoundaryFaceType BFType = BoundaryFaceElementBase<RealType>::UNDEF_BOUNDARY ) :
    BoundaryFaceElementBase<RealType> ( X, Y, Z, Level, Type, Index, BFType ) {}

  BoundaryFaceElement ( Element & el, BoundaryFaceType BFType = BoundaryFaceElementBase<RealType>::UNDEF_BOUNDARY ) :
    BoundaryFaceElementBase<RealType> ( el, BFType )  {}

  BoundaryFaceElement ( ) :
    BoundaryFaceElementBase<RealType> ( ) {}

  //! Returns the outer normal to the boundary face of the boundary element.
  aol::Vec2<RealType> getNormal ( ) const {
    static const aol::Vec2<RealType> Lookup[4] = { aol::Vec2<RealType> ( -1., 0. ),
                                                   aol::Vec2<RealType> ( 1., 0. ),
                                                   aol::Vec2<RealType> ( 0., -1. ),
                                                   aol::Vec2<RealType> ( 0., 1. ) };
    return Lookup[ static_cast<short>( this->_bfType ) ];
  }

  //! The boundary face is a 1D object, but the element itself is a 2D object.
  //! This function gets a 1D reference coordinate and returns the corresponding 2D reference coordinate.
  aol::Vec2<RealType> getRefCoordOnElementFromRefCoordOnFace ( const RealType RefCoord1D ) const {
    const RealType rcEntries[3] = { 0., 1., RefCoord1D };
    static const short Lookup[4][3] = { { 0, 2 }, { 1, 2 }, { 2, 0 }, { 2, 1 } };
    aol::Vec2<RealType> refCoord ( rcEntries [ Lookup [ static_cast<short>( this->_bfType ) ][ 0 ]],
                                   rcEntries [ Lookup [ static_cast<short>( this->_bfType ) ][ 1 ]] );
    return refCoord;
  }

  aol::Vec2<RealType> getRefCoordOnElementFromRefCoordOnFace ( aol::Vec<1, RealType> RefCoord1D ) const {
    const RealType rcEntries[3] = { 0., 1., RefCoord1D[0] };
    static const short Lookup[4][3] = { { 0, 2 }, { 1, 2 }, { 2, 0 }, { 2, 1 } };
    aol::Vec2<RealType> refCoord ( rcEntries [ Lookup [ static_cast<short>( this->_bfType ) ][ 0 ]],
                                   rcEntries [ Lookup [ static_cast<short>( this->_bfType ) ][ 1 ]] );
    return refCoord;
  }

  //! Returns the coordinates of the ith boundary node on the boundary face of the boundary element.
  aol::Vec2<RealType> getBoundaryNode ( const short i ) const {
    if ( i < 0 || i > 2 )
      throw aol::Exception ( "Tried to access boundary node i for i not in [0,1].", __FILE__, __LINE__ );
    static const aol::Vec2<RealType> Lookup[4][2] = { { aol::Vec2<RealType> ( 0., 0. ), aol::Vec2<RealType> ( 0., 1. ) },
                                                      { aol::Vec2<RealType> ( 1., 0. ), aol::Vec2<RealType> ( 1., 1. ) },
                                                      { aol::Vec2<RealType> ( 0., 0. ), aol::Vec2<RealType> ( 1., 0. ) },
                                                      { aol::Vec2<RealType> ( 0., 1. ), aol::Vec2<RealType> ( 1., 1. ) } };
    aol::Vec2<RealType> node ( Lookup[ static_cast<short>( this->_bfType ) ][ i ] );
    return node;
  }

  //! Converts the local nodes index on the boundary face to the local node index on the whole element.
  short getLocalNodeOnBoundary ( const short i ) const {
    if ( i < 0 || i > 2 )
      throw aol::Exception ( "Tried to access boundary node i for i not in [0,1].", __FILE__, __LINE__ );
    static const short Lookup[4][2] = { { 0, 2 }, { 1, 3 }, { 0, 1 }, { 2, 3 } };
    return Lookup[ static_cast<short>( this->_bfType ) ][ i ];
  }
};

}

#endif
