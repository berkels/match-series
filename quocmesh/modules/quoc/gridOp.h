#ifndef __GRIDOP_H
#define __GRIDOP_H

#include <op.h>
#include <gridBase.h>
#include <scalarArray.h>

namespace qc {


/**
 * Sets all values specified by the iterator to a specified value.
 * \author Droske
 */
template <typename RealType, typename IteratorType, typename GridType = qc::GridDefinition>
class SetValuesOp : public aol::Op<aol::Vector<RealType> > {
public:
  SetValuesOp ( const GridType /*qc::GridDefinition*/ &Grid )
      : _grid ( Grid ), _val ( static_cast<RealType> ( 0. ) ) {}

  void setValue ( RealType Val ) {
    _val = Val;
  }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> tmp ( Dest.size() );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( Arg.size() != Dest.size() || Arg.size() != _grid.getNumberOfNodes() ) {
      throw aol::Exception ( "SetValuesOp: size of vectors should be equal to number of nodes of grid.\n", __FILE__, __LINE__ );
    }

    // If the given vectors do not point to the same address, copy Arg to Dest.
    if ( &Dest != &Arg ) {
      Dest = Arg;
    }

    IteratorType it;
    for ( it = _grid.begin(); it != _grid.end(); ++it ) {
      Dest[ it->y() * _grid.getWidth() + it->x() ] = _val;
    }
  }
protected:
  const /*qc::GridDefinition*/GridType &_grid;
  RealType _val;
};


// class provides an apply that operates only on specified elements. Therefore in the
// derived class an iterator has to be implemented, which traverses the elements, that
// should not be included in the apply.
template <typename RealType, typename Imp, typename VecType = aol::Vector<RealType> >
class RestrictValuesOpInterface : public aol::Op<VecType>  {
protected:
  const aol::Op<VecType> &_op;
public:
  RestrictValuesOpInterface ( const aol::Op<VecType> &Op )
      : _op ( Op ) {}

  vector<int>::const_iterator begin() const {
    return asImp().begin();
  }

  vector<int>::const_iterator end() const {
    return asImp().end();
  }

  int dofNumber ( const CoordType &c ) const {
    return asImp().dofNumber ( c );
  }

  void applyAdd ( const VecType &Arg, VecType &Dest ) const {
    VecType tmp ( Dest );
    apply ( Arg, tmp );
    Dest += tmp;
  }


  // apply only applies the op on elements, that aren't reached by the iterator of
  // the mask.
  void apply ( const VecType &Arg, VecType &Dest ) const {
    if ( static_cast<int> ( Arg.size() ) !=  static_cast<int> ( Dest.size() ) ) {
      throw aol::Exception ( "RestrictValuesOp: size of vectors should be equal to number of nodes of grid.\n", __FILE__, __LINE__ );
    }

    VecType tmp ( Arg );
    tmp = Arg;
    vector<int>::const_iterator it;

    // set all nodes, that are reached by the mask-iterator, to zero.
    for ( it = begin(); it != end(); ++it ) {
      int index = *it;
      tmp[ index ] =  static_cast<RealType> ( 0. );
    }

    _op.apply ( tmp, Dest );

    // copy the original values to the mask-elements
    for ( it = begin(); it != end(); ++it ) {
      int index = *it;
      Dest[ index ] = Arg[ index ];
    }
  }

protected:
  Imp& asImp() {
    return static_cast<Imp&> ( *this );
  }
  const Imp& asImp()  const {
    return static_cast<const Imp&> ( *this );
  }

};

/**
 * performs an apply of a restricted version of the given operator.
 * the IteratorType specifies the iterator walking through the nodes to be restricted.
 * \author Droske
 */
template <typename RealType>
class RestrictValuesOpDef : public RestrictValuesOpInterface<RealType, RestrictValuesOpDef<RealType> > {
  vector<int> _bnd;
public:
  RestrictValuesOpDef ( const aol::Op<aol::Vector<RealType> > &Op, const qc::GridDefinition &Grid )
      : RestrictValuesOpInterface<RealType, RestrictValuesOpDef<RealType> > ( Op ), _grid ( Grid ) {

    qc::GridDefinition::OldFullBoundaryNodeIterator it;
    for ( it = _grid.begin(); it != _grid.end(); ++it ) {
      _bnd.push_back ( dofNumber ( *it ) );
    }
  }

  virtual ~RestrictValuesOpDef() {}

  vector<int>::const_iterator begin() const {
    return _bnd.begin();
  }

  vector<int>::const_iterator end() const {
    return _bnd.end();
  }

  int dofNumber ( const CoordType &c ) const {
    return ( c.z() * _grid.getWidth()  + c.y() ) * _grid.getWidth() + c.x();
  }

protected:
  const qc::GridDefinition &_grid;
};


// *********************************************************************************
// RestrictValuesFromVectorOp
// restricts an operator using a vector as a mask
// mask != 0 => calculate here
// mask == 0 => touch this region and die!
// *********************************************************************************

template <typename RealType, typename VecType>
class RestrictValuesFromVectorOp : public qc::RestrictValuesOpInterface < RealType,
      RestrictValuesFromVectorOp<RealType, VecType>, VecType > {
protected:
  const qc::GridDefinition &_grid;
  vector<int> _bnd;
public:

  // constructor
  explicit RestrictValuesFromVectorOp ( const aol::Op<VecType> &Op, const aol::Vector<RealType> &maskVec,
                                        const qc::GridDefinition &grid )
      : qc::RestrictValuesOpInterface<RealType, RestrictValuesFromVectorOp<RealType, VecType>, VecType > ( Op ),
      _grid ( grid ) {
    for ( int i = 0; i < maskVec.size(); i++ )
      if ( maskVec[i] == 0 ) _bnd.push_back ( i );
  }

virtual ~RestrictValuesFromVectorOp() {}

  vector<int>::const_iterator begin() const {
    return _bnd.begin();
  }

  vector<int>::const_iterator end() const {
    return _bnd.end();
  }

  int dofNumber ( const CoordType &c ) const {
    return ( c.z() * _grid.getWidth()  + c.y() ) * _grid.getWidth() + c.x();
  }

  // test the iterator: Elements, that are used in the computation are colored
  // black, elements, that are deleted are white
  void testIt ( qc::ScalarArray<RealType, qc::QC_2D> &test ) const {
    test.setAll ( 0. );
    vector<int>::const_iterator it;

    for ( it = begin(); it != end(); ++it ) {
      int index = *it;
      test[ index ] =  static_cast<RealType> ( 1. );
    }
    test.save ( "test.bz2", 9 );
  }


};

}

#endif
