#ifndef __PERIODICBC_H
#define __PERIODICBC_H

#include <arrayExtensions.h>
#include <gridBase.h>
#include <iterators.h>
#include <matrix.h>
#include <vectorExtensions.h>

namespace qc {

/** Basis class for handling of periodic boundary conditions.
 *  Pure virtual functions need to be overloaded to represent appropriate periodicity
 *  \author Schwen
 */
template< typename GridType, typename DataType, qc::Dimension Dim >
class PeriodicityHandlerBase {
protected:
  const GridType &_grid;
  qc::FastILexMapper< Dim > _map;
  bool _quietMode;

public:
  PeriodicityHandlerBase ( const GridType &Grid ) : _grid ( Grid ), _map ( Grid ), _quietMode ( false ) {
    if ( _grid.getDimOfWorld() != Dim )
      throw aol::Exception ( "PeriodicityHandlerBase: dimension of grid does not match Dim", __FILE__, __LINE__ );
  }

  virtual ~PeriodicityHandlerBase ( ) {
  }

  // default copy constructor, assignment operator and destructor are not implemented explicitly.

  //! check whether a node is not a DOF but a periodic copy of another DOF
  virtual bool isPeriodicNode ( const qc::CoordType &/*pos*/ ) const = 0;

  //! for periodic nodes, returns the position where DOF is actually stored, for all other nodes, returns the position
  virtual qc::CoordType facingPresentNode ( const qc::CoordType &/*pos*/ ) const = 0;

  //! remove rows and columns corresponding to periodic nodes, adding entries to facing present positions, writing diagonal entry (generic slow implementation, overload on derived classes)
  virtual void periodicallyCollapseMatrix ( aol::Matrix<DataType> &matrix, const DataType diagValue = aol::NumberTrait<DataType>::one ) const {
    aol::ProgressBar<> pb ( "Slowly Collapsing matrix" );
    if ( !_quietMode ) pb.start ( _grid.getNumberOfBoundaryNodes() );

    for ( qc::RectangularBoundaryIterator< Dim > bit ( this->_grid ); bit.notAtEnd(); ++bit ) {
      const qc::CoordType pos ( *bit ), fpos ( facingPresentNode ( pos ) );

      if ( pos != fpos ) { // need to do something
        const int gP = _map ( pos ), gF = _map ( fpos );

        // eliminate col
        for ( int i = 0; i < matrix.getNumRows(); ++i ) {
          matrix.add ( i, gF, matrix.get ( i, gP ) );
          matrix.set ( i, gP, aol::NumberTrait<DataType>::zero );
        }

        // eliminate row // could be done via makeRowEntries
        for ( int j = 0; j < matrix.getNumCols(); ++j ) {
          matrix.add ( gF, j, matrix.get ( gP, j ) );
          matrix.set ( gP, j, aol::NumberTrait<DataType>::zero );
        }

        matrix.set ( gP, gP, diagValue );
      }
      if ( !_quietMode ) pb++;
    }
    if ( !_quietMode ) pb.finish();
  }

  //! periodic collapsing of Block matrix, writing one to diagonal for diagonal blocks and zero for off-diagonal blocks
  template< typename BlockMatrixType >
  void periodicallyCollapseBlockMatrix ( BlockMatrixType &blockMatrix ) const {
    for ( int i = 0; i < blockMatrix.getNumRows(); ++i )
      for ( int j = 0; j < blockMatrix.getNumCols(); ++j )
        periodicallyCollapseMatrix ( blockMatrix.getReference ( i, j ), ( i == j ? aol::NumberTrait<DataType>::one : aol::NumberTrait<DataType>::zero ) );
  }

  void restrictPeriodicBC ( aol::Vector<DataType> &Arg ) const {
    for ( qc::RectangularBoundaryIterator< Dim > bit ( _grid ); bit.notAtEnd();++bit ) { // Boundary??
      const qc::CoordType pos ( *bit ), fpos ( facingPresentNode ( pos ) );
      const int gP = _map ( pos );
      if ( pos != fpos ) {
        Arg[gP] = aol::NumberTrait<DataType>::zero;
      }
    }
  }

  void collapsePeriodicBC ( aol::Vector<DataType> &Arg ) const {
    for ( qc::RectangularIterator< Dim > bit ( _grid ); bit.notAtEnd();++bit ) {         // not boundary??
      const qc::CoordType pos ( *bit ), fpos ( facingPresentNode ( pos ) );
      const int gP = _map ( pos ), fPP = _map ( fpos );
      if ( gP != fPP ) {
        Arg[fPP] += Arg[gP];
        Arg[gP] = aol::NumberTrait<DataType>::zero;
      }
    }
  }

  virtual void extendPeriodicBC ( aol::Vector<DataType> &Arg ) const {
    for ( qc::RectangularBoundaryIterator<qc::QC_3D> bit ( _grid ); bit.notAtEnd(); ++bit ) {
      const qc::CoordType pos ( *bit ), fpos ( facingPresentNode ( pos ) );
      const int gP = _map ( pos );
      if ( pos != fpos ) {
        Arg[gP] = Arg [ _map ( fpos ) ];
      }
    }
  }


  void restrictPeriodicBC ( aol::MultiVector<DataType> &multiVector ) const {
    for ( int i = 0; i < multiVector.numComponents(); ++i )
      restrictPeriodicBC ( multiVector[i] );
  }

  void collapsePeriodicBC ( aol::MultiVector<DataType> &multiVector ) const {
    for ( int i = 0; i < multiVector.numComponents(); ++i )
      collapsePeriodicBC ( multiVector[i] );
  }

  void extendPeriodicBC ( aol::MultiVector<DataType> &multiVector ) const {
    for ( int i = 0; i < multiVector.numComponents(); ++i )
      extendPeriodicBC ( multiVector[i] );
  }

  void setQuietMode ( bool beQuiet )  {
    _quietMode = beQuiet;
  }

};


/** Periodicity in all space directions on Quoc grids
 *  \author Schwen
 */
template< typename GridType, typename DataType, qc::Dimension Dim >
class QuocPeriodicityHandler : public PeriodicityHandlerBase< GridType, DataType, Dim > {
public:
  QuocPeriodicityHandler ( const GridType &Grid ) : PeriodicityHandlerBase< GridType, DataType, Dim > ( Grid ) {
  }

  // default copy constructor, assignment operator and destructor are correct.

  virtual bool isPeriodicNode ( const qc::CoordType &pos ) const {
    switch ( Dim ) {
    case qc::QC_2D:
      return ( ( pos[0] == ( this->_grid.getNumX() - 1 ) ) || ( pos[1] == ( this->_grid.getNumY() - 1 ) ) );
    case qc::QC_3D:
      return ( ( pos[0] == ( this->_grid.getNumX() - 1 ) ) || ( pos[1] == ( this->_grid.getNumY() - 1 ) ) || ( pos[2] == ( this->_grid.getNumZ() - 1 ) ) );
    default:
      throw aol::Exception ( "QuocPeriodicityHandler::isPeriodicNode: unsupported dimension", __FILE__, __LINE__ );
      return false;
    }
  }

  virtual qc::CoordType facingPresentNode ( const qc::CoordType &pos ) const {
    switch( Dim ) {
    case qc::QC_2D:
    {
  qc::CoordType ret ( pos ), sz ( this->_grid.getNumX() - 1, this->_grid.getNumY() - 1);
  for ( int i = 0; i < 3; ++i )
    if ( ret[i] == sz[i] )
      ret[i] = 0;
  return ( ret );
    }
    case qc::QC_3D:
    {
      qc::CoordType ret ( pos ), sz ( this->_grid.getNumX() - 1, this->_grid.getNumY() - 1, this->_grid.getNumZ() - 1 );
      for ( int i = 0; i < 3; ++i )
          if ( ret[i] == sz[i] )
            ret[i] = 0;
      return ( ret );
    }
    default:
  throw aol::Exception ( "QuocPeriodictyHandler::facingPresentNode: unsupported dimension", __FILE__, __LINE__ );
  return pos;
      }
  }

  // periodicallyCollapseMatrix could be implemented much faster than in generic case, note difference between 2D and 3D.

};


/** Periodicity in just one space direction on Quoc grids
*  \author Geihe
*/
template< typename GridType, typename DataType, qc::Dimension Dim >
class QuocOneDirectionPeriodicityHandler : public PeriodicityHandlerBase< GridType, DataType, Dim > {
  protected:
    const int _dir;

  public:
    //! Construtor. direction for periodicity given by dir. dir=0 means X, dir=1 means Y, dir=2 means Z
    QuocOneDirectionPeriodicityHandler ( const GridType &Grid, const int dir = 0 )
    : PeriodicityHandlerBase< GridType, DataType, Dim > ( Grid ),
      _dir ( dir )
    {
      if ( ( dir < 0 ) || ( dir > 2 ) )                   throw aol::Exception ( "QuocOneDirectionPeriodicityHandler: invalid direction", __FILE__, __LINE__ );
      if ( ( Dim == qc::QC_2D ) && ( dir == 2 ) )         throw aol::Exception ( "QuocOneDirectionPeriodicityHandler: invalid direction", __FILE__, __LINE__ );
      if ( ( Dim != qc::QC_2D ) && ( Dim != qc::QC_3D ) ) throw aol::Exception ( "QuocOneDirectionPeriodicityHandler: unsupported dimension", __FILE__, __LINE__ );
    }

    // default copy constructor, assignment operator and destructor are correct.

    virtual bool isPeriodicNode ( const qc::CoordType &pos ) const {
      switch ( _dir ) {
        case 0: return ( pos[0] == ( this->_grid.getNumX() - 1 ) );
        case 1: return ( pos[1] == ( this->_grid.getNumY() - 1 ) );
        case 2: return ( pos[2] == ( this->_grid.getNumZ() - 1 ) );
        default: throw aol::Exception ( "QuocOneDirectionPeriodicityHandler: invalid direction", __FILE__, __LINE__ );
      }
      return false;
    }

    virtual qc::CoordType facingPresentNode ( const qc::CoordType &pos ) const {
      qc::CoordType ret ( pos );
      switch ( _dir ) {
        case 0: { if ( ret[0] == this->_grid.getNumX() - 1 ) ret[0] = 0; break; }
        case 1: { if ( ret[1] == this->_grid.getNumY() - 1 ) ret[1] = 0; break; }
        case 2: { if ( ret[2] == this->_grid.getNumZ() - 1 ) ret[2] = 0; break; }
        default: throw aol::Exception ( "QuocOneDirectionPeriodicityHandler: invalid direction", __FILE__, __LINE__ );
      }
      return ( ret );
    }
};

}

#endif
