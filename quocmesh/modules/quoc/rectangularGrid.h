#ifndef __RECTANGULARGRID_H
#define __RECTANGULARGRID_H

#include <aol.h>
#include <baseFunctionSet.h>
#include <gridBase.h>
#include <indexMapper.h>
#include <iterators.h>

namespace qc {

/** Class for uniform grids whose extensions do not have to be a power of 2. Moreover the grid can have
 *  different extent in the various coordinate directions
 */

template <qc::Dimension Dim>
class RectangularGrid : public GridStructure {
public:
  typedef RectangularGrid<Dim>          Self;
  typedef Self                          CubicGridType;

protected:
  double _h;
  qc::FastILexMapper< Dim > _indexMapper;

protected:
  void initializeIterators( ) {
    begin_it.getCurrentPosition().set ( 0, 0, 0 );

    if ( Dim == qc::QC_1D ) {
      end_it.getCurrentPosition().set ( 0, 1, 0 );
    } else if ( Dim == qc::QC_2D ) {
      end_it.getCurrentPosition().set ( 0, 0, 1 );
    } else if ( Dim == qc::QC_3D ) {
      end_it.getCurrentPosition().set ( 0, 0, this->getNumZ() - 1 );
    } else {
      throw aol::Exception ( "unsupported dimension", __FILE__, __LINE__ );
    }
  }

public:
  explicit RectangularGrid ( const RectangularGrid &other )
    : GridStructure ( other ),
      _h ( other._h ),
      begin_it ( other.getSize() ), end_it ( other.getSize() ) {
    _indexMapper.resize ( qc::CoordType ( this->getNumX(), this->getNumY(), this->getNumZ() ) );
    initializeIterators( );
  }

  explicit RectangularGrid ( const aol::Vec3<int> &size )
    : GridStructure ( size, Dim ),
      _h ( 1. / aol::Max ( size[0] - 1, size[1] - 1 , size[2] - 1 ) ),
      begin_it ( size ),
      end_it ( size ) {
    _indexMapper.resize ( qc::CoordType ( this->getNumX(), this->getNumY(), this->getNumZ() ) );
    initializeIterators( );
  }

  explicit RectangularGrid ( const qc::GridSize<Dim> &gridSize )
    : GridStructure ( aol::Vec3<int> ( gridSize.getNumX(), ( Dim > QC_1D ? gridSize.getNumY() : 1 ), ( Dim > QC_2D ? gridSize.getNumZ() : 1 ) ), Dim ),
      _h ( 1. / aol::Max ( static_cast<int>( gridSize.getNumX() - 1 ), ( Dim > QC_1D ? gridSize.getNumY() - 1 : 1 ), ( Dim > QC_2D ? gridSize.getNumZ() - 1 : 1 ) ) ),
      begin_it ( aol::Vec3<int> ( gridSize.getNumX(), ( Dim > QC_1D ? gridSize.getNumY() : 1 ), ( Dim > QC_2D ? gridSize.getNumZ() : 1 ) ) ),
      end_it ( aol::Vec3<int> ( gridSize.getNumX(), ( Dim > QC_1D ? gridSize.getNumY() : 1 ), ( Dim > QC_2D ? gridSize.getNumZ() : 1 ) ) ) {
    _indexMapper.resize ( qc::CoordType ( gridSize.getNumX(), ( Dim > QC_1D ? gridSize.getNumY() : 1 ), ( Dim > QC_2D ? gridSize.getNumZ() : 1 ) ) );
    initializeIterators( );
  }

  /**
   * This constructor mimics the behavior of the usual GridDefinition constructor.
   * Used for compatibility reasons.
   *
   * \author Berkels
   */
  explicit RectangularGrid ( const int GridDepth, const Dimension DimOfWorld = qc::QC_2D )
    : GridStructure ( aol::Vec3<int>( (( 1 << GridDepth ) + 1), (( 1 << GridDepth ) + 1), (Dim == QC_2D) ? 1 : (( 1 << GridDepth ) + 1) ), Dim ),
      _h ( 1. / ( aol::Max ( aol::Max ( this->getNumX() - 1, this->getNumY() - 1 ), this->getNumZ() - 1 ) ) ),
      begin_it ( this->getSize() ), // ???????
      end_it ( this->getSize() ) {
    if ( DimOfWorld != Dim )
      throw aol::Exception ( "DimOfWorld != Dim", __FILE__, __LINE__ );
    _indexMapper.resize ( qc::CoordType ( this->getNumX(), this->getNumY(), this->getNumZ() ) );
    initializeIterators( );
  }

  int getNumberOfElements() const {
    if(Dim == QC_3D) {
      return (this->getNumX()-1)*(this->getNumY()-1)*(this->getNumZ()-1);
    } else {
      return (this->getNumX()-1)*(this->getNumY()-1);
    }
  }

  const Self & getCubicGrid() const {
    return *this;
  }

  const Self & getFullGrid() const {
    return *this;
  }

  int getGridDepth() const {
    throw aol::Exception ( "RectangularGrid does not support getGridDepth()!", __FILE__, __LINE__ );
  }

  int getElementIndex ( const Element & /*El*/ ) const {
    throw aol::Exception ( "RectangularGrid does not support getElementIndex()!", __FILE__, __LINE__ );
  }

  // this iterator iterates over the elements of this grid
  struct OldAllElementIterator {
    typedef OldAllElementIterator _Self;
    typedef Element  IteratedType;
    typedef _Self    BeginType;
    typedef _Self    EndType;
    aol::Vec<3, int>  _size;

    OldAllElementIterator() : _size() {}

    explicit OldAllElementIterator ( const aol::Vec<3, int> &size ) : _size ( size )  { }

    qc::Element &getCurrentPosition() {
      return _cur;
    }

    const qc::Element &getCurrentPosition() const {
      return _cur;
    }

    OldAllElementIterator &operator= ( const OldAllElementIterator& Other ) {
      _size = Other._size;
      _cur = Other._cur;
      return *this;
    }

    inline OldAllElementIterator& operator++() {
      _cur.xref() ++;
      if ( _cur.x() == ( _size[0] - 1 ) ) {
        _cur.xref() = 0;
        _cur.yref() ++;
        if ( _cur.y() == ( _size[1] - 1 ) ) {
          _cur.yref() = 0;
          _cur.zref() ++;
        }
      }
      return *this;
    }

    inline OldAllElementIterator operator++ ( int ) {
      OldAllElementIterator copy ( *this );
      ++ ( *this );
      return copy;
    }

    qc::Element& operator*() {
      return _cur;
    }

    const qc::Element& operator*() const {
      return _cur;
    }

    const qc::Element* operator->() const {
      return &_cur;
    }

    inline bool operator!= ( const _Self& Other ) const {
      return ( _cur != Other._cur );
    }

    inline bool operator== ( const _Self& Other ) const {
      return ( _cur == Other._cur );
    }

  protected:
    qc::Element _cur;
    // End of internal class OldAllElementIterator
  };

  OldAllElementIterator begin_it;
  OldAllElementIterator end_it;

  double H() const {
    return _h;
  }

  void setH ( double h ) {
    _h = h;
  }

  //! return width = number of nodes in x direction
  int getWidth() const {
    return ( getNumX() );
  }

  //! return height = number of nodes in y direction
  int getHeight() const {
    return ( getNumY() );
  }

  //! return depth = number of nodes in z direction
  int getDepth() const {
    return ( getNumZ() );
  }

  const qc::FastILexMapper<Dim>& getIndexMapperRef() const {
    return ( _indexMapper );
  }

  //! return grid depth for initializing an Element iterator, returning -1 for rectangularGrids and to be overloaded on derived classes.
  virtual short getElementGridDepth ( ) const {
    return ( -1 );
  }

  //! return the volume fraction of the brick that is represented by this RectangularGrid relative to the unit cube
  double getVolumeFractionOf01Cube ( ) const {
    switch ( Dim ) {
    case QC_2D:
      return ( ( getNumX() - 1 ) * ( getNumY() - 1 ) * aol::Sqr ( _h ) );
    case QC_3D:
      return ( ( getNumX() - 1 ) * ( getNumY() - 1 ) * ( getNumZ() - 1 ) * aol::Cub ( _h ) );
    default:
      throw aol::UnimplementedCodeException ( "qc::RectangularGrid::getVolumeFractionOf01Cube not implemented for dimension != 2, 3", __FILE__, __LINE__ );
    }
  }

  /** Iterator that iterates over all elements in the grid
   *  \author Schwen
   */
  class FullElementIterator : public qc::RectangularIteratorBase< Dim, qc::Element > {
  public:
    explicit FullElementIterator ( const RectangularGrid<Dim> &grid ) : qc::RectangularIteratorBase< Dim, qc::Element > ( qc::Element(), qc::Element() ) {
      this->_lower.set ( 0, 0, 0, grid.getElementGridDepth() );
      this->_current.set ( 0, 0, 0, grid.getElementGridDepth() );
      this->_upper.set ( grid.getSize()[0] - 1, grid.getSize()[1] - 1, grid.getSize()[2] - 1, grid.getElementGridDepth() );
    }

    //! prefix increment operator.
    const qc::Element& operator++ ( ) {
      this->increment();
      return ( this->_current );
    }
  };

  /** Iterator that iterates over all nodes in the grid
   *  \author Schwen
   */
  class FullNodeIterator : public qc::RectangularIterator< Dim > {
  public:
    explicit FullNodeIterator ( const qc::RectangularGrid<Dim> &grid ) : qc::RectangularIterator<Dim> ( qc::CoordType(), qc::CoordType ( grid.getSize() ) ) {
    }
  };

  /** Iterator that iterates over all boundary nodes in the grid
   *  \warning currently uses slow but simple boundary node iterator
   *  \author Schwen
   */
  class FullBoundaryNodeIterator : public qc::RectangularBoundaryIterator< Dim > {
  public:
    explicit FullBoundaryNodeIterator ( const qc::RectangularGrid<Dim> &grid ) : qc::RectangularBoundaryIterator<Dim> ( qc::CoordType(), qc::CoordType( grid.getSize() ) ) {
    }
  };


};


/** Possible replacement for qc::GridDefinition
 *  \todo think about relation to RectangularGrid, move _indexMapper there?
 *  \todo think about placement of code of basis classes (better here ...)
 */
template <qc::Dimension Dim>
class CubicGrid : public RectangularGrid<Dim> {
protected:
  int _depth;

public:
  explicit CubicGrid ( const int Depth ) : RectangularGrid<Dim> ( Depth, Dim ), _depth ( Depth ) {
  }

  //! numX, numY and numZ are the same
  inline int getNumXYZ ( ) const {
    return ( this->getNumX() );
  }

  int getGridDepth() const {
    return ( _depth );
  }

  int getNumberOfBoundaryNodes() const {
    if ( Dim == qc::QC_2D ) {
      return ( 4 * getNumXYZ() - 4 );
    } else if ( Dim == qc::QC_3D ) {
      return ( 6 * ( getNumXYZ() - 1 ) * ( getNumXYZ() - 1 ) + 2 );
    } else {
      throw aol::Exception ( "unsupported dimension", __FILE__, __LINE__ );
      return ( - 1 );
    }
  }

  virtual short getElementGridDepth ( ) const {
    return ( _depth );
  }

};

} // end of namespace qc.

#endif
