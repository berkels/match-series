#ifndef __GRIDBASE_H
#define __GRIDBASE_H

#include <op.h>
#include <quoc.h>
#include <gridSize.h>
#include <qmElement.h>

namespace qc {

template <Dimension dim> class ElementInfoTrait { };
template <> class ElementInfoTrait<QC_2D> {
public:
  enum { numNodesPerElement = 4 };
  enum { numElementsAdjacentToNode = 4 };
};

template <> class ElementInfoTrait<QC_3D> {
public:
  enum { numNodesPerElement = 8 };
  enum { numElementsAdjacentToNode = 8 };
};


struct NodeOffsetsStruct {
  aol::Vec3<int> _offsets[8];
  NodeOffsetsStruct( ) {
    _offsets[0].set ( 0, 0, 0 );
    _offsets[1].set ( 1, 0, 0 );
    _offsets[2].set ( 0, 1, 0 );
    _offsets[3].set ( 1, 1, 0 );

    _offsets[4].set ( 0, 0, 1 );
    _offsets[5].set ( 1, 0, 1 );
    _offsets[6].set ( 0, 1, 1 );
    _offsets[7].set ( 1, 1, 1 );
  }
};

extern NodeOffsetsStruct NodeOffsets;


//! traverses all nodes of a given element
template <Dimension Dim>
class ElementNodeIterator {

public:
  ElementNodeIterator ( const Element &el, int idx = 0 )
      : _el ( el ), _curIndex ( idx ) {
    fill();
  }

  bool operator== ( const ElementNodeIterator &it ) const { return ( it._curIndex == _curIndex ); }

  bool operator!= ( const ElementNodeIterator &it ) const { return ( it._curIndex != _curIndex ); }

  CoordType& operator*() { return _cur; }

  CoordType* operator->() { return &_cur; }

  inline ElementNodeIterator& operator++() {
    ++_curIndex;
    fill( );
    return *this;
  }
protected:

  void fill( ) {
    if ( _curIndex < ElementInfoTrait<Dim>::numNodesPerElement ) {
      aol::Vec3<int> &o = NodeOffsets._offsets[_curIndex];
      for ( int c = 0; c < 3; c++ ) {
        _cur[c] = _el[c] + o[c];
      }
    }
  }

  const Element &_el;
  int _curIndex;
  CoordType _cur;
};


//! Some iterators can compare themselves not to each other, but only to an EndElement. This is a dummy class.
class EndElement { };

/**
 * This class collects some very basic grid properties (currently of the grid
 * classes GridDefinition and RectangularGrid). This improves the interchangability
 * of these grid types.
 *
 * \author Berkels
 */
class GridStructure {
public:
  typedef Element ElementType;
  static const aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE;

  /**
   * Iterates over all nodes in the grid.
   *
   * \todo Does the same as qc::GridDefinition::OldFullNodeIterator. Can the latter be deleted?
   */
  struct OldAllNodeIterator {
    typedef OldAllNodeIterator Self;
    aol::Vec<3, int>  _size;

    OldAllNodeIterator() : _size() {}

    explicit OldAllNodeIterator ( const aol::Vec<3, int> &size ) : _size ( size )  { }

    qc::Element &getCurrentPosition() {
      return _cur;
    }

    const qc::Element &getCurrentPosition() const {
      return _cur;
    }

    OldAllNodeIterator &operator= ( const OldAllNodeIterator& Other ) {
      _size = Other._size;
      _cur = Other._cur;
      return *this;
    }

    inline OldAllNodeIterator& operator++() {
      _cur.xref() ++;
      if ( _cur.x() == ( _size[0] ) ) {
        _cur.xref() = 0;
        _cur.yref() ++;
        if ( _cur.y() == ( _size[1] ) ) {
          _cur.yref() = 0;
          _cur.zref() ++;
        }
      }
      return *this;
    }

    inline OldAllNodeIterator operator++ ( int ) {
      OldAllNodeIterator copy ( *this );
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

    inline bool operator!= ( const Self& Other ) const {
      return ( _cur != Other._cur );
    }

    inline bool operator== ( const Self& Other ) const {
      return ( _cur == Other._cur );
    }

protected:
    qc::Element _cur;
  }; // End of internal class OldAllNodeIterator

  /**
   * Iterates over all boundary nodes in the grid.
   *
   * \note Implementation heavily based on qc::GridDefinition::OldFullBoundaryNodeIterator, but unlike the latter
   *       it supports non-quadratic grids and uses the "new" iterator syntax.
   *
   * \todo Deleted qc::GridDefinition::OldFullBoundaryNodeIterator?
   */
  class AllBoundaryNodeIterator {
  public:
    typedef AllBoundaryNodeIterator Self;
    typedef CoordType        IteratedType;

    AllBoundaryNodeIterator ( const qc::GridStructure &Grid )
      : _cur ( 0, 0, 0 ),
        _size ( Grid.getNumX(), Grid.getNumY(), ( Grid.getDimOfWorld() == QC_3D ) ? Grid.getNumZ() : 1 ),
        done ( false ),
        vert ( false ),
        topbottom ( false ),
        frontback ( true ),
        rightleft ( true ) {

      if ( Grid.getDimOfWorld() == QC_1D )
        throw aol::Exception ( "qc::GridStructure::AllBoundaryNodeIterator not implemented in 1D", __FILE__, __LINE__ );

      if ( ( aol::Min( _size[0], _size[1] ) < 3 ) || ( ( Grid.getDimOfWorld() == QC_3D ) && ( _size[2] < 3 ) ) )
        throw aol::Exception ( "qc::GridStructure::AllBoundaryNodeIterator does not work for grids with less than 3 nodes in any direction", __FILE__, __LINE__ );
    }

    const IteratedType& operator*() const {
      return _cur;
    }

    const IteratedType* operator->() const {
      return &_cur;
    }

    Self & operator++ ( ) {
      if ( _size[2] == 1 ) {
        if ( !vert ) {                       // 2D-Iterator-BEGIN...
          _cur.xref() ++;
          if ( _cur.x() == _size[0] ) {
            if ( _cur.y() == _size[1] - 1 ) {
              vert = true;
              _cur.xref() = 0;
              _cur.yref() = 1;
            } else {
              _cur.xref() = 0;
              _cur.yref() = _size[1] - 1;
            }
          }
        } else {
          _cur.yref() ++;
        }
        while ( vert && !done && _cur.y() == _size[1] - 1 ) {
          if ( _cur.x() == 0 ) {
            _cur.xref() = _size[0] - 1;
            _cur.yref() = 1;
          } else {
            done = true;
          }
        }
      }                                    // 2D-Iterator-END
      else {

        if ( !topbottom ) {                // 3D-Iterator-BEGIN...
          _cur.xref() ++;
          if ( _cur.xref() == _size[0] ) {
            _cur.yref() ++;
            _cur.xref() = 0;
            if ( _cur.yref() == _size[1] ) {
              _cur.yref() = 0;
              if ( _cur.zref() == 0 ) {
                _cur.zref() = _size[2] - 1;
              } else {
                topbottom = true;
                frontback = false;
                _cur.xref() = -1;
                _cur.yref() = 0;
                _cur.zref() = 1;
              }
            }
          }
        }

        if ( !frontback ) {
          _cur.xref() ++;
          if ( _cur.xref() == _size[0] ) {
            _cur.zref() ++;
            _cur.xref() = 0;
            if ( _cur.zref() == _size[2] - 1 ) {
              _cur.zref() = 1;
              if ( _cur.yref() == 0 ) {
                _cur.yref() = _size[1] - 1;
              } else {
                frontback = true;
                rightleft = false;
                _cur.xref() = 0;
                _cur.yref() = 0;
                _cur.zref() = 1;
              }
            }
          }
        }

        if ( !rightleft ) {
          _cur.yref() ++;
          if ( _cur.yref() == _size[1] - 1 ) {
            _cur.zref() ++;
            _cur.yref() = 1;
            if ( _cur.zref() == _size[2] - 1 ) {
              _cur.zref() = 1;
              if ( _cur.xref() == 0 ) {
                _cur.xref() = _size[0] - 1;
              } else {
                rightleft = false;
                done = true;
              }
            }
          }
        }                                     // 3D-Iterator-END
      }
      return *this;
    }

    bool atEnd () const {
      return done;
    }

    bool notAtEnd () const {
      return !done;
    }

  protected:
    IteratedType _cur;
    const aol::Vec3<short> _size;
    bool done, vert, topbottom, frontback, rightleft;
  }; // End of internal class AllBoundaryNodeIterator

private:
  const Dimension _dimOfWorld;
  const aol::Vec3<int> _size;

public:
  OldAllNodeIterator _nBeginIt;
  OldAllNodeIterator _nEndIt;

  explicit GridStructure ( const GridStructure &Grid )
      : _dimOfWorld ( Grid._dimOfWorld ),
      _size ( Grid._size ),
      _nBeginIt ( Grid._nBeginIt ),
      _nEndIt ( Grid._nEndIt ) {
  }

  explicit GridStructure ( const aol::Vec3<int> &Size, const Dimension DimOfWorld )
      : _dimOfWorld ( DimOfWorld ),
      _size ( Size ),
      _nBeginIt ( Size ),
      _nEndIt ( Size ) {
    _nBeginIt.getCurrentPosition().set ( 0, 0, 0, 0 );
    // Due to the design of OldAllNodeIterator::operator++(), _nEndIt is the same for 1D and 2D
    if ( ( _dimOfWorld == qc::QC_1D ) || ( _dimOfWorld == qc::QC_2D ) ) {
      _nEndIt.getCurrentPosition().set ( 0, 0, 1, 0 );
    } else if ( _dimOfWorld == qc::QC_3D ) {
      _nEndIt.getCurrentPosition().set ( 0, 0, getNumZ(), 0 );
    } else {
      throw aol::Exception ( "unsupported dimension" );
    }
  }

  template <qc::Dimension DimOfWorld>
  explicit GridStructure ( const GridSize<DimOfWorld> & gridSize )
      : _dimOfWorld ( DimOfWorld ),
      _size ( aol::Vec3<int> ( gridSize.getNumX(),
                               DimOfWorld > QC_1D ? gridSize.getNumY() : 1,
                               DimOfWorld > QC_2D ? gridSize.getNumZ() : 1 ) ),
      _nBeginIt ( _size ),
      _nEndIt ( _size ) {
    _nBeginIt.getCurrentPosition().set ( 0, 0, 0, 0 );
    if ( DimOfWorld == qc::QC_2D ) {
      _nEndIt.getCurrentPosition().set ( 0, 0, 1, 0 );
    } else if ( DimOfWorld == qc::QC_3D ) {
      _nEndIt.getCurrentPosition().set ( 0, 0, getNumZ(), 0 );
    } else {
      throw aol::Exception ( "unsupported dimension" );
    }
  }

  virtual ~GridStructure() {}

  inline const aol::Vec3<int>& getSize() const {
    return _size;
  }

  inline int getNumX() const {
    return _size[0];
  }

  inline int getNumY() const {
    return _size[1];
  }

  inline int getNumZ() const {
    return _size[2];
  }

  //! Returns 1, 2 or 3.
  inline Dimension getDimOfWorld() const {
    return _dimOfWorld;
  }

  /** Return whether error estimator exists or not
   */
  virtual bool isAdaptive() const {
    return false;
  }

  //! \todo Rename, document, arguments
  //! Checks whether a node depends on parent nodes or not.
  virtual int checkForHangingNode ( const Element &, int ) const {
    return -1;
  }

  /** Returns number of nodes in the grid
   * @todo for adaptive case
   */
  int getNumberOfNodes() const {
    if ( isAdaptive() ) {
      static bool adaptiveNumNodesWarningPrinted = false;
      if ( ! adaptiveNumNodesWarningPrinted ) {
        cerr << "qc::GridStructure::getNumberOfNodes returns number of nodes for full grid in case of adaptive grids" << endl;
        adaptiveNumNodesWarningPrinted = true;
      }
    }

    switch ( this->getDimOfWorld() ) {
      case QC_1D :
        return ( getNumX() );
      case QC_2D :
        return ( getNumX() * getNumY() );
      case QC_3D :
        return ( getNumX() * getNumY() * getNumZ() );
      default:
        throw aol::Exception ( "qc::GridStructure::getNumberOfNodes: unsupported dimension", __FILE__, __LINE__ );
        return -1;
    }
  }

  /** Returns number of dofs (=number of nodes for non-adaptive grids) in the grid
   * @todo for adaptive case
   */
  int getNumberOfDofs() const {
    if ( isAdaptive() )
      throw aol::Exception ( "getNumberOfDofs not implemented for adaptive grids!", __FILE__, __LINE__ );
    return getNumberOfNodes();
  }

  /** Returns number of boundary nodes in the grid
   * @todo for adaptive case
   */
  int getNumberOfBoundaryNodes() const {
    if ( isAdaptive() ) {
      throw aol::Exception ( "qc::GridStructure::getNumberOfBoundaryNodes not implemented for adaptive grids", __FILE__, __LINE__ );
      return -1;
    }

    // # boundary nodes = # nodes - # inner nodes
    switch ( this->getDimOfWorld() ) {
      case QC_1D :
        return ( 2 );
      case QC_2D :
        return ( getNumX() * getNumY() - ( getNumX() - 2 ) * ( getNumY() - 2  ) );
      case QC_3D :
        return ( getNumX() * getNumY() * getNumZ() - ( getNumX() - 2 ) * ( getNumY() - 2 ) * ( getNumZ() - 2 )  );
      default:
        throw aol::Exception ( "unsupported dimension", __FILE__, __LINE__ );
        return -1;
    }
  }

  //! checks whether pos lies inside [(0,0,0), size]
  inline bool isAdmissibleNode ( const qc::CoordType & pos ) const {
    return (    pos[0] >= 0 && pos[0] < _size[0]
                && pos[1] >= 0 && pos[1] < _size[1]
                && pos[2] >= 0 && pos[2] < _size[2] );
  }

};

class GridDefinition : public GridStructure {
public:
  typedef GridDefinition          Self;

  typedef Self                    BeginIterType;
  typedef EndElement              EndIterType;
  typedef Self                    CubicGridType;
  typedef Self                    FullGridType;
  typedef Element                 ElementType;


  const EndElement _end_element;

  int  *h;         /**< h[i] giving element sizes on level i */
  double *h_x;     /**< Computational width of the elements on levels */
  double *h_y;     /**< Computational height of the elements on levels */
  double *h_z;     /**< Computational depth of the elements on levels */
  double *h_sqr;   /**< In 2D: h_sqr = h_x * h_y, in 3D: NULL */
  double *h_cub;   /**< In 3D: h_cub = h_x * h_y * h_z, in 2D: NULL */

// private:
//   GridDefinition ( const GridDefinition&  ) : _end_element() {};

public:
  explicit GridDefinition ( const GridDefinition &Grid  )
      : GridStructure ( Grid ),
      _end_element(),
      begin_it ( ( 1 << Grid.getGridDepth() ) + 1 ),
      end_it ( ( 1 << Grid.getGridDepth() ) + 1 ),
      gridDepth ( Grid.getGridDepth() ),
      width ( ( 1 << gridDepth ) + 1 ) {
    init ( );
  }

  explicit GridDefinition ( const int GridDepth, const Dimension DimOfWorld = qc::QC_2D ) :
      GridStructure ( aol::Vec3<int> ( ( ( 1 << GridDepth ) + 1 ),
                                       ( DimOfWorld == QC_1D ) ? 1 : ( ( 1 << GridDepth ) + 1 ),
                                       ( DimOfWorld == QC_3D ) ? ( ( 1 << GridDepth ) + 1 ) : 1 ),
                      DimOfWorld ),
      _end_element(),
      begin_it ( ( 1 << GridDepth ) + 1 ),
      end_it ( ( 1 << GridDepth ) + 1 ) ,
      gridDepth ( GridDepth ),
      width ( ( 1 << gridDepth ) + 1 ) {

    init ( );
  }

  explicit GridDefinition ( const aol::Vec3< int > &size )
      : GridStructure ( size, size[2] == 1 ? ( size[1] == 1 ? QC_1D : QC_2D ) : QC_3D ),
      _end_element(),
      begin_it ( ( 1 << qc::logBaseTwo ( size[0] - 1 ) ) + 1 ),
      end_it ( ( 1 << qc::logBaseTwo ( size[0] - 1 ) ) + 1 ),
      gridDepth ( qc::logBaseTwo ( size[0] - 1 ) ),
      width ( ( 1 << gridDepth ) + 1 ) {

    if ( ( size[0] != size[1] ) || ( ( this->getDimOfWorld() == QC_3D ) && ( size[0] != size[2] ) ) ) {
      cerr << "\nsize[0]=" << size[0] << "  size[1]=" << size[1] << "  size[2]=" << size[2] << "   dimOfWorld=" << this->getDimOfWorld() << endl;
      throw aol::Exception ( "not possible to use GridDefinition for a nonquadratic domain, try RectangularGrid!", __FILE__, __LINE__ );
    }
    if ( size[0] != ( 1 << gridDepth ) + 1 )
      throw aol::Exception ( "Your domain is too small or too large (the appropriate width is 2^(gridDepth) + 1)!", __FILE__, __LINE__ );
    init ( );
  }

  template <Dimension Dim>
  explicit GridDefinition ( const GridSize<Dim> & size )
      : GridStructure ( size ),
      _end_element(),
      begin_it ( ( 1 << qc::logBaseTwo ( size[0] - 1 ) ) + 1 ),
      end_it ( ( 1 << qc::logBaseTwo ( size[0] - 1 ) ) + 1 ),
      gridDepth ( qc::logBaseTwo ( size[0] - 1 ) ),
      width ( ( 1 << gridDepth ) + 1 ) {

    if ( ( size[0] != size[1] ) || ( ( this->getDimOfWorld() == QC_3D ) && ( size[0] != size[2] ) ) ) {
      cerr << "\nsize[0]=" << size[0] << "  size[1]=" << size[1] << "  size[2]=" << size[2] << "   dimOfWorld=" << this->getDimOfWorld() << endl;
      throw aol::Exception ( "not possible to use GridDefinition for a nonquadratic domain, try RectangularGrid!", __FILE__, __LINE__ );
    }
    if ( size[0] != ( 1 << gridDepth ) + 1 )
      throw aol::Exception ( "Your domain is too small or too large (the appropriate width is 2^(gridDepth) + 1)!", __FILE__, __LINE__ );
    init ( );
  }

  void init ( ) {
    h   = new int[gridDepth+1];
    h_x = new double[gridDepth+1];
    h_y = new double[gridDepth+1];
    h_z = new double[gridDepth+1];

    if ( this->getDimOfWorld() == QC_3D ) {
      h_sqr = NULL;
      h_cub = new double[gridDepth+1];
    } else if ( this->getDimOfWorld() == QC_2D ) {
      h_sqr = new double[gridDepth+1];
      h_cub = NULL;
    } else {
      h_sqr = NULL;
      h_cub = NULL;
    }

    for ( int i = 0; i <= gridDepth; i++ ) {
      h[i] = ( 1 << ( gridDepth - i ) );

      h_x[i] = h_y[i] = h[i] / static_cast<double> ( width - 1 );

      if ( this->getDimOfWorld() == QC_3D ) {
        h_z[i] = h_y[i];
        h_cub[i] = h_x[i] * h_y[i] * h_z[i];
      } else if ( this->getDimOfWorld() == QC_2D ) {
        h_z[i] = 0;
        h_sqr[i] = h_x[i] * h_y[i];
      } else {
        h_z[i] = 0;
        h_y[i] = 0;
      }
    }

    begin_it.getCurrentPosition().set ( 0, 0, 0, gridDepth );
    if ( this->getDimOfWorld() == QC_1D ) {
      end_it.getCurrentPosition().set ( 0, 1, 0, gridDepth );
    } else if ( this->getDimOfWorld() == QC_2D ) {
      end_it.getCurrentPosition().set ( 0, 0, 1, gridDepth );
    } else if ( this->getDimOfWorld() == QC_3D ) {
      end_it.getCurrentPosition().set ( 0, 0, width - 1, gridDepth );
    }


    if ( isAdaptive() ) {
      uniformGridIndexOffset = NULL;
    } else {
      int n;
      switch ( this->getDimOfWorld() ) {
        case QC_1D:
          uniformGridIndexOffset = new int [3];
          uniformGridIndexOffset[0] =     - 1 ;
          uniformGridIndexOffset[1] =       0 ;
          uniformGridIndexOffset[2] =       1 ;
          break;
        case QC_2D:
          n = ( 1 << gridDepth ) + 1;
          uniformGridIndexOffset = new int[9];
          uniformGridIndexOffset[0] = - n - 1 ;
          uniformGridIndexOffset[1] = - n     ;
          uniformGridIndexOffset[2] = - n + 1 ;
          uniformGridIndexOffset[3] =     - 1 ;
          uniformGridIndexOffset[4] =       0 ;
          uniformGridIndexOffset[5] =       1 ;
          uniformGridIndexOffset[6] =   n - 1 ;
          uniformGridIndexOffset[7] =   n     ;
          uniformGridIndexOffset[8] =   n + 1 ;
          break;
        case QC_3D:
          n = ( 1 << gridDepth ) + 1;
          uniformGridIndexOffset = new int[27];
          uniformGridIndexOffset[ 0] = - n * n - n - 1 ;
          uniformGridIndexOffset[ 1] = - n * n - n     ;
          uniformGridIndexOffset[ 2] = - n * n - n + 1 ;
          uniformGridIndexOffset[ 3] = - n * n     - 1 ;
          uniformGridIndexOffset[ 4] = - n * n         ;
          uniformGridIndexOffset[ 5] = - n * n     + 1 ;
          uniformGridIndexOffset[ 6] = - n * n + n - 1 ;
          uniformGridIndexOffset[ 7] = - n * n + n     ;
          uniformGridIndexOffset[ 8] = - n * n + n + 1 ;
          uniformGridIndexOffset[ 9] =       - n - 1 ;
          uniformGridIndexOffset[10] =       - n     ;
          uniformGridIndexOffset[11] =       - n + 1 ;
          uniformGridIndexOffset[12] =           - 1 ;
          uniformGridIndexOffset[13] =             0 ;
          uniformGridIndexOffset[14] =             1 ;
          uniformGridIndexOffset[15] =         n - 1 ;
          uniformGridIndexOffset[16] =         n     ;
          uniformGridIndexOffset[17] =         n + 1 ;
          uniformGridIndexOffset[18] =   n * n - n - 1 ;
          uniformGridIndexOffset[19] =   n * n - n     ;
          uniformGridIndexOffset[20] =   n * n - n + 1 ;
          uniformGridIndexOffset[21] =   n * n     - 1 ;
          uniformGridIndexOffset[22] =   n * n         ;
          uniformGridIndexOffset[23] =   n * n     + 1 ;
          uniformGridIndexOffset[24] =   n * n + n - 1 ;
          uniformGridIndexOffset[25] =   n * n + n     ;
          uniformGridIndexOffset[26] =   n * n + n + 1 ;
          break;
        default:
          uniformGridIndexOffset = NULL;
      }
    }
  }

  virtual ~GridDefinition( ) {
    if ( h != NULL ) delete[] h;
    if ( h_x != NULL ) delete[] h_x;
    if ( h_y != NULL ) delete[] h_y;
    if ( h_z != NULL ) delete[] h_z;
    if ( h_sqr != NULL ) delete[] h_sqr;
    if ( h_cub != NULL ) delete[] h_cub;
    if ( uniformGridIndexOffset != NULL ) delete[] uniformGridIndexOffset;
  }

private:
  GridDefinition& operator= ( const GridDefinition& ) {
    throw aol::UnimplementedCodeException ( "qc::GridDefinition::operator= not implemented", __FILE__, __LINE__ );
    return ( *this );
  }

public:
  double H() const {
    return aol::Max ( aol::Max ( h_x[ gridDepth ], h_y[ gridDepth ] ), h_z[ gridDepth ] );
  }

  const EndElement &end( ) const {
    return _end_element;
  }

  const GridDefinition &begin() const {
    return *this;
  }


  //! Returns element that contains a given point without checking
  //! whether the point is inside [0,1]^3
  void fastGetElementByPoint ( aol::Vec3<double> point, Element &el ) const {
    el.xref() = static_cast<short> ( point[0] / H() );
    el.yref() = static_cast<short> ( point[1] / H() );
    el.zref() = static_cast<short> ( point[2] / H() );
  }

  //! Returns element that contains a given point and checks
  //! whether the point is inside [0,1]^3 or not
  void getElementByPoint ( aol::Vec3<double> point, Element &el ) const {
    if ( point[0] >= 0 && point[0] <= 1 &&
         point[1] >= 0 && point[1] <= 1 &&
         point[2] >= 0 && point[2] <= 1 )  {
      el.xref() = static_cast<short> ( point[0] / H() );
      el.yref() = static_cast<short> ( point[1] / H() );
      el.zref() = static_cast<short> ( point[2] / H() );
    } else throw aol::Exception ( "getElement: point is NOT inside [0,1]^3" );
  }


  // --------------------- iterators -------------------------------------------



  template <qc::Dimension Dim> static
  ElementNodeIterator<Dim> elementNodeBegin ( const qc::Element &el ) {
    return ElementNodeIterator<Dim> ( el, 0 );
  }

  template <qc::Dimension Dim> static
  ElementNodeIterator<Dim> elementNodeEnd ( const qc::Element &el ) {
    return ElementNodeIterator<Dim> ( el, ElementInfoTrait<Dim>::numNodesPerElement );
  }


  //! this is an element with an additional Vec2 which gives you the nodes of the vertex
  //! that is at the boundary.
  struct ElementOnBoundary : public Element {
    ElementOnBoundary ( short X, short Y, short Z = 0, short Level = -1, unsigned char Type = 9, int Index = 0  )
        : Element ( X, Y, Z, Level, Type, Index ) {
    }

    aol::Vec2<short> _locNode;      // stores the nodes of the boundary-vertex
  };

  //! iterator traverses the boundary elements of a 2d-grid (i.e. the four sides of the unit-square)
  class OldBoundaryHyperfaceIterator2D {
  public:
    typedef OldBoundaryHyperfaceIterator2D         Self;
    typedef ElementOnBoundary IteratedType;
    typedef Self             BeginType;
    typedef Self             EndType;

    OldBoundaryHyperfaceIterator2D ( const GridDefinition &Grid )
        : _cur ( 0, 0, 0, Grid.getGridDepth() ),
        _num ( Grid.getWidth() - 1 ), done ( false ), vert ( false ) {


      if ( Grid.getDimOfWorld() != QC_2D ) {
        throw aol::Exception ( "only implemented for 2D", __FILE__, __LINE__ );
      }
      _cur._locNode.set ( 0, 1 );
    }

    IteratedType & operator*() { return _cur; }

    IteratedType * operator->() { return &_cur; }

    inline Self & operator++() {
      if ( !vert ) {                       // 2D-Iterator-BEGIN...
        _cur.xref() ++;
        if ( _cur.x() == _num ) {
          if ( _cur.y() == _num - 1 ) {
            vert = true;
            _cur.xref() = 0;
            _cur.yref() = 0;
            _cur._locNode.set ( 0, 2 );
          } else {
            _cur.xref() = 0;
            _cur.yref() = _num - 1;
            _cur._locNode.set ( 2, 3 );
          }
        }
      } else {
        _cur.yref() ++;
      }
      while ( vert && !done && _cur.y() == _num ) {
        if ( _cur.x() == 0 ) {
          _cur.xref() = _num - 1;
          _cur.yref() = 0;
          _cur._locNode.set ( 1, 3 );
        } else {
          done = true;
        }
      }
      return *this;
    }

    operator bool() const {
      return !done;
    }

  protected:
    IteratedType _cur;
    short  _num;
    bool done, vert;
  };


  const OldBoundaryHyperfaceIterator2D ebegin() const {
    return OldBoundaryHyperfaceIterator2D ( *this );
  }

  struct OldAllElementIterator {
    typedef OldAllElementIterator Self;
    typedef Element  IteratedType;
    typedef Self    BeginType;
    typedef Self    EndType;

    int _width;

    OldAllElementIterator( ) : _width ( -1 )  {}

    OldAllElementIterator ( int Width ) {
      _width = Width - 1;
      _cur.setZero( );
    }

    IteratedType & getCurrentPosition( ) {
      return _cur;
    }

    OldAllElementIterator &operator= ( const BeginType & Other ) {
      _width = Other._width;
      _cur = Other._cur;

      return *this;
    }

    inline OldAllElementIterator& operator++() {
      _cur.xref() ++;
      if ( _cur.x() == _width ) {
        _cur.xref() = 0;
        _cur.yref() ++;
        if ( _cur.y() == _width ) {
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

    IteratedType & operator*() {
      return _cur;
    }

    const IteratedType & operator*() const {
      return _cur;
    }

    const IteratedType * operator->() const {
      return &_cur;
    }

    inline bool operator!= ( const Self & Other ) const {
      return ( _cur != Other._cur );
    }

    inline bool operator== ( const Self & Other ) const {
      return ( _cur == Other._cur );
    }

protected:
    Element _cur;
  };



  OldAllElementIterator begin_it;
  OldAllElementIterator end_it;

  // ----------------- end iterators -------------------------------------------

  //! Returns grid depth.
  int getGridDepth() const {
    return gridDepth;
  }

  /** Returns number of nodes of the grid in one dimension
   */
  int getWidth() const {
    return getNumX();
  }
  int getHeight() const {
    return getNumY();
  }

  /** Returns number of elements in the grid
   * @todo for adaptive case
   */
  int getNumberOfElements() const {
    if ( !isAdaptive() ) return ( this->getDimOfWorld() == QC_2D ) ? aol::Sqr ( width - 1 ) : ( ( this->getDimOfWorld() == QC_1D ) ? width - 1 : aol::Cub ( width - 1 ) );
    else {
      throw aol::Exception ( "getNumberOfElements not implemented for adaptive grids" );
      return -1;
    }
  }

  /** Returns index offsets in case of uniform grid
   */
  int getUniformGridIndexOffset ( int i ) const {
    if ( uniformGridIndexOffset ) {
      return uniformGridIndexOffset[i];
    } else {
      throw aol::Exception ( "not a uniform grid", __FILE__, __LINE__ );
    };
  }

  /** Return the grid coordinate for a given global index
   */
  void indexToCoordGlobal ( int index, aol::Vec3<int> &coord ) const {
    // TP: Please note: The jumps over the case labels is intended. Please do not add breaks. Then it won't work any more!
    switch ( this->getDimOfWorld() ) {
      case QC_3D: {
        const int f = aol::Sqr ( width );
        coord[2] = index / f;
        index = index % f;
      }
      /* fall-thru */
      case QC_2D: {
        const int f = width;
        coord[1] = index / f;
        index = index % f;
      }
      /* fall-thru */
      case QC_1D:
        coord[0] = index;
        break; // ok here
      default:
          throw aol::Exception ( "unsupported dimension", __FILE__, __LINE__ );
    }
  }

  //! returns global index of the given element.
  //! TODO: replace implementation of this function
  //!       by call to an ElementIndexMapper (to be written).
  //!       At the moment, the configurator's index mapper
  //!       has this ability, but it element to index mapping
  //!       should belong to the grid.
  int getElementIndex ( const Element & El ) const {
    return El.x() + this->getNumX() * El.y() + this->getNumX() * this->getNumY() * El.z();
  }

  //! Checks whether there exists a node in the grid of given Level,
  //! with coordinates Coords.
  int nodeInGrid ( const aol::Vec3<short> &Coords, int Level ) const {
    int elSize = 1 << ( gridDepth - Level );
    switch ( this->getDimOfWorld() ) {
      case QC_1D:
        return ( ( ( Coords.x() % elSize ) == 0 ) );
        break;
      case QC_2D:
        return ( ( ( Coords.x() % elSize ) == 0 ) &&
                 ( ( Coords.y() % elSize ) == 0 ) );
        break;
      case QC_3D:
        return ( ( ( Coords.x() % elSize ) == 0 ) &&
                 ( ( Coords.y() % elSize ) == 0 ) &&
                 ( ( Coords.z() % elSize ) == 0 ) );
        break;
      default:
          throw aol::Exception ( "unsupported dimension", __FILE__, __LINE__ );
    }
    return 0;
  }

  //! Checks whether the passed Element is valid.
  //! \todo Adpaptive grids?
  int elementInGrid ( const Element &El ) const {
    if ( !nodeInGrid ( El, El.level() ) ) return 0;
    int elSize = 1 << ( gridDepth - El.level() );

    // Preusser corrected an error here: this->width -1 was replaced by
    // this->width
    switch ( this->getDimOfWorld() ) {
      case QC_1D:
        return ( El.x() >= 0 &&
                 El.x() + elSize < this->width );
        break;
      case QC_2D:
        return ( El.x() >= 0 &&
                 El.x() + elSize < this->width &&
                 El.y() >= 0 &&
                 El.y() + elSize < this->width );
        break;
      case QC_3D:
        return ( El.x() >= 0 &&
                 El.x() + elSize < this->width &&
                 El.y() >= 0 &&
                 El.y() + elSize < this->width &&
                 El.z() >= 0 &&
                 El.z() + elSize < this->width  );
        break;
      default:
        throw aol::Exception ( "unsupported dimension", __FILE__, __LINE__ );
    }
    return 0;
  }

  class OldFullNodeIterator {
  public:
    typedef OldFullNodeIterator Self;
    typedef CoordType        IteratedType;
    typedef GridDefinition   BeginType;
    typedef EndElement       EndType;

    OldFullNodeIterator( ) : _numX ( -1 ), _numY ( -1 ), _numZ ( -1 ) {}

    Self & operator= ( const BeginType & Grid ) {
      _numX = Grid.getNumX();
      _numY = Grid.getNumY();
      _numZ = Grid.getNumZ();

      _cur.setZero( );
      done = false;
      return *this;
    }

    IteratedType& operator*() {
      return _cur;
    }

    IteratedType* operator->() {
      return &_cur;
    }

    inline Self & operator++ ( ) {
      _cur.xref() ++;
      if ( _cur.x() == _numX ) {
        _cur.xref() = 0;
        _cur.yref() ++;
        if ( _cur.y() == _numY ) {
          _cur.yref() = 0;
          _cur.zref() ++;
          if ( _cur.z() == _numZ ) {
            done = true;
          }
        }
      }
      return *this;
    }

    int index () {
      return ( _numY * _cur.z() + _cur.y() ) * _numX + _cur.x();
    }

    bool operator== ( const EndType & ) const {
      return done;
    }

    bool operator!= ( const EndType & ) const {
      return !done;
    }

    bool atEnd () const {
      return done;
    }

  protected:
    IteratedType _cur;
    short _numX, _numY, _numZ;
    bool done;
  };

  class OldFullBoundaryNodeIterator {
  public:
    typedef OldFullBoundaryNodeIterator Self;
    typedef CoordType        IteratedType;
    typedef GridDefinition   BeginType;
    typedef EndElement       EndType;

    OldFullBoundaryNodeIterator( ) : _cur ( -1, -1, -1 ), _width ( -1 ), _numZ ( -1 ), done ( false ), vert ( false ), obenunten ( false ), vornehinten ( false ), rechtslinks ( false ) {} // make sure all members are initialized

    Self &operator= ( const BeginType & Grid ) {
      _width = Grid.getWidth();

      if ( Grid.getGridDepth() < 1 ) {
        throw aol::Exception ( "qc::GridDefinition::OldFullBoundaryNodeIterator does not work for grid depth < 1", __FILE__, __LINE__ );
      }

      if ( Grid.getDimOfWorld() == QC_2D ) {
        _numZ = 1;
      } else if ( Grid.getDimOfWorld() == QC_3D ) {
        _numZ = _width;
      };

      _cur.setZero( );
      done = false;
      vert = false;
      obenunten = false;
      vornehinten = true;
      rechtslinks = true;
      return *this;
    }

    IteratedType & operator*() {
      return _cur;
    }

    IteratedType * operator->() {
      return &_cur;
    }

    inline Self & operator++ ( ) {
      if ( _numZ == 1 ) {
        if ( !vert ) {                       // 2D-Iterator-BEGIN...
          _cur.xref() ++;
          if ( _cur.x() == _width ) {
            if ( _cur.y() == _width - 1 ) {
              vert = true;
              _cur.xref() = 0;
              _cur.yref() = 1;
            } else {
              _cur.xref() = 0;
              _cur.yref() = _width - 1;
            }
          }
        } else {
          _cur.yref() ++;
        }
        while ( vert && !done && _cur.y() == _width - 1 ) {
          if ( _cur.x() == 0 ) {
            _cur.xref() = _width - 1;
            _cur.yref() = 1;
          } else {
            done = true;
          }
        }
      }                                    // 2D-Iterator-END
      else {

        if ( !obenunten ) {                // 3D-Iterator-BEGIN...
          _cur.xref() ++;
          if ( _cur.xref() == _width ) {
            _cur.yref() ++;
            _cur.xref() = 0;
            if ( _cur.yref() == _width ) {
              _cur.yref() = 0;
              if ( _cur.zref() == 0 ) {
                _cur.zref() = _width - 1;
              } else {
                obenunten = true;
                vornehinten = false;
                _cur.xref() = -1;
                _cur.yref() = 0;
                _cur.zref() = 1;
              }
            }
          }
        }

        if ( !vornehinten ) {
          _cur.xref() ++;
          if ( _cur.xref() == _width ) {
            _cur.zref() ++;
            _cur.xref() = 0;
            if ( _cur.zref() == _width - 1 ) {
              _cur.zref() = 1;
              if ( _cur.yref() == 0 ) {
                _cur.yref() = _width - 1;
              } else {
                vornehinten = true;
                rechtslinks = false;
                _cur.xref() = 0;
                _cur.yref() = 0;
                _cur.zref() = 1;
              }
            }
          }
        }

        if ( !rechtslinks ) {
          _cur.yref() ++;
          if ( _cur.yref() == _width - 1 ) {
            _cur.zref() ++;
            _cur.yref() = 1;
            if ( _cur.zref() == _width - 1 ) {
              _cur.zref() = 1;
              if ( _cur.xref() == 0 ) {
                _cur.xref() = _width - 1;
              } else {
                rechtslinks = false;
                done = true;
              }
            }
          }
        }                                     // 3D-Iterator-END
      }



      return *this;
    }

    bool operator== ( const EndType & ) const {
      return done;
    }

    bool operator!= ( const EndType & ) const {
      return !done;
    }

  protected:
    IteratedType _cur;
    short _width, _numZ;
    bool done, vert, obenunten, vornehinten, rechtslinks;
  };

  //! Iterator over all elements in a quoc grid
  class OldFullElementIterator {
  public:
    typedef OldFullElementIterator Self;
    typedef Element             IteratedType;
    typedef GridDefinition      BeginType;
    typedef EndElement          EndType;

    OldFullElementIterator( ) : _width ( -1 ) {}

    OldFullElementIterator ( const GridDefinition & Grid ) {
      *this = Grid;
    }

    Self &operator= ( const GridDefinition & Grid ) {
      _width = Grid.getWidth() - 1;

      _numZ = ( Grid.getDimOfWorld() == QC_3D ) ? _width : 1;
      _numY = ( Grid.getDimOfWorld() == QC_1D ) ? 1 : _width;

      _cur.set ( 0, 0, 0, Grid.getGridDepth( ) );
      done = false;
      return *this;
    }

    IteratedType & operator*() {
      return _cur;
    }

    IteratedType * operator->() {
      return &_cur;
    }

    inline OldFullElementIterator& operator++ ( ) {
      _cur.xref() ++;
      if ( _cur.x() == _width ) {
        _cur.xref() = 0;
        _cur.yref() ++;
        if ( _cur.y() == _numY ) {
          _cur.yref() = 0;
          _cur.zref() ++;
          if ( _cur.z() == _numZ ) {
            done = true;
          }
        }
      }
      return *this;
    }

    bool operator== ( const EndType & ) const {
      return done;
    }

    bool operator!= ( const EndType & ) const {
      return !done;
    }

    bool atEnd () const {
      return done;
    }

  protected:
    short _numY;
    short _numZ;
    IteratedType _cur;
    short _width;
    bool done;
  };

  //! traverses all elements adjacent to an element
  template <Dimension Dim>
  class FullNodeNeighborIterator {
  public:
    typedef FullNodeNeighborIterator Self;
    typedef Element                  IteratedType;
    typedef EndElement               EndType;

    FullNodeNeighborIterator( ) : _done ( false ), _width ( -1 ), _curIndex ( 0 ) {}

    //! specify also the level of elements to be traversed. -1 defaults to same maxlevel
    void start ( const GridDefinition &Grid, const CoordType &Node, short ElementLevel = -1 ) {
      _curIndex = -1;
      _width = Grid.getWidth( ) - 1;
      _done = false;

      if ( ElementLevel == -1 ) {
        ElementLevel = Grid.getGridDepth();
      }
      short s = 1 << ( Grid.getGridDepth() - ElementLevel );

      short x = Node.x(), y = Node.y(), z = Node.z();
      _nbEls[ 0 ].set ( x  , y  , z  , ElementLevel, 9 );
      _nbEls[ 1 ].set ( x - s, y  , z  , ElementLevel, 9 );
      _nbEls[ 2 ].set ( x  , y - s, z  , ElementLevel, 9 );
      _nbEls[ 3 ].set ( x - s, y - s, z  , ElementLevel, 9 );
      if ( Dim == QC_3D ) {
        _nbEls[ 4 ].set ( x  , y  , z - s, ElementLevel, 9 );
        _nbEls[ 5 ].set ( x - s, y  , z - s, ElementLevel, 9 );
        _nbEls[ 6 ].set ( x  , y - s, z - s, ElementLevel, 9 );
        _nbEls[ 7 ].set ( x - s, y - s, z - s, ElementLevel, 9 );
      }

      operator++ ( );
    }

    bool operator== ( const EndType & ) const {
      return _done;
    }

    bool operator!= ( const EndType & ) const {
      return !_done;
    }

    IteratedType& operator*() {
      return _cur;
    }

    IteratedType* operator->() {
      return &_cur;
    }

    inline FullNodeNeighborIterator<Dim>& operator++ ( ) {
      if ( _curIndex == ElementInfoTrait<Dim>::numElementsAdjacentToNode - 1 ) {
        _done = true;
        return *this;
      }
      _cur = _nbEls[ ++_curIndex ];
      while ( _curIndex < ElementInfoTrait<Dim>::numElementsAdjacentToNode - 1 && ( _cur.x() < 0 || _cur.y() < 0 || _cur.z() < 0 ||
                                                                                    _cur.x() >= _width || _cur.y() >= _width || _cur.z() >= _width ) ) {
        _cur = _nbEls[ ++_curIndex ];
        if ( _curIndex == ElementInfoTrait<Dim>::numElementsAdjacentToNode - 1 ) {
          _done = true;
          break;
        }
      }
      return *this;
    }

  protected:
    qc::Element _nbEls[ ElementInfoTrait<Dim>::numElementsAdjacentToNode ];
    bool _done;
    int _width;
    int _curIndex;
    IteratedType _cur;
  };

  typedef FullNodeNeighborIterator<qc::QC_3D> FullNodeNeighborIterator3D;
  typedef FullNodeNeighborIterator<qc::QC_2D> FullNodeNeighborIterator2D;



  class NodeNeighborHoodIterator2D {
  public:
    typedef NodeNeighborHoodIterator2D Self;
    typedef CoordType                  IteratedType;
    typedef EndElement                 EndType;

    NodeNeighborHoodIterator2D( ) : _done ( false ), _width ( -1 ) {}

    void begin ( const GridDefinition &Grid, CoordType &Node, int PatchSize = 5 ) {
      _width = Grid.getWidth( ) - 1;
      _done = false;
      short s = ( PatchSize - 1 ) >> 1;

      _x1 = aol::Max ( 0, Node.x() - s );
      _x2 = aol::Min ( _width, Node.x() + s );

      _y1 = aol::Max ( 0, Node.y() - s );
      _y2 = aol::Min ( _width, Node.y() + s );

      _cur.set ( _x1, _y1, 0 );
    }

    bool operator== ( const EndType & ) const {
      return _done;
    }

    bool operator!= ( const EndType & ) const {
      return !_done;
    }

    IteratedType& operator*() {
      return _cur;
    }

    IteratedType* operator->() {
      return &_cur;
    }

    inline Self& operator++ ( ) {

      if ( _cur.x() < _x2 ) {
        _cur.xref() ++;
      } else if ( _cur.y() < _y2 ) {
        _cur.xref() = _x1;
        _cur.yref() ++;
      } else {
        _done = true;
      }

      return *this;
    }


  protected:
    bool _done;
    int _width;
    short _x1, _x2, _y1, _y2;
    IteratedType _cur;
  };


  //! for a grid with periodic boundary condition
  class OldFullPeriodicElementIterator {
  public:
    typedef OldFullPeriodicElementIterator Self;
    typedef Element                    IteratedType;
    typedef GridDefinition             BeginType;
    typedef EndElement                 EndType;

    OldFullPeriodicElementIterator() : _width ( -1 ) {}

    Self & operator= ( const BeginType & Grid ) {
      _width = Grid.getWidth();
      if ( Grid.getDimOfWorld() == QC_2D ) {
        _numZ = 1;
      } else {
        _numZ = _width;
      }
      _cur.set ( 0, 0, 0, Grid.getGridDepth() );
      done = false;
      return *this;
    }

    IteratedType & operator*() {
      return _cur;
    }

    IteratedType * operator->() {
      return &_cur;
    }

    inline Self & operator++ ( ) {
      _cur.xref() ++;
      if ( _cur.xref() == _width ) {
        _cur.xref() = 0;
        _cur.yref() ++;
        if ( _cur.yref() == _width ) {
          _cur.yref() = 0;
          _cur.zref() ++;
          if ( _cur.zref() == _numZ ) {
            done = true;
          }
        }
      }
      return *this;
    }

    bool operator== ( const EndType & ) const {
      return done;
    }

    bool operator!= ( const EndType & ) const {
      return !done;
    }

  protected:
    short _numZ;
    IteratedType _cur;
    short _width;
    bool done;
  };


  //! periodic boundary condition in x-direction
  class OldXPeriodicElementIterator {
  public:
    typedef OldXPeriodicElementIterator Self;
    typedef Element                    IteratedType;
    typedef GridDefinition             BeginType;
    typedef EndElement                 EndType;

    OldXPeriodicElementIterator() : _width ( -1 ) {}

    Self &operator= ( const BeginType & Grid ) {
      _width = Grid.getWidth() - 1;
      if ( Grid.getDimOfWorld() == QC_2D ) {
        _numZ = 1;
      } else {
        _numZ = _width;
      }
      _cur.set ( 0, 0, 0, Grid.getGridDepth() );
      done = false;
      return *this;
    }

    IteratedType & operator*() {
      return _cur;
    }

    IteratedType * operator->() {
      return &_cur;
    }

    inline Self & operator++ ( ) {
      _cur.xref() ++;
      if ( _cur.xref() == _width + 1 ) {
        _cur.xref() = 0;
        _cur.yref() ++;
        if ( _cur.yref() == _width ) {
          _cur.yref() = 0;
          _cur.zref() ++;
          if ( _cur.zref() == _numZ ) {
            done = true;
          }
        }
      }
      return *this;
    }

    bool operator== ( const EndType & ) const {
      return done;
    }

    bool operator!= ( const EndType & ) const {
      return !done;
    }

  protected:
    short _numZ;
    IteratedType _cur;
    short _width;
    bool done;
  };


  //! iterator for L-shape
  class OldLShapeElementIterator {
  public:
    typedef OldLShapeElementIterator Self;
    typedef Element               IteratedType;
    typedef GridDefinition        BeginType;
    typedef EndElement            EndType;

    OldLShapeElementIterator( ) : _width ( -1 ) {}

    Self &operator= ( const BeginType & Grid ) {
      _width = Grid.getWidth() - 1;

      if ( Grid.getDimOfWorld() == QC_2D ) {
        _numZ = 1;
      } else {
        _numZ = _width;
      }

      _cur.set ( 0, 0, 0, Grid.getGridDepth( ) );
      done = false;
      return *this;
    }

    IteratedType & operator*() {
      return _cur;
    }

    IteratedType * operator->() {
      return &_cur;
    }

    inline Self & operator++ ( ) {
      _cur.xref() ++;
      if ( _cur.x() == _width ) {
        _cur.yref() ++;
        if ( _cur.yref() < _width / 2 ) {
          _cur.xref() = 0;
        } else {
          _cur.xref() = _width / 2;
        }
        if ( _cur.y() == _width ) {
          _cur.yref() = 0;
          _cur.zref() ++;
          if ( _cur.z() == _numZ ) {
            done = true;
          }
        }
      }
      return *this;
    }

    bool operator== ( const EndType & ) const {
      return done;
    }

    bool operator!= ( const EndType & ) const {
      return !done;
    }

  protected:
    short _numZ;
    IteratedType _cur;
    short _width;
    bool done;
  };

  //! \brief Iterator for all elements adjacent to a node
  //! \author Toelkes
  //! \warning Only 2D version implemented!
  //! \note First version, inefficient
  class NodeElementIterator {
    typedef GridDefinition GridType;

    static int _origNodeLocalIndex[5];
  public:
    //! \brief Constructor
    //! \param grid The grid containing the node and elements to be iterated
    //! \param level The grid's level
    //! \param coord The coordinate of the node that should be iterated around
    NodeElementIterator ( const GridType &grid, const int level, const CoordType coord )
  : _grid ( grid ), _level ( level ), _currentCoord ( coord ), _currentPos ( NE ) {
      _current.set ( _currentCoord[0], _currentCoord[1], _currentCoord[2], _level );

      if ( _currentCoord[0] == grid.getNumX () - 1 && _currentCoord[1] == 0  ) {
        --_current[0];
        _currentPos = NW;
      }
      else if ( _currentCoord[0] == 0 && _currentCoord[1] == grid.getNumY () - 1 ) {
        --_current[1];
        _currentPos = SE;
      }
      else if ( _currentCoord[1] == grid.getNumY () - 1 ) {
        --_current[0];
        --_current[1];
        _currentPos = SW;
      }
      else if ( _currentCoord[0] == grid.getNumX () - 1 ) {
        --_current[0];
        _currentPos = NW;
      }
    }

    const NodeElementIterator& operator++ () {
      switch ( _currentPos ) {
      case NE:
        if ( _current[0] != 0 ) {
          // Set NW
          --_current[0];
          _currentPos = NW;
        }
        else if ( _current[1] != 0 ) {
          --_current[1];
          _currentPos = SE;
        }
        break;
      case NW:
        if ( _current[1] != 0 ) {
          // Set SW
          --_current[1];
          _currentPos = SW;
        }
        else {
          _currentPos = END;
        }
        break;
      case SW:
        if ( _current[0] < _grid.getNumX () - 2 ) {
          // Set SE
          ++_current[0];
          _currentPos = SE;
        }
        else {
          _currentPos = END;
        }
        break;
      case SE:
        // END
        _currentPos = END;
        break;
      case END:
        throw aol::Exception ( "Iterator should not be incremented more than 4 times!", __FILE__, __LINE__ );
        break;
      default:
        throw aol::Exception ( "Undefined position!", __FILE__, __LINE__ );
      }

      return *this;
    }

    inline bool atEnd () const {
      return ( _currentPos == END );
    }

    inline bool notAtEnd () const {
      return !atEnd ();
    }

    const qc::Element& operator* () const {
      return this->_current;
    }

    const qc::Element* operator-> () const {
      return &( this->_current );
    }

    //! \brief Returns the local index of the node that is being iterated around (the node given in the constructor).
    int getCentralNodeIndex () const {
      return this->_origNodeLocalIndex[this->_currentPos];
    }

  protected:
    enum Position {
      NE = 0,
      NW = 1,
      SW = 2,
      SE = 3,
      END = 4
    };
    const GridType &_grid;
    const int _level;
    Element _current;
    CoordType _currentCoord;
    Position _currentPos;
  };


protected:
  int           gridDepth;  /**< Maximum depth of grid */
  const int     width;      /**< Number of nodes in one direction */

  int *         uniformGridIndexOffset; /**< Index offsets in case of uniform grid */
};



// --------------------------------------------------------------
//! class Ball for traversing all elements in an epsilon ball
//! with a flexible epicenter
// --------------------------------------------------------------

class Ball {

protected:
  double _epsilon;                   // the radius of the ball
  aol::Vec3<double> _center;         // the epicenter of the ball
  Element _centerEl;                 // the element that contains the epicenter

  const GridDefinition &_grid;
  bool _initialized;                 // is the ball already calculated
  double _h;                         // the width in each direction


public:
  // list of offsets of the elements that are in the epsilon-ball
  // offset means (x,y,z)-offset from the epicenter
  vector< aol::Vec3<short> > _offsets;


  //! Method to calculate the offsets of the elements in the ball to the center.
  //! They are saved in the list _offsets!
  //! if someone needs this very fast, it could be optimized by exploiting the
  //! symmetry of the ball ...

  void calcOffsets ( const GridDefinition &grid ) {

    // get the lowest, leftest, frontest (english for runaways...) element
    aol::Vec3<double> shiftetCenter ( _center );
    for ( int i = 0; i < 3; i++ ) shiftetCenter[i] -= _epsilon;
    Element startEl;

    grid.getElementByPoint ( shiftetCenter, startEl );

    int xStart = startEl.x();
    int yStart = startEl.y();
    int zStart = startEl.z();

    int n = 2 * ( static_cast<int> ( _epsilon / _h ) + 1 );   // diameter of the ball in #elements
    bool inside;                                         // element inside the ball?
    aol::Vec3<double> pointCoords, temp;                 // the actual point
    aol::Vec3<short>  listElement;                       // element to be added to the list


//     cerr<<"\n\nIterator-Ausgaben:\nn: "<<n<<"\ncenter-coords: "<<_center<<endl;
//     cerr<<"Center-el: "<<_centerEl<<"\nshiftetCenter-coords: "<<shiftetCenter<<endl;
//     cerr<<"startEl: "<<startEl<<endl<<endl;
//

    // now traverse the cube around the center and check, whether the elements
    // touch the ball or not
    for ( int x = xStart; x < xStart + n; x++ )
      for ( int y = yStart; y < yStart + n; y++ )
        for ( int z = zStart; z < zStart + n; z++ ) {
          // get the world coordinates of node # 0 of the actual element
          pointCoords[0] = x * _h;
          pointCoords[1] = y * _h;
          pointCoords[2] = z * _h;

          inside = false;

          // now check all 8 nodes of the element
          for ( int i = 0; i < 2; i++ ) {
            for ( int j = 0; j < 2; j++ ) {
              for ( int k = 0; k < 2; k++ ) {

                // vector between actual point and center
                temp = _center - pointCoords;
                if ( temp.norm() <= _epsilon ) {    // distance <= eps? => inside the ball

                  if ( !inside ) {                  // first time => take it into the list
                    inside = true;
                    listElement[0] = x - _centerEl.x();
                    listElement[1] = y - _centerEl.y();
                    listElement[2] = z - _centerEl.z();
                    _offsets.push_back ( listElement );
                  }

                }
                pointCoords[0] += _h;
              }
              pointCoords[1] += _h;
              pointCoords[0] = x * _h;
            }
            pointCoords[2] += _h;
            pointCoords[0] = x * _h;
            pointCoords[1] = y * _h;
          }

        }  // end for z
  }  // end method


  // ---------------- constructors -------------------------------------

  Ball ( const GridDefinition &grid, double epsilon ) :
      _epsilon ( epsilon ), _grid ( grid )  {

    if ( epsilon > 0.5 ) throw aol::Exception ( "iteratorBall: only epsilon<=0.5 implemented yet!" );

    _h   = grid.H();

    for ( int i = 0; i < 3; i++ ) _center[i] = 0.5;      // pre-initializiation
    grid.getElementByPoint ( _center, _centerEl );       // element that contains the epicenter

    int n = 2 * ( static_cast<int> ( epsilon / _h ) + 1 );    // diameter of the ball in #elements
    _offsets.reserve ( n*n*n );                          // reserve memory for the whole cube
    // => enough at any rate
    calcOffsets ( grid );                                // calculate the list of offsets
  }

  // ------------------ methods ----------------------------------------

  // for testing:
  void printList ( void ) {
    cerr << "\n\nContent of the offset-list: \n";
    for ( unsigned int i = 0; i < _offsets.size(); i++ )
      cerr << _offsets[i] << "; ";
    cerr << endl << endl;
  }




  // ---------------- the iterator in this class -------------------------

  class iteratorBall {
  public:
    typedef iteratorBall Self;

    Element _centerEl;                    // the element that contains the epicenter of the ball
    Element _cur;
    vector< aol::Vec3<short> >::iterator _offsetsIt;
    vector< aol::Vec3<short> >::iterator _offsetsItEnd;
    int _N;                               // size of one side of the cube


    // ----------------------- constructors --------------------------------------------

    iteratorBall( ) : _centerEl(), _offsetsIt() {
      _cur.setZero();
    }


    iteratorBall ( vector< aol::Vec3<short> >::iterator it,
                   vector< aol::Vec3<short> >::iterator itEnd,  Element centerEl, int N ) :
        _centerEl ( centerEl ), _offsetsIt ( it ), _offsetsItEnd ( itEnd ), _N ( N - 1 ) {

      // search the first element that is inside [0,1]^3

      bool counterFlag = 0;      // flag for the case of it = itEnd (then don't increase _offsetsIt)
      short xk, yk, zk;

      do {
        if ( counterFlag ) _offsetsIt++;
        counterFlag = true;

        if ( _offsetsIt != _offsetsItEnd ) {
          xk = _centerEl.x() + ( *_offsetsIt ).x();
          yk = _centerEl.y() + ( *_offsetsIt ).y();
          zk = _centerEl.z() + ( *_offsetsIt ).z();
        } else {
          xk = yk = zk = -1;
//         cerr<<"("<<xk<<","<<yk<<","<<zk<<"), ";
        }
      } while ( ! ( xk >= 0 && xk < _N && yk >= 0 && yk < _N && zk >= 0 && zk < _N )
                && _offsetsIt != _offsetsItEnd );

      _cur.xref() = xk;
      _cur.yref() = yk;
      _cur.zref() = zk;

    }


    // ------------------------- methods -----------------------------------------------

    Element& getCurrentPosition( ) {
      return _cur;
    }

    iteratorBall &operator= ( const iteratorBall& Other ) {
      _centerEl       = Other._centerEl;
      _offsetsIt      = Other._offsetsIt;
      _offsetsItEnd   = Other._offsetsItEnd;
      _cur = Other._cur;
      _N = Other._N;

      return *this;
    }


    inline iteratorBall& operator++ ( int ) {
      short xk, yk, zk;

      do {
        _offsetsIt++;

        if ( _offsetsIt != _offsetsItEnd ) {
          xk = _centerEl.x() + ( *_offsetsIt ).x();
          yk = _centerEl.y() + ( *_offsetsIt ).y();
          zk = _centerEl.z() + ( *_offsetsIt ).z();
        } else {
          xk = yk = zk = -1;
          //cerr<<"("<<xk<<","<<yk<<","<<zk<<"), ";
        }
      } while ( ! ( xk >= 0 && xk < _N && yk >= 0 && yk < _N && zk >= 0 && zk < _N )
                && _offsetsIt != _offsetsItEnd );

      _cur.xref() = xk;
      _cur.yref() = yk;
      _cur.zref() = zk;

      return *this;
    }


    Element& operator*() {
      return _cur;
    }

    const Element& operator*() const {
      return _cur;
    }

    const Element* operator->() const {
      return &_cur;
    }

    inline bool operator!= ( const Self& Other ) const {
      if ( _centerEl != Other._centerEl )
        throw aol::Exception ( "Ball-Iterator: Vergleich nur sinnvoll bei identischen Zentren!" );
      return ( _offsetsIt != Other._offsetsIt );
    }

    inline bool operator== ( const Self& Other ) const {
      if ( _centerEl != Other._centerEl )
        throw aol::Exception ( "Ball-Iterator: Vergleich nur sinnvoll bei identischen Zentren!" );
      return ( _offsetsIt == Other._offsetsIt );
    }

  };


  iteratorBall begin ( aol::Vec3<double> center ) {
    Element centerEl;
    _grid.getElementByPoint ( center, centerEl );
    return iteratorBall ( _offsets.begin(), _offsets.end(), centerEl, _grid.getWidth() );
  }

  iteratorBall end ( aol::Vec3<double> center ) {
    Element centerEl;
    _grid.getElementByPoint ( center, centerEl );
    return iteratorBall ( _offsets.end(), _offsets.end(), centerEl, _grid.getWidth() );
  }
};


} // end namespace

#endif

