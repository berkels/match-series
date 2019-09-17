#ifndef __ITERATORS_H
#define __ITERATORS_H

#include <quoc.h>
#include <arrayExtensions.h>

namespace qc {

//! Basis class for RectangularIterator and RectangularBoundaryIterator in 2D, do not use directly.
template< qc::Dimension Dim, typename IteratedType = qc::CoordType >
class RectangularIteratorBase {
protected:
  IteratedType _lower, _upper;
  IteratedType _current;

  inline void increment ( ) {   // has to be implemented differently for different dimension
    switch ( Dim ) {
    case QC_1D: {
      ++ ( _current[0] );
      break;
    }
    case QC_2D: {
      ++ ( _current[0] );
      if ( _current[0] == _upper[0] ) {
        _current[0] = _lower[0];
        ++ ( _current[1] );
      }
      break;
    }
    case QC_3D: {
      ++ ( _current[0] );
      if ( _current[0] == _upper[0] ) {
        _current[0] = _lower[0];
        ++ ( _current[1] );
        if ( _current[1] == _upper[1] ) {
          _current[1] = _lower[1];
          ++ ( _current[2] );
        }
      }
      break;
    }
    default:
      throw aol::Exception ( "RectangularIteratorBase::increment(): illegal dimension", __FILE__, __LINE__ );
    }
  }

public:
  //! Default constructor
  RectangularIteratorBase ( ) : _lower (), _upper (), _current ( _lower ) { }

  RectangularIteratorBase ( const IteratedType &Lower, const IteratedType &Upper ) : _lower ( Lower ), _upper ( Upper ), _current ( _lower ) {
    for ( int i = 0; i < Dim; ++i )
      if ( _lower[i] >= _upper[i] ) // then the iterated set is empty, i.e. we may set _current such that atEnd returns true
        _current[Dim-1] = _upper[Dim-1];
  }

  RectangularIteratorBase ( const typename aol::VecDimTrait<int, Dim>::VecType &Lower, const typename aol::VecDimTrait<int, Dim>::VecType &Upper ) : _lower ( ), _upper ( ), _current ( ) {
    for ( short int i = 0; i < Dim; ++i ) {
      _lower[i] = Lower[i];
      _upper[i] = Upper[i];
      if ( Lower[i] != _lower[i] || Upper[i] != _upper[i] )
        throw aol::Exception ( "RectangularIteratorBase: integer to short conversion produced overflow", __FILE__, __LINE__ );
    }
    _current = _lower;
    for ( int i = 0; i < Dim; ++i )
      if ( _lower[i] >= _upper[i] ) // then the iterated set is empty, i.e. we may set _current such that atEnd returns true
        _current[Dim-1] = _upper[Dim-1];
  }

  template < typename Structure >
  explicit RectangularIteratorBase ( const Structure &struc ) : _lower ( ), _upper ( ), _current ( ) {
    restart ( struc );
  }

  template < typename ContainedType, typename ContainedArrayType >
  explicit RectangularIteratorBase ( const RectangularContainer< ContainedType, ContainedArrayType, Dim >  &rCont ) : _lower ( ), _upper (  ), _current ( ) {
    for ( short int i = 0; i < Dim; ++i ) {
      _lower[i] = rCont.getLower()[i];
      _upper[i] = rCont.getUpper()[i];
    }
    _current = _lower;
    for ( int i = 0; i < Dim; ++i )
      if ( _lower[i] >= _upper[i] ) // then the iterated set is empty, i.e. we may set _current such that atEnd returns true
        _current[Dim-1] = _upper[Dim-1];
  }

  inline bool atEnd ( ) const {
    return ( _current[ Dim-1 ] == _upper[ Dim-1 ] );
  }

  inline bool notAtEnd ( ) const {
    return ( !atEnd() );
  }

  const IteratedType& operator* ( ) const {
    return ( _current );
  }

  const IteratedType* operator-> ( ) const {
    return ( &_current );
  }

  template < typename Structure >
  void restart ( const Structure &struc ) {
      for ( short int i = 0; i < Dim; ++i ) {
      _lower[i] = 0;
      _upper[i] = struc.getSize()[i];
    }
    _current = _lower;
    for ( int i = 0; i < Dim; ++i )
      if ( _lower[i] >= _upper[i] ) // then the iterated set is empty, i.e. we may set _current such that atEnd returns true
        _current[Dim-1] = _upper[Dim-1];
  }

};


/** This iterator iterates over all points in a 2D or 3D brick of indices
 *  \author Schwen
 */
template< qc::Dimension Dim, typename IteratedType = qc::CoordType >
class RectangularIterator : public RectangularIteratorBase< Dim, IteratedType  > {
public:
  //! Constructor setting up brick iterator for brick [Lower, upper)
  RectangularIterator ( const typename aol::VecDimTrait<int, Dim>::VecType &Lower, const typename aol::VecDimTrait<int, Dim>::VecType &Upper ) : RectangularIteratorBase<Dim, IteratedType> ( Lower, Upper ) { }

  //! Constructor setting up brick iterator for brick [Lower, upper)
  RectangularIterator ( const IteratedType &Lower, const IteratedType &Upper ) : RectangularIteratorBase<Dim, IteratedType> ( Lower, Upper ) { }

  //! Constructor setting up brick iterator for Structure that supports getNum{X,Y}
  template < typename Structure >
  explicit RectangularIterator ( const Structure &struc ) : RectangularIteratorBase< Dim, IteratedType > ( struc ) { }

  //! prefix increment operator.
  const IteratedType& operator++ ( ) {
    this->increment();
    return ( this->_current );
  }
};


/** Class for template specialization */
template< qc::Dimension Dim>
class LocalLInfBoxIterator;

/** 2D Iterator over an L infinity box centered at a given point and limited to inside bounding brick [ Point - radius, Point + radius ]
 *  \author Mevenkamp
 */
template<>
class LocalLInfBoxIterator<qc::QC_2D> : public qc::RectangularIterator<qc::QC_2D> {
public:
  LocalLInfBoxIterator ( const qc::CoordType &Center, const short radius, const aol::Vec3<int> &min, const aol::Vec3<int> &max )
  : qc::RectangularIterator<qc::QC_2D> ( qc::CoordType ( static_cast<short> ( aol::Max ( Center[0] - radius     , min[0] ) ),
                                                         static_cast<short> ( aol::Max ( Center[1] - radius     , min[1] ) ),
                                                         0 ),
                                         qc::CoordType ( static_cast<short> ( aol::Min ( Center[0] + radius + 1 , max[0] ) ),
                                                         static_cast<short> ( aol::Min ( Center[1] + radius + 1 , max[1] ) ),
                                                         0 ) ) {
  }
};

/** Iterator over an L infinity box centered at a given point and limited to inside bounding brick [ Point - radius, Point + radius ]
 *  \author Schwen
 */
template<>
class LocalLInfBoxIterator<qc::QC_3D> : public qc::RectangularIterator<qc::QC_3D> {
public:
  LocalLInfBoxIterator ( const qc::CoordType &Center, const short radius, const aol::Vec3<int> &min, const aol::Vec3<int> &max )
  : qc::RectangularIterator<qc::QC_3D> ( qc::CoordType ( static_cast<short> ( aol::Max ( Center[0] - radius     , min[0] ) ),
                                                         static_cast<short> ( aol::Max ( Center[1] - radius     , min[1] ) ),
                                                         static_cast<short> ( aol::Max ( Center[2] - radius     , min[2] ) ) ),
                                         qc::CoordType ( static_cast<short> ( aol::Min ( Center[0] + radius + 1 , max[0] ) ),
                                                         static_cast<short> ( aol::Min ( Center[1] + radius + 1 , max[1] ) ),
                                                         static_cast<short> ( aol::Min ( Center[2] + radius + 1 , max[2] ) ) ) ) {
  }
};


/** Class for template specialization */
template< qc::Dimension Dim>
class LocalOnesidedLInfBoxIterator;

// 2D version similar, but not implemented yet

/** Iterator over an L infinity box starting at a given point and limited to inside bounding brick [ Point, Point + boxSize )
 *  \author Schwen
 */
template<>
class LocalOnesidedLInfBoxIterator<qc::QC_3D> : public qc::RectangularIterator<qc::QC_3D> {
public:
  LocalOnesidedLInfBoxIterator ( const qc::CoordType &Point, const short boxSize, const aol::Vec3<int> &max )
  : qc::RectangularIterator<qc::QC_3D> ( Point,
                                         qc::CoordType ( static_cast<short> ( aol::Min ( Point[0] + boxSize , max[0] ) ),
                                                         static_cast<short> ( aol::Min ( Point[1] + boxSize , max[1] ) ),
                                                         static_cast<short> ( aol::Min ( Point[2] + boxSize , max[2] ) ) ) ) {
  }
};


/** This iterator iterates over all boundary points of a 2D or 3D brick of indices.
 *  \warning It wastes time by iterating over all points and checking whether it is a boundary point, but this is fast and easy to implement correctly.
 *  \todo replace by nicer version of loop in qc::GridDefinition::OldFullBoundaryNodeIterator
 *  \author Schwen
 */
template< qc::Dimension Dim, typename IteratedType = qc::CoordType >
class RectangularBoundaryIterator: public RectangularIteratorBase< Dim, IteratedType > {
  static bool speedWarningPrinted;

  void printSpeedWarningIfNecessary ( ) const {
#ifdef _OPENMP
#pragma omp critical (BoundaryIteratorSpeedWarning)
#endif
    {
      if ( !speedWarningPrinted ) {
        cerr << aol::color::red << "\nWarning: RectangularBoundaryIterator is slow.\n\n" << aol::color::reset;
        speedWarningPrinted = true;
      }
    }
  }

public:
  //! Constructor setting up brick iterator for brick [Lower, upper)
  RectangularBoundaryIterator ( const typename aol::VecDimTrait<int, Dim>::VecType &Lower, const typename aol::VecDimTrait<int, Dim>::VecType &Upper ) : RectangularIteratorBase< Dim, IteratedType > ( Lower, Upper ) {
    printSpeedWarningIfNecessary();
  }

  //! Constructor setting up brick iterator for brick [Lower, upper)
  RectangularBoundaryIterator ( const qc::CoordType &Lower, const qc::CoordType &Upper ) : RectangularIteratorBase< Dim, IteratedType > ( Lower, Upper ) {
    printSpeedWarningIfNecessary();
  }

  //! Constructor setting up brick iterator for Structure that supports getNum{X,Y}, e. g. qc::GridStructure, qc::Array, ...
  template < typename Structure >
  explicit RectangularBoundaryIterator ( const Structure &struc ) : RectangularIteratorBase< Dim, IteratedType > ( struc ) {
    printSpeedWarningIfNecessary();
  }

  //! prefix increment operator.
  const qc::CoordType& operator++ ( ) {
    this->increment();
    switch ( Dim ) {
    case QC_2D: {
      while ( ( this->_current[0] > this->_lower[0]     ) &&
              ( this->_current[0] < this->_upper[0] - 1 ) &&
              ( this->_current[1] > this->_lower[1]     ) &&
              ( this->_current[1] < this->_upper[1] - 1 ) ) // inner node
        this->increment();
      break;
    }
    case QC_3D: {
      while ( ( this->_current[0] > this->_lower[0]     ) &&
              ( this->_current[0] < this->_upper[0] - 1 ) &&
              ( this->_current[1] > this->_lower[1]     ) &&
              ( this->_current[1] < this->_upper[1] - 1 ) &&
              ( this->_current[2] > this->_lower[2]     ) &&
              ( this->_current[2] < this->_upper[2] - 1 ) ) // inner node
        this->increment();
      break;
    }
    default:
      throw aol::Exception ( "RectangularBoundaryIterator::operator++: illegal dimension", __FILE__, __LINE__ );
    }
    return ( this->_current );
  }
};

template< qc::Dimension Dim, typename IteratedType >
bool RectangularBoundaryIterator< Dim, IteratedType >::speedWarningPrinted = false;


/** Class for template specialization */
template< qc::Dimension Dim>
class LocalLInfBoxBoundaryIterator;

// 2D version similar, but not implemented yet

/** Iterator over the boundary of an L infinity box centered at a given point and limited to inside bounding brick [ Point - radius, Point + radius ]
 *  \author Schwen
 */
template<>
class LocalLInfBoxBoundaryIterator<qc::QC_3D> : public qc::RectangularBoundaryIterator<qc::QC_3D> {
public:
  LocalLInfBoxBoundaryIterator ( const qc::CoordType &Center, const short radius, const aol::Vec3<int> &min, const aol::Vec3<int> &max )
  : qc::RectangularBoundaryIterator<qc::QC_3D> ( qc::CoordType ( static_cast<short> ( aol::Max ( Center[0] - radius     , min[0] ) ),
                                                                 static_cast<short> ( aol::Max ( Center[1] - radius     , min[1] ) ),
                                                                 static_cast<short> ( aol::Max ( Center[2] - radius     , min[2] ) ) ),
                                                 qc::CoordType ( static_cast<short> ( aol::Min ( Center[0] + radius + 1 , max[0] ) ),
                                                                 static_cast<short> ( aol::Min ( Center[1] + radius + 1 , max[1] ) ),
                                                                 static_cast<short> ( aol::Min ( Center[2] + radius + 1 , max[2] ) ) ) ) {
  }
};


/**
 * This iterator iterates over all boundary elements of a 3D grid .
 * It is possible that this class moves to another file.
 *  \author Teusner
 */
template <typename ConfiguratorType>
class RectangularBoundaryFaceElementIterator :
  public qc::RectangularIteratorBase<ConfiguratorType::Dim, BoundaryFaceElement<typename ConfiguratorType::RealType,ConfiguratorType::Dim> > {

public:
  typedef BoundaryFaceElement<typename ConfiguratorType::RealType,ConfiguratorType::Dim> BFElementType;
  typedef BoundaryFaceElement<typename ConfiguratorType::RealType,ConfiguratorType::Dim> IteratedType;

  explicit RectangularBoundaryFaceElementIterator ( const typename ConfiguratorType::InitType &Initializer ) :
  qc::RectangularIteratorBase< ConfiguratorType::Dim, BFElementType > ( BFElementType(), BFElementType() ) {
      int gridDepth = -1;
      try 
      {
          gridDepth = Initializer.getGridDepth();
      } catch ( aol::Exception &ex ) {
          // just let gridDepth be -1 and suppress the exception message.
          ex.consume();
          // cout << " setting griddepth to -1 " << endl;
      }

      this->_lower.set ( 0, 0, 0, gridDepth, BFElementType::X_LOWER_BOUNDARY );
      this->_current.set ( 0, 0, 0, gridDepth, BFElementType::X_LOWER_BOUNDARY );
      if ( ConfiguratorType::Dim == qc::QC_3D )
          this->_upper.set ( Initializer.getNumX() - 1, Initializer.getNumY() - 1, Initializer.getNumZ() - 1, gridDepth, BFElementType::Z_UPPER_BOUNDARY );
      else
          this->_upper.set ( Initializer.getNumX() - 1, Initializer.getNumY() - 1, 1, gridDepth, BFElementType::Y_UPPER_BOUNDARY );
  }

  //! Iterates over the boundary faces in the order X_LOWER, X_UPPER, Y_LOWER, Y_UPPER, Z_LOWER, Z_UPPER.
  //! On each side of the grid, elements are ordered in lexicographic order.
  const BFElementType& operator++ ( ) {
    switch( this->_current.getBoundaryFaceType() ) {

    case BFElementType::X_LOWER_BOUNDARY:
      ++ ( this->_current[1] );
      if ( this->_current[1] >= this->_upper[1] ) {
        this->_current[1] = this->_lower[1];
        ++ ( this->_current[2] );
      }
      if ( this->_current[2] >= this->_upper[2] )
        this->_current.set ( this->_upper[0] - 1, this->_lower[1], this->_lower[2], this->_current.level(), BFElementType::X_UPPER_BOUNDARY );
      break;

    case BFElementType::X_UPPER_BOUNDARY:
      ++ ( this->_current[1] );
      if ( this->_current[1] >= this->_upper[1] ) {
        this->_current[1] = this->_lower[1];
        ++ ( this->_current[2] );
      }
      if ( this->_current[2] >= this->_upper[2] )
        this->_current.set ( this->_lower[0], this->_lower[1], this->_lower[2], this->_current.level(), BFElementType::Y_LOWER_BOUNDARY );
      break;

    case BFElementType::Y_LOWER_BOUNDARY:
      ++ ( this->_current[0] );
      if ( this->_current[0] >= this->_upper[0] ) {
        this->_current[0] = this->_lower[0];
        ++ ( this->_current[2] );
      }
      if ( this->_current[2] >= this->_upper[2] )
        this->_current.set ( this->_lower[0], this->_upper[1] - 1, this->_lower[2], this->_current.level(), BFElementType::Y_UPPER_BOUNDARY );
      break;

    case BFElementType::Y_UPPER_BOUNDARY:
      ++ ( this->_current[0] );
      if ( this->_current[0] >= this->_upper[0] ) {
        this->_current[0] = this->_lower[0];
        ++ ( this->_current[2] );
      }
      if ( this->_current[2] >= this->_upper[2] ) {
        if ( ConfiguratorType::Dim == qc::QC_3D ) // in 2D, we have arrived at the last face here
          this->_current.set ( this->_lower[0], this->_lower[1], this->_lower[2], this->_current.level(), BFElementType::Z_LOWER_BOUNDARY );
        else
          this->_current = this->_upper;
      }
      break;

    case BFElementType::Z_LOWER_BOUNDARY:
      ++ ( this->_current[0] );
      if ( this->_current[0] >= this->_upper[0] ) {
        this->_current[0] = this->_lower[0];
        ++ ( this->_current[1] );
      }
      if ( this->_current[1] >= this->_upper[1] )
        this->_current.set ( this->_lower[0], this->_lower[1], this->_upper[2] - 1, this->_current.level(), BFElementType::Z_UPPER_BOUNDARY );
      break;

    case BFElementType::Z_UPPER_BOUNDARY:
      ++ ( this->_current[0] );
      if ( this->_current[0] >= this->_upper[0] ) {
        this->_current[0] = this->_lower[0];
        ++ ( this->_current[1] );
      }
      if ( this->_current[1] >= this->_upper[1] )
        this->_current[2] = this->_upper[2]; // in 3D, we have arrived at the last face here
      break;

    default:
      this->_current = this->_upper;
    };

    return ( this->_current );
  }

private:
  using qc::RectangularIteratorBase<ConfiguratorType::Dim, BFElementType >::increment;
};

  

template< qc::Dimension Dim>
class LocalSpiralIterator;
  
/** Iterator over a counter-clockwise spiral (in square snake form), starting at a given point
*  \author Mevenkamp
*
*  Original implementation suggested at stack overflow:
*  http://stackoverflow.com/questions/3706219/algorithm-for-iterating-over-an-outward-spiral-on-a-discrete-2d-grid-from-the-or
*/
template<>
class LocalSpiralIterator<qc::QC_2D> {
protected:
  int _leg, _layer;
  aol::Vec2<short> _cur;
public:
  LocalSpiralIterator ( const aol::Vec2<short> &X ) : _leg ( 0 ), _layer ( 1 ), _cur ( X ) { }
  
  aol::Vec2<short>& operator* ( ) {
    return _cur;
  }
  
  LocalSpiralIterator& operator++ ( ) {
    switch ( _leg ) {
      case 0: ++_cur[0]; if ( _cur[0] == _layer ) ++_leg; break;
      case 1: ++_cur[1]; if ( _cur[1] == _layer ) ++_leg; break;
      case 2: --_cur[0]; if ( -_cur[0] == _layer ) ++_leg; break;
      case 3: --_cur[1]; if ( -_cur[1] == _layer ) { _leg = 0; ++_layer; } break;
      default: throw aol::Exception ( "Illegal leg", __FILE__, __LINE__ ); break;
    }
    
    return *this;
  }
};
  

} // end namespace

#endif
