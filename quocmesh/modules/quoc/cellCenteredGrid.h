#ifndef __CELLCENTEREDGRID_H
#define __CELLCENTEREDGRID_H

#include <rectangularGrid.h>
#include <configurators.h>
#include <linearSmoothOp.h>

namespace qc {

/**
 * \author Berkels
 */
template <qc::Dimension Dim>
class CellCenteredCubicGrid : public qc::RectangularGrid<Dim> {
  const int _depth;

  void properSizeOrDie ( ) {
    // Ensure that the grid is "quadratic".
    qc::GridSize<Dim> ( this->getSize() ).quadraticOrDie();

    // Ensure that the width is 2^d.
    if ( this->getNumX() != ( 1 << _depth ) )
      throw aol::Exception ( aol::strprintf ( "qc::CellCenteredCubicGrid: The width of your domain (%d) is not \"2^d\"!", this->getNumX() ).c_str(), __FILE__, __LINE__ );
  }
public:
  // The second argument is necessary to be compatible with the qc::GridDefinition constructor syntax.
  explicit CellCenteredCubicGrid ( const int Depth, const qc::Dimension InputDim = Dim )
    : qc::RectangularGrid<Dim> ( aol::Vec3<int>( ( 1 << Depth ), (Dim == qc::QC_1D) ? 1 : ( 1 << Depth ), (Dim == qc::QC_3D) ? ( 1 << Depth ) : 1 ) ),
      _depth ( Depth ) {
    if ( Dim != InputDim ) {
      throw aol::Exception ( "qc::CellCenteredCubicGrid: Incompatible dimensions given", __FILE__, __LINE__ );
    }
  }

  explicit CellCenteredCubicGrid ( const aol::Vec3<int> &Size )
    : qc::RectangularGrid<Dim> ( Size ),
      _depth ( qc::logBaseTwo ( Size[0] ) ) {
    properSizeOrDie();
  }

  explicit CellCenteredCubicGrid ( const qc::GridSize<Dim> &Size )
    : qc::RectangularGrid<Dim> ( Size ),
      _depth ( qc::logBaseTwo ( Size[0] ) ) {
    properSizeOrDie();
  }

  //! numX, numY and numZ are the same
  inline int getNumXYZ ( ) const {
    return ( this->getNumX() );
  }

  int getGridDepth() const {
    return ( _depth );
  }
};

/**
 * \author Berkels
 */
template <qc::Dimension Dim>
class GeneralCellCenteredCubicGrid : public qc::RectangularGrid<Dim> {
  int _depth;
  aol::Vec<Dim, int> _baseSize;

  void init ( const aol::Vec3<int> &Size ) {

    aol::Vec<Dim, int> depths;
    for ( int i = 0; i < Dim; ++i ) {
      depths[i] = aol::getFactorExponent ( Size[i], 2 );
#ifdef VERBOSE
      cerr << Size[i] << " = 2^" << depths[i] << "*" << Size[i] / aol::Pow(2,depths[i] ) << endl;
#endif
    }

    _depth = depths.getMinValue();
    for ( int i = 0; i < Dim; ++i ) {
      _baseSize[i] = Size[i] / aol::Pow ( 2, _depth );
#ifdef VERBOSE
      cerr << Size[i] << " = 2^" << _depth << "*" << _baseSize[i] << endl;
#endif
    }

    if ( _depth == 0 )
      throw aol::Exception ( "At least one of the input dimensions is not divisible by 2", __FILE__, __LINE__ );
  }

public:
  // The second argument is necessary to be compatible with the qc::GridDefinition constructor syntax.
  explicit GeneralCellCenteredCubicGrid ( const int Depth, const qc::Dimension InputDim = Dim )
    : qc::RectangularGrid<Dim> ( aol::Vec3<int>( ( 1 << Depth ), (Dim == qc::QC_1D) ? 1 : ( 1 << Depth ), (Dim == qc::QC_3D) ? ( 1 << Depth ) : 1 ) ),
      _depth ( Depth ) {
    _baseSize.setAll ( 1 );
    if ( Dim != InputDim ) {
      throw aol::Exception ( "qc::CellCenteredCubicGrid: Incompatible dimensions given", __FILE__, __LINE__ );
    }
  }

  explicit GeneralCellCenteredCubicGrid ( const aol::Vec3<int> &Size )
    : qc::RectangularGrid<Dim> ( Size ) {
    init ( Size );
  }

  int getGridDepth() const {
    return ( _depth );
  }

  const aol::Vec<Dim, int> &getBaseSize ( ) const {
    return _baseSize;
  }
};

/**
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
class CellCenteredProlongRestrictOpBase : public aol::BiOp<aol::Vector<RealType> > {
protected:
  const RectangularGrid<Dim> &_coarseGrid;
  const RectangularGrid<Dim> &_fineGrid;
  const int _factor;
  const qc::RectangularGrid<Dim> _localGrid;
public:
  CellCenteredProlongRestrictOpBase ( const RectangularGrid<Dim> &Coarse,
                                      const RectangularGrid<Dim> &Fine )
    : _coarseGrid ( Coarse ),
      _fineGrid ( Fine ),
      _factor ( _fineGrid.getNumX() / _coarseGrid.getNumX() ),
      _localGrid ( qc::GridSize<Dim> ( typename aol::VecDimTrait<int, Dim>::VecType ( this->_factor ) ) ) {
    const typename aol::VecDimTrait<int, Dim>::VecType fineSize = qc::GridSize<Dim> ( _fineGrid ).getSizeAsVecDim();
    const typename aol::VecDimTrait<int, Dim>::VecType coarseSize = qc::GridSize<Dim> ( _coarseGrid ).getSizeAsVecDim();

    if ( _factor < 2 )
      throw aol::Exception ( "CellCenteredProlongRestrictOpBase: Only makes sense for scaling factors bigger than 1", __FILE__, __LINE__ );

    if ( fineSize != _factor * coarseSize ) {
      cerr << "Incompatible grid sizes! Fine size = " << fineSize << ", coarse size = " << coarseSize << endl;
      throw aol::Exception ( "CellCenteredProlongRestrictOpBase: Incompatible grid levels", __FILE__, __LINE__ );
    }
  }

  int getFactor ( ) const {
    return _factor;
  }

  const qc::RectangularGrid<Dim> &getLocalGrid ( ) const {
    return _localGrid;
  }
private:
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException ( "qc::CellCenteredProlongRestrictOpBase::applyAdd not implemented.", __FILE__, __LINE__ );
  }
};

/**
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
class CellCenteredRestrictOp : public CellCenteredProlongRestrictOpBase<RealType, Dim> {

public:
  CellCenteredRestrictOp ( const RectangularGrid<Dim> &Coarse,
                           const RectangularGrid<Dim> &Fine )
   : CellCenteredProlongRestrictOpBase<RealType, Dim> ( Coarse, Fine ) { }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    const RealType scaling = aol::ZOTrait<RealType>::one / aol::Pow ( this->_factor, Dim );

    // iterate over the coarse grid nodes
    for ( qc::RectangularIterator<Dim> cit ( this->_coarseGrid ); cit.notAtEnd(); ++cit ) {
      RealType value = 0;
      // iterate over the fine grid nodes corresponding to the current coarse grid node.
      for ( qc::RectangularIterator<Dim> fit ( this->_localGrid ); fit.notAtEnd(); ++fit ) {
        qc::CoordType fineCoord ( *fit );
        for ( int i = 0; i < Dim; ++i )
          fineCoord[i] += this->_factor * (*cit)[i];
        value += Arg [ this->_fineGrid.getIndexMapperRef().getGlobalIndex( fineCoord ) ];
      }
      Dest [ this->_coarseGrid.getIndexMapperRef().getGlobalIndex( *cit ) ] = scaling * value;
    }
  }

  using aol::BiOp< aol::Vector<RealType> >::apply;
};

/**
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
class CellCenteredProlongOp : public CellCenteredProlongRestrictOpBase<RealType, Dim> {
public:
  CellCenteredProlongOp ( const RectangularGrid<Dim> &Coarse,
                          const RectangularGrid<Dim> &Fine )
   : CellCenteredProlongRestrictOpBase<RealType, Dim> ( Coarse, Fine ) { }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    // iterate over the coarse grid nodes
    for ( qc::RectangularIterator<Dim> cit ( this->_coarseGrid ); cit.notAtEnd(); ++cit ) {
      const RealType value = Arg [ this->_coarseGrid.getIndexMapperRef().getGlobalIndex( *cit ) ];
      // iterate over the fine grid nodes corresponding to the current coarse grid node.
      for ( qc::RectangularIterator<Dim> fit ( this->_localGrid ); fit.notAtEnd(); ++fit ) {
        qc::CoordType fineCoord ( *fit );
        for ( int i = 0; i < Dim; ++i )
          fineCoord[i] += this->_factor * (*cit)[i];
        Dest [ this->_fineGrid.getIndexMapperRef().getGlobalIndex( fineCoord ) ] = value;
      }
    }
  }

  using aol::BiOp< aol::Vector<RealType> >::apply;
};

/**
 * \author Berkels
 */
template <typename _RealType, qc::Dimension Dim>
class CellCenteredGridTrait {
public:
  typedef _RealType RealType;
  typedef qc::CellCenteredCubicGrid<Dim> GridType;
  typedef qc::CellCenteredProlongOp<RealType, Dim> ProlongOpType;
  typedef qc::CellCenteredRestrictOp<RealType, Dim> RestrictOpType;
  typedef qc::MultilevelArray<RealType, qc::ScalarArray<RealType, Dim>, ProlongOpType, RestrictOpType, GridType> MultilevelArrayType;
};

/**
 * \author Berkels
 */
template <typename RealType, typename GridType>
class MultilevelArrayTrait {};

template <typename RealType>
class MultilevelArrayTrait<RealType, qc::CellCenteredCubicGrid<qc::QC_1D> > {
public:
  typedef qc::CellCenteredProlongOp<RealType, qc::QC_1D> ProlongOpType;
  typedef qc::CellCenteredRestrictOp<RealType, qc::QC_1D> RestrictOpType;
  typedef qc::MultilevelArray<RealType, qc::ScalarArray<RealType, qc::QC_1D>, ProlongOpType, RestrictOpType, qc::CellCenteredCubicGrid<qc::QC_1D> > MultilevelArrayType;
  typedef qc::MultiDimMultilevelArray<RealType, qc::ScalarArray<RealType, qc::QC_1D>, ProlongOpType, RestrictOpType, qc::CellCenteredCubicGrid<qc::QC_1D> > MultiDimMultilevelArrayType;
  typedef qc::CellCenteredGridTrait<RealType, qc::QC_1D> GridTraitType;

  //! \todo Implement an optimized 1D SmoothOp
  typedef qc::RectangularGridConfigurator<RealType, qc::QC_1D, aol::GaussQuadrature<RealType,qc::QC_1D,3>, qc::CellCenteredCubicGrid<qc::QC_1D> > LinSmoothConfType;
  typedef qc::GeneralLinearSmoothOp<LinSmoothConfType> LinSmoothType;
};
  
template <typename RealType>
class MultilevelArrayTrait<RealType, qc::CellCenteredCubicGrid<qc::QC_2D> > {
public:
  typedef qc::CellCenteredProlongOp<RealType, qc::QC_2D> ProlongOpType;
  typedef qc::CellCenteredRestrictOp<RealType, qc::QC_2D> RestrictOpType;
  typedef qc::MultilevelArray<RealType, qc::ScalarArray<RealType, qc::QC_2D>, ProlongOpType, RestrictOpType, qc::CellCenteredCubicGrid<qc::QC_2D> > MultilevelArrayType;
  typedef qc::MultiDimMultilevelArray<RealType, qc::ScalarArray<RealType, qc::QC_2D>, ProlongOpType, RestrictOpType, qc::CellCenteredCubicGrid<qc::QC_2D> > MultiDimMultilevelArrayType;
  typedef qc::CellCenteredGridTrait<RealType, qc::QC_2D> GridTraitType;
  typedef qc::LinearSmoothOp<RealType, GridTraitType> LinSmoothType;
};

template <typename RealType>
class MultilevelArrayTrait<RealType, qc::CellCenteredCubicGrid<qc::QC_3D> > {
public:
  typedef qc::CellCenteredProlongOp<RealType, qc::QC_3D> ProlongOpType;
  typedef qc::CellCenteredRestrictOp<RealType, qc::QC_3D> RestrictOpType;
  typedef qc::MultilevelArray<RealType, qc::ScalarArray<RealType, qc::QC_3D>, ProlongOpType, RestrictOpType, qc::CellCenteredCubicGrid<qc::QC_3D> > MultilevelArrayType;
  typedef qc::MultiDimMultilevelArray<RealType, qc::ScalarArray<RealType, qc::QC_3D>, ProlongOpType, RestrictOpType, qc::CellCenteredCubicGrid<qc::QC_3D> > MultiDimMultilevelArrayType;
  typedef qc::CellCenteredGridTrait<RealType, qc::QC_3D> GridTraitType;
  typedef qc::LinearSmoothOp<RealType, GridTraitType> LinSmoothType;
};

template <typename RealType>
class MultilevelArrayTrait<RealType, qc::GridDefinition> {
public:
  typedef qc::ProlongOp<RealType> ProlongOpType;
  typedef qc::RestrictOp<RealType, qc::STD_QUOC_RESTRICT> RestrictOpType;
  typedef qc::MultilevelArray<RealType> MultilevelArrayType;
  typedef qc::MultiDimMultilevelArray<RealType> MultiDimMultilevelArrayType;
  typedef qc::DyadicGridTrait<RealType> GridTraitType;
  typedef qc::LinearSmoothOp<RealType, GridTraitType> LinSmoothType;
};

// Dummy trait specialization that allows to use LinSmoothType for qc::RectangularGrid
// in case the size of the grid is compatible with qc::CellCenteredCubicGrid.
template <typename RealType>
class MultilevelArrayTrait<RealType, qc::RectangularGrid<qc::QC_2D> > {
public:
  typedef qc::LinearSmoothOp<RealType, qc::CellCenteredGridTrait<RealType, qc::QC_2D> > LinSmoothType;
};

} // namespace qc

#endif // __CELLCENTEREDGRID_H
