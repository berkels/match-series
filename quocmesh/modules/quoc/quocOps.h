#ifndef __QUOCOPS_H
#define __QUOCOPS_H

#include <rectangularGrid.h>
#include <configurators.h>

namespace qc {

/**
 * Mass Op for multi-linear Finite Elements on a regular quadrilateral grid.
 * Needs almost no memory and should be noticeably faster to apply than even a
 * qc::FastUniformGridMatrix since unnecessary memory access is avoided.
 *
 * \author Berkels
 */
template <class RealType, Dimension Dim>
class UniformGridMassOp {};

template <class RealType>
class UniformGridMassOp<RealType, qc::QC_1D> : public aol::Op<aol::Vector<RealType> > {
  const int numX;
  const RealType _scale;
  const RealType _weightFull;
  const RealType _weightHalf;
  const RealType _weightQuarter;
public:
  template <typename InitType>
  UniformGridMassOp ( const InitType &InitObj, aol::OperatorType /*OpType*/ = aol::ASSEMBLED)
  : numX ( InitObj.getNumX() ),
  _scale ( (2*qc::RectangularGrid<qc::QC_1D> ( qc::GridSize<qc::QC_1D> ( InitObj ) ).H() )/3 ),
  _weightFull ( _scale * 1 ),
  _weightHalf ( _scale * 0.5 ),
  _weightQuarter ( _scale * 0.25 ) {
  }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    Dest [0] += ( _weightHalf * Arg[0] + _weightQuarter * Arg [1] );
    for ( int i = 1; i < numX-1; ++i )
      Dest[i] += ( _weightQuarter * Arg[i-1] + _weightFull* Arg[i] + _weightQuarter * Arg [i+1] );
    Dest [numX-1] += ( _weightQuarter * Arg[numX-2] + _weightHalf * Arg [numX-1] );
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    Dest [0] = ( _weightHalf * Arg[0] + _weightQuarter * Arg [1] );
    for ( int i = 1; i < numX-1; ++i )
      Dest[i] = ( _weightQuarter * Arg[i-1] + _weightFull* Arg[i] + _weightQuarter * Arg [i+1] );
    Dest [numX-1] = ( _weightQuarter * Arg[numX-2] + _weightHalf * Arg [numX-1] );
  }
};

template <class RealType>
class UniformGridMassOp<RealType, qc::QC_2D> : public aol::Op<aol::Vector<RealType> > {
  const int _numX;
  const int _numY;
  const RealType _scale;
  const RealType _weightFull;
  const RealType _weightHalf;
  const RealType _weightQuarter;
  const RealType _weightEighth;
  const RealType _weightSixteenth;
public:
  template <typename InitType>
  UniformGridMassOp ( const InitType &InitObj, aol::OperatorType /*OpType*/ = aol::ASSEMBLED )
  : _numX ( InitObj.getNumX() ),
  _numY ( InitObj.getNumY() ),
  _scale ( aol::Sqr( 2*qc::RectangularGrid<qc::QC_2D> ( qc::GridSize<qc::QC_2D> ( InitObj ) ).H() / 3 ) ),
  _weightFull ( _scale * 1 ),
  _weightHalf ( _scale * 0.5 ),
  _weightQuarter ( _scale * 0.25 ),
  _weightEighth ( _scale * 0.125 ),
  _weightSixteenth ( _scale * 0.0625 ) {}

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    // Bottom left node
    Dest [0] += ( _weightQuarter * Arg[0]
                 + _weightEighth * ( Arg[1] + Arg[_numX] )
                 + _weightSixteenth * Arg [_numX+1] );
    // Bottom right node
    Dest [_numX-1] += ( _weightQuarter * Arg[_numX-1]
                       + _weightEighth * ( Arg[_numX-2] + Arg[2*_numX-1] )
                       + _weightSixteenth * Arg [2*_numX-2] );
    // Top left node
    {
      const int index = qc::ILexCombine2 ( 0, _numY-1, _numX );
      Dest [index] += ( _weightQuarter * Arg[index]
                       + _weightEighth * ( Arg[index+1] + Arg[index-_numX] )
                       + _weightSixteenth * Arg [index-_numX+1] );
    }
    // Top right node
    {
      const int index = qc::ILexCombine2 ( _numX-1, _numY-1, _numX );
      Dest [index] += ( _weightQuarter * Arg[index]
                       + _weightEighth * ( Arg[index-1] + Arg[index-_numX] )
                       + _weightSixteenth * Arg [index-_numX-1] );
    }
    for ( int i = 1; i < _numX-1; ++i ) {
      // Bottom row
      {
        const int index = qc::ILexCombine2 ( i, 0, _numX );
        Dest[index] += ( _weightHalf * Arg [index]
                        + _weightEighth * ( Arg [index-1]+Arg[index+1] )
                        + _weightSixteenth * ( Arg[index + _numX - 1] + Arg[index + _numX + 1] )
                        + _weightQuarter * Arg[index + _numX ] );
      }
      // Top row
      {
        const int index = qc::ILexCombine2 ( i, _numY-1, _numX );
        Dest[index] += ( _weightHalf * Arg [index]
                        + _weightEighth * ( Arg [index-1]+Arg[index+1] )
                        + _weightSixteenth * ( Arg[index - _numX - 1] + Arg[index - _numX + 1] )
                        + _weightQuarter * Arg[index - _numX ] );
      }
    }

    for ( int j = 1; j < _numY-1; ++j ) {
      // Left column
      {
        const int index = qc::ILexCombine2 ( 0, j, _numX );
        Dest[index] += ( _weightHalf * Arg [index]
                        + _weightEighth * ( Arg [index-_numX]+Arg[index+_numX] )
                        + _weightSixteenth * ( Arg[index + _numX + 1] + Arg[index - _numX + 1] )
                        + _weightQuarter * Arg[index + 1 ] );
      }
      // Right column
      {
        const int index = qc::ILexCombine2 ( _numX-1, j, _numX );
        Dest[index] += ( _weightHalf * Arg [index]
                        + _weightEighth * ( Arg [index-_numX]+Arg[index+_numX] )
                        + _weightSixteenth * ( Arg[index + _numX - 1] + Arg[index - _numX - 1] )
                        + _weightQuarter * Arg[index - 1 ] );
      }
    }

    // Inner nodes
    for ( int j = 1; j < _numY-1; ++j ) {
      for ( int i = 1; i < _numX-1; ++i ) {
        const int index = qc::ILexCombine2 ( i, j, _numX );
        Dest[index] += ( _weightQuarter * ( Arg[index - _numX] + Arg[index-1] + Arg [index+1] + Arg[index + _numX] )
                        + _weightSixteenth * ( Arg[index - _numX - 1] + Arg[index - _numX + 1] + Arg [index+_numX-1] + Arg[index + _numX + 1] )
                        + _weightFull * Arg[index] );
      }
    }
  }
};

template <class RealType>
class UniformGridMassOp<RealType, qc::QC_3D> : public aol::Op<aol::Vector<RealType> > {
  typedef qc::RectangularGridConfigurator<RealType, qc::QC_3D, aol::GaussQuadrature<RealType,qc::QC_3D,3> > ConfiguratorType;
  const typename ConfiguratorType::InitType _grid;
  const aol::MassOp<ConfiguratorType> _massOp;
public:
  template <typename InitType>
  UniformGridMassOp ( const InitType &InitObj, aol::OperatorType /*OpType*/ = aol::ASSEMBLED )
    : _grid ( qc::GridSize<qc::QC_3D>::createFrom ( InitObj ) ),
      _massOp ( _grid, aol::ASSEMBLED ) {
    static bool performanceWarningPrinted = false;
    if ( ! performanceWarningPrinted ) {
      cerr << aol::color::red << "qc::UniformGridMassOp is not implemented in 3D yet. Using aol::MassOp as fallback." << aol::color::reset << endl;
      performanceWarningPrinted = true;
    }
  }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    _massOp.applyAdd ( Arg, Dest );
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    _massOp.apply ( Arg, Dest );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class PixelWiseToRowWiseOp : public aol::Op<aol::Vector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Op<aol::Vector<RealType>, aol::Scalar<RealType> > &_pixelWiseOp;
  const aol::Op<aol::Vector<RealType> > &_pixelWiseOpDerivative;
public:
  PixelWiseToRowWiseOp ( const typename ConfiguratorType::InitType &Grid,
                         const aol::Op<aol::Vector<RealType>, aol::Scalar<RealType> > &PixelWiseOp,
                         const aol::Op<aol::Vector<RealType> > &PixelWiseOpDerivative )
    : _grid ( Grid ),
      _pixelWiseOp ( PixelWiseOp ),
      _pixelWiseOpDerivative ( PixelWiseOpDerivative ) {}
  
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    qc::ScalarArray<RealType, ConfiguratorType::Dim> argArray ( _grid );
    argArray.setInXDirectionFrom ( Arg );
    _pixelWiseOp.applyAdd ( argArray, Dest );
  }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    qc::ScalarArray<RealType, ConfiguratorType::Dim> argArray ( _grid );
    argArray.setInXDirectionFrom ( Arg );
    _pixelWiseOp.apply ( argArray, Dest );
  }
  
  void applyDerivative ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    qc::ScalarArray<RealType, ConfiguratorType::Dim> argArray ( _grid );
    qc::ScalarArray<RealType, ConfiguratorType::Dim> destArray ( _grid );
    argArray.setInXDirectionFrom ( Arg );
    _pixelWiseOpDerivative.apply( argArray, destArray );
    destArray.sumInXDirectionTo ( Dest );
  }
};
  
} // end of namespace qc.

#endif // __QUOCOPS_H
