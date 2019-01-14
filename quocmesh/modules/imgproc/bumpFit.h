#ifndef __BUMPFIT_H
#define __BUMPFIT_H


#include <indexMapper.h>
#include <iterators.h>
#include <regression.h>
#include <linearSmoothOp.h>

namespace im {


/**
 * \brief \f$ h \exp\left( \frac{- 1}{ 2 ( 1 - r^2 ) } \left( \left( \frac{x_1-p_1}{\sigma_x} \right)^2 + \left( \frac{x_2-p_2}{\sigma_y} \right)^2 - \frac{2r}{ \sigma_x \sigma_y} (x_1-p_1 ) (x_2-p_2) \right) \right) + o\f$
 * \author Mevenkamp, Berkels
 */
template <typename RealType, bool OffsetIsDOF = true>
class AsymmetricGaussianBumpFunction {
public:
  static const qc::Dimension Dim = qc::QC_2D;
  static const int NumberOfParameters = 6 + OffsetIsDOF;
  static const int NumberOfShareableParameters = 3;

private:
  const aol::Vec2<RealType> _center;
  const RealType _height;
  const RealType _sigmaX, _sigmaY;
  const RealType _rotation;
  const RealType _offset;
  const RealType _c1, _c2;

public:
  AsymmetricGaussianBumpFunction (  const aol::Vec2<RealType> &Center,
                                    const RealType Height,
                                    const RealType SigmaX, const RealType SigmaY,
                                    const RealType Rotation,
                                    const RealType Offset )
   : _center ( Center ),
     _height ( Height ),
     _sigmaX ( SigmaX ), _sigmaY ( SigmaY ),
     _rotation ( Rotation ),
     _offset ( Offset ),
     _c1 ( - 1 / ( 2 * ( 1 - aol::Sqr<RealType> ( _rotation ) ) ) ),
     _c2 ( 2 * _rotation / ( _sigmaX * _sigmaY ) ) {}

  AsymmetricGaussianBumpFunction ( const aol::Vec<Dim, RealType> &Center,
                                   const RealType Height,
                                   const RealType Sigma )
   : _center ( Center ),
     _height ( Height ),
     _sigmaX ( Sigma ), _sigmaY ( Sigma ),
     _rotation ( 0 ),
     _offset ( 0 ),
     _c1 ( - 1 / ( 2 * ( 1 - aol::Sqr<RealType> ( _rotation ) ) ) ),
     _c2 ( 2 * _rotation / ( _sigmaX * _sigmaY ) ) {}

  AsymmetricGaussianBumpFunction ( const aol::Vector<RealType> &Parameters )
   : _center ( Parameters.getData() ),
     _height ( Parameters[2] ),
     _sigmaX ( Parameters[3] ), _sigmaY ( Parameters[4] ),
     _rotation ( Parameters[5] ),
     _offset ( OffsetIsDOF ? Parameters[6] : 0 ),
     _c1 ( - 1 / ( 2 * ( 1 - aol::Sqr<RealType> ( _rotation ) ) ) ),
     _c2 ( 2 * _rotation / ( _sigmaX * _sigmaY ) ) { }

  static void initializeParameterVec ( const aol::Vec<Dim, RealType> &Center,
                                       const RealType Height,
                                       const RealType Sigma,
                                       aol::Vector<RealType> &Parameters ) {
    Parameters[0] = Center[0];
    Parameters[1] = Center[1];
    Parameters[2] = Height;
    Parameters[3] = Sigma;
    Parameters[4] = Sigma;
    Parameters[5] = 0;
    if ( OffsetIsDOF )
      Parameters[6] = 0;
  }

  template <typename VectorType>
  static void initializeSharedParameters ( const RealType Sigma, VectorType &SharedParams ) {
    SharedParams[0] = Sigma;
    SharedParams[1] = Sigma;
    SharedParams[2] = 0;
  }

  RealType evaluate ( const aol::Vec<2, RealType> &X ) const {
    aol::Vec<2, RealType> pos ( X );
    pos -= _center;

    return _height * exp ( _c1 * ( aol::Sqr<RealType> ( pos[0] / _sigmaX ) + aol::Sqr<RealType> ( pos[1] / _sigmaY ) - _c2 * pos[0] * pos[1] ) ) + _offset;
  }

  void evaluateSpatialGradient( const aol::Vec<Dim, RealType> &X, aol::Vec<Dim, RealType> &Grad ) const {
    aol::Vec<Dim, RealType> pos ( X );
    pos -= _center;
    //throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__);
    //Grad.setZero();
    const RealType f = evaluate ( X ) - _offset;

    Grad[0] = - f * _c1 * ( - 2 * pos[0] / aol::Sqr<RealType> ( _sigmaX ) + _c2 * pos[1] ); // Grad[7]
    Grad[1] = - f * _c1 * ( - 2 * pos[1] / aol::Sqr<RealType> ( _sigmaY ) + _c2 * pos[0] ); // Grad[8]
  }
  
  void evaluateSpatialHessian( const aol::Vec<Dim, RealType> &X, aol::Mat<Dim, Dim, RealType> &Hess ) const {
    aol::Vec<Dim, RealType> pos ( X );
    pos -= _center;
    //throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__);
    //Hess.setZero();
    const RealType f = evaluate ( X ) - _offset;
    Hess[0][0] = 2 * f * _c1 / _sigmaX / _sigmaX + _c1 * _c1 * f * aol::Sqr<RealType> ( 2 * pos[0] / _sigmaX / _sigmaX - _c2 * pos[1] );  // Grad[9]
    Hess[0][1] = - f * _c1 * _c2 + f * _c1 * _c1 * ( 2 * pos[0] / _sigmaX / _sigmaX - _c2 * pos[1] ) * (2 * pos[1] / _sigmaY / _sigmaY - _c2 * pos[0] );  // Grad[10]
    Hess[1][0] = Hess[0][1];
    Hess[1][1] = 2 * f * _c1 / _sigmaY / _sigmaY + f * _c1 * _c1 * aol::Sqr<RealType> ( 2 * pos[1] / _sigmaY / _sigmaY - _c2 * pos[0] );  // Grad[11]
  }
  
  template <typename VecType>
  void evaluateParameterGradient ( const aol::Vec<2, RealType> &X, VecType &Grad ) const {
    aol::Vec<2, RealType> pos ( X );
    pos -= _center;
    
    const RealType f = evaluate ( X ) - _offset;
    Grad[0] = f * _c1 * ( - 2 * pos[0] / aol::Sqr<RealType> ( _sigmaX ) + _c2 * pos[1] );
    Grad[1] = f * _c1 * ( - 2 * pos[1] / aol::Sqr<RealType> ( _sigmaY ) + _c2 * pos[0] );
    Grad[2] = f / _height;
    Grad[3] = f * _c1 * ( - 2 * aol::Sqr<RealType> ( pos[0] / _sigmaX ) / _sigmaX + _c2 / _sigmaX * pos[0] * pos[1] );
    Grad[4] = f * _c1 * ( - 2 * aol::Sqr<RealType> ( pos[1] / _sigmaY ) / _sigmaY + _c2 / _sigmaY * pos[0] * pos[1] );
    Grad[5] = f * ( pos[0] * pos[1] / ( _sigmaX * _sigmaY * ( 1 - aol::Sqr<RealType> ( _rotation ) ) )
              - _rotation / aol::Sqr<RealType> ( ( 1 - aol::Sqr<RealType> ( _rotation ) ) ) * ( aol::Sqr<RealType> ( pos[0] / _sigmaX )
              + aol::Sqr<RealType> ( pos[1] / _sigmaY ) - _c2 * pos[0] * pos[1] ) );
    if ( OffsetIsDOF )
      Grad[6] = 1;
  }

  void evaluateParameterHessian( const aol::Vec<Dim, RealType> &X, aol::Mat<NumberOfParameters,NumberOfParameters,RealType> &Hess ) const {
    aol::Vec<Dim, RealType> pos ( X );
    pos -= _center;
    //throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__);
    //Hess.setZero();
    const RealType f = evaluate ( X ) - _offset;

    Hess[0][0] = 2 * f * _c1 / _sigmaX / _sigmaX + _c1 * _c1 * f * aol::Sqr<RealType> ( 2 * pos[0] / _sigmaX / _sigmaX - _c2 * pos[1] );

    Hess[0][1] = - f * _c1 * _c2 + f * _c1 * _c1 * ( 2 * pos[0] / _sigmaX / _sigmaX - _c2 * pos[1] ) * (2 * pos[1] / _sigmaY / _sigmaY - _c2 * pos[0] );

    Hess[1][1] = 2 * f * _c1 / _sigmaY / _sigmaY + f * _c1 * _c1 * aol::Sqr<RealType> ( 2 * pos[1] / _sigmaY / _sigmaY - _c2 * pos[0] );

    Hess[0][2] = - f * _c1 / _height * ( 2 * pos[0] / _sigmaX / _sigmaX - _c2 * pos[1] );

    Hess[1][2] = f * _c1 / _height * ( - 2 * pos[1] / _sigmaY / _sigmaY  + _c2 * pos[0] );

    Hess[2][2] = 0 ; 

    Hess[0][3] = - f * _c1 * ( - 4 * pos[0] / _sigmaX / _sigmaX / _sigmaX + _c2 * pos[1] / _sigmaX ) - f * _c1 * _c1 * ( 2 * pos[0] / _sigmaX / _sigmaX - _c2 * pos[1] ) * ( - 2 * aol::Sqr<RealType> ( pos[0] ) / _sigmaX / _sigmaX / _sigmaX + _c2 * pos[0] * pos[1] / _sigmaX );

    Hess[1][3] = - f * _c1 * _c2 * pos[0] / _sigmaX - f * _c1 * _c1 * ( - _c2 * pos[0] + 2 * pos[1] / _sigmaY / _sigmaY ) * ( - 2 * aol::Sqr<RealType> ( pos[0] ) / _sigmaX / _sigmaX / _sigmaX + _c2 * pos[0] * pos[1] / _sigmaX );

    Hess[2][3] = f * _c1 / _height * ( - 2 * aol::Sqr<RealType> ( pos[0] ) / _sigmaX / _sigmaX / _sigmaX + _c2 * pos[0] * pos[1] / _sigmaX );

    Hess[3][3] = f * _c1 / _sigmaX / _sigmaX * ( 6 * aol::Sqr<RealType> ( pos[0] / _sigmaX ) - 2 * _c2 * pos[0] * pos[1] ) + f * _c1 * _c1 / _sigmaX / _sigmaX * aol::Sqr<RealType> ( - 2 * aol::Sqr<RealType> ( pos[0] / _sigmaX ) + _c2 * pos[0] * pos[1] );

    Hess[0][4] = - f * _c1 * _c2 * pos[1] / _sigmaY - f * _c1 * _c1 * ( 2 * pos[0] / _sigmaX / _sigmaX - _c2 * pos[1] ) * ( _c2 * pos[0] * pos[1] / _sigmaY - 2 * aol::Sqr<RealType> ( pos[1] / _sigmaY ) / _sigmaY );

    Hess[1][4] = - f * _c1 * ( _c2 * pos[0] / _sigmaY - 4 * pos[1] / _sigmaY / _sigmaY / _sigmaY ) - f * _c1 * _c1 * ( - _c2 * pos[0] + 2 * pos[1] / _sigmaY / _sigmaY ) * ( _c2 * pos[0] * pos[1] / _sigmaY - 2 * aol::Sqr<RealType> ( pos[1] ) / _sigmaY / _sigmaY / _sigmaY) ;

    Hess[2][4] = f * _c1 / _height * ( _c2 * pos[0] * pos[1] / _sigmaY - 2 * aol::Sqr<RealType> ( pos[1] ) / _sigmaY / _sigmaY / _sigmaY );

    Hess[3][4] = - f * _c1 * _c2 * pos[0] * pos[1] / _sigmaX / _sigmaY + f * _c1 * _c1 * ( - 2 * aol::Sqr<RealType> ( pos[0] ) / _sigmaX / _sigmaX / _sigmaX + _c2 * pos[0] * pos[1] / _sigmaX ) * ( _c2 * pos[0] * pos[1] / _sigmaY - 2 * aol::Sqr<RealType> ( pos[1] ) / _sigmaY / _sigmaY / _sigmaY );

    Hess[4][4] = f * _c1 * ( - 2 * _c2 * pos[0] * pos[1] / _sigmaY / _sigmaY + 6 * aol::Sqr<RealType> ( pos[1] / _sigmaY / _sigmaY ) ) + f * _c1 * _c1 * aol::Sqr<RealType> ( _c2 * pos[0] * pos[1] / _sigmaY / _sigmaY - 2 * aol::Sqr<RealType> ( pos[1] ) / _sigmaY / _sigmaY / _sigmaY );

    Hess[0][5] = - f * _c1 * ( - 2 * pos[1] / _sigmaX / _sigmaY - 4 * _c1 * _rotation * ( 2 * pos[0] / _sigmaX / _sigmaX - _c2 * pos[1] ) ) - f * _c1 * ( 2 * pos[0] / _sigmaX / _sigmaX - _c2 * pos[1] ) * ( - 2 * _c1 * pos[0] * pos[1] / _sigmaX / _sigmaY - 4 * _c1 * _c1 * _rotation * ( aol::Sqr<RealType> ( pos[0] / _sigmaX ) - _c2 * pos[0] * pos[1] + aol::Sqr<RealType> ( pos[1] / _sigmaY ) ) );

    Hess[1][5] = - f * _c1 * ( - 2 * pos[0] / _sigmaX / _sigmaY - 4 * _c1 * _rotation * ( 2 * pos[1] / _sigmaY / _sigmaY - _c2 * pos[0] ) ) - f * _c1 * ( 2 * pos[1] / _sigmaY / _sigmaY - _c2 * pos[0] ) * ( - 2 * _c1 * pos[0] * pos[1] / _sigmaX / _sigmaY - 4 * _c1 * _c1 * _rotation * ( aol::Sqr<RealType> ( pos[0] / _sigmaX ) - _c2 * pos[0] * pos[1] + aol::Sqr<RealType> ( pos[1] / _sigmaY ) ) );

    Hess[2][5] = f / _height * ( - 2 * _c1 * pos[0] * pos[1] / _sigmaX / _sigmaY - 4 * _c1 * _c1 * _rotation * ( aol::Sqr<RealType> ( pos[0] / _sigmaX ) - _c2 * pos[0] * pos[1] + aol::Sqr<RealType> ( pos[1] / _sigmaY ) ) );

    Hess[3][5] = f * ( 2 * _c1 * pos[0] * pos[1] / _sigmaX / _sigmaY / _sigmaX - 4 * _c1 * _c1 * _rotation * ( - 2 * aol::Sqr<RealType> ( pos[0] ) / _sigmaX / _sigmaX / _sigmaX + _c2 * pos[0] * pos[1] / _sigmaX ) ) + f * _c1 * ( - 2 * aol::Sqr<RealType> ( pos[0] ) / _sigmaX / _sigmaX / _sigmaX + _c2 * pos[0] * pos[1] / _sigmaX ) * ( - 2 * _c1 * pos[0] * pos[1] / _sigmaX / _sigmaY - 4 * _c1 * _c1 * _rotation * ( aol::Sqr<RealType> ( pos[0] / _sigmaX ) - _c2 * pos[0] * pos[1] + aol::Sqr<RealType> ( pos[1] / _sigmaY ) ) );

    Hess[4][5] = f * ( 2 * _c1 * pos[0] * pos[1] / _sigmaX / _sigmaY / _sigmaY - 4 * _c1 * _c1 * _rotation * ( - 2 * aol::Sqr<RealType> ( pos[1] ) / _sigmaY / _sigmaY / _sigmaY + _c2 * pos[0] * pos[1] / _sigmaY ) ) + f * _c1 * ( - 2 * aol::Sqr<RealType> ( pos[1] ) / _sigmaY / _sigmaY / _sigmaY + _c2 * pos[0] * pos[1] / _sigmaY ) * ( - 2 * _c1 * pos[0] * pos[1] / _sigmaX / _sigmaY - 4 * _c1 * _c1 * _rotation * ( aol::Sqr<RealType> ( pos[0] / _sigmaX ) - _c2 * pos[0] * pos[1] + aol::Sqr<RealType> ( pos[1] / _sigmaY ) ) );

    Hess[5][5] = f * ( 8 * _c1 * _c1 * _c2 * pos[0] * pos[1] - 32 * _c1 * _c1 * _c1 * _rotation * _rotation * ( aol::Sqr<RealType> ( pos[0] / _sigmaX ) - _c2 * pos[0] * pos[1] + aol::Sqr<RealType> ( pos[1] / _sigmaY ) ) - 4 * _c1 * _c1 * ( aol::Sqr<RealType> ( pos[0] / _sigmaX ) - _c2 * pos[0] * pos[1] + aol::Sqr<RealType> ( pos[1] / _sigmaY ) ) ) + f * aol::Sqr<RealType> ( - 2 * _c1 * pos[0] * pos[1] / _sigmaX / _sigmaY - 4 * _c1 * _c1 * _rotation * ( aol::Sqr<RealType> ( pos[0] / _sigmaX ) - _c2 * pos[0] * pos[1] + aol::Sqr<RealType> ( pos[1] / _sigmaY ) ) );  

    Hess[1][0] = Hess[0][1];
    Hess[2][0] = Hess[0][2];
    Hess[3][0] = Hess[0][3];
    Hess[4][0] = Hess[0][4];
    Hess[5][0] = Hess[0][5];
    Hess[2][1] = Hess[1][2];
    Hess[3][1] = Hess[1][3];
    Hess[4][1] = Hess[1][4];
    Hess[5][1] = Hess[1][5];
    Hess[3][2] = Hess[2][3];
    Hess[4][2] = Hess[2][4];
    Hess[5][2] = Hess[2][5];
    Hess[4][3] = Hess[3][4];
    Hess[5][3] = Hess[3][5];
    Hess[5][4] = Hess[4][5];
  }
  
  void evaluateMixedHessian( const aol::Vec<Dim, RealType> &X, aol::Mat<NumberOfParameters,Dim,RealType> &Hess ) const {
    aol::Vec<Dim, RealType> pos ( X );
    pos -= _center;
    //throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__);
    Hess.setZero();
    const RealType f = evaluate ( X ) - _offset;

    Hess[0][0] = - 2 * f * _c1 / _sigmaX / _sigmaX - _c1 * _c1 * f * aol::Sqr<RealType> ( 2 * pos[0] / _sigmaX / _sigmaX - _c2 * pos[1] );

    Hess[0][1] = f * _c1 * _c2 - f * _c1 * _c1 * ( 2 * pos[0] / _sigmaX / _sigmaX - _c2 * pos[1] ) * (2 * pos[1] / _sigmaY / _sigmaY - _c2 * pos[0] );

    Hess[1][0] = Hess[0][1];

    Hess[1][1] = - 2 * f * _c1 / _sigmaY / _sigmaY - f * _c1 * _c1 * aol::Sqr<RealType> ( 2 * pos[1] / _sigmaY / _sigmaY - _c2 * pos[0] );

    Hess[2][0] = f * _c1 / _height * ( 2 * pos[0] / _sigmaX / _sigmaX - _c2 * pos[1] );

    Hess[2][1] = f * _c1 / _height * ( 2 * pos[1] / _sigmaY / _sigmaY - _c2 * pos[0] );

    Hess[3][0] = f * _c1 * ( - 4 * pos[0] / _sigmaX / _sigmaX / _sigmaX + _c2 * pos[1] / _sigmaX ) + f * _c1 * _c1 * ( 2 * pos[0] / _sigmaX / _sigmaX - _c2 * pos[1] ) * ( - 2 * aol::Sqr<RealType> ( pos[0] ) / _sigmaX / _sigmaX / _sigmaX + _c2 * pos[0] * pos[1] / _sigmaX );

    Hess[3][1] = f * _c1 * _c2 * pos[0] / _sigmaX + f * _c1 * _c1 * ( - _c2 * pos[0] + 2 * pos[1] / _sigmaY / _sigmaY ) * ( - 2 * aol::Sqr<RealType> ( pos[0] ) / _sigmaX / _sigmaX / _sigmaX + _c2 * pos[0] * pos[1] / _sigmaX );

    Hess[4][0] = f * _c1 * _c2 * pos[1] / _sigmaY + f * _c1 * _c1 * ( 2 * pos[0] / _sigmaX / _sigmaX - _c2 * pos[1] ) * ( _c2 * pos[0] * pos[1] / _sigmaY - 2 * aol::Sqr<RealType> ( pos[1] / _sigmaY ) / _sigmaY );

    Hess[4][1] = f * _c1 * ( _c2 * pos[0] / _sigmaY - 4 * pos[1] / _sigmaY / _sigmaY / _sigmaY ) + 
    f * _c1 * _c1 * ( - _c2 * pos[0] + 2 * pos[1] / _sigmaY / _sigmaY ) * ( _c2 * pos[0] * pos[1] / _sigmaY - 2 * aol::Sqr<RealType> ( pos[1] ) / _sigmaY / _sigmaY / _sigmaY) ;

    Hess[5][0] = f * _c1 * ( - 2 * pos[1] / _sigmaX / _sigmaY - 4 * _c1 * _rotation * ( 2 * pos[0] / _sigmaX / _sigmaX - _c2 * pos[1] ) ) + 
    f * _c1 * ( 2 * pos[0] / _sigmaX / _sigmaX - _c2 * pos[1] ) * ( - 2 * _c1 * pos[0] * pos[1] / _sigmaX / _sigmaY - 4 * _c1 * _c1 * _rotation * ( aol::Sqr<RealType> ( pos[0] / _sigmaX ) - _c2 * pos[0] * pos[1] + aol::Sqr<RealType> ( pos[1] / _sigmaY ) ) );

    Hess[5][1] = f * _c1 * ( - 2 * pos[0] / _sigmaX / _sigmaY - 4 * _c1 * _rotation * ( 2 * pos[1] / _sigmaY / _sigmaY - _c2 * pos[0] ) ) + 
    f * _c1 * ( 2 * pos[1] / _sigmaY / _sigmaY - _c2 * pos[0] ) * ( - 2 * _c1 * pos[0] * pos[1] / _sigmaX / _sigmaY - 4 * _c1 * _c1 * _rotation * ( aol::Sqr<RealType> ( pos[0] / _sigmaX ) - _c2 * pos[0] * pos[1] + aol::Sqr<RealType> ( pos[1] / _sigmaY ) ) );

  }

  const aol::Vec2<RealType>& getCenter ( ) const {
    return _center;
  }

  RealType estimateSupportRadius ( const RealType MinValue ) const {
    if ( aol::appeqAbsolute ( _offset, aol::ZTrait<RealType>::zero ) )
      return sqrt ( - 2 * aol::Sqr ( aol::Max ( _sigmaX, _sigmaY ) ) * ( 1 + _rotation ) * log ( ( MinValue - _offset ) / _height ) );
    else
      throw aol::ParameterRangeException( "estimateSupportRadius only makes sense if the offset is zero.", __FILE__, __LINE__ );
  }
};


template <typename RealType>
class SingleAsymmetricBumpFitTargetFunctional : public aol::Op<aol::Vector<RealType> > {
  const qc::ScalarArray<RealType, qc::QC_2D> &_data;
  const qc::FastILexMapper<qc::QC_2D> _mapper;
public:
  SingleAsymmetricBumpFitTargetFunctional ( const qc::ScalarArray<RealType, qc::QC_2D> &Data )
    : _data ( Data ), _mapper ( Data ) {}

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    AsymmetricGaussianBumpFunction<RealType> bump ( Arg );
    for ( qc::RectangularIterator<qc::QC_2D> it ( _data ); it.notAtEnd ( ) ; ++it ) {
      aol::Vec<qc::QC_2D, RealType> pos;
      pos.readFromBuffer ( (*it).getData ( ) );
      const RealType g = _data.get ( *it );
      Dest[_mapper.getGlobalIndex ( *it )] += ( !aol::isNaN ( g ) ? bump.evaluate ( pos ) - g : 0 );
    }
  }
};


template <typename RealType>
class SingleAsymmetricBumpFitTargetJacobian : public aol::Op<aol::Vector<RealType>, aol::FullMatrix<RealType> > {
  const qc::ScalarArray<RealType, qc::QC_2D> &_data;
  const qc::FastILexMapper<qc::QC_2D> _mapper;
public:
  SingleAsymmetricBumpFitTargetJacobian ( const qc::ScalarArray<RealType, qc::QC_2D> &Data )
    : _data ( Data ), _mapper ( Data ) { }

  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::FullMatrix<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__);
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::FullMatrix<RealType> &Dest ) const {
    AsymmetricGaussianBumpFunction<RealType> bump ( Arg );
    aol::Vector<RealType> gradient ( bump.NumberOfParameters );
    int i = 0;
    for ( qc::RectangularIterator<qc::QC_2D> it ( _data ); it.notAtEnd ( ) ; ++it ) {
      aol::Vec<qc::QC_2D, RealType> pos;
      pos.readFromBuffer ( (*it).getData ( ) );
      if ( aol::isFinite ( _data.get ( *it ) ) )
        bump.evaluateParameterGradient ( pos, gradient );
      else {
        gradient.setZero ( );
        gradient[bump.NumberOfParameters - 1] = 1;
      }
      for ( int j=0; j<bump.NumberOfParameters ; ++j )
        Dest.set ( i, j, gradient[j] );
      ++i;
    }
  }
};


template <typename RealType>
class AsymmetricGaussianDoubleBumpFunction {
public:
  static const int NumberOfParameters = 13;
  
private:
  const aol::Vec2<RealType> _center1, _center2;
  const RealType _height1, _height2;
  const RealType _sigmaX1, _sigmaY1, _sigmaX2, _sigmaY2;
  const RealType _rotation1, _rotation2;
  const RealType _offset;
  const RealType _c11, _c21, _c12, _c22;
  
public:
  AsymmetricGaussianDoubleBumpFunction ( const aol::Vec2<RealType> &Center1, const aol::Vec2<RealType> &Center2,
                                         const RealType Height1, const RealType Height2,
                                         const RealType SigmaX1, const RealType SigmaY1, const RealType SigmaX2, const RealType SigmaY2,
                                         const RealType Rotation1, const RealType Rotation2,
                                         const RealType Offset )
  : _center1 ( Center1 ), _center2 ( Center2 ),
    _height1 ( Height1 ), _height2 ( Height2 ),
    _sigmaX1 ( SigmaX1 ), _sigmaY1 ( SigmaY1 ), _sigmaX2 ( SigmaX2 ), _sigmaY2 ( SigmaY2 ),
    _rotation1 ( Rotation1 ), _rotation2 ( Rotation2 ),
    _offset ( Offset ),
    _c11 ( - 1 / ( 2 * ( 1 - aol::Sqr<RealType> ( _rotation2 ) ) ) ),
    _c21 ( 2 * _rotation1 / ( _sigmaX1 * _sigmaY1 ) ),
    _c12 ( - 1 / ( 2 * ( 1 - aol::Sqr<RealType> ( _rotation2 ) ) ) ),
    _c22 ( 2 * _rotation2 / ( _sigmaX2 * _sigmaY2 ) ) { }
  
  AsymmetricGaussianDoubleBumpFunction ( const aol::Vector<RealType> &Parameters )
    : _center1 ( Parameters[0], Parameters[1] ), _center2 ( Parameters[2], Parameters[3] ),
      _height1 ( Parameters[4] ), _height2 ( Parameters[5] ),
      _sigmaX1 ( Parameters[6] ), _sigmaY1 ( Parameters[7] ), _sigmaX2 ( Parameters[8] ), _sigmaY2 ( Parameters[9] ),
      _rotation1 ( Parameters[10] ), _rotation2 ( Parameters[11] ),
      _offset ( Parameters[12] ),
      _c11 ( - 1 / ( 2 * ( 1 - pow ( _rotation1, 2 ) ) ) ),
      _c21 ( 2 * _rotation1 / ( _sigmaX1 * _sigmaY1 ) ),
      _c12 ( - 1 / ( 2 * ( 1 - pow ( _rotation2, 2 ) ) ) ),
      _c22 ( 2 * _rotation2 / ( _sigmaX2 * _sigmaY2 ) ) { }
  
  RealType evaluate ( const aol::Vec<2, RealType> &X ) const {
    aol::Vec<2, RealType> pos1 ( X ), pos2 ( X );
    pos1 -= _center1;
    pos2 -= _center2;
    
    return _height1 * exp ( _c11 * ( aol::Sqr<RealType> ( pos1[0] / _sigmaX1 )
             + aol::Sqr<RealType> ( pos1[1] / _sigmaY1 ) - _c21 * pos1[0] * pos1[1] ) )
           + _height2 * exp ( _c12 * ( aol::Sqr<RealType> ( pos2[0] / _sigmaX2 )
             + aol::Sqr<RealType> ( pos2[1] / _sigmaY2 ) - _c22 * pos2[0] * pos2[1] ) )
           + _offset;
  }
  
  void evaluateParameterGradient ( const aol::Vec<2, RealType> &X, aol::Vector<RealType> &Grad ) const {
    aol::Vec<2, RealType> pos1 ( X ), pos2 ( X );
    pos1 -= _center1; pos2 -= _center2;
    
    const RealType c11 = aol::Sqr<RealType> ( pos1[0] / _sigmaX1 ), c21 = aol::Sqr<RealType> ( pos1[1] / _sigmaY1 ),
                   c12 = aol::Sqr<RealType> ( pos2[0] / _sigmaX2 ), c22 = aol::Sqr<RealType> ( pos2[1] / _sigmaY2 ),
                   f1 = _height1 * exp ( _c11 * ( pow ( pos1[0] / _sigmaX1, 2 ) + pow ( pos1[1] / _sigmaY1, 2 ) - _c21 * pos1[0] * pos1[1] ) ),
                   f2 = _height2 * exp ( _c12 * ( pow ( pos2[0] / _sigmaX2, 2 ) + pow ( pos2[1] / _sigmaY2, 2 ) - _c22 * pos2[0] * pos2[1] ) );
    Grad[0] = f1 * _c11 * ( -2 * pos1[0] / aol::Sqr<RealType> ( _sigmaX1 ) + _c21 * pos1[1] );
    Grad[1] = f1 * _c11 * ( -2 * pos1[1] / aol::Sqr<RealType> ( _sigmaY1 ) + _c21 * pos1[0] );
    Grad[2] = f2 * _c12 * ( -2 * pos2[0] / aol::Sqr<RealType> ( _sigmaX2 ) + _c22 * pos2[1] );
    Grad[3] = f2 * _c12 * ( -2 * pos2[1] / aol::Sqr<RealType> ( _sigmaY2 ) + _c22 * pos2[0] );
    Grad[4] = f1 / _height1;
    Grad[5] = f2 / _height2;
    Grad[6] = f1 * _c11 * ( -2 * c11 / _sigmaX1 + _c21 / _sigmaX1 * pos1[0] * pos1[1] );
    Grad[7] = f1 * _c11 * ( -2 * c21 / _sigmaY1 + _c21 / _sigmaY1 * pos1[0] * pos1[1] );
    Grad[8] = f2 * _c12 * ( -2 * c12 / _sigmaX2 + _c22 / _sigmaX2 * pos2[0] * pos2[1] );
    Grad[9] = f2 * _c12 * ( -2 * c22 / _sigmaY2 + _c22 / _sigmaY2 * pos2[0] * pos2[1] );
    Grad[10] = f1 * ( pos1[0] * pos1[1] / ( _sigmaX1 * _sigmaY1 * ( 1 - aol::Sqr<RealType> ( _rotation1 ) ) )
               - _rotation1 / aol::Sqr<RealType> ( ( 1 - aol::Sqr<RealType> ( _rotation1 ) ) ) * ( c11 + c21 - _c21 * pos1[0] * pos1[1] ) );
    Grad[11] = f2 * ( pos2[0] * pos2[1] / ( _sigmaX2 * _sigmaY2 * ( 1 - aol::Sqr<RealType> ( _rotation2 ) ) )
               - _rotation2 / aol::Sqr<RealType> ( ( 1 - aol::Sqr<RealType> ( _rotation2 ) ) ) * ( c12 + c22 - _c22 * pos2[0] * pos2[1] ) );
    Grad[12] = 1;
  }
};


template <typename RealType>
class SingleAsymmetricDoubleBumpFitTargetFunctional : public aol::Op<aol::Vector<RealType> > {
  const qc::ScalarArray<RealType, qc::QC_2D> &_data;
  const qc::FastILexMapper<qc::QC_2D> _mapper;
public:
  SingleAsymmetricDoubleBumpFitTargetFunctional ( const qc::ScalarArray<RealType, qc::QC_2D> &Data )
    : _data ( Data ), _mapper ( Data ) {}
  
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    AsymmetricGaussianDoubleBumpFunction<RealType> bump ( Arg );
    for ( qc::RectangularIterator<qc::QC_2D> it ( _data ); it.notAtEnd ( ) ; ++it ) {
      aol::Vec<qc::QC_2D, RealType> pos;
      pos.readFromBuffer ( (*it).getData ( ) );
      const RealType g = _data.get ( *it );
      Dest[_mapper.getGlobalIndex ( *it )] += ( !aol::isNaN ( g ) ? bump.evaluate ( pos ) - g : 0 );
    }
  }
};


template <typename RealType>
class SingleAsymmetricDoubleBumpFitTargetJacobian : public aol::Op<aol::Vector<RealType>, aol::FullMatrix<RealType> > {
  const qc::ScalarArray<RealType, qc::QC_2D> &_data;
  const qc::FastILexMapper<qc::QC_2D> _mapper;
public:
  SingleAsymmetricDoubleBumpFitTargetJacobian ( const qc::ScalarArray<RealType, qc::QC_2D> &Data )
  : _data ( Data ), _mapper ( Data ) { }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::FullMatrix<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__);
  }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::FullMatrix<RealType> &Dest ) const {
    AsymmetricGaussianDoubleBumpFunction<RealType> bump ( Arg );
    aol::Vector<RealType> gradient ( bump.NumberOfParameters );
    int i = 0;
    for ( qc::RectangularIterator<qc::QC_2D> it ( _data ); it.notAtEnd ( ) ; ++it ) {
      aol::Vec<qc::QC_2D, RealType> pos;
      pos.readFromBuffer ( (*it).getData ( ) );
      if ( aol::isFinite ( _data.get ( *it ) ) )
        bump.evaluateParameterGradient ( pos, gradient );
      else {
        gradient.setZero ( );
        gradient[bump.NumberOfParameters - 1] = 1;
      }
      for ( int j=0; j<bump.NumberOfParameters ; ++j )
        Dest.set ( i, j, gradient[j] );
      ++i;
    }
  }
};


//void testAsymmetricGaussianBumpFunction ( int argc, char **argv ) {
//  aol::ParameterParser parser ( argc, argv, "bumpFunctionTest.par" );
//  aol::Vector<RealType> parameters;
//  parser.getRealVec( "parameters", parameters );
//
//  aol::Vector<RealType> positionTmp;
//  parser.getRealVec ( "position", positionTmp );
//  aol::Vec2<RealType> x ( positionTmp[0], positionTmp[1] );
//
//  AsymmetricGaussianBumpFunction<RealType> bump ( parameters );
//  aol::Vector<RealType> grad ( 6 );
//  bump.evaluateParameterGradient ( x, grad );
//
//  std::cerr << bump.evaluate ( x ) << std::endl;
//  std::cerr << grad << std::endl;
//}
  
  
/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class BumpFitCenterEstimator {
  typedef typename ConfiguratorType::RealType RealType;
  const RealType _smoothSigma;
  const int _eraseInfinityRadius;
  const RealType _minCenterValue;
  const int _maxNumCenters;
  const int _minBoundaryDistance;
public:
  BumpFitCenterEstimator ( const RealType SmoothSigma, const int EraseInfinityRadius, const RealType MinCenterValue = 0.5, const int MaxNumCenters = 0 )
  : _smoothSigma ( SmoothSigma ),
  _eraseInfinityRadius ( EraseInfinityRadius ),
  _minCenterValue ( MinCenterValue ),
  _maxNumCenters ( MaxNumCenters ),
  _minBoundaryDistance ( 0 ) {}

  BumpFitCenterEstimator ( const aol::ParameterParser &Parser, const typename ConfiguratorType::ArrayType &Image )
  : _smoothSigma ( Parser.getReal<RealType> ( "iniGuessSmoothSigma" ) ),
  _eraseInfinityRadius ( Parser.getInt ( "iniGuessEraseInfinityRadius" ) ),
  _minCenterValue ( Image.getMinValue() + Parser.getReal<RealType> ( "iniGuessMinCenterIntensityPercentage" ) * ( Image.getMaxValue() - Image.getMinValue() ) ),
  _maxNumCenters ( Parser.getInt ( "iniGuessMaxNumCenters" ) ),
  _minBoundaryDistance ( Parser.getIntOrDefault ( "iniGuessMinBoundaryDistance", 0 ) ) {}

  void guessCenters ( const typename ConfiguratorType::ArrayType &Image, aol::MultiVector<RealType> &GuessedCenters ) {
    GuessedCenters.reallocate ( 0, 0 );

    typename ConfiguratorType::ArrayType image ( Image );

    if ( image.checkForNANsAndINFs () )
      throw aol::InconsistentDataException( "Input data contains non-finite values.", __FILE__, __LINE__);

    if ( _smoothSigma > 0 ) {
      typename ConfiguratorType::InitType grid ( qc::GridSize<ConfiguratorType::Dim>::createFrom ( image ) );
      qc::GeneralLinearSmoothOp<ConfiguratorType> smoother ( grid, _smoothSigma * grid.H() );
      smoother.applySingle ( image );
    }

    if ( _minBoundaryDistance > 0 ) {
      image.setBlock ( aol::Vec2<int> ( 0, 0 ), aol::Vec2<int> ( image.getNumX(), _minBoundaryDistance ), 0 );
      image.setBlock ( aol::Vec2<int> ( 0, 0 ), aol::Vec2<int> ( _minBoundaryDistance, image.getNumY() ), 0 );
      image.setBlock ( aol::Vec2<int> ( image.getNumX() - _minBoundaryDistance, 0 ), aol::Vec2<int> ( image.getNumX(), image.getNumY() ), 0 );
      image.setBlock ( aol::Vec2<int> ( 0, image.getNumY() - _minBoundaryDistance ), aol::Vec2<int> ( image.getNumX(), image.getNumY() ), 0 );
    }

    while ( ( image.getMaxIndexAndValue().second > _minCenterValue )
           && ( ( _maxNumCenters == 0 ) || ( GuessedCenters.numComponents() < _maxNumCenters ) ) ) {
      const int index = image.getMaxIndexAndValue().first;
      const qc::CoordType pos ( index % image.getNumX(), index / image.getNumX() );

      const int centerNum = GuessedCenters.numComponents() + 1;
      GuessedCenters.resize ( centerNum, 2 );
      for ( int i = 0; i < 2; ++i )
        GuessedCenters[centerNum -1][i] = pos[i];

      for ( int i = aol::Max ( 0, pos[0] - _eraseInfinityRadius ); i < aol::Min ( image.getNumX(), pos[0]+_eraseInfinityRadius ); ++i )
        for ( int j = aol::Max ( 0, pos[1] - _eraseInfinityRadius ); j < aol::Min ( image.getNumY(), pos[1]+_eraseInfinityRadius ); ++j )
          image.set( i, j, 0 );
    }
  }
};


} // end namespace


#endif
