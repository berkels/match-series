#ifndef __AMBROSIOTORTORELLI_H
#define __AMBROSIOTORTORELLI_H

#include <FEOpInterface.h>
#include <deformations.h>

namespace qc {

/**
 * \brief \f$ \left(\int_\Omega v(\phi(x))^2 \nabla\varphi_i(x)\cdot\nabla\varphi_j(x)dx\right)_{ij} \f$.
 * \ingroup MatrixFEOp
 * \todo rename to qc::SquaredDeformedWeightStiffOp
 */
template <typename ConfiguratorType>
class ATStiffV2DeformOp : public aol::FELinScalarWeightedStiffInterface<ConfiguratorType, ATStiffV2DeformOp<ConfiguratorType> > {
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrV;
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType, ConfiguratorType::Dim> _discrPhi;
public:
  typedef typename ConfiguratorType::RealType RealType;

  ATStiffV2DeformOp ( const typename ConfiguratorType::InitType &Initializer,
                      const aol::Vector<RealType> &VDofs,
                      const aol::MultiVector<RealType> &PhiDofs,
                      aol::OperatorType OpType = aol::ONTHEFLY )
      : aol::FELinScalarWeightedStiffInterface<ConfiguratorType, ATStiffV2DeformOp<ConfiguratorType> > ( Initializer, OpType ),
      _discrV ( Initializer, VDofs ),
      _discrPhi ( Initializer, PhiDofs ) {}
  inline RealType getCoeff ( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::DomVecType& RefCoord ) const {
    typename ConfiguratorType::VecType offset, transformed_local_coord;
    qc::Element transformed_el;
    _discrPhi.evaluateAtQuadPoint ( El, QuadPoint, offset );
    if ( !qc::transformCoord<ConfiguratorType> ( this->_grid, El, RefCoord, offset, transformed_el, transformed_local_coord ) ) {
      return 0;
    }
    return aol::Sqr ( _discrV.evaluate ( transformed_el, transformed_local_coord ) );
  }
};


/**
 * \brief \f$ \left(\int_\Omega w(x) v(\phi(x))^2 \nabla\varphi_i(x)\cdot\nabla\varphi_j(x)dx\right)_{ij} \f$.
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType>
class ATWeightedStiffV2DeformOp : public aol::FELinScalarWeightedStiffInterface<ConfiguratorType, ATWeightedStiffV2DeformOp<ConfiguratorType> > {
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrWeight;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrV;
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType, ConfiguratorType::Dim> _discrPhi;
public:
  typedef typename ConfiguratorType::RealType RealType;

  ATWeightedStiffV2DeformOp ( const typename ConfiguratorType::InitType &Initializer,
                              const aol::Vector<RealType> &WeightDofs,
                              const aol::Vector<RealType> &VDofs,
                              const aol::MultiVector<RealType> &PhiDofs,
                              aol::OperatorType OpType = aol::ONTHEFLY )
      : aol::FELinScalarWeightedStiffInterface<ConfiguratorType, ATWeightedStiffV2DeformOp<ConfiguratorType> > ( Initializer, OpType ),
      _discrWeight ( Initializer, WeightDofs ),
      _discrV ( Initializer, VDofs ),
      _discrPhi ( Initializer, PhiDofs ) {}

  inline RealType getCoeff ( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::DomVecType& RefCoord ) const {
    const RealType weight = _discrWeight.evaluateAtQuadPoint ( El, QuadPoint );
    if ( weight != 0. ) {
      typename ConfiguratorType::VecType offset, transformed_local_coord;
      qc::Element transformed_el;
      _discrPhi.evaluateAtQuadPoint ( El, QuadPoint, offset );
      if ( !qc::transformCoord<ConfiguratorType> ( this->_grid, El, RefCoord, offset, transformed_el, transformed_local_coord ) ) {
        return 0;
      }
      return weight * aol::Sqr ( _discrV.evaluate ( transformed_el, transformed_local_coord ) );
    } else return 0.;
  }
};


//! Deprecated, use aol::SquaredWeightMassOp instead.
template <typename ConfiguratorType>
class ATMassV2Op : public aol::FELinScalarWeightedMassInterface<ConfiguratorType, ATMassV2Op<ConfiguratorType> > {
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrV;
public:
  typedef typename ConfiguratorType::RealType RealType;

  ATMassV2Op ( const typename ConfiguratorType::InitType &Initializer,
               const aol::Vector<RealType> &VDofs,
               aol::OperatorType OpType = aol::ONTHEFLY )
      : aol::FELinScalarWeightedMassInterface<ConfiguratorType, ATMassV2Op<ConfiguratorType> > ( Initializer, OpType ), _discrV ( Initializer, VDofs ) {}

  inline RealType getCoeff ( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::DomVecType& /*RefCoord*/ ) const {
    return aol::Sqr ( _discrV.evaluateAtQuadPoint ( El, QuadPoint ) );
  }
};

//! Deprecated, use aol::SquaredWeightStiffOp instead.
template <typename ConfiguratorType>
class ATStiffV2Op : public aol::FELinScalarWeightedStiffInterface<ConfiguratorType, ATStiffV2Op<ConfiguratorType> > {
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrV;
public:
  typedef typename ConfiguratorType::RealType RealType;

  ATStiffV2Op ( const typename ConfiguratorType::InitType &Initializer,
                const aol::Vector<RealType> &VDofs,
                aol::OperatorType OpType = aol::ONTHEFLY )
      : aol::FELinScalarWeightedStiffInterface<ConfiguratorType, ATStiffV2Op<ConfiguratorType> > ( Initializer, OpType ), _discrV ( Initializer, VDofs ) {}

  inline RealType getCoeff ( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::DomVecType& /*RefCoord*/ ) const {
    return aol::Sqr ( _discrV.evaluateAtQuadPoint ( El, QuadPoint ) );
  }
};

/**
 * \brief \f$f(x)=\sqrt{(\lambda_1x_1)^2 +(\lambda_2x_2)^2+ \delta^2)}\f$
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template<typename ConfiguratorType>
class ATAnisotropy2Norm {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  const RealType _deltaSqr;
  aol::Vec2<RealType> _lambda;
public:
  ATAnisotropy2Norm ( const RealType Delta, const RealType lambda_x, const RealType lambda_y ) : _deltaSqr ( pow ( Delta, 2 ) ) {
    if ( ConfiguratorType::Dim != 2 )
      throw aol::Exception ( "ATAnisotropy2Norm only implemented for d=2 yet", __FILE__, __LINE__ );
    _lambda[0] = lambda_x;
    _lambda[1] = lambda_y;
  }
~ATAnisotropy2Norm() {}
  RealType evaluate ( const typename ConfiguratorType::VecType &x ) const {
    return ( sqrt ( pow ( _lambda[0]*x[0], 2 ) + pow ( _lambda[1]*x[1], 2 ) + _deltaSqr ) );
  }
  void evaluateGradient ( const typename ConfiguratorType::VecType &x, typename ConfiguratorType::VecType &Gradient ) const {
    /*
        if( (x[0] == 0.) && (x[1] == 0.) ){
          Gradient[0] = 0.;
          Gradient[1] = 0.;
        } else {
    */
    RealType tmp = evaluate ( x );
    Gradient[0] = ( _lambda[0] * x[0] ) / tmp;
    Gradient[1] = ( _lambda[1] * x[1] ) / tmp;
//    }
  }
  void evaluateHessian ( const typename ConfiguratorType::VecType &x, typename ConfiguratorType::MatType &Hessian ) const {
    /*
        if( (x[0] == 0.) && (x[1] == 0.) ){
          Hessian[0][0] = 0.;
          Hessian[1][1] = 0.;
          Hessian[0][1] = 0.;
          Hessian[1][0] = 0.;
        } else {
    */
    RealType tmp = evaluate ( x );
    Hessian[0][0] = _lambda[0] / tmp - ( pow ( _lambda[0] * x[0], 2 ) ) / ( pow ( tmp, 3 ) );
    Hessian[1][1] = _lambda[1] / tmp - ( pow ( _lambda[1] * x[1], 2 ) ) / ( pow ( tmp, 3 ) );
    Hessian[0][1] = ( -1. * _lambda[0] * _lambda[1] * x[0] * x[1] ) / ( pow ( tmp, 3 ) );
    Hessian[1][0] = Hessian[0][1];
//    }
  }
private:
  ATAnisotropy2Norm ( const ATAnisotropy2Norm &Anisotropy ) : _deltaSqr ( Anisotropy._deltaSqr ), _lambda ( Anisotropy._lambda ) {}
};

/**
 * \brief \f$f(x)=\sum_{i=1}^d \sqrt{x_i^2 + \delta^2}\f$
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template<typename ConfiguratorType>
class ATAnisotropy1Norm {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  const RealType _deltaSqr;
  const int _dim;
public:
  ATAnisotropy1Norm ( const RealType Delta ) : _deltaSqr ( pow ( Delta, 2 ) ), _dim ( ConfiguratorType::Dim ) {}
  ~ATAnisotropy1Norm() {}
  RealType evaluate ( const typename ConfiguratorType::VecType &x ) const {
    RealType value = 0.;
    for ( int i = 0; i < _dim; i++ ) {
      value += sqrt ( pow ( x[i], 2 ) + _deltaSqr );
    }
    return value;
  }
  void evaluateGradient ( const typename ConfiguratorType::VecType &x, typename ConfiguratorType::VecType &Gradient ) const {
    for ( int i = 0; i < _dim; i++ ) {
      RealType tmp = sqrt ( pow ( x[i], 2 ) + _deltaSqr );
      Gradient[i] = x[i] / tmp;
    }
  }
  void evaluateHessian ( const typename ConfiguratorType::VecType &x, typename ConfiguratorType::MatType &Hessian ) const {
    for ( int i = 0; i < _dim; i++ ) {
      RealType tmp = sqrt ( pow ( x[i], 2 ) + _deltaSqr );
      Hessian[i][i] = 1. / tmp - ( pow ( x[i], 2 ) ) / ( pow ( tmp, 3 ) );
      for ( int j = 0; j < i; j++ ) {
        Hessian[i][j] = 0.;
        Hessian[j][i] = 0.;
      }
    }
  }
private:
  ATAnisotropy1Norm ( const ATAnisotropy1Norm &Anisotropy ) : _deltaSqr ( Anisotropy._deltaSqr ) {}
};

/**
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template<typename ConfiguratorType>
class ATAnisotropyOktaeder {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  const int _numberOfDirections;
  vector<typename ConfiguratorType::VecType> _directions;
  RealType _scalingConstant;
  const RealType _deltaSqr;
public:
  ATAnisotropyOktaeder ( const RealType Delta ) : _numberOfDirections ( 4 ), _directions ( _numberOfDirections ), _deltaSqr ( pow ( Delta, 2 ) ) {
    if ( ConfiguratorType::Dim != 2 )
      throw aol::Exception ( "ATAnisotropyOktaeder only implemented for d=2 yet", __FILE__, __LINE__ );
    _directions[0][0] = 1.;
    _directions[1][1] = 1.;
    _directions[2][0] = 1. / sqrt ( 2. );
    _directions[2][1] = _directions[2][0];
    _directions[3][0] = -1. * _directions[2][0];
    _directions[3][1] = _directions[2][0];

    typename ConfiguratorType::VecType x;
    x[0] = 1.;
    _scalingConstant = 0.;
    for ( int i = 0; i < _numberOfDirections; i++ )
      _scalingConstant += fabs ( _directions[i] * x );
  }
~ATAnisotropyOktaeder() {}
  RealType evaluate ( const typename ConfiguratorType::VecType &x ) const {
    RealType value = 0.;
    for ( int i = 0; i < _numberOfDirections; i++ )
      value += sqrt ( pow ( ( _directions[i] * x ), 2 ) + _deltaSqr );
    return value / _scalingConstant;
  }
  void evaluateGradient ( const typename ConfiguratorType::VecType &x, typename ConfiguratorType::VecType &Gradient ) const {
    for ( int i = 0; i < 2; i++ ) {
      Gradient[i] = 0.;
      for ( int k = 0; k < _numberOfDirections; k++ ) {
        Gradient[i] += _directions[k][i] * ( x * _directions[k] ) / sqrt ( pow ( ( _directions[k] * x ), 2 ) + _deltaSqr );
      }
      Gradient[i] /= _scalingConstant;
    }
  }
  void evaluateHessian ( const typename ConfiguratorType::VecType &x, typename ConfiguratorType::MatType &Hessian ) const {
    RealType tmp;
    for ( int i = 0; i < 2; i++ ) {
      for ( int j = 0; j < 2; j++ ) {
        Hessian[i][j] = 0.;
        for ( int k = 0; k < _numberOfDirections; k++ ) {
          tmp = sqrt ( pow ( ( _directions[k] * x ), 2 ) + _deltaSqr );
          Hessian[i][j] += _directions[k][i] * _directions[k][j] * ( 1. / tmp - pow ( ( x * _directions[k] ), 2 ) / pow ( tmp, 3 ) );
        }
        Hessian[i][j] /= _scalingConstant;
      }
    }
  }
};


template <typename ConfiguratorType, typename AnisotropyType>
class ATAnisoGradSqrOp : public aol::FENonlinOpInterface< ConfiguratorType, ATAnisoGradSqrOp<ConfiguratorType, AnisotropyType> > {
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrV;
  const AnisotropyType &_anisotropy;
public:
  ATAnisoGradSqrOp ( const typename ConfiguratorType::InitType &Initializer,
                     const aol::Vector<RealType> &VDofs,
                     const AnisotropyType &Anisotropy )
      : aol::FENonlinOpInterface< ConfiguratorType, ATAnisoGradSqrOp<ConfiguratorType, AnisotropyType> > ( Initializer ),
      _discrV ( Initializer, VDofs ),
      _anisotropy ( Anisotropy ) {}

  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &/*DiscFunc*/,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/,
                         typename ConfiguratorType::RealType &NL ) const {
    typename ConfiguratorType::VecType gradV;
    _discrV.evaluateGradientAtQuadPoint ( El, QuadPoint, gradV );
    RealType AnisotropyValue = _anisotropy.evaluate ( gradV );
    NL = pow ( AnisotropyValue, 2 );
  }
};

template <typename ConfiguratorType, typename AnisotropyType>
class ATAnisotropyStiffOp : public aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, ATAnisotropyStiffOp<ConfiguratorType, AnisotropyType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrV;
  const AnisotropyType &_anisotropy;
  const int _dim;
public:
  ATAnisotropyStiffOp ( const typename ConfiguratorType::InitType &Initializer,
                        const aol::Vector<RealType> &VDofs,
                        const AnisotropyType &Anisotropy,
                        aol::OperatorType OpType = aol::ONTHEFLY )
      : aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, ATAnisotropyStiffOp<ConfiguratorType, AnisotropyType> > ( Initializer, OpType ),
      _discrV ( Initializer, VDofs ),
      _anisotropy ( Anisotropy ),
      _dim ( ConfiguratorType::Dim ) {}

  inline void getCoeffMatrix ( const qc::Element &El, int QuadPoint,
                               const typename ConfiguratorType::DomVecType& /*RefCoord*/, typename ConfiguratorType::MatType &Matrix ) const {
    typename ConfiguratorType::VecType gradV;
    _discrV.evaluateGradientAtQuadPoint ( El, QuadPoint, gradV );
    RealType AnisotropyValue = _anisotropy.evaluate ( gradV );
    typename ConfiguratorType::VecType AnisotropyGradient;
    _anisotropy.evaluateGradient ( gradV, AnisotropyGradient );
    typename ConfiguratorType::MatType AnisotropyHessian;
    _anisotropy.evaluateHessian ( gradV, AnisotropyHessian );

    for ( int i = 0; i < _dim; i++ ) {
      for ( int j = 0; j < _dim; j++ ) {
        Matrix[i][j] = AnisotropyGradient[i] * AnisotropyGradient[j] + AnisotropyValue * AnisotropyHessian[i][j];
      }
    }
  }
};

template <typename ConfiguratorType, typename AnisotropyType>
class ATAnisotropyStiffOpForU : public aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, ATAnisotropyStiffOpForU<ConfiguratorType, AnisotropyType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrU;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrV;
  const AnisotropyType &_anisotropy;
  const int _dim;
public:
  ATAnisotropyStiffOpForU ( const typename ConfiguratorType::InitType &Initializer,
                            const aol::Vector<RealType> &UDofs,
                            const aol::Vector<RealType> &VDofs,
                            const AnisotropyType &Anisotropy,
                            aol::OperatorType OpType = aol::ONTHEFLY )
      : aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, ATAnisotropyStiffOpForU<ConfiguratorType, AnisotropyType> > ( Initializer, OpType ),
      _discrU ( Initializer, UDofs ),
      _discrV ( Initializer, VDofs ),
      _anisotropy ( Anisotropy ),
      _dim ( ConfiguratorType::Dim ) {}

  inline void getCoeffMatrix ( const qc::Element &El, int QuadPoint,
                               const typename ConfiguratorType::DomVecType& /*RefCoord*/, typename ConfiguratorType::MatType &Matrix ) const {
    typename ConfiguratorType::VecType gradU;
    _discrU.evaluateGradientAtQuadPoint ( El, QuadPoint, gradU );
    RealType AnisotropyValue = _anisotropy.evaluate ( gradU );
    typename ConfiguratorType::VecType AnisotropyGradient;
    _anisotropy.evaluateGradient ( gradU, AnisotropyGradient );
    typename ConfiguratorType::MatType AnisotropyHessian;
    _anisotropy.evaluateHessian ( gradU, AnisotropyHessian );

    for ( int i = 0; i < _dim; i++ ) {
      for ( int j = 0; j < _dim; j++ ) {
        Matrix[i][j] = pow ( _discrV.evaluateAtQuadPoint ( El, QuadPoint ), 2 ) * ( AnisotropyGradient[i] * AnisotropyGradient[j] + AnisotropyValue * AnisotropyHessian[i][j] );
      }
    }
  }
};

template <typename ConfiguratorType, typename AnisotropyType>
class ATAnisotropyDiffOp : public aol::FENonlinDiffOpInterface< ConfiguratorType, ATAnisotropyDiffOp<ConfiguratorType, AnisotropyType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrV;
  const AnisotropyType &_anisotropy;
  const int _dim;
public:
  ATAnisotropyDiffOp ( const typename ConfiguratorType::InitType &Initializer,
                       const aol::Vector<RealType> &VDofs,
                       const AnisotropyType &Anisotropy )
      : aol::FENonlinDiffOpInterface< ConfiguratorType, ATAnisotropyDiffOp<ConfiguratorType, AnisotropyType> > ( Initializer ),
      _discrV ( Initializer, VDofs ),
      _anisotropy ( Anisotropy ),
      _dim ( ConfiguratorType::Dim ) {}


  inline void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &/*DiscFunc*/,
                                const typename ConfiguratorType::ElementType &El,
                                int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/,
                                typename ConfiguratorType::VecType &NL ) const {
    typename ConfiguratorType::VecType gradV;
    _discrV.evaluateGradientAtQuadPoint ( El, QuadPoint, gradV );
    RealType AnisotropyValue = _anisotropy.evaluate ( gradV );
    typename ConfiguratorType::VecType AnisotropyGradient;
    _anisotropy.evaluateGradient ( gradV, AnisotropyGradient );

    for ( int i = 0; i < _dim; i++ ) {
      NL[i] = AnisotropyValue * AnisotropyGradient[i];
    }
  }
};

template <typename ConfiguratorType, typename AnisotropyType>
class ATAnisotropyDiffOpForU : public aol::FENonlinDiffOpInterface< ConfiguratorType, ATAnisotropyDiffOpForU<ConfiguratorType, AnisotropyType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrV;
  const AnisotropyType &_anisotropy;
  const int _dim;
public:
  ATAnisotropyDiffOpForU ( const typename ConfiguratorType::InitType &Initializer,
                           const aol::Vector<RealType> &VDofs,
                           const AnisotropyType &Anisotropy )
      : aol::FENonlinDiffOpInterface< ConfiguratorType, ATAnisotropyDiffOpForU<ConfiguratorType, AnisotropyType> > ( Initializer ),
      _discrV ( Initializer, VDofs ),
      _anisotropy ( Anisotropy ),
      _dim ( ConfiguratorType::Dim ) {}


  inline void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                                const typename ConfiguratorType::ElementType &El,
                                int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/,
                                typename ConfiguratorType::VecType &NL ) const {
    typename ConfiguratorType::VecType gradU;
    DiscFunc.evaluateGradientAtQuadPoint ( El, QuadPoint, gradU );
    RealType AnisotropyValue = _anisotropy.evaluate ( gradU );
    typename ConfiguratorType::VecType AnisotropyGradient;
    _anisotropy.evaluateGradient ( gradU, AnisotropyGradient );

    for ( int i = 0; i < _dim; i++ ) {
      NL[i] = pow ( _discrV.evaluateAtQuadPoint ( El, QuadPoint ), 2 ) * AnisotropyValue * AnisotropyGradient[i];
    }
  }
};

//! Deprecated, use aol::SquaredDiffWeightMassOp
template <typename ConfiguratorType>
class ATMassGradU2Op : public aol::FELinScalarWeightedMassInterface<ConfiguratorType, ATMassGradU2Op<ConfiguratorType> > {
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrU;
public:
  typedef typename ConfiguratorType::RealType RealType;

  ATMassGradU2Op ( const typename ConfiguratorType::InitType &Initializer,
                   const aol::Vector<RealType> &UDofs,
                   aol::OperatorType OpType = aol::ONTHEFLY )
      : aol::FELinScalarWeightedMassInterface<ConfiguratorType, ATMassGradU2Op<ConfiguratorType> > ( Initializer, OpType ),
      _discrU ( Initializer, UDofs ) {}

  inline RealType getCoeff ( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::DomVecType& /*RefCoord*/ ) const {
    typename ConfiguratorType::VecType grad;
    _discrU.evaluateGradientAtQuadPoint ( El, QuadPoint, grad );
    return grad.normSqr();
  }
};

template <typename ConfiguratorType, typename AnisotropyType>
class ATMassAnisoU2Op : public aol::FELinScalarWeightedMassInterface<ConfiguratorType, ATMassAnisoU2Op<ConfiguratorType, AnisotropyType> > {
protected:
  const AnisotropyType &_anisotropy;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrU;
public:
  typedef typename ConfiguratorType::RealType RealType;

  ATMassAnisoU2Op ( const typename ConfiguratorType::InitType &Initializer,
                    const aol::Vector<RealType> &UDofs,
                    const AnisotropyType &Anisotropy,
                    aol::OperatorType OpType = aol::ONTHEFLY )
      : aol::FELinScalarWeightedMassInterface<ConfiguratorType, ATMassAnisoU2Op<ConfiguratorType, AnisotropyType> > ( Initializer, OpType ),
      _anisotropy ( Anisotropy ),
      _discrU ( Initializer, UDofs ) {}

  inline RealType getCoeff ( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::DomVecType& /*RefCoord*/ ) const {
    typename ConfiguratorType::VecType grad;
    _discrU.evaluateGradientAtQuadPoint ( El, QuadPoint, grad );
    return pow ( _anisotropy.evaluate ( grad ), 2 );
  }
};

/*
 * \todo Check if the calculation of phi inverse can be speed up somehow
 *       or make it possible to supply phi inverse in the constructor.
 */
template <typename ConfiguratorType>
class ATMassGradU2DeformOp : public aol::FELinScalarWeightedMassInterface<ConfiguratorType, ATMassGradU2DeformOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrU;
  aol::MultiVector<RealType> _discrPhiVector;
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType, ConfiguratorType::Dim> _discrPhiInv;
  typename qc::BitArray<ConfiguratorType::Dim> _valuesSet;
public:
  ATMassGradU2DeformOp ( const typename ConfiguratorType::InitType &Initializer,
                         const aol::Vector<RealType> &UDofs,
                         const aol::MultiVector<RealType> &PhiDofs,
                         aol::OperatorType OpType = aol::ONTHEFLY )
    : aol::FELinScalarWeightedMassInterface<ConfiguratorType, ATMassGradU2DeformOp<ConfiguratorType> > ( Initializer, OpType ),
      _discrU( Initializer, UDofs ),
      _discrPhiVector( ConfiguratorType::Dim, UDofs.size() ),
      _discrPhiInv ( Initializer, _discrPhiVector ) {
    aol::MultiVector<RealType> identity( Initializer );
    qc::DataGenerator<ConfiguratorType>( Initializer ).generateIdentity( identity );
    qc::TransformFunction<RealType, ConfiguratorType::Dim> transformOp( Initializer );
    transformOp.setDeformation( PhiDofs );
    transformOp.transform( identity, _discrPhiVector, _valuesSet );
  }

  inline RealType getCoeff ( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::DomVecType& /*RefCoord*/ ) const {
    // if $\phi^{-1}$ is not defined at this position, return 0
    if ( !_valuesSet.elementTrue( El ) )
      return 0.;

    const RealType gradSqr = _discrU.evaluateAtQuadPoint ( El, QuadPoint );
    typename ConfiguratorType::MatType mat;
    _discrPhiInv.evaluateGradientAtQuadPoint( El, QuadPoint, mat );

    if ( !aol::isFinite( gradSqr) ) {
      cerr << "gradSqr is nan!\n";
    }
    RealType detInverse = fabs ( mat.det() );
    if( detInverse > 10. ){
      detInverse = 10.;
    }
    return gradSqr * detInverse;
  }
};

//! \f$ \int_\Omega \lVert\nabla u_{R}\rVert^2(v\circ\phi) \nabla(v\circ\phi )\vartheta \f$, where \f$ u_{R} \f$ = _image
template <typename ConfiguratorType, qc::Dimension Dim>
class ATDeformationGradient
      : public aol::FENonlinVectorOpInterface<ConfiguratorType, Dim, Dim, ATDeformationGradient<ConfiguratorType, Dim> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _v;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _image;
  const qc::GridDefinition &_grid;
protected:
public:
  ATDeformationGradient ( const qc::GridDefinition &Grid,
                          const aol::Vector<RealType> &VDofs,
                          const aol::Vector<RealType> &ImageDofs )
      : aol::FENonlinVectorOpInterface<ConfiguratorType, Dim, Dim, ATDeformationGradient<ConfiguratorType, Dim> > ( Grid ),
      _v ( Grid, VDofs ),
      _image ( Grid, ImageDofs ),
      _grid ( Grid ) {}


  void getNonlinearity ( aol::auto_container<Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord,
                         aol::Vec<Dim, typename ConfiguratorType::RealType> &NL ) const {
    typename ConfiguratorType::VecType transformed_local_coord;
    qc::Element transformed_el;
    if ( !qc::transformCoord<ConfiguratorType> ( _grid, DiscFuncs, El, QuadPoint, RefCoord, transformed_el, transformed_local_coord ) ) {
      NL.setZero();
      return;
    }

    typename ConfiguratorType::VecType gradImg;
    _image.evaluateGradientAtQuadPoint ( El, QuadPoint, gradImg );
    _v.evaluateGradient ( transformed_el, transformed_local_coord, NL );

    NL *= gradImg.normSqr() * _v.evaluate ( transformed_el, transformed_local_coord );
  }
};

//! \f$ \int_\Omega c_{R}\lVert\nabla u_{R}\rVert^2((v+1-c_{R})\circ\phi) \nabla(v+1-c_{R}\circ\phi )\vartheta \f$, where \f$ u_{R} = \f$ _image
template <typename ConfiguratorType, qc::Dimension Dim>
class ATWeightedDeformationGradient
      : public aol::FENonlinVectorOpInterface<ConfiguratorType, Dim, Dim, ATWeightedDeformationGradient<ConfiguratorType, Dim> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrWeight;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _v;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _image;
  const qc::GridDefinition &_grid;
protected:
public:
  ATWeightedDeformationGradient ( const qc::GridDefinition &Grid,
                                  const aol::Vector<RealType> &WeightDofs,
                                  const aol::Vector<RealType> &VDofs,
                                  const aol::Vector<RealType> &ImageDofs )
      : aol::FENonlinVectorOpInterface<ConfiguratorType, Dim, Dim, ATWeightedDeformationGradient<ConfiguratorType, Dim> > ( Grid ),
      _discrWeight ( Grid, WeightDofs ),
      _v ( Grid, VDofs ),
      _image ( Grid, ImageDofs ),
      _grid ( Grid ) {}

  void getNonlinearity ( aol::auto_container<Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord,
                         aol::Vec<Dim, typename ConfiguratorType::RealType> &NL ) const {
    const RealType weight = _discrWeight.evaluateAtQuadPoint ( El, QuadPoint );
    if ( weight != 0. ) {
      typename ConfiguratorType::VecType transformed_local_coord;
      qc::Element transformed_el;
      if ( !qc::transformCoord<ConfiguratorType> ( _grid, DiscFuncs, El, QuadPoint, RefCoord, transformed_el, transformed_local_coord ) ) {
        NL.setZero();
        return;
      }

      typename ConfiguratorType::VecType gradImg;
      _image.evaluateGradientAtQuadPoint ( El, QuadPoint, gradImg );
      _v.evaluateGradient ( transformed_el, transformed_local_coord, NL );

      NL *= weight * gradImg.normSqr() * _v.evaluate ( transformed_el, transformed_local_coord );
    } else {
      NL.setZero();
      return;
    }
  }
};

template <typename ConfiguratorType>
class ATGradSqrOp : public aol::FENonlinOpInterface< ConfiguratorType, ATGradSqrOp<ConfiguratorType> > {
public:
  ATGradSqrOp ( const typename ConfiguratorType::InitType &Initializer )
      : aol::FENonlinOpInterface< ConfiguratorType, ATGradSqrOp<ConfiguratorType> > ( Initializer ) {}

  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/,
                         typename ConfiguratorType::RealType &NL ) const {
    typename ConfiguratorType::VecType grad;
    DiscFunc.evaluateGradientAtQuadPoint ( El, QuadPoint, grad );
    NL = grad.normSqr();
  }
};

//! Computes \f$ \alpha\int (u-u_0)^2 + \beta\int v^2|\nabla u|^2 +... \f$, where arg[0]=u, arg[1]=v.
template <typename ConfiguratorType>
class ATEnergyOp
      : public aol::FENonlinIntegrationVectorInterface<ConfiguratorType, ATEnergyOp<ConfiguratorType>, 2 > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const RealType _alpha;
  const RealType _beta;
  const RealType _nu;
  const RealType _epsilon;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrU0;
public:
  ATEnergyOp ( const typename ConfiguratorType::InitType &Initializer,
               const aol::Vector<RealType> &U0,
               const RealType Alpha,
               const RealType Beta,
               const RealType Nu,
               const RealType Epsilon )
      : aol::FENonlinIntegrationVectorInterface < ConfiguratorType,
      ATEnergyOp<ConfiguratorType>, 2 > ( Initializer ),
      _alpha ( Alpha ),
      _beta ( Beta ),
      _nu ( Nu ),
      _epsilon ( Epsilon ),
      _discrU0 ( Initializer, U0 ) {}

  RealType evaluateIntegrand ( const aol::auto_container<2, aol::DiscreteFunctionDefault<ConfiguratorType> > &discrFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    const RealType u0 = _discrU0.evaluateAtQuadPoint ( El, QuadPoint );
    const RealType u = discrFuncs[0].evaluateAtQuadPoint ( El, QuadPoint );
    const RealType v = discrFuncs[1].evaluateAtQuadPoint ( El, QuadPoint );
    RealType integrand = 0.5 * _alpha * aol::Sqr ( u - u0 );

    typename ConfiguratorType::VecType tmp;
    discrFuncs[0].evaluateGradientAtQuadPoint ( El, QuadPoint, tmp );
    integrand += 0.5 * _beta * aol::Sqr ( v ) * tmp.normSqr();

    discrFuncs[1].evaluateGradientAtQuadPoint ( El, QuadPoint, tmp );
    integrand += 0.5 * _nu * ( _epsilon * tmp.normSqr() + aol::Sqr ( v - 1 ) / ( 4. * _epsilon ) );

    return integrand;
  }
};

/**
 * \brief Computes \f$ \frac{1}{2}\int (\psi(\phi(x))-x)^2 \mathrm{d}x \f$, where arg = \f$ \phi \f$ and
 *        \f$ \psi \f$ is supplied in the constructor.
 *
 * \author Han, Berkels
 */
template <typename ConfiguratorType>
class ConsistencyEnergyOp
      : public aol::FENonlinIntegrationVectorInterface < ConfiguratorType,
      ConsistencyEnergyOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  aol::auto_container< ConfiguratorType::Dim, const aol::DiscreteFunctionDefault<ConfiguratorType> > _discrPhiInverse;
protected:
public:
  ConsistencyEnergyOp ( const typename ConfiguratorType::InitType &Grid,
                        const aol::MultiVector<RealType> &PhiInverseDofs )
      : aol::FENonlinIntegrationVectorInterface < ConfiguratorType,
      ConsistencyEnergyOp<ConfiguratorType> > ( Grid ) {
    for ( int c = 0; c < PhiInverseDofs.numComponents(); c++ ) {
      aol::DiscreteFunctionDefault<ConfiguratorType> temp ( this->getConfigurator(), PhiInverseDofs[c] );
      _discrPhiInverse.set_copy ( c, temp );
    }
  };

  RealType evaluateIntegrand ( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &discrFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    typename ConfiguratorType::VecType offset, transformed_local_coord;
    qc::Element transformed_el;
    RealType energy = 0.;

    // calculate \Abs{g(h(x))-x}^2
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      offset[i] = discrFuncs[i].evaluateAtQuadPoint ( El, QuadPoint );
    }

    qc::transformAndClipCoord<ConfiguratorType> ( this->getConfigurator(), El, RefCoord, offset, transformed_el, transformed_local_coord );

    // Since _discrPhiInverse only contains the displacement, we have to omit "-x",
    // but add the offset (comes from evaluating x at the deformed position), i.e.
    // g( h(x) ) - x = ( x + u(x) ) + v( x + u(x) ) - x
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      energy += aol::Sqr ( _discrPhiInverse[i].evaluate ( transformed_el, transformed_local_coord ) + offset[i] );
    }

    return 0.5 * energy;
  }
};

/**
 * \brief Computes \f$ \left(\int (D\psi)(\phi(x))^T(\psi(\phi(x))-x)\cdot\varphi_i(x) \mathrm{d}x\right)_i \f$, where arg = \f$ \phi \f$ and
 *        \f$ \psi \f$ is supplied in the constructor.
 *
 * \author Han, Berkels
 */
template <typename ConfiguratorType>
class VariationOfConsistencyEnergyOp:
public aol::FENonlinVectorOpInterface < ConfiguratorType,
                                        ConfiguratorType::Dim,
                                        ConfiguratorType::Dim,
                                        VariationOfConsistencyEnergyOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType, ConfiguratorType::Dim> _known;
public:
  VariationOfConsistencyEnergyOp ( const typename ConfiguratorType::InitType &Grid,
                                   const aol::MultiVector<RealType> &KnownDofs )
  : aol::FENonlinVectorOpInterface<ConfiguratorType,
    ConfiguratorType::Dim,
    ConfiguratorType::Dim,
    VariationOfConsistencyEnergyOp<ConfiguratorType> > ( Grid ),
    _known( Grid, KnownDofs ) {
  }

  void getNonlinearity ( aol::auto_container<ConfiguratorType::Dim,aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint,
                         const typename ConfiguratorType::VecType &RefCoord,
                         aol::Vec<ConfiguratorType::Dim,typename ConfiguratorType::RealType> &NL ) const {
    typename ConfiguratorType::VecType offset, transformed_local_coord;
    qc::Element transformed_el;

    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      offset[i] = DiscFuncs[i].evaluateAtQuadPoint( El, QuadPoint );
    }

    aol::Vec<ConfiguratorType::Dim, bool> coordinateWithinLimits;
    qc::transformAndClipCoord<ConfiguratorType> ( this->getConfigurator(), El, RefCoord, offset, transformed_el, transformed_local_coord, coordinateWithinLimits );

    typename ConfiguratorType::MatType Jacobi_known;
    typename ConfiguratorType::VecType VecKnown;

    //compute [g\circ h](x)-x
    // Since _known only contains the displacement, we have to omit "-x",
    // but add the offset (comes from evaluating x at the deformed position), i.e.
    // g( h(x) ) - x = ( x + u(x) ) + v( x + u(x) ) - x
    _known.evaluate( transformed_el, transformed_local_coord, VecKnown);
    VecKnown += offset;

    //compute D[g](h(x))
    _known.evaluateGradient(transformed_el, transformed_local_coord, Jacobi_known );
    for ( int i = 0; i < ConfiguratorType::Dim; i ++ ) {
      Jacobi_known[i][i] += 1.;
    }
      
    //NL = (g((h(x))-x)D[g](h(x))
    Jacobi_known.multTransposed ( VecKnown, NL );

    for ( int i = 0; i < ConfiguratorType::Dim; ++i ) {
      if ( coordinateWithinLimits[i] == false )
        NL[i] = 0;
    }
  }
};
  
/**
 * \brief Computes \f$ \left(\int (\phi(\psi(x))-x)\cdot\varphi_i(\psi(x)) \mathrm{d}x\right)_i \f$, where arg = \f$ \phi \f$ and
 *        \f$ \psi \f$ is supplied in the constructor. If a weight is specified in the constructor, the integrand is multiplied
 *        by \f$ w^2(x) \f$.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class VariationOfConsistencyEnergyOpWRTOuterDef : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_grid;
  aol::DeleteFlagPointer<aol::Op<aol::Vector<RealType> > > _pDefMassOp;
  aol::DiagonalBlockOp<RealType> _defMassBlockOp;
  aol::MultiVector<RealType> _rightDefMassPsi;
public:
  VariationOfConsistencyEnergyOpWRTOuterDef ( const typename ConfiguratorType::InitType &Grid,
                                              const aol::MultiVector<RealType> &PsiDofs )
    : _grid ( Grid ),
      _pDefMassOp ( new qc::DeformMassOp<ConfiguratorType> ( _grid, PsiDofs, aol::ONTHEFLY, qc::BOTH ), true ),
      _defMassBlockOp ( *_pDefMassOp ),
      _rightDefMassPsi ( _grid ) {
    qc::DeformMassOp<ConfiguratorType> rightDefMassOp ( _grid, PsiDofs, aol::ONTHEFLY, qc::RIGHT );
    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      rightDefMassOp.apply ( PsiDofs[i], _rightDefMassPsi[i] );
  }

  VariationOfConsistencyEnergyOpWRTOuterDef ( const typename ConfiguratorType::InitType &Grid,
                                              const aol::MultiVector<RealType> &PsiDofs,
                                              const aol::Vector<RealType> &Weight )
    : _grid ( Grid ),
      _pDefMassOp ( new qc::SquaredWeightDeformMassOp<ConfiguratorType> ( _grid, PsiDofs, Weight, aol::ONTHEFLY, qc::BOTH ), true ),
      _defMassBlockOp ( *_pDefMassOp ),
      _rightDefMassPsi ( _grid ) {
    qc::SquaredWeightDeformMassOp<ConfiguratorType> rightDefMassOp ( _grid, PsiDofs, Weight, aol::ONTHEFLY, qc::RIGHT );
    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      rightDefMassOp.apply ( PsiDofs[i], _rightDefMassPsi[i] );
  }

  virtual void applyAdd( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    _defMassBlockOp.applyAdd ( MArg, MDest );
    // Since _psi only contains the displacement, we have to omit "-x",
    // but take care of the displacement part of psi
    MDest += _rightDefMassPsi;
  }
};

/**
 * Joint edge detection and denoising with the Ambrosio-Tortorelli model.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class ATSegmentation {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

protected:
  //! controls the phase field transition width
  RealType _epsilon;
  //! weight of the phase field length term
  const RealType _alpha;
  //! weight of the fidelity term, i.e. \f$ \frac{1}{2}\int_\Omega (u-u_0)^2 dx \f$
  const RealType _beta;

  const typename ConfiguratorType::InitType _grid;

  //! Input image
  const ArrayType _u0;

  //! Temporary storage vector.
  mutable ArrayType _rhs;

  const aol::MassOp<ConfiguratorType> _mass;
  const aol::StiffOp<ConfiguratorType> _stiff;

  mutable typename ConfiguratorType::MatrixType _mat;
  aol::SSORPreconditioner<aol::Vector<RealType>, typename ConfiguratorType::MatrixType> _ssor;
  aol::PCGInverse<aol::Vector<RealType> > _solver;

  const bool _cacheMassAndStiffMat;

public:
  ATSegmentation ( const RealType Epsilon,
                   const RealType Alpha,
                   const RealType Beta,
                   const qc::Array<RealType> &U0,
                   const bool CacheMassAndStiffMat = true )
  : _epsilon(Epsilon),
    _alpha(Alpha),
    _beta(Beta),
    _grid( qc::GridSize<ConfiguratorType::Dim>::createFrom ( U0 ) ),
    _u0( U0, aol::FLAT_COPY ),
    _rhs( _grid ),
    _mass (_grid ),
    _stiff (_grid),
    _mat( _grid ),
    _ssor( _mat ),
    _solver( _mat, _ssor, 1e-10, 1000 ),
    _cacheMassAndStiffMat ( CacheMassAndStiffMat )
  {
    _solver.setQuietMode( false );
    _solver.setStopping( aol::STOPPING_ABSOLUTE );
  }

  virtual ~ATSegmentation () {}

  //! Update the image while keeping the phase field fixed.
  virtual void updateU ( aol::Vector<RealType> &U,
                         const aol::Vector<RealType> &V,
                         const int InnerIterations = 0,
                         const int /*OuterIterations*/ = 0 ) const {
    cerr << "Solving for u at iteration " << InnerIterations << endl;
    cerr << "Assembling matrices ";
    // \beta \int (u \phi) +  \int (v^2 \nabla u \cdot \nabla \phi) = \beta \int (u_0 \phi)
    qc::ATStiffV2Op<ConfiguratorType> atStiffV2( _grid, V );

    _mat.setZero();

    // mat = beta * M
    if ( _cacheMassAndStiffMat )
      _mat.addMultiple ( _mass.getMatrix(), _beta );
    else {
      _mass.assembleAddMatrix( _mat );
      _mat *= ( _beta );
    }

    atStiffV2.assembleAddMatrix( _mat ); // mat += L[v^2]
    cerr << "done\n";
    // rhs could be cached to get a minor speed increase
    _mass.apply( _u0, _rhs ); //rhs = beta * M u_0
    _rhs *= _beta;

    _solver.apply( _rhs, U );
  }

  //! Update the phase field while keeping the image fixed.
  void updateV( const aol::Vector<RealType> &U,
                aol::Vector<RealType> &V,
                const int InnerIterations = 0,
                const int /*OuterIterations*/ = 0 ) const {
    // \int (\|\nabla u \|^2 v \phi) + \alpha \epsilon \int (\nabla v \cdot \nabla \phi)
    // + \frac{\alpha}{4 \epsilon} \int ( v \phi) = \frac{\alpha}{4 \epsilon} \int \phi
    cerr << "Solving for v at iteration " << InnerIterations << endl;
    cerr << "Assembling matrices ";
    qc::ATMassGradU2Op<ConfiguratorType> atMassGradU2( _grid, U );

    _mat.setZero();

    //mat = alpha*epsilon*L
    if ( _cacheMassAndStiffMat )
      _mat.addMultiple ( _stiff.getMatrix(), _alpha * _epsilon );
    else {
      _stiff.assembleAddMatrix( _mat );
      _mat *= ( _alpha * _epsilon );
    }

    atMassGradU2.assembleAddMatrix( _mat ); //mat += M[GradU2]

    _mat *= ( (4. * _epsilon) / _alpha); // mat += (alpha/(4*epsilon))*M
    _mass.assembleAddMatrix( _mat );
    _mat *= ( _alpha / ( 4. * _epsilon) );
    cerr << "done\n";

    ArrayType tmp( _grid );
    tmp.setAll( _alpha / ( _epsilon * 4. ) ); // rhs = M tmp
    _mass.apply( tmp, _rhs );

    _solver.apply( _rhs, V );
  }

  void segment( ArrayType &U, ArrayType &V, const int MaxIterations = 10 ) const {
    aol::Vector<RealType> diffU ( U );
    aol::Vector<RealType> diffV ( V );

    for ( int iter = 0; iter < MaxIterations; ++iter ) {
      updateU ( U, V, iter );
      updateV ( U, V, iter );
      diffU -= U;
      diffV -= V;
      if ( aol::appeqAbsolute ( diffU.norm() + diffV.norm(), aol::ZTrait<RealType>::zero ) )
        break;
      diffU = U;
      diffV = V;
    }
  }
};

} // namespace qc

#endif // __AMBROSIOTORTORELLI_H
