#ifndef __ANISOTROPIES_H
#define __ANISOTROPIES_H

#include <aol.h>
#include <polarCoords.h>
#include <scalarArray.h>

namespace qc {

//! (Hopefully) efficient graph interface for 3d-anisotropies to be used in the 2d-graph-case.
/*! \author Nemitz */
template <typename RealType, typename Imp>
class Anisotropy3dGraphInterface {
public:
  Anisotropy3dGraphInterface(  )  { _graphFlag = false; }

  //! call this method if you want to use the graph-methods.
  void setGraphFlag() { _graphFlag = true; }
  bool isGraphCase() const { return( _graphFlag ); }

  //! This method is for the graph-case. Here the normal is (-grad u, 1).
  //! Notice: To use this method you have to set the graph-flag before!!
  RealType gammaNorm ( const aol::Vec2<RealType> &z ) const {
    if ( !_graphFlag )
      throw aol::Exception ( "\nYou have called the gammaNorm-Method with a 2d-vector.\nSet the graph-flag before!\n", __FILE__, __LINE__ );
    aol::Vec3<RealType> z3 ( z[0], z[1], -1 );
    return this->asImp().gammaNorm ( z3 );
  }

  //! This method is for the graph-case. Here the normal is (-grad u, 1).
  //! Notice: To use this method you have to set the graph-flag before!!
  void gammaFirstDerivative ( const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    if ( !_graphFlag ) {
      throw aol::Exception ( "\nYou have called the gammaFirstDerivative-Method with a 2d-vector.\nSet the graph-flag before!\n", __FILE__, __LINE__ );
    }
    aol::Vec3<RealType> z3 ( z[0], z[1], -1 );
    aol::Vec3<RealType> v3;
    this->asImp().gammaFirstDerivative ( z3, v3 );
    v[0] =  v3[0];     // copy the first two arguments of v3
    v[1] =  v3[1];
  }

  void gammaFirstDerivative ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    gammaFirstDerivative( z, v );
  }
  void gammaFirstDerivative ( const qc::Element &, const int &, const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    gammaFirstDerivative( z, v );
  }


  void implicitPart ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z,
                      aol::Matrix22<RealType> &mat ) const {
    implicitPart ( z, mat );
  }
  void implicitPart ( const qc::Element &, const int &, const aol::Vec2<RealType> &z,
                      aol::Matrix22<RealType> &mat ) const {
                      implicitPart ( z, mat );
  }

  void explicitPart ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z,
                      aol::Matrix22<RealType> &mat ) const {
    explicitPart ( z, mat );
  }

  //! second derivative, implicit part of the p-regularized willmore flow (graph case)
  void implicitPart ( const aol::Vec2<RealType> &z, aol::Matrix22<RealType> &mat ) const {
    aol::Vec3<RealType> z3 ( z[0], z[1], -1 );
    aol::Matrix33<RealType> mat33;
    this->asImp().implicitPart( z3, mat33 );
    mat[0][0] = mat33[0][0];
    mat[0][1] = mat33[0][1];
    mat[1][0] = mat33[1][0];
    mat[1][1] = mat33[1][1];
  }

  //! second derivative, explicit part of the p-regularized willmore flow (graph case)
  void explicitPart ( const aol::Vec2<RealType> &z, aol::Matrix22<RealType> &mat ) const {
    aol::Vec3<RealType> z3 ( z[0], z[1], -1 );
    aol::Matrix33<RealType> mat33;
    this->asImp().implicitPart( z3, mat33 );
    mat[0][0] = mat33[0][0];
    mat[0][1] = mat33[0][1];
    mat[1][0] = mat33[1][0];
    mat[1][1] = mat33[1][1];
  }



protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

  bool  _graphFlag;             //! set this to use the gamma-methods with a 2d-vector-argument (graph-case)
};



//! General class, that applies a linear transformation (2x2-matrix) to a given anisotropy.
template <typename RealType, typename AnisoType>
class LinearTransformedAnisotropy2d {
private:
  AnisoType &_anisotropy;
  aol::Matrix22<RealType> _mat;              //! the linear transformation matrix
  aol::Matrix22<RealType> _matTransposed;

public:
  LinearTransformedAnisotropy2d ( AnisoType &Anisotropy, aol::Matrix22<RealType> Mat ) :
      _anisotropy ( Anisotropy ), _mat ( Mat ), _matTransposed ( Mat ) { _matTransposed.transpose(); }


  //! Evaluating the function itself, i.e. the regularization of
  //! \f$ \gamma(z) = \sqrt{z_1^2 + \epsilon^2} + \sqrt{z_1^2 + \epsilon^2} \f$.
  RealType gammaNorm ( const aol::Vec2<RealType> &z ) const {
    aol::Vec2<RealType> v;
    _mat.mult ( z, v );
    return ( _anisotropy.gammaNorm ( v ) );
  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    aol::Vec2<RealType> temp1 ( 0., 0. );
    aol::Vec2<RealType> temp2 ( 0., 0. );
    _mat.mult ( z, temp1 );
    _anisotropy.gammaFirstDerivative ( temp1, temp2 );
    _matTransposed.mult ( temp2, v );
//     _anisotropy.gammaFirstDerivative ( temp1, v );      // HACK
  }
};



//! Class, that generates a 3d-anisotropy by rotating a 2d-anisotropy
//! around the z-axis. For the derivatives, the 2d-anisotropy has to
//! provide partial derivatives wrt r and phi in the methods
//! partialDerivativesWRTPhiR( const z, pd& )
template <typename RealType, typename AnisoType>
class Rotated3dAnisotropy : public Anisotropy3dGraphInterface<RealType, Rotated3dAnisotropy<RealType, AnisoType> > {
public:
  using Anisotropy3dGraphInterface<RealType, Rotated3dAnisotropy<RealType, AnisoType> >::gammaNorm;
  using Anisotropy3dGraphInterface<RealType, Rotated3dAnisotropy<RealType, AnisoType> >::gammaFirstDerivative;

private:
  const AnisoType &_anisotropy;

  const RealType _delta;     //! Regularization-parameter for the radius
  RealType _deltaSqr;        //! uh, don't know ...

public:
  Rotated3dAnisotropy ( const AnisoType &Anisotropy, const RealType delta ) :
      Anisotropy3dGraphInterface<RealType, Rotated3dAnisotropy<RealType, AnisoType> > ( ),
      _anisotropy ( Anisotropy ),
      _delta ( delta ) { _deltaSqr = _delta * _delta; }

  //! Evaluating the function itself
  RealType gammaNorm ( const aol::Vec3<RealType> &z ) const {
    aol::Vec2<RealType> v ( sqrt ( z[0]*z[0] + z[1]*z[1] + _deltaSqr ), z[2] );
    return ( _anisotropy.gammaNorm ( v ) );
  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec3<RealType> &, const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v )
  const {    gammaFirstDerivative ( z, v );  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    RealType r   = sqrt ( z[0]*z[0] + z[1]*z[1] + _deltaSqr );

    // compute the partial derivatives of the 2d anisotropy
    aol::Vec2<RealType> z2Arg( r, z[2] );
    aol::Vec2<RealType> pdGamma;

    _anisotropy.gammaFirstDerivative ( z2Arg, pdGamma );

    RealType drdx = z[0] / r;       //! partial derivatives of r(x,y,z)
    RealType drdy = z[1] / r;

    v[0] = pdGamma[0] * drdx;
    v[1] = pdGamma[0] * drdy;
    v[2] = pdGamma[1];
  }
};


//! Class for blending two different anisotropies dependent on the image value, e.g. value>0.5 => use
//! anisotropy1 else use anisotropy2 and use a linear blending in an interval between.
template <typename ConfiguratorType, typename AnisoType1, typename AnisoType2>
class LinearBlendedAnisotropy2dGraph {

public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;

private:
  AnisoType1 &_anisotropy1;       // for image values < _xMin
  AnisoType2 &_anisotropy2;       // for image values > _xMax
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discFunc;
  ConfiguratorType _config;
  qc::GridDefinition _grid;
  RealType _xMin;
  RealType _xMax;
  RealType _delta;
  RealType _deltaSqr;
  mutable qc::ScalarArray<RealType, qc::QC_2D> _testImg;

public:
  LinearBlendedAnisotropy2dGraph ( AnisoType1 &Anisotropy1, AnisoType2 &Anisotropy2,
                              const typename ConfiguratorType::InitType &Grid, RealType XMin, RealType XMax, RealType delta ) :
      _anisotropy1 ( Anisotropy1 ),
      _anisotropy2 ( Anisotropy2 ),
      _discFunc( NULL ),
      _config( Grid ), _grid( Grid ),
      _xMin( XMin ), _xMax( XMax ),
      _delta ( delta ), _deltaSqr ( delta*delta ),
      _testImg( Grid ) { _testImg.setAll( -1. );
  }


  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( _discFunc )
      delete _discFunc;
    _discFunc = new aol::DiscreteFunctionDefault<ConfiguratorType>( this->_config.getInitializer(), Image );
  }

  void saveImage() const {
    _testImg.save( "BlendingTest.bz2", 9 );
  }

  //! this is the blending function.
  RealType eta( const RealType x ) const {
    if ( x < _xMin ) return 1.;
    if ( x > _xMax ) return 0.;
    return ( _xMax - x ) / ( _xMax - _xMin );
  }


  //! Evaluating the function itself, i.e. the blending
  //! \f$ \gamma(z) = \eta \gamma_1(z) + (1-\eta)\gamma_2(z) \f$.
  RealType gammaNorm ( const ElementType &El, const int QuadPoint, const aol::Vec2<RealType> &z ) const {
    if (!_discFunc) throw aol::Exception ( "\nSet Image first!\n", __FILE__, __LINE__ );
    aol::Vec3<RealType> z3 ( z[0], z[1], -1. );
    aol::Vec3<RealType> minusZ3( -z[0], -z[1], 1. );

    // get the weight of the blending function
    RealType x = _discFunc->evaluateAtQuadPoint( El, QuadPoint );
    RealType etaX = eta(x);

    // here is the linear blending
    if ( etaX == 0. ) return _anisotropy2.gammaNorm( z3 );
    else if ( etaX == 1. ) return _anisotropy1.gammaNorm( minusZ3 );
    else return ( 1. - etaX ) * _anisotropy2.gammaNorm( z3 ) + etaX * _anisotropy1.gammaNorm( minusZ3 );
  }


  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const ElementType &El, const int QuadPoint,
                              const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    aol::Vec3<RealType> z3 ( z[0], z[1], -1. );
    aol::Vec3<RealType> minusZ3( -z[0], -z[1], 1. );

    // get the weight of the blending function
    RealType x = _discFunc->evaluateAtQuadPoint( El, QuadPoint );
    RealType etaX = eta(x);

    aol::Vec3<RealType> resultAniso1, resultAniso2, v3;
    if ( etaX == 0. ) {
      _anisotropy2.gammaFirstDerivative( z3, resultAniso2 );
      resultAniso1.setZero();
    } else if ( etaX == 1. ) {
      _anisotropy1.gammaFirstDerivative( minusZ3, resultAniso1 );
      resultAniso2.setZero();
    } else {
      _anisotropy1.gammaFirstDerivative( minusZ3, resultAniso1 );
      _anisotropy2.gammaFirstDerivative( z3, resultAniso2 );
    }

    // here is the linear blending
    v3 = ( 1. - eta(x) ) * resultAniso2 - eta(x) * resultAniso1;

    v[0] = v3[0];     // copy the first two arguments of v3 (negative sign from the inner derivative)
    v[1] = v3[1];
  }


  //! this method computes \f$ gamma_zz \f$, it's called implicitPart just for technical reasons.
  void implicitPart ( const ElementType &El, const int QuadPoint,
                      const aol::Vec2<RealType> &z, aol::Matrix22<RealType> &mat ) const {
    aol::Vec3<RealType> z3 ( z[0], z[1], -1. );
    aol::Vec3<RealType> minusZ3( -z[0], -z[1], 1. );

    // get the weight of the blending function
    RealType x = _discFunc->evaluateAtQuadPoint( El, QuadPoint );
    RealType etaX = eta(x);

    aol::Matrix33<RealType> resultAniso1, resultAniso2;
    if ( etaX == 0. ) {
      _anisotropy2.implicitPart( z3, resultAniso2 );
      resultAniso1.setZero();
    } else if ( etaX == 1. ) {
      _anisotropy1.implicitPart( minusZ3, resultAniso1 );
      resultAniso2.setZero();
    } else {
      _anisotropy1.implicitPart( minusZ3, resultAniso1 );
      _anisotropy2.implicitPart( z3, resultAniso2 );
    }

    resultAniso2 *= ( 1. - eta(x) );                        // weight the two anisotropies
    resultAniso1 *= eta(x);

    mat[0][0] = resultAniso1[0][0] + resultAniso2[0][0];    // linear blending
    mat[0][1] = resultAniso1[0][1] + resultAniso2[0][1];
    mat[1][0] = resultAniso1[1][0] + resultAniso2[1][0];
    mat[1][1] = resultAniso1[1][1] + resultAniso2[1][1];
  }

};

//! Class for blending two different anisotropies dependent on the image value, e.g. value>0.5 => use
//! anisotropy1 else use anisotropy2 and use a linear blending in an interval between.
template <typename ConfiguratorType, typename AnisoType1, typename AnisoType2>
class LinearBlended3dAnisotropy {

public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;

private:
  AnisoType1 &_anisotropy1;       // for image values < _xMin
  AnisoType2 &_anisotropy2;       // for image values > _xMax
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discFunc;    // either use an image for local blending values
  RealType _eta;                                                      // or one fixed value
  bool _blendingValueIsSet;                                           // flag whether user has defined at least one
  ConfiguratorType _config;
  qc::GridDefinition _grid;
  RealType _xMin;
  RealType _xMax;
  RealType _delta;
  RealType _deltaSqr;
  mutable qc::ScalarArray<RealType, qc::QC_3D> _testImg;

public:
  LinearBlended3dAnisotropy ( AnisoType1 &Anisotropy1, AnisoType2 &Anisotropy2,
                              const typename ConfiguratorType::InitType &Grid, RealType XMin, RealType XMax, RealType delta ) :
      _anisotropy1 ( Anisotropy1 ),
      _anisotropy2 ( Anisotropy2 ),
      _discFunc( NULL ),
      _config( Grid ), _grid( Grid ),
      _xMin( XMin ), _xMax( XMax ),
      _delta ( delta ), _deltaSqr ( delta*delta ),
      _testImg( Grid ) { _testImg.setAll( -1. );
  }


  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( _discFunc )
      delete _discFunc;
    _discFunc = new aol::DiscreteFunctionDefault<ConfiguratorType>( this->_config.getInitializer(), Image );
    this->reset();
  }

  void setBlendingValue ( const RealType Eta ) {
    _eta = Eta;
    _blendingValueIsSet = true;
  }

  void saveImage() const {
    _testImg.save( "BlendingTest.bz2", 9 );
  }

  //! this is the blending function.
  RealType eta( const RealType x ) const {
    if ( x < _xMin ) return 1.;
    if ( x > _xMax ) return 0.;
    return ( _xMax - x ) / ( _xMax - _xMin );
  }

#if 0
  // this method cannot compile (const method trying to set non-mutable member)
  //! Evaluating the function itself, i.e. the blending
  //! \f$ \gamma(z) = \eta \gamma_1(z) + (1-\eta)\gamma_2(z) \f$.
  RealType gammaNorm ( const ElementType &El, const int QuadPoint, const aol::Vec3<RealType> &z3 ) const {
    if (!_discFunc) throw aol::Exception ( "\nLinearBlended3dAnisotropy::Set Image first!\n", __FILE__, __LINE__ );

    // get the weight of the blending function
    RealType x = _discFunc->evaluateAtQuadPoint( El, QuadPoint );
    _eta = eta(x);

    // this flag is here, because otherwise it might happen, that a blending-image is defined, but
    // only operators that call the gammaNorm-function without Element and QuadPoint arguments
    // are used. Then _eta would be undefined.
    _blendingValueIsSet = true;

    // call the linear blending
    gammaNorm( z3 );
  }
#endif

  RealType gammaNorm ( const aol::Vec3<RealType> &z3 ) const {
    if ( !_blendingValueIsSet ) throw aol::Exception ( "\nLinearBlended3dAnisotropy: Set Blending-value first!\n", __FILE__, __LINE__ );
    aol::Vec3<RealType> minusZ3( -z3[0], -z3[1], -z3[2] );

    // here is the linear blending
    if ( _eta == 0. ) return _anisotropy2.gammaNorm( z3 );
    else if ( _eta == 1. ) return _anisotropy1.gammaNorm( minusZ3 );
    else return ( 1. - _eta ) * _anisotropy2.gammaNorm( z3 ) + _eta * _anisotropy1.gammaNorm( minusZ3 );
  }


#if 0
  // this method cannot compile (const method trying to set non-mutable member)
  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const ElementType &El, const int QuadPoint,
                              const aol::Vec3<RealType> &z3, aol::Vec3<RealType> &v ) const {
    // get the weight of the blending function
    RealType x = _discFunc->evaluateAtQuadPoint( El, QuadPoint );
    _eta = eta(x);

    _blendingValueIsSet = true;

    gammaFirstDerivative( z3, v );
  }
#endif

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec3<RealType> &z3, aol::Vec3<RealType> &v ) const {
    if ( !_blendingValueIsSet ) throw aol::Exception ( "\nLinearBlended3dAnisotropy: Set Blending-value first!\n", __FILE__, __LINE__ );
    aol::Vec3<RealType> minusZ3( -z3[0], -z3[1], -z3[2] );

    aol::Vec3<RealType> resultAniso1, resultAniso2;
    if ( _eta == 0. ) {
      _anisotropy2.gammaFirstDerivative( z3, resultAniso2 );
      resultAniso1.setZero();
    } else if ( _eta == 1. ) {
      _anisotropy1.gammaFirstDerivative( minusZ3, resultAniso1 );
      resultAniso2.setZero();
    } else {
      _anisotropy1.gammaFirstDerivative( minusZ3, resultAniso1 );
      _anisotropy2.gammaFirstDerivative( z3, resultAniso2 );
    }

    // here is the linear blending
    v = ( 1. - _eta ) * resultAniso2 - _eta * resultAniso1;
  }

};


//! Class for the simple L1-norm-anisotropy in 2d.
template <typename RealType>
class L1Norm2d {
private:
  RealType _epsilon;
  RealType _epsSqr;

public:
  L1Norm2d ( RealType epsilon ) : _epsilon ( epsilon ), _epsSqr ( epsilon*epsilon ) { }


  //! Evaluating the function itself, i.e. the regularization of
  //! \f$ \gamma(z) = \sqrt{z_0^2 + \epsilon^2} + \sqrt{z_1^2 + \epsilon^2} \f$.
  RealType gammaNorm ( const aol::Vec2<RealType> &z ) const {
    return ( sqrt ( z[0]*z[0] + _epsSqr ) + sqrt ( z[1]*z[1] + _epsSqr ) );
  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    RealType norm = gammaNorm ( z );
    v[0] = z[0] / norm;
    v[1] = z[1] / norm;
  }
};


//! Class for the simple parallelogram-anisotropy in 2d.
template <typename RealType>
class Parallelogram2d {
private:
  RealType _epsilon;
  RealType _epsSqr;
  RealType _alpha;        // angle of x-axis with one side
  RealType _beta;         // angle of x-axis with the other side
  aol::Vec2<RealType> p;  // vector to one corner
  aol::Vec2<RealType> q;  // vector to another corner (not -p)

public:
  Parallelogram2d ( RealType Epsilon, RealType Alpha, RealType Beta )
    : _epsilon ( Epsilon ),
      _epsSqr ( Epsilon*Epsilon ),
      _alpha ( Alpha ),
      _beta ( Beta )
  {
    p[0] = cos(_alpha)*cos(_beta) + sin(_alpha);    // the vector p
    p[1] = cos(_alpha)*sin(_beta);
    q[0] = cos(_alpha) - sin(_alpha)*cos(_beta);    // the vector q
    q[1] = - sin(_alpha)*sin(_beta);
  }

  //! Evaluating the function itself
  RealType gammaNorm ( const aol::Vec2<RealType> &z ) const {

    const RealType temp1 =  aol::Sqr( p*z );    // p*z
    const RealType temp2 =  aol::Sqr( q*z );    // q*z

    return sqrt(temp1+_epsSqr) + sqrt(temp2+_epsSqr);
  }

};


//! Class for the simple L1-norm-anisotropy in 2d.
template <typename RealType>
class LInfNorm2d {
private:
  RealType _epsilon;
  RealType _epsSqr;

public:
  LInfNorm2d ( RealType epsilon ) : _epsilon ( epsilon ), _epsSqr ( epsilon*epsilon ) { }


  //! Evaluating the function itself, i.e. the regularization of
  //! \f$ \gamma(z) = \sqrt{z_0^2 + \epsilon^2} + \sqrt{z_1^2 + \epsilon^2} \f$.
  RealType gammaNorm ( const aol::Vec2<RealType> &z ) const {
    return ( 0.5 * ( sqrt ( aol::Sqr( z[0] - z[1] ) + _epsSqr ) + sqrt ( aol::Sqr( z[0] + z[1] ) + _epsSqr ) ) );
  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    v[0] = 0.5 * ( aol::signum1at0( z[0] - z[1] ) + aol::signum1at0( z[0] + z[1] ) );
    v[1] = 0.5 * ( aol::signum1at0( z[0] + z[1] ) - aol::signum1at0( z[0] - z[1] ) );
  }
};

//! Class for the simple LInf-norm-anisotropy in 3d:
//! \f$ \gamma(z) = \max \{ |z_1|,|z_2|,|z_3| \} \f$. While the norm itself is simple,
//! the regularized norm is really awful and is therefore not given in TeX here!
template <typename RealType>
class LInfNorm3d : public Anisotropy3dGraphInterface<RealType, LInfNorm3d<RealType> > {
public:
  using Anisotropy3dGraphInterface<RealType, LInfNorm3d<RealType> >::gammaNorm;
  using Anisotropy3dGraphInterface<RealType, LInfNorm3d<RealType> >::gammaFirstDerivative;
  using Anisotropy3dGraphInterface<RealType, LInfNorm3d<RealType> >::implicitPart;

private:
  RealType _delta;
  RealType _deltaSqr;

public:
  LInfNorm3d ( RealType delta )  :
    Anisotropy3dGraphInterface<RealType, LInfNorm3d<RealType> > ( ), _delta ( delta ), _deltaSqr ( delta*delta ) { }

  //! Evaluating the function itself, i.e. the regularization of
  //! \f$ \gamma(z) = \max \{ |z_1|,|z_2|,|z_3| \} \f$. Code is generated by Maple and
  //! thus not very readable, but should be correct (otherwise blame Mapele ;-)
  RealType gammaNorm ( const aol::Vec3<RealType> &z ) const {
      RealType t1 = z[0]*z[0];            RealType t3 = 2.0*z[0]*z[1];
      RealType t4 = z[1]*z[1];            RealType t5 = _deltaSqr;
      RealType t7 = sqrt(t1+t3+t4+t5);    RealType t9 = sqrt(t1-t3+t4+t5);
      RealType t10 = 2.0*z[0];            RealType t12 = aol::Sqr( t7+t9+t10 );
      RealType t14 = sqrt(t12+t5);        RealType t17 = aol::Sqr( t7+t9-t10 );
      RealType t19 = sqrt(t17+t5);

      return ( 0.25*t14+0.25*t19 );
  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec3<RealType> &, const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
      RealType t1 = z[0]*z[0];             RealType t3 = 2.0*z[0]*z[1];
      RealType t4 = z[1]*z[1];             RealType t5 = _deltaSqr;
      RealType t7 = sqrt(t1+t3+t4+t5);     RealType t9 = sqrt(t1-t3+t4+t5);
      RealType t10 = 2.0*z[2];             RealType t11 = t7+t9+t10;
      RealType t12 = t11*t11;              RealType t14 = sqrt(t12+t5);
      RealType t15 = 1/t14;                RealType t16 = t15*t11;
      RealType t19 = 2.0/t7*(z[0]+z[1]);   RealType t20 = 1/t9;
      RealType t21 = z[0]-z[1];            RealType t23 = t19+2.0*t20*t21;
      RealType t26 = t7+t9-t10;            RealType t27 = t26*t26;
      RealType t29 = sqrt(t27+t5);         RealType t30 = 1/t29;
      RealType t31 = t30*t26;              RealType t36 = t19-2.0*t20*t21;
      RealType t42 = 4.0*t7;               RealType t43 = 4.0*t9;
      RealType t44 = 8.0*z[2];

      v[0] = 0.125*t16*t23+0.125*t31*t23;
      v[1] = 0.125*t16*t36+0.125*t31*t36;
      v[2] = 0.125*t15*(t42+t43+t44)+0.125*t30*(-t42-t43+t44);
  }


  //! this method computes \f$ gamma_zz \f$, it's called implicitPart just for technical reasons.
  void implicitPart ( const aol::Vec3<RealType> &z, aol::Matrix33<RealType> &mat ) const {
      RealType t1 = z[0]*z[0];             RealType t3 = 2.0*z[0]*z[1];
      RealType t4 = z[1]*z[1];             RealType t5 = _deltaSqr;
      RealType t6 = t1+t3+t4+t5;           RealType t7 = sqrt(t6);
      RealType t8 = t1-t3+t4+t5;           RealType t9 = sqrt(t8);
      RealType t10 = 2.0*z[2];             RealType t11 = t7+t9+t10;
      RealType t12 = t11*t11;              RealType t13 = t12+t5;
      RealType t14 = sqrt(t13);            RealType t16 = 1/t14/t13;
      RealType t17 = t16*t12;              RealType t18 = 1/t7;
      RealType t19 = z[0]+z[1];            RealType t20 = 2.0*t18*t19;
      RealType t21 = 1/t9;                 RealType t22 = z[0]-z[1];
      RealType t24 = t20+2.0*t21*t22;      RealType t25 = t24*t24/4.0;
      RealType t28 = 1/t14;                RealType t31 = t28*t11;
      RealType t36 = 1/t7/t6*t19*t19;      RealType t38 = 1/t9/t8;
      RealType t42 = -t36+t18-t38*t22*t22+t21;
      RealType t44 = 0.25*t31*t42;         RealType t45 = t7+t9-t10;
      RealType t46 = t45*t45;              RealType t47 = t46+t5;
      RealType t48 = sqrt(t47);            RealType t50 = 1/t48/t47;
      RealType t51 = t50*t46;              RealType t54 = 1/t48;
      RealType t57 = t54*t45;              RealType t59 = 0.25*t57*t42;
      RealType t62 = t20-2.0*t21*t22;      RealType t63 = t24*t62/4.0;
      RealType t66 = t28*t62/2.0;          RealType t72 = -t36+t18+t38*t22*t22-t21;
      RealType t77 = t54*t62/2.0;
      RealType t82 = -0.25*t17*t63+0.125*t66*t24+0.25*t31*t72-0.25*t51*t63+0.125*t77*t24+0.25*t57*t72;
      RealType t83 = t16*t11;              RealType t84 = 4.0*t7;
      RealType t85 = 4.0*t9;               RealType t86 = 8.0*z[2];
      RealType t87 = t84+t85+t86;          RealType t90 = 0.625E-1*t83*t24*t87;
      RealType t93 = t50*t45;              RealType t94 = -t84-t85+t86;
      RealType t97 = 0.625E-1*t93*t24*t94; RealType t101 = t62*t62/4.0;
      RealType t113 = 0.625E-1*t83*t62*t87;
      RealType t117 = 0.625E-1*t93*t62*t94;
      RealType t130 = t87*t87;             RealType t134 = t94*t94;

      mat[0][0] = -0.25*t17*t25+0.25*t28*t25+t44-0.25*t51*t25+0.25*t54*t25+t59;
      mat[0][1] = t82;
      mat[0][2] = -t90+0.25*t28*t24-t97-0.25*t54*t24;
      mat[1][0] = t82;
      mat[1][1] = -0.25*t17*t101+0.25*t28*t101+t44-0.25*t51*t101+0.25*t54*t101+t59;
      mat[1][2] = -t113+0.5*t66-t117-0.5*t77;
      mat[2][0] = -t90+0.25*t28*t24-t97-0.25*t54*t24;
      mat[2][1] = -t113+0.25*t28*t62-t117-0.25*t54*t62;
      mat[2][2] = -0.625E-1*t16*t130+0.1E1*t28-0.625E-1*t50*t134+0.1E1*t54;

#ifdef DEBUG
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        if ( aol::isNaN( mat[i][j] ) ) cerr << "\n******* InfNorm3D: Found NAN in HessMat!! ******* ";
#endif

  }

};



//! Class for the simple L1-norm-anisotropy in 3d (Wulff shape is a cube).
template <typename RealType>
class L1Norm3d : public Anisotropy3dGraphInterface<RealType, L1Norm3d<RealType> > {
public:
  using Anisotropy3dGraphInterface<RealType, L1Norm3d<RealType> >::gammaNorm;
  using Anisotropy3dGraphInterface<RealType, L1Norm3d<RealType> >::gammaFirstDerivative;
  using Anisotropy3dGraphInterface<RealType, L1Norm3d<RealType> >::implicitPart;

private:
  RealType _epsilon;
  RealType _epsSqr;

public:
  L1Norm3d ( RealType epsilon ) :
      Anisotropy3dGraphInterface<RealType, L1Norm3d<RealType> > ( ), _epsilon ( epsilon ), _epsSqr ( epsilon*epsilon ) { }

  //! Evaluating the function itself, i.e. the regularization of
  //! \f$ \gamma(z) = \sqrt{z_0^2 + \epsilon^2} + \sqrt{z_1^2 + \epsilon^2} + \sqrt{z_2^2 + \epsilon^2} \f$.
  RealType gammaNorm ( const aol::Vec3<RealType> &z ) const {
    return ( sqrt ( z[0]*z[0] + _epsSqr ) + sqrt ( z[1]*z[1] + _epsSqr ) + sqrt ( z[2]*z[2] + _epsSqr ) );
  }

  //! evaluating the first derivative: \f$ \gamma_z = \frac{z}{|z|_{\varepsilon}} \f$
  void gammaFirstDerivative ( const aol::Vec3<RealType> &, const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  //! evaluating the first derivative: \f$ \gamma_z = \frac{z}{|z|_{\varepsilon}} \f$
  void gammaFirstDerivative ( const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    v[0] = z[0] / sqrt ( z[0]*z[0] + _epsSqr );
    v[1] = z[1] / sqrt ( z[1]*z[1] + _epsSqr );
    v[2] = z[2] / sqrt ( z[2]*z[2] + _epsSqr );
  }

  //! this method computes \f$ gamma_zz \f$, it's called implicitPart just for technical reasons.
  void implicitPart ( const aol::Vec3<RealType> &z, aol::Matrix33<RealType> &mat ) const {
    mat[0][0] = _epsSqr / pow( z[0]*z[0] + _epsSqr, 1.5 );
    mat[1][1] = _epsSqr / pow( z[1]*z[1] + _epsSqr, 1.5 );
    mat[2][2] = _epsSqr / pow( z[2]*z[2] + _epsSqr, 1.5 );
    mat[0][1] = mat[0][2] = mat[1][0] = mat[1][2] = mat[2][0] = mat[2][1] = 0.;

#ifdef DEBUG
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        if ( aol::isNaN( mat[i][j] ) ) cerr << "\n******* Found NAN in HessMat!! ******* ";
#endif

  }

};

//! Class for the double cone anisotropy in 3d (Frank diagram is a cylinder aligned to the z-axis).
template <typename RealType>
class DoubleCone3d : public Anisotropy3dGraphInterface<RealType, DoubleCone3d<RealType> > {
public:
  using Anisotropy3dGraphInterface<RealType, DoubleCone3d<RealType> >::gammaNorm;
  using Anisotropy3dGraphInterface<RealType, DoubleCone3d<RealType> >::gammaFirstDerivative;
  using Anisotropy3dGraphInterface<RealType, DoubleCone3d<RealType> >::implicitPart;

private:
  RealType _delta;      //! Regularization parameter
  RealType _deltaSqr;

public:
  DoubleCone3d ( RealType delta ) :
    Anisotropy3dGraphInterface<RealType, DoubleCone3d<RealType> > ( ), _delta ( delta ), _deltaSqr ( delta*delta ) { }

  //! Evaluating the function itself, i.e. the regularization of
  //! \f$ \gamma(z) = \sqrt{z_0^2 + \epsilon^2} + \sqrt{z_1^2 + \epsilon^2} + \sqrt{z_2^2 + \epsilon^2} \f$.
  RealType gammaNorm ( const aol::Vec3<RealType> &z ) const {
    return ( 0.5 * ( sqrt( aol::Sqr( sqrt ( z[0]*z[0] + z[1]*z[1] ) + z[2] ) + _deltaSqr )
                   + sqrt( aol::Sqr( sqrt ( z[0]*z[0] + z[1]*z[1] ) - z[2] ) + _deltaSqr ) ) );
  }

  //! evaluating the first derivative: \f$ \gamma_z = \frac{z}{|z|_{\varepsilon}} \f$
  void gammaFirstDerivative ( const aol::Vec3<RealType> &, const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  //! evaluating the first derivative: \f$ \gamma_z = \frac{z}{|z|_{\varepsilon}} \f$
  void gammaFirstDerivative ( const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    const RealType t1 = z[0]*z[0];
    const RealType t2 = z[1]*z[1];
    const RealType t3 = _delta*_delta;
    const RealType t5 = sqrt(t1+t2+t3);
    const RealType t6 = t5+z[2];
    const RealType t7 = t6*t6;
    const RealType t9 = sqrt(t7+t3);
    const RealType t10 = 1/t9;
    const RealType t11 = t10*t6;
    const RealType t12 = 1/t5;
    const RealType t13 = t12*z[0];
    const RealType t16 = t5-z[2];
    const RealType t17 = t16*t16;
    const RealType t19 = sqrt(t17+t3);
    const RealType t20 = 1/t19;
    const RealType t21 = t20*t16;
    const RealType t25 = t12*z[1];
    v[0] = 0.5*t11*t13+0.5*t21*t13;
    v[1] = 0.5*t11*t25+0.5*t21*t25;
    v[2] = 0.5*t10*t6-0.5*t20*t16;
  }

  //! this method computes \f$ gamma_zz \f$, it's called implicitPart just for technical reasons.
  void implicitPart ( const aol::Vec3<RealType> &z, aol::Matrix33<RealType> &mat ) const {
    // maple part
    RealType x = z[0]; RealType y = z[1];

    const RealType t1 = x*x;
    const RealType t2 = y*y;
    const RealType t3 = _delta*_delta;
    const RealType t4 = t1+t2+t3;
    const RealType t5 = sqrt(t4);
    const RealType t6 = t5+z[2];
    const RealType t7 = t6*t6;
    const RealType t8 = t7+t3;
    const RealType t9 = sqrt(t8);
    const RealType t11 = 1/t9/t8;
    const RealType t12 = t11*t7;
    const RealType t13 = 1/t4;
    const RealType t14 = t13*t1;
    const RealType t17 = 1/t9;
    const RealType t18 = t17*t13;
    const RealType t21 = t17*t6;
    const RealType t23 = 1/t5/t4;
    const RealType t24 = t23*t1;
    const RealType t27 = 1/t5;
    const RealType t29 = 0.5*t21*t27;
    const RealType t30 = t5-z[2];
    const RealType t31 = t30*t30;
    const RealType t32 = t31+t3;
    const RealType t33 = sqrt(t32);
    const RealType t35 = 1/t33/t32;
    const RealType t36 = t35*t31;
    const RealType t39 = 1/t33;
    const RealType t40 = t39*t13;
    const RealType t43 = t39*t30;
    const RealType t47 = 0.5*t43*t27;
    const RealType t50 = t13*x*y;
    const RealType t53 = x*y;
    const RealType t57 = t23*x*y;
    const RealType t66 = -0.5*t12*t50+0.5*t18*t53-0.5*t21*t57-0.5*t36*t50+0.5*t40*t53-0.5*t43*t57;
    const RealType t68 = t27*x;
    const RealType t71 = 0.5*t11*t6*t6*t68;
    const RealType t72 = t17*t27;
    const RealType t73 = t72*x;
    const RealType t78 = -0.5*t35*t30*t30*t68;
    const RealType t79 = t39*t27;
    const RealType t80 = t79*x;
    const RealType t83 = t13*t2;
    const RealType t88 = t23*t2;
    const RealType t102 = 0.5*t11*t6*t6*t27*y;
    const RealType t103 = t72*y;
    const RealType t109 = -0.5*t35*t30*t30*t27*y;
    const RealType t110 = t79*y;
    mat[0][0] = -0.5*t12*t14+0.5*t18*t1-0.5*t21*t24+t29-0.5*t36*t14+0.5*t40*t1-0.5*t43*t24+t47;
    mat[0][1] = t66;
    mat[0][2] = -t71+0.5*t73-t78-0.5*t80;
    mat[1][0] = t66;
    mat[1][1] = -0.5*t12*t83+0.5*t18*t2-0.5*t21*t88+t29-0.5*t36*t83+0.5*t40*t2-0.5*t43*t88+t47;
    mat[1][2] = -t102+0.5*t103-t109-0.5*t110;
    mat[2][0] = -t71+0.5*t73-t78-0.5*t80;
    mat[2][1] = -t102+0.5*t103-t109-0.5*t110;
    mat[2][2] = -0.5*t11*t6*t6+0.5*t17-0.5*t35*t30*t30+0.5*t39;
  }


};

//! Class for a hexagon in 2d.
template <typename RealType>
class Hexagon2d {

private:
  RealType _delta;          //! Regularization parameter
  RealType _deltaSqr;
  RealType _alpha;          //! the angle of one hexagon-site and the x-axis (or y-axis)
  aol::Vec2<RealType> _Vm;  //! V1-V2
  aol::Vec2<RealType> _Vp;  //! V1+V2
  aol::Vec2<RealType> _V3;  //! +-V1,+-V2,+-V3 are the corners of the hexagon in the r-z-plane (r=sqrt(x^2+y^2))


public:
  Hexagon2d ( RealType delta, RealType alpha ) : _delta ( delta ), _deltaSqr ( delta*delta ), _alpha( alpha )
  {
//     _Vm[0] = 0.5;                         // regelmaessiges Hexagon, alpha=60 degrees
//     _Vm[1] = -0.8660254037844386467326525391730029923565;
//     _Vp[0] = 1.5;
//     _Vp[1] = 0.8660254037844386467326525391730029923565;
//     _V3[0] = -0.5;
//     _V3[1] = 0.8660254037844386467326525391730029923565;

//     _Vm[0] = cos( _alpha );        // gleichseitiges Hexagon
//     _Vm[1] = -sin( _alpha );
//     _Vp[0] = 1. + cos( _alpha );
//     _Vp[1] = sin( _alpha );
//     _V3[0] = -0.5;
//     _V3[1] = sin( _alpha );

    _Vm[0] = cos( _alpha ) / sin ( _alpha );  // Die Normalen auf die oberen Seitenflaechen nehmen den gleichen Wert an
    _Vm[1] = -1.;
    _Vp[0] = ( 2. - cos( _alpha ) ) / sin( _alpha );
    _Vp[1] = 1.;
    _V3[0] = (cos( _alpha ) - 1.) / sin( _alpha );
    _V3[1] = 1.;
  }

  //! Evaluating the function itself, i.e. the regularization of
  //! \f$ \gamma(z) = \sqrt{z_0^2 + \epsilon^2} + \sqrt{z_1^2 + \epsilon^2} + \sqrt{z_2^2 + \epsilon^2} \f$.
  RealType gammaNorm ( const aol::Vec2<RealType> &z ) const {
    const RealType t1 = _Vm[0]*_Vm[0];
    const RealType t2 = z[0]*z[0];
    const RealType t8 = _Vm[1]*_Vm[1];
    const RealType t9 = z[1]*z[1];
    const RealType t11 = _delta*_delta;
    const RealType t13 = sqrt(t1*t2+2.0*_Vm[0]*z[0]*_Vm[1]*z[1]+t8*t9+t11);
    const RealType t14 = 0.5*t13;
    const RealType t15 = _Vp[0]*_Vp[0];
    const RealType t21 = _Vp[1]*_Vp[1];
    const RealType t24 = sqrt(t15*t2+2.0*_Vp[0]*z[0]*_Vp[1]*z[1]+t21*t9+t11);
    const RealType t25 = 0.5*t24;
    const RealType t26 = _V3[0]*z[0];
    const RealType t27 = _V3[1]*z[1];
    const RealType t29 = pow(t14+t25-t26-t27,2.0);
    const RealType t31 = sqrt(t29+t11);
    const RealType t34 = pow(t14+t25+t26+t27,2.0);
    const RealType t36 = sqrt(t34+t11);

    return( 0.5*t31+0.5*t36 );
  }

  //! evaluating the first derivative: \f$ \gamma_z = \frac{z}{|z|_{\varepsilon}} \f$
  void gammaFirstDerivative ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  //! evaluating the first derivative: \f$ \gamma_z = \frac{z}{|z|_{\varepsilon}} \f$
  void gammaFirstDerivative ( const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    const RealType t1 = _Vm[0]*_Vm[0];
    const RealType t2 = z[0]*z[0];
    const RealType t4 = _Vm[0]*z[0];
    const RealType t8 = _Vm[1]*_Vm[1];
    const RealType t9 = z[1]*z[1];
    const RealType t11 = _delta*_delta;
    const RealType t13 = sqrt(t1*t2+2.0*t4*_Vm[1]*z[1]+t8*t9+t11);
    const RealType t14 = 0.5*t13;
    const RealType t15 = _Vp[0]*_Vp[0];
    const RealType t17 = _Vp[0]*z[0];
    const RealType t21 = _Vp[1]*_Vp[1];
    const RealType t24 = sqrt(t15*t2+2.0*t17*_Vp[1]*z[1]+t21*t9+t11);
    const RealType t25 = 0.5*t24;
    const RealType t26 = _V3[0]*z[0];
    const RealType t27 = _V3[1]*z[1];
    const RealType t28 = t14+t25-t26-t27;
    const RealType t29 = t28*t28;
    const RealType t31 = sqrt(t29+t11);
    const RealType t33 = 1/t31*t28;
    const RealType t34 = 1/t13;
    const RealType t40 = 0.5*t34*(t1*z[0]+_Vm[0]*_Vm[1]*z[1]);
    const RealType t41 = 1/t24;
    const RealType t47 = 0.5*t41*(t15*z[0]+_Vp[0]*_Vp[1]*z[1]);
    const RealType t51 = t14+t25+t26+t27;
    const RealType t52 = t51*t51;
    const RealType t54 = sqrt(t52+t11);
    const RealType t56 = 1/t54*t51;
    const RealType t65 = 0.5*t34*(t4*_Vm[1]+t8*z[1]);
    const RealType t70 = 0.5*t41*(t17*_Vp[1]+t21*z[1]);

    v[0] = 0.5*t33*(t40+t47-_V3[0])+0.5*t56*(t40+t47+_V3[0]);
    v[1] = 0.5*t33*(t65+t70-_V3[1])+0.5*t56*(t65+t70+_V3[1]);
  }

};


//! Class for an rotated hexagon in 3d.
template <typename RealType>
class RotatedHexagon3d : public Anisotropy3dGraphInterface<RealType, RotatedHexagon3d<RealType> > {

public:
  using Anisotropy3dGraphInterface<RealType, RotatedHexagon3d<RealType> >::gammaNorm;
  using Anisotropy3dGraphInterface<RealType, RotatedHexagon3d<RealType> >::gammaFirstDerivative;
  using Anisotropy3dGraphInterface<RealType, RotatedHexagon3d<RealType> >::implicitPart;

private:
  RealType _delta;          //! Regularization parameter
  RealType _deltaSqr;
  RealType _alpha;          //! the angle of one hexagon-site and the x-axis (or y-axis)
  aol::Vec2<RealType> _Vm;  //! V1-V2
  aol::Vec2<RealType> _Vp;  //! V1+V2
  aol::Vec2<RealType> _V3;  //! +-V1,+-V2,+-V3 are the corners of the hexagon in the r-z-plane (r=sqrt(x^2+y^2))


public:
  RotatedHexagon3d ( RealType delta, RealType alpha ) :
    Anisotropy3dGraphInterface<RealType, RotatedHexagon3d<RealType> > ( ), _delta ( delta ), _deltaSqr ( delta*delta ), _alpha( alpha )
  {
/*    _Vm[0] = cos( _alpha );        // gleichseitiges Hexagon
    _Vm[1] = -sin( _alpha );
    _Vp[0] = 1. + cos( _alpha );
    _Vp[1] = sin( _alpha );
    _V3[0] = -0.5;
    _V3[1] = sin( _alpha ); */

    _Vm[0] = cos( _alpha ) / sin ( _alpha );  // Die Normalen auf die oberen Seitenflaechen nehmen den gleichen Wert an
    _Vm[1] = -1.;
    _Vp[0] = ( 2. - cos( _alpha ) ) / sin( _alpha );
    _Vp[1] = 1.;
    _V3[0] = (cos( _alpha ) - 1.) / sin( _alpha );
    _V3[1] = 1.;

//     _Vm[0] = 1.;        // zeichnerisch vorgegebene Punkte, so dass die Parallele zur r-Achse erheblich billiger
//     _Vm[1] = -1.;       // ist als die Seitenwaende. Die Seitenwaende haben zur r-Achse einen Winkel von 45 Grad.
//     _Vp[0] = 5.;
//     _Vp[1] = 1.;
//     _V3[0] = -2.;
//     _V3[1] = 1.;

//     _Vm[0] = 2.5;        // zeichnerisch vorgegebene Punkte, so dass die Parallele zur r-Achse erheblich teurer
//     _Vm[1] = -2.5;       // ist als die Seitenwaende. Die Seitenwaende haben zur r-Achse einen Winkel von 45 Grad.
//     _Vp[0] = 3.5;
//     _Vp[1] = 2.5;
//     _V3[0] = -0.5;
//     _V3[1] = 2.5;

  }

  //! Evaluating the function itself, i.e. the regularization of
  //! \f$ \gamma(z) = \sqrt{z_0^2 + \epsilon^2} + \sqrt{z_1^2 + \epsilon^2} + \sqrt{z_2^2 + \epsilon^2} \f$.
  RealType gammaNorm ( const aol::Vec3<RealType> &z ) const {
    const RealType t1 = z[0]*z[0];
    const RealType t2 = z[1]*z[1];
    const RealType t3 = _delta*_delta;
    const RealType t5 = sqrt(t1+t2+t3);
    const RealType t9 = aol::Sqr(_Vm[0]*t5+_Vm[1]*z[2]);
    const RealType t11 = sqrt(t9+t3);
    const RealType t12 = 0.5*t11;
    const RealType t16 = aol::Sqr(_Vp[0]*t5+_Vp[1]*z[2]);
    const RealType t18 = sqrt(t16+t3);
    const RealType t19 = 0.5*t18;
    const RealType t20 = _V3[0]*t5;
    const RealType t21 = _V3[1]*z[2];
    const RealType t23 = aol::Sqr(t12+t19-t20-t21);
    const RealType t25 = sqrt(t23+t3);
    const RealType t28 = aol::Sqr(t12+t19+t20+t21);
    const RealType t30 = sqrt(t28+t3);
    return( 0.5*t25+0.5*t30 );
  }

  //! evaluating the first derivative: \f$ \gamma_z = \frac{z}{|z|_{\varepsilon}} \f$
  void gammaFirstDerivative ( const aol::Vec3<RealType> &, const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  //! evaluating the first derivative: \f$ \gamma_z = \frac{z}{|z|_{\varepsilon}} \f$
  void gammaFirstDerivative ( const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    const RealType t1 = z[0]*z[0];
    const RealType t2 = z[1]*z[1];
    const RealType t3 = _delta*_delta;
    const RealType t5 = sqrt(t1+t2+t3);
    const RealType t8 = _Vm[0]*t5+_Vm[1]*z[2];
    const RealType t9 = t8*t8;
    const RealType t11 = sqrt(t9+t3);
    const RealType t12 = 0.5*t11;
    const RealType t15 = _Vp[0]*t5+_Vp[1]*z[2];
    const RealType t16 = t15*t15;
    const RealType t18 = sqrt(t16+t3);
    const RealType t19 = 0.5*t18;
    const RealType t20 = _V3[0]*t5;
    const RealType t21 = _V3[1]*z[2];
    const RealType t22 = t12+t19-t20-t21;
    const RealType t23 = t22*t22;
    const RealType t25 = sqrt(t23+t3);
    const RealType t27 = 1/t25*t22;
    const RealType t29 = 1/t11*t8;
    const RealType t30 = 1/t5;
    const RealType t31 = _Vm[0]*t30;
    const RealType t34 = 0.5*t29*t31*z[0];
    const RealType t36 = 1/t18*t15;
    const RealType t37 = _Vp[0]*t30;
    const RealType t40 = 0.5*t36*t37*z[0];
    const RealType t41 = _V3[0]*t30;
    const RealType t42 = t41*z[0];
    const RealType t46 = t12+t19+t20+t21;
    const RealType t47 = t46*t46;
    const RealType t49 = sqrt(t47+t3);
    const RealType t51 = 1/t49*t46;
    const RealType t58 = 0.5*t29*t31*z[1];
    const RealType t61 = 0.5*t36*t37*z[1];
    const RealType t62 = t41*z[1];
    const RealType t71 = 0.5*t29*_Vm[1];
    const RealType t73 = 0.5*t36*_Vp[1];
    v[0] = 0.5*t27*(t34+t40-t42)+0.5*t51*(t34+t40+t42);
    v[1] = 0.5*t27*(t58+t61-t62)+0.5*t51*(t58+t61+t62);
    v[2] = 0.5*t27*(t71+t73-_V3[1])+0.5*t51*(t71+t73+_V3[1]);

//       for (int j=0; j<3; j++)
//         cerr << v[j] << ", ";
//         if ( aol::isNaN( v[j] ) ) cerr << "\n******* Found NAN in Grad!! ******* ";
  }

   //! this method computes \f$ gamma_zz \f$, it's called implicitPart just for technical reasons.
  void implicitPart ( const aol::Vec3<RealType> &z, aol::Matrix33<RealType> &mat ) const {

    const RealType t1 = z[0]*z[0];              const RealType t2 = z[1]*z[1];                const RealType t3 = _delta*_delta;
    const RealType t4 = t1+t2+t3;               const RealType t5 = sqrt(t4);                 const RealType t8 = _Vm[0]*t5+_Vm[1]*z[2];
    const RealType t9 = t8*t8;                  const RealType t10 = t9+t3;                   const RealType t11 = sqrt(t10);
    const RealType t12 = 0.5*t11;               const RealType t15 = _Vp[0]*t5+_Vp[1]*z[2];   const RealType t16 = t15*t15;
    const RealType t17 = t16+t3;                const RealType t18 = sqrt(t17);               const RealType t19 = 0.5*t18;
    const RealType t20 = _V3[0]*t5;             const RealType t21 = _V3[1]*z[2];             const RealType t22 = t12+t19-t20-t21;
    const RealType t23 = t22*t22;               const RealType t24 = t23+t3;                  const RealType t25 = sqrt(t24);
    const RealType t28 = 1/t25/t24*t23;         const RealType t29 = 1/t11;                   const RealType t30 = t29*t8;
    const RealType t31 = 1/t5;                  const RealType t32 = _Vm[0]*t31;              const RealType t33 = t32*z[0];
    const RealType t35 = 0.5*t30*t33;           const RealType t36 = 1/t18;                   const RealType t37 = t36*t15;
    const RealType t38 = _Vp[0]*t31;            const RealType t39 = t38*z[0];                const RealType t41 = 0.5*t37*t39;
    const RealType t42 = _V3[0]*t31;            const RealType t43 = t42*z[0];                const RealType t44 = t35+t41-t43;
    const RealType t45 = t44*t44;               const RealType t48 = 1/t25;                   const RealType t51 = t48*t22;
    const RealType t54 = 1/t11/t10*t9;          const RealType t55 = _Vm[0]*_Vm[0];           const RealType t56 = 1/t4;
    const RealType t57 = t55*t56;               const RealType t60 = 0.5*t54*t57*t1;          const RealType t61 = t29*t55;
    const RealType t62 = t56*t1;                const RealType t64 = 0.5*t61*t62;             const RealType t66 = 1/t5/t4;
    const RealType t67 = _Vm[0]*t66;            const RealType t70 = 0.5*t30*t67*t1;          const RealType t72 = 0.5*t30*t32;
    const RealType t75 = 1/t18/t17*t16;         const RealType t76 = _Vp[0]*_Vp[0];           const RealType t77 = t76*t56;
    const RealType t80 = 0.5*t75*t77*t1;        const RealType t81 = t36*t76;                 const RealType t83 = 0.5*t81*t62;
    const RealType t84 = _Vp[0]*t66;            const RealType t87 = 0.5*t37*t84*t1;          const RealType t89 = 0.5*t37*t38;
    const RealType t90 = _V3[0]*t66;            const RealType t91 = t90*t1;                  const RealType t95 = t12+t19+t20+t21;
    const RealType t96 = t95*t95;               const RealType t97 = t96+t3;                  const RealType t98 = sqrt(t97);
    const RealType t101 = 1/t98/t97*t96;        const RealType t102 = t35+t41+t43;            const RealType t103 = t102*t102;
    const RealType t106 = 1/t98;                const RealType t109 = t106*t95;               const RealType t114 = t32*z[1];
    const RealType t116 = 0.5*t30*t114;         const RealType t117 = t38*z[1];               const RealType t119 = 0.5*t37*t117;
    const RealType t120 = t42*z[1];             const RealType t121 = t116+t119-t120;         const RealType t130 = t56*z[0]*z[1];
    const RealType t132 = 0.5*t54*t55*t130;     const RealType t134 = 0.5*t61*t130;           const RealType t137 = t66*z[0]*z[1];
    const RealType t139 = 0.5*t30*_Vm[0]*t137;  const RealType t142 = 0.5*t75*t76*t130;       const RealType t144 = 0.5*t81*t130;
    const RealType t147 = 0.5*t37*_Vp[0]*t137;  const RealType t149 = t90*z[0]*z[1];          const RealType t153 = t116+t119+t120;
    const RealType t163 = -0.5*t28*t44*t121+0.5*t48*t121*t44+0.5*t51*(-t132+t134-t139-t142+t144-t147+t149)-0.5*t101*t102*t153
                          +0.5*t106*t153*t102+0.5*t109*(-t132+t134-t139-t142+t144-t147-t149);
    const RealType t165 = 0.5*t30*_Vm[1];       const RealType t167 = 0.5*t37*_Vp[1];         const RealType t168 = t165+t167-_V3[1];
    const RealType t172 = t48*t168;             const RealType t175 = t54*_Vm[0];             const RealType t176 = t31*z[0];
    const RealType t180 = t29*_Vm[1];           const RealType t183 = t75*_Vp[0];             const RealType t187 = t36*_Vp[1];
    const RealType t190 = -0.5*t175*t176*_Vm[1]+0.5*t180*t33-0.5*t183*t176*_Vp[1]+0.5*t187*t39;
    const RealType t193 = t165+t167+_V3[1];     const RealType t197 = t106*t193;
    const RealType t202 = -0.5*t28*t44*t168+0.5*t172*t44+0.5*t51*t190-0.5*t101*t102*t193+0.5*t197*t102+0.5*t109*t190;
    const RealType t203 = t121*t121;            const RealType t210 = 0.5*t54*t57*t2;         const RealType t211 = t56*t2;
    const RealType t213 = 0.5*t61*t211;         const RealType t216 = 0.5*t30*t67*t2;         const RealType t219 = 0.5*t75*t77*t2;
    const RealType t221 = 0.5*t81*t211;         const RealType t224 = 0.5*t37*t84*t2;         const RealType t225 = t90*t2;
    const RealType t229 = t153*t153;            const RealType t243 = t31*z[1];
    const RealType t254 = -0.5*t175*t243*_Vm[1]+0.5*t180*t114-0.5*t183*t243*_Vp[1]+0.5*t187*t117;
    const RealType t264 = -0.5*t28*t121*t168+0.5*t172*t121+0.5*t51*t254-0.5*t101*t153*t193+0.5*t197*t153+0.5*t109*t254;
    const RealType t265 = t168*t168;            const RealType t270 = _Vm[1]*_Vm[1];          const RealType t275 = _Vp[1]*_Vp[1];
    const RealType t280 = -0.5*t54*t270+0.5*t29*t270-0.5*t75*t275+0.5*t36*t275;               const RealType t283 = t193*t193;

    mat[0][0] = -0.5*t28*t45+0.5*t48*t45+0.5*t51*(-t60+t64-t70+t72-t80+t83-t87+t89+t91-t42)-0.5*t101*t103
                +0.5*t106*t103+0.5*t109*(-t60+t64-t70+t72-t80+t83-t87+t89-t91+t42);
    mat[0][1] = t163;
    mat[0][2] = t202;
    mat[1][0] = t163;
    mat[1][1] = -0.5*t28*t203+0.5*t48*t203+0.5*t51*(-t210+t213-t216+t72-t219+t221-t224+t89+t225-t42)-0.5*t101*t229
                +0.5*t106*t229+0.5*t109*(-t210+t213-t216+t72-t219+t221-t224+t89-t225+t42);
    mat[1][2] = t264;
    mat[2][0] = t202;
    mat[2][1] = t264;
    mat[2][2] = -0.5*t28*t265+0.5*t48*t265+0.5*t51*t280-0.5*t101*t283+0.5*t106*t283+0.5*t109*t280;

    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        if ( aol::isNaN( mat[i][j] ) ) cerr << "\n******* Found NAN in HessMat!! ******* ";

  }

};

/** class for computing several terms for the anisotropy
    \f$\tilde{\gamma}(z) = \left( \gamma^p(z) + \epsilon^p \right)^{\frac{1}{p}}\f$
    with \f$\gamma(z) = \left( z_1^q + z_2^q \right)^{\frac{1}{q}}\f$.
    */
template <typename RealType>
class LpAnisotropy {
private:
  RealType _p, _q, _epsilon;
public:
  LpAnisotropy ( RealType P, RealType Q, RealType Epsilon ) : _p ( P ), _q ( Q ), _epsilon ( Epsilon ) { }

  RealType p() const { return _p; }

  void implicitPart ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z,
                      aol::Matrix22<RealType> &mat ) const {
    implicitPart ( z, mat );
  }

  void explicitPart ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z,
                      aol::Matrix22<RealType> &mat ) const {
    explicitPart ( z, mat );
  }

  /** factor of the qc::WillmoreStiffImpOp:
    *  \f$\frac{\left[ \gamma^p(z) \right]_{zz} }{p \|z\|^{p-1}_{\gamma,\epsilon,p}}\f$
    */
  void implicitPart ( const aol::Vec2<RealType> &z, aol::Matrix22<RealType> &mat ) const {
    aol::Vec2<RealType> Z ( aol::Abs ( z[0] ), aol::Abs ( z[1] ) );
    RealType a = pow ( Z[0], _q ) + pow ( Z[1], _q );

    if ( a != 0 ) {
      RealType factor1 = 0.;
      RealType factor2 = 0.;

      if ( _p != _q ) {
        factor1 = ( _p - _q ) * pow ( a, _p / _q - 2. );
        factor2 = pow ( a, _p / _q - 1. ) * ( _q - 1. );
      }

      mat[0][0] = factor1 * pow ( Z[0], 2 * ( _q - 1. ) ) + factor2 * pow ( Z[0], _q - 2. );
      mat[1][1] = factor1 * pow ( Z[1], 2 * ( _q - 1. ) ) + factor2 * pow ( Z[1], _q - 2. );

      mat[1][0] = mat[0][1] = factor1 * pow ( Z[0] * Z[1], _q - 1. ) * aol::signum1at0 ( z[0] ) * aol::signum1at0 ( z[1] );

      // Now mat = 0.5 * gamma^2_zz
      mat /= pow ( pow ( a, _p / _q ) + pow ( _epsilon, _p ), ( _p - 1. ) / _p ) ;
    } else {
      mat[0][0] = mat[1][1] = mat[1][0] = mat[1][0] = 0.;
    }
  }

  /** factor of the qc::WillmoreStiffExpOp:
    * \f$\frac{p-1}{p^2} \frac{\gamma^p(z)_z \otimes \gamma^p(z)_z}{\|z\|^{2p-1}_{\gamma,\epsilon,p}}\f$
    */
  // I've removed a little bug from Marc (at least I think so) and optimized the code a little bit.
  // For the latest version of Marc's computation see revision 3738.
  void explicitPart ( const aol::Vec2<RealType> &z, aol::Matrix22<RealType> &mat ) const {
    aol::Vec2<RealType> Z ( aol::Abs ( z[0] ), aol::Abs ( z[1] ) );
    RealType a = pow ( Z[0], _q ) + pow ( Z[1], _q );

    if ( a != 0 ) {
      mat[0][0] = pow ( Z[0], 2. * ( _q - 1. ) );
      mat[1][1] = pow ( Z[1], 2. * ( _q - 1. ) );

      mat[1][0] = mat[0][1] = pow ( Z[0], _q - 1. ) * pow ( Z[1], _q - 1. ) * aol::signum1at0 ( z[0] ) * aol::signum1at0 ( z[1] );

      RealType gammaFactor = pow ( a, 2. * ( _p / _q - 1. ) );
      RealType normSqr = pow ( pow ( a, _p / _q ) + pow ( _epsilon, _p ), ( 2 * _p - 1. ) / _p );

      mat *= gammaFactor * ( _p - 1. ) / normSqr;
    } else {
      mat[0][0] = mat[1][1] = mat[1][0] = mat[1][0] = 0.;
    }
  }

  //! evaluating the function itself \f$\gamma(z) = \left( z_1^q + z_2^q \right)^{\frac{1}{q}}\f$
  void gammaFirstDerivative ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  /** method computes:
   * \f$ \frac{ \gamma^{p-q}(z) \left( |z_i|^{q-1} sgn(z_i) \right)_i }{ \| z \|^{p-1}_{\gamma,\epsilon,p} }  \f$
   */
  void gammaFirstDerivative ( const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    aol::Vec2<RealType> Z ( aol::Abs ( z[0] ), aol::Abs ( z[1] ) );
    RealType a = pow ( Z[0], _q ) + pow ( Z[1], _q );
    if ( a != 0 ) {
      v[0] = pow ( Z[0], _q - 1. ) * aol::signum1at0 ( z[0] );
      v[1] = pow ( Z[1], _q - 1. ) * aol::signum1at0 ( z[1] );

      v *= pow ( a, _p / _q - 1. ) *  pow ( pow ( a, _p / _q ) + pow ( _epsilon, _p ), ( 1. - _p ) / _p );
    } else {
      v[0] = v[1] = 0.;
    }
  }

  //! evaluating the function itself \f$\gamma(z) = \left( z_1^q + z_2^q \right)^{\frac{1}{q}}\f$
  RealType gammaNorm ( const aol::Vec2<RealType> &z ) const {
    return pow ( pow ( aol::Abs ( z[0] ), _q ) + pow ( aol::Abs ( z[1] ), _q ) + pow ( _epsilon, _q ), 1. / _q );
  }
};



/* *******************************************************************************
 * gamma = Lp-Norm in 3d
 * ******************************************************************************* */
/*! class for computing several terms for the anisotropy
    \f$\tilde{\gamma}(z) = \left( \gamma^p(z) + \epsilon^p \right)^{\frac{1}{p}}\f$
    with \f$\gamma(z) = \left( z_1^q + z_2^q \right)^{\frac{1}{q}}\f$.
    */
template <typename RealType>
class LpAnisotropy3d {
private:
  RealType _p, _q, _epsilon;
public:
  LpAnisotropy3d ( RealType P, RealType Q, RealType Epsilon ) : _p ( P ), _q ( Q ), _epsilon ( Epsilon ) {
    if ( Q < 1 ) cerr << aol::color::red << "\n\nWARNING: Q is smaller 1, there might be singularities!!\n";
  }

RealType p() const { return _p; }

  void implicitPart ( const aol::Vec3<RealType> &, const aol::Vec3<RealType> &z,
                      aol::Matrix33<RealType> &mat ) const {
    implicitPart ( z, mat );
  }

  void explicitPart ( const aol::Vec3<RealType> &, const aol::Vec3<RealType> &z,
                      aol::Matrix33<RealType> &mat ) const {
    explicitPart ( z, mat );
  }

  /** factor of the qc::WillmoreStiffImpOp:
    *  \f$\frac{\left[ \gamma^p(z) \right]_{zz} }{p \|z\|^{p-1}_{\gamma,\epsilon,p}}\f$
    */
  void implicitPart ( const aol::Vec3<RealType> &z, aol::Matrix33<RealType> &mat ) const {
    aol::Vec3<RealType> Z ( aol::Abs ( z[0] ), aol::Abs ( z[1] ), aol::Abs ( z[2] ) );
    RealType a = pow ( Z[0], _q ) + pow ( Z[1], _q ) + pow ( Z[2], _q );

    if ( a != 0 ) {
      RealType factor2 = 0.;
      RealType factor1 = 0.;
      if ( _p != _q ) {
        factor1 = ( _p - _q ) * pow ( a, _p / _q - 2. );
        factor2 = ( _q - 1. ) * pow ( a, _p / _q - 1. );
      }

      mat[0][0] = factor1 * pow ( Z[0], 2 * ( _q - 1. ) ) + factor2 * pow ( Z[0], _q - 2. );
      mat[1][1] = factor1 * pow ( Z[1], 2 * ( _q - 1. ) ) + factor2 * pow ( Z[1], _q - 2. );
      mat[2][2] = factor1 * pow ( Z[2], 2 * ( _q - 1. ) ) + factor2 * pow ( Z[2], _q - 2. );

      mat[1][0] = mat[0][1] = factor1 * pow ( Z[0] * Z[1], _q - 1. ) * aol::signum1at0 ( z[0] ) * aol::signum1at0 ( z[1] );
      mat[2][0] = mat[0][2] = factor1 * pow ( Z[0] * Z[2], _q - 1. ) * aol::signum1at0 ( z[0] ) * aol::signum1at0 ( z[2] );
      mat[2][1] = mat[1][2] = factor1 * pow ( Z[2] * Z[1], _q - 1. ) * aol::signum1at0 ( z[2] ) * aol::signum1at0 ( z[1] );

      // Now mat = 0.5 * gamma^2_zz
      mat /= pow ( pow ( a, _p / _q ) + pow ( _epsilon, _p ), ( _p - 1. ) / _p ) ;
    } else {
      mat[0][0] = mat[0][1] = mat[0][2] = 0.;
      mat[1][0] = mat[1][1] = mat[1][2] = 0.;
      mat[2][0] = mat[2][1] = mat[2][2] = 0.;
    }
  }

  /** factor of the qc::WillmoreStiffExpOp:
    * \f$\frac{p-1}{p^2} \frac{\gamma^p(z)_z \otimes \gamma^p(z)_z}{\|z\|^{2p-1}_{\gamma,\epsilon,p}}\f$
    */
  void explicitPart ( const aol::Vec3<RealType> &z, aol::Matrix33<RealType> &mat ) const {
    aol::Vec3<RealType> Z ( aol::Abs ( z[0] ), aol::Abs ( z[1] ), aol::Abs ( z[2] ) );
    RealType a = pow ( Z[0], _q ) + pow ( Z[1], _q ) + pow ( Z[2], _q );

    if ( a != 0 ) {
      mat[0][0] = pow ( Z[0], 2. * ( _q - 1. ) );
      mat[1][1] = pow ( Z[1], 2. * ( _q - 1. ) );
      mat[2][2] = pow ( Z[2], 2. * ( _q - 1. ) );

      mat[1][0] = mat[0][1] = pow ( Z[0], _q - 1. ) * pow ( Z[1], _q - 1. ) * aol::signum1at0 ( z[0] ) * aol::signum1at0 ( z[1] );
      mat[2][0] = mat[0][2] = pow ( Z[0], _q - 1. ) * pow ( Z[2], _q - 1. ) * aol::signum1at0 ( z[0] ) * aol::signum1at0 ( z[2] );
      mat[1][2] = mat[2][1] = pow ( Z[2], _q - 1. ) * pow ( Z[1], _q - 1. ) * aol::signum1at0 ( z[2] ) * aol::signum1at0 ( z[1] );

      RealType gammaFactor = pow ( a, 2. * ( _p / _q - 1. ) );
      RealType normSqr = pow ( pow ( a, _p / _q ) + pow ( _epsilon, _p ), ( 2 * _p - 1. ) / _p );

      mat *= gammaFactor * ( _p - 1. ) / normSqr;
    } else {
      mat[0][0] = mat[0][1] = mat[0][2] = 0.;
      mat[1][0] = mat[1][1] = mat[1][2] = 0.;
      mat[2][0] = mat[2][1] = mat[2][2] = 0.;
    }
  }


  //! evaluating the function itself \f$\gamma(z) = \left( z_1^q + z_2^q \right)^{\frac{1}{q}}\f$
  RealType gamma ( const aol::Vec3<RealType> &z ) const {
    return pow ( pow ( aol::Abs ( z[0] ), _q ) + pow ( aol::Abs ( z[1] ), _q )
                 + pow ( aol::Abs ( z[2] ), _q ) + pow ( _epsilon, _q ),  1 / _q );
  }


  void gammaFirstDerivative ( const aol::Vec3<RealType> &, const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  /** method computes:
   * \f$ \frac{ \gamma^{p-q}(z) \left( |z_i|^{q-1} sgn(z_i) \right)_i }{ \| z \|^{p-1}_{\gamma,\epsilon,p} }  \f$
     */
  void gammaFirstDerivative ( const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    aol::Vec3<RealType> Z ( aol::Abs ( z[0] ), aol::Abs ( z[1] ), aol::Abs ( z[2] ) );
    RealType a = pow ( Z[0], _q ) + pow ( Z[1], _q ) + pow ( Z[2], _q );

    if ( a != 0 ) {
      v[0] = pow ( Z[0], _q - 1. ) * aol::signum1at0 ( z[0] );
      v[1] = pow ( Z[1], _q - 1. ) * aol::signum1at0 ( z[1] );
      v[2] = pow ( Z[2], _q - 1. ) * aol::signum1at0 ( z[2] );

      v *= pow ( a, _p / _q - 1. ) *  pow ( pow ( a, _p / _q ) + pow ( _epsilon, _p ), ( 1. - _p ) / _p );
    } else {
      v[0] = v[1] = v[2] = 0.;
    }
  }

  //! evaluating the function itself \f$\gamma(z) = \left( z_1^q + z_2^q \right)^{\frac{1}{q}}\f$
  RealType gammaNorm ( const aol::Vec3<RealType> &z ) const {
    return pow ( pow ( aol::Abs ( z[0] ), _q ) + pow ( aol::Abs ( z[1] ), _q )
                 + pow ( aol::Abs ( z[2] ), _q ) + pow ( _epsilon, _q ), 1. / _q );
  }
};



//! Class for a pedestal-like Wulffshape in 2D.
template <typename RealType>
class Pedestal2d {
private:
  RealType _alpha;                  //! angle of the walls
  RealType _alphaHalf;              //! 0.5 * alpha (often needed)
  RealType _deltaAbs;               //! Regularization-parameter for the absolute-value
  RealType _deltaAbsSqr;            //! Square of Delta (very complicated...)
  RealType _deltaRad;               //! Regularization-parameter for the radius
  RealType _deltaRadSqr;
  RealType _c1, _c2;                //! two constants for the computation (Franck-Diagram-Ansatz)
  RealType _sinAlpha, _cosAlpha;    //! two more constants for the computation (WulffShape-Ansatz)

  RealType _phi;                    //! the angle of z in polar coordinates
  RealType _r;                      //! |z|   [ z = r (cos phi, sin phi) ]

public:
  Pedestal2d ( RealType alpha, RealType deltaAbs, RealType deltaRad ) :
      _alpha ( alpha ), _deltaAbs ( deltaAbs ), _deltaRad ( deltaRad ) {
    _alphaHalf = 0.5 * _alpha;
    _deltaAbsSqr = _deltaAbs * _deltaAbs;
    _deltaRadSqr = _deltaRad * _deltaRad;
    _c1 = cos ( _alphaHalf );
    _c2 = sin ( _alphaHalf );
    _sinAlpha = sin ( _alpha );
    _cosAlpha = cos ( _alpha );
  }

  //! Evaluating the function itself, i.e. the regularization of
  //! \f$ \gamma(z) = \frac12 |z| max\left\{ \frac{\sin(\frac{\alpha}{2}+\varphi)}{\cos \frac{\alpha}{2}}, \frac{\cos(\frac{\alpha}{2}+\varphi)}{\sin \frac{\alpha}{2}} \right\}  \f$.
  RealType gammaNormFromDistFct ( const aol::Vec2<RealType> &z ) const {
    double r = sqrt ( z.normSqr() + _deltaRadSqr );
    double phi = asin ( z[1] / r );    // in [-pi/2, pi/2]

    double CSin = sin ( _alphaHalf + phi );
    double CCos = cos ( _alphaHalf + phi );

    double term1 = sqrt ( aol::Sqr ( CSin / _c1 + CCos / _c2 ) + _deltaAbsSqr );
    double term2 = sqrt ( aol::Sqr ( CSin / _c1 - CCos / _c2 ) + _deltaAbsSqr );
    return 0.5 * r * ( term1 + term2 );
  }

  RealType gammaNormFromWulffShape ( const aol::Vec2<RealType> &z ) const {
    double r = sqrt ( z.normSqr() + _deltaRadSqr );
    double phi = asin ( z[1] / r );    // in [-pi/2, pi/2]

    double CSin = sin ( phi );
    double CCos = cos ( phi );

    double term1 = sqrt ( aol::Sqr ( CCos / _sinAlpha ) + _deltaAbsSqr );
    double term2 = sqrt ( aol::Sqr ( CSin - CCos * _cosAlpha / _sinAlpha ) + _deltaAbsSqr );

    return r * ( term1 + term2 );
  }

  RealType gammaNorm ( const aol::Vec2<RealType> &z ) const {
    return gammaNormFromWulffShape ( z );
  }

  RealType evaluate ( const aol::Vec2<RealType> &z ) const {
    return gammaNorm ( z );
  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    double r       = z.normSqr();
    double r_delta = sqrt ( r + _deltaRadSqr );
    double phi     = asin ( z[1] / r_delta );         //! in [-pi/2, pi/2]

    double CSin = sin ( phi );
    double CCos = cos ( phi );

    //! \f$ \frac{d}{dr} \gamma \f$
    double dGammaR = gammaNorm ( z );
    dGammaR *= r / ( r_delta * r_delta * r_delta );   // here is already 1/r_delta from the whole formula incorporated

    double C = _cosAlpha / _sinAlpha;
    double S = sqrt ( r_delta * r_delta - z[1] * z[1] );

    //! \f$ \frac{d}{d \varphi} \gamma \f$
    double dGammaPhi = ( CSin - C * CCos ) * ( CCos + C * CSin ) / sqrt ( aol::Sqr ( CSin - C * CCos ) + _deltaAbsSqr );
    dGammaPhi -= CSin * CCos / ( _sinAlpha * _sinAlpha * sqrt ( C * C + _deltaAbsSqr ) );
    //! Some further computing
    dGammaPhi /= r_delta;

    v[0] = z[0] * ( dGammaR - dGammaPhi * z[1] / S );
    v[1] = z[1] * dGammaR + dGammaPhi * S;

  }
};

//! Class for a pedestal-like Wulffshape.
template <typename RealType>
class Pedestal3d : public Anisotropy3dGraphInterface<RealType, Pedestal3d<RealType> > {

public:
  using Anisotropy3dGraphInterface<RealType, Pedestal3d<RealType> >::gammaNorm;
  using Anisotropy3dGraphInterface<RealType, Pedestal3d<RealType> >::gammaFirstDerivative;

private:
  RealType _alpha;        //! angle of the walls
  RealType _alphaHalf;    //! 0.5 * alpha (often needed)
  RealType _deltaAbs;     //! Regularization-parameter for the absolute-value
  RealType _deltaAbsSqr;  //! Square of Delta (very complicated...)
  RealType _deltaRad;     //! Regularization-parameter for the radius
  RealType _deltaRadSqr;
  RealType _c1, _c2;      //! two constants for the computation
  RealType _sinAlpha, _cosAlpha;    //! two more constants for the computation (WulffShape-Ansatz)

public:
  Pedestal3d ( RealType alpha, RealType deltaAbs, RealType deltaRad ) :
      Anisotropy3dGraphInterface<RealType, Pedestal3d<RealType> > ( ),
      _alpha ( alpha ), _deltaAbs ( deltaAbs ), _deltaRad ( deltaRad ) {
    _alphaHalf = 0.5 * _alpha;
    _deltaAbsSqr = _deltaAbs * _deltaAbs;
    _deltaRadSqr = _deltaRad * _deltaRad;
    _c1 = cos ( _alphaHalf );
    _c2 = sin ( _alphaHalf );
    _sinAlpha = sin ( _alpha );
    _cosAlpha = cos ( _alpha );
  }

  //! Evaluating the function itself, i.e. the regularization of
  //! \f$ \gamma(z) = \frac12 |z| max\left\{ \frac{\sin(\frac{\alpha}{2}+\varphi)}{\cos \frac{\alpha}{2}}, \frac{\cos(\frac{\alpha}{2}+\varphi)}{\sin \frac{\alpha}{2}} \right\}  \f$.
  RealType gammaNormFromFranckDiagram ( const aol::Vec3<RealType> &z ) const {
    double r   = sqrt ( z.normSqr() + _deltaRadSqr );
    double phi = asin ( z[2] / r );    //! in [-pi/2, pi/2]

    double CSin = sin ( _alphaHalf + phi );
    double CCos = cos ( _alphaHalf + phi );

    double term1 = sqrt ( aol::Sqr ( CSin / _c1 + CCos / _c2 ) + _deltaAbsSqr );
    double term2 = sqrt ( aol::Sqr ( CSin / _c1 - CCos / _c2 ) + _deltaAbsSqr );
    return 0.5 * r * ( term1 + term2 );
  }

  RealType gammaNormFromWulffShape ( const aol::Vec3<RealType> &z ) const {
    double r   = sqrt ( z.normSqr() + _deltaRadSqr );
    double phi = asin ( z[2] / r );    //! in [-pi/2, pi/2]

    double CSin = sin ( phi );
    double CCos = cos ( phi );

    double term1 = sqrt ( aol::Sqr ( CCos / _sinAlpha ) + _deltaAbsSqr );
    double term2 = sqrt ( aol::Sqr ( CSin - CCos * _cosAlpha / _sinAlpha ) + _deltaAbsSqr );
    return r * ( term1 + term2 );
  }

  RealType gammaNorm ( const aol::Vec3<RealType> &z ) const {
    return gammaNormFromWulffShape ( z );
  }

  //! evaluating the first derivative: (too long for TeXing)
  void gammaFirstDerivativeFromFranckDiagram ( const aol::Vec3<RealType> &z,
                                               aol::Vec3<RealType> &v ) const {
    double r   = sqrt ( z.normSqr() + _deltaRadSqr );
    double phi = 0.;
    if ( r != 0. ) phi = asin ( z[2] / r );  //! in [-pi/2, pi/2]

    double CSin = sin ( _alphaHalf + phi );
    double CCos = cos ( _alphaHalf + phi );

    //! \f$ \frac{d}{dr} \gamma \f$
    double dGammaR = gammaNorm ( z );
    if ( r == 0 )  dGammaR = 0.;
    else         dGammaR /= r * r;    // one 1/r is already from the derivative of r

    //! \f$ \frac{d}{d \varphi} \gamma \f$
    double dGammaPhi = ( CSin / _c1 + CCos / _c2 ) * ( CCos / _c1 - CSin / _c2 ) / sqrt ( aol::Sqr ( CSin / _c1 + CCos / _c2 ) + _deltaAbsSqr );
    dGammaPhi += ( CSin / _c1 - CCos / _c2 ) * ( CCos / _c1 + CSin / _c2 ) / sqrt ( aol::Sqr ( CSin / _c1 - CCos / _c2 ) + _deltaAbsSqr );
    dGammaPhi *= 0.5 * r;
    //! further summing-up
    dGammaPhi /= r * r * r * sqrt ( 1. - z[1] * z[1] / ( r * r ) );


    v[0] = z[0] * ( dGammaR - z[1] * dGammaPhi );
    v[1] = z[1] * dGammaR + ( r * r - z[1] * z[1] ) * dGammaPhi;
    v[2] = z[2] * ( dGammaR - z[1] * dGammaPhi );
  }


  //! evaluating the first derivative: (too long for TeXing)
  void gammaFirstDerivativeFromWulffShape ( const aol::Vec3<RealType> &z,
                                            aol::Vec3<RealType> &v ) const {
    double r   = sqrt ( z.normSqr() + _deltaRadSqr );
    double phi = asin ( z[2] / r );    //! in [-pi/2, pi/2]

    double CSin = sin ( phi );
    double CCos = cos ( phi );

    //! \f$ \frac{d}{dr} \gamma \f$
    double dGammaR = gammaNorm ( z );
    dGammaR /= r * r;

    double C = _cosAlpha / _sinAlpha;
    double S = sqrt ( r * r - z[2] * z[2] );

    //! \f$ \frac{d}{d \varphi} \gamma \f$
    double dGammaPhi = ( CSin - C * CCos ) * ( CCos + C * CSin ) / sqrt ( aol::Sqr ( CSin - C * CCos ) + _deltaAbsSqr );
    dGammaPhi -= CSin * CCos / ( _sinAlpha * _sinAlpha * sqrt ( C * C + _deltaAbsSqr ) );
    //! Some further computing
    dGammaPhi /= r;

    v[0] = z[0] * ( dGammaR - dGammaPhi * z[2] / S );
    v[1] = z[1] * ( dGammaR - dGammaPhi * z[2] / S );
    v[2] = z[2] * dGammaR + dGammaPhi * S;
  }

  //! evaluating the first derivative: (too long for TeXing)
  void gammaFirstDerivative ( const aol::Vec3<RealType> &, const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }


  void gammaFirstDerivative ( const aol::Vec3<RealType> &z,
                              aol::Vec3<RealType> &v ) const {
    gammaFirstDerivativeFromWulffShape ( z, v );
  }

};



//! Class for a pentagon-like Wulffshape in 2D.
template <typename RealType>
class Pentagon2d {
private:
  RealType _deltaAbs;     //! Regularization-parameter for the absolute-value
  RealType _deltaAbsSqr;  //! Square of Delta (very complicated...)
  RealType _deltaRad;     //! Regularization-parameter for the radius
  RealType _deltaRadSqr;

public:
  Pentagon2d ( RealType deltaAbs, RealType deltaRad ) :
      _deltaAbs ( deltaAbs ), _deltaRad ( deltaRad ) {
    _deltaAbsSqr = _deltaAbs * _deltaAbs;
    _deltaRadSqr = _deltaRad * _deltaRad;
  }


  //! Evaluating the function itself, i.e. the regularization of
  //! \f$ \gamma(z) = r_{\delta} \left\{ 1 + \sqrt{\sin^2(\frac52(\phi + \frac{\pi}{2})) + \delta^2 } \right\}  \f$.
  RealType gammaNorm ( const aol::Vec2<RealType> &z ) const {
    double r   = sqrt ( z.normSqr() + _deltaRadSqr );
    double phi = asin ( z[1] / r );    //! in [-pi/2, pi/2]

    double CSin = aol::Sqr ( sin ( 2.5 * ( phi + 1.570796326499999999979674536465523715378 ) ) );
    return r * ( 1. + sqrt ( CSin + _deltaAbsSqr ) );
  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    double r   = sqrt ( z.normSqr() + _deltaRadSqr );
    double phi = asin ( z[1] / r );    //! in [-pi/2, pi/2]

    double CSin = sin ( 2.5 * ( phi + 1.570796326499999999979674536465523715378 ) );
    double CCos = cos ( 2.5 * ( phi + 1.570796326499999999979674536465523715378 ) );

    double root  = sqrt ( CSin * CSin + _deltaAbsSqr );
    double term1 = ( 1. + root ) / r;
    double term2 = CSin * CCos / ( root * r );

    v[0] = z[0] * term1 - 2.5 * term2 * z[0] * z[1] / sqrt ( r * r - z[1] * z[1] );
    v[1] = z[1] * term1 + 2.5 * term2 * sqrt ( r * r - z[1] * z[1] );
  }
};


//! Class for a rotated pentagon-like Wulffshape in 3D.
template <typename RealType>
class Pentagon3d : public Anisotropy3dGraphInterface<RealType, Pentagon3d<RealType> > {
public:
  using Anisotropy3dGraphInterface<RealType, Pentagon3d<RealType> >::gammaNorm;
  using Anisotropy3dGraphInterface<RealType, Pentagon3d<RealType> >::gammaFirstDerivative;

private:
  RealType _deltaAbs;     //! Regularization-parameter for the absolute-value
  RealType _deltaAbsSqr;  //! Square of Delta (very complicated...)
  RealType _deltaRad;     //! Regularization-parameter for the radius
  RealType _deltaRadSqr;

public:
  Pentagon3d ( RealType deltaAbs, RealType deltaRad ) :
      Anisotropy3dGraphInterface<RealType, Pentagon3d<RealType> > ( ),
      _deltaAbs ( deltaAbs ), _deltaRad ( deltaRad ) {
    _deltaAbsSqr = _deltaAbs * _deltaAbs;
    _deltaRadSqr = _deltaRad * _deltaRad;
  }

  //! Evaluating the function itself, i.e. the regularization of
  //! \f$ \gamma(z) = r_{\delta} \left\{ 1 + \sqrt{\sin^2(\frac52(\phi + \frac{\pi}{2})) + \delta^2 } \right\}  \f$.
  RealType gammaNorm ( const aol::Vec3<RealType> &z ) const {
    RealType r   = sqrt ( z.normSqr() + _deltaRadSqr );
    RealType phi = asin ( z[2] / r );    //! in [-pi/2, pi/2]

    RealType CSin = aol::Sqr ( sin ( 2.5 * ( phi + aol::NumberTrait<RealType>::pi / 2. ) ) );
    return r * ( 1. + sqrt ( CSin + _deltaAbsSqr ) );
  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec3<RealType> &, const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    RealType r   = sqrt ( z.normSqr() + _deltaRadSqr );
    RealType phi = asin ( z[2] / r );    //! in [-pi/2, pi/2]

    RealType CSin = sin ( 2.5 * ( phi + aol::NumberTrait<RealType>::pi / 2. ) );
    RealType CCos = cos ( 2.5 * ( phi + aol::NumberTrait<RealType>::pi / 2. ) );

    RealType root1 = sqrt ( CSin * CSin + _deltaAbsSqr );
    RealType root2 = sqrt ( r * r - z[2] * z[2] );
    RealType term1 = ( 1. + root1 ) / r;
    RealType term2 = CSin * CCos / ( root1 * r );

    v[0] = z[0] * term1 - 2.5 * term2 * z[0] * z[2] / root2;
    v[1] = z[1] * term1 - 2.5 * term2 * z[1] * z[2] / root2;
    v[2] = z[2] * term1 + 2.5 * term2 * root2;
  }
};


//! Class for a rotated polygonal Wulffshape in 3D generated by the function
//! \f$ \gamma(r,\varphi,\vartheta) = r ( a + b| \sin(c(\varphi+d)) | ) \f$.
//! Regularized version: \f$ \gamma(r,\varphi,\vartheta) = r ( a + b \sqrt( \sin^2(c(\varphi+d)) + \delta^2) ) \f$
//! Special cases:
//! \f$ a=1, b=3, c=\frac{3}{2}, d=\frac{\pi}{2}: a rotated triangle                 \f$
//! \f$ a=1, b=3, c=2,           d=0:             a rotated square (diamond like)    \f$
//! \f$ a=1, b=1, c=\frac{5}{2}, d=\frac{\pi}{2}: a rotated pentagon                 \f$
//! \f$ a=1, b=1, c=3,           d=\frac{\pi}{2}: a rotated hexagon                  \f$
//! \f$ a=1, b=1, c=\frac{7}{2}, d=\frac{\pi}{2}: a rotated heptagon                 \f$
//! \f$ a=1, b=1, c=4,           d=0:             a rotated octagon                  \f$
template <typename RealType>
class SinePolygonRotated3d : public Anisotropy3dGraphInterface<RealType, SinePolygonRotated3d<RealType> > {
public:
  using Anisotropy3dGraphInterface<RealType, SinePolygonRotated3d<RealType> >::gammaNorm;
  using Anisotropy3dGraphInterface<RealType, SinePolygonRotated3d<RealType> >::gammaFirstDerivative;

private:
  const RealType _a;            //! the parameters for specifying the polygon
  const RealType _b;
  const RealType _c;
  const RealType _d;

  const RealType _deltaAbs;     //! Regularization-parameter for the absolute-value
  const RealType _deltaRad;     //! Regularization-parameter for the radius
  RealType _deltaAbsSqr;        //! Square of DeltaAbs (very complicated...)
  RealType _deltaRadSqr;        //! uh, don't know ...

public:
  SinePolygonRotated3d ( RealType a, RealType b, RealType c, RealType d, RealType deltaAbs, RealType deltaRad ) :
      Anisotropy3dGraphInterface<RealType, SinePolygonRotated3d<RealType> > ( ),
      _a ( a ), _b ( b ), _c ( c ), _d ( d ),
      _deltaAbs ( deltaAbs ), _deltaRad ( deltaRad ) {
    _deltaAbsSqr = _deltaAbs * _deltaAbs;
    _deltaRadSqr = _deltaRad * _deltaRad;
  }

  //! Evaluating the function itself, i.e. the regularization of
  //! \f$ \gamma(z) = r_{\delta} \left\{ a + b \sqrt{\sin^2( c (\phi + d)) + \delta^2 } \right\}  \f$.
  RealType gammaNorm ( const aol::Vec3<RealType> &z ) const {
    RealType r   = sqrt ( z.normSqr() + _deltaRadSqr );
    RealType phi = asin ( z[2] / r );    //! in [-pi/2, pi/2]

    RealType CSin = aol::Sqr ( sin ( _c * ( phi + _d ) ) );
    return r * ( _a + _b*sqrt ( CSin + _deltaAbsSqr ) );
  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec3<RealType> &, const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  //! evaluating the first derivative: \f$ ... \f$
  void gammaFirstDerivative ( const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    RealType r   = sqrt ( z.normSqr() + _deltaRadSqr );
    RealType phi = asin ( z[2] / r );    //! in [-pi/2, pi/2]

    RealType CSin = sin ( _c * ( phi + _d ) );
    RealType CCos = cos ( _c * ( phi + _d ) );

    RealType root1 = sqrt ( CSin * CSin + _deltaAbsSqr );
    RealType root2 = sqrt ( r * r - z[2] * z[2] );
    RealType term1 = ( _a + _b * root1 ) / r;
    RealType term2 = CSin * CCos / ( root1 * r );

    v[0] = z[0] * term1 - _b * _c * term2 * z[0] * z[2] / root2;
    v[1] = z[1] * term1 - _b * _c * term2 * z[1] * z[2] / root2;
    v[2] = z[2] * term1 + _b * _c * term2 * root2;
  }
};





/* *******************************************************************************
 * gamma = Lp-Norm in 3d (rotated)
 * ******************************************************************************* */

template <typename RealType>
class RotatedLpAnisotropy3d {
private:
  RealType _p, _q, _epsilon;
  aol::Vec3<RealType> _v;        // direction, on which the x-axxis lies
  aol::Matrix33<RealType> _B;    // matrix to perform a basis-transformation
  aol::Matrix33<RealType> _BT;   // the transposed transformation-matrix
public:
  RotatedLpAnisotropy3d ( RealType P, RealType Q, RealType Epsilon ) : _p ( P ), _q ( Q ), _epsilon ( Epsilon ) {
    if ( Q < 1 ) cerr << aol::color::red << "\n\nWARNING: Q is smaller 1, there might be singularities!!\n";
    _v[0] = 1.; _v[1] = 0.; _v[2] = 0.;
    calcB();
  }

  RotatedLpAnisotropy3d ( RealType P, RealType Q, RealType Epsilon, aol::Vec3<RealType> V ) :
      _p ( P ), _q ( Q ), _epsilon ( Epsilon ),  _v ( V ) {
    if ( Q < 1 ) cerr << aol::color::red << "\n\nWARNING: Q is smaller 1, there might be singularities!!\n";
    calcB();
  }


  // calc the matrix B for transforming the coordinate-system for the x-axxis
  // lying on v
  void calcB() {
    // the first row is the direction itself
    _B.setRow ( 0, _v );
    // the other directions must be orthogonal to _v, but then, they can be
    // chosen arbitrary, first choose y orthogonal to _v
    // TODO: arranging by size to avoid taking the smallest values (=> stability)
    aol::Vec3<RealType> y ( 0, 0, 0 );
    if ( _v[1] != 0 ) {
      if ( _v[2] != 0 ) {
        y[1] = -_v[2];
        y[2] = _v[1];
        y /= y.norm();
      } else y[2] = 1;
    } else y[1] = 1;

    _B.setRow ( 1, y );

    // now z = vector-product of v and y, both normed => z is normed too
    aol::Vec3<RealType> z ( _v.crossProduct ( y ) );
    _B.setRow ( 2, z );

    // finally save the transposed matrix
    _BT = _B;
    _BT.transpose();
  }


  void printB() {
    cerr << "Drehmatrix: \n";
    cerr << _BT;
  }

  // set the rotate-direction v
  void setRotateDirection ( const aol::Vec3<RealType> &v ) {
    // only anisotropic, if v is not 0
    // if v is not 0, it is very improbable, that v is equal to _v
    // (that means = v from the timestep before)
    // => it's not worthwile to check this

    if ( v[0] != 0 || v[1] != 0 || v[2] != 0 ) {
      _v = v / v.norm();
      calcB();
    } else cerr << "\n\n'anisotropy:************* VECTOR FOR ROTATING THE ANISOTROPY IS 0 !!!! ************** \n\n";
  }

  // set the rotate-direction v
  void setRotateDirection ( RealType v0, RealType v1, RealType v2 ) {
    if ( v0 != 0 || v1 != 0 || v2 != 0 ) {
      _v[0] = v0; _v[1] = v1; _v[2] = v2;
      _v /= _v.norm();
      calcB();
    } else cerr << "\n\n'anisotropy:************* VECTOR FOR ROTATING THE ANISOTROPY IS 0 !!!! ************** \n\n";
  }



RealType p() const { return _p; }


  // evaluating the function itself
  RealType gamma ( const aol::Vec3<RealType> &arg ) const {
    // first apply B on the argument => rotation
    aol::Vec3<RealType> z;
    _B.mult ( arg, z );

    return pow ( pow ( aol::Abs ( z[0] ), _q ) + pow ( aol::Abs ( z[1] ), _q )
                 + pow ( aol::Abs ( z[2] ), _q ) + pow ( _epsilon, _q ),  1 / _q );
  }


  void gammaFirstDerivative ( const aol::Vec3<RealType> &, const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  void gammaFirstDerivative ( const aol::Vec3<RealType> &arg, aol::Vec3<RealType> &v ) const {
    // first apply B on the argument => rotation
    aol::Vec3<RealType> z;
    aol::Vec3<RealType> temp;
    _B.mult ( arg, z );

    aol::Vec3<RealType> Z ( aol::Abs ( z[0] ), aol::Abs ( z[1] ), aol::Abs ( z[2] ) );
    RealType a = pow ( Z[0], _q ) + pow ( Z[1], _q ) + pow ( Z[2], _q ) + pow ( _epsilon, _q );

    // It has to be q>1, because for z[i]=0 there appears a singularity
    temp[0] = pow ( Z[0], _q - 1. ) * aol::signum1at0 ( z[0] );
    temp[1] = pow ( Z[1], _q - 1. ) * aol::signum1at0 ( z[1] );
    temp[2] = pow ( Z[2], _q - 1. ) * aol::signum1at0 ( z[2] );
    temp *= pow ( a, _p / _q - 1. ) *  pow ( pow ( a, _p / _q ) + pow ( _epsilon, _p ), ( 1. - _p ) / _p );

    // due to the chain rule we have to multiply the rotation-matrix once again
//     cerr << "gamma_z( " << z << ") = " << temp << endl;
    _BT.mult ( temp, v );
  }

  RealType gammaNorm ( const aol::Vec3<RealType> &arg ) const {
    // first apply B on the argument => rotation
    aol::Vec3<RealType> z;
    _B.mult ( arg, z );

    return pow ( pow ( aol::Abs ( z[0] ), _q ) + pow ( aol::Abs ( z[1] ), _q )
                 + pow ( aol::Abs ( z[2] ), _q ) + pow ( _epsilon, _q ), 1. / _q );
  }
};




template <typename RealType>
class RegMaxAnisotropy {
private:
  RealType _epsilon;
public:
  RegMaxAnisotropy ( RealType eps ) : _epsilon ( eps ) { }

  RealType p() const { return this->_p; }

  void implicitPart ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z, aol::Matrix22<RealType> &mat ) const {
    implicitPart ( z, mat );
  }

  void explicitPart ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z, aol::Matrix22<RealType> &mat ) const {
    explicitPart ( z, mat );
  }

  void implicitPart ( const aol::Vec2<RealType> &/*z*/, aol::Matrix22<RealType> &/*mat*/ ) const {
    cerr << "implicit part";
  }

  void explicitPart ( const aol::Vec2<RealType> &/*z*/, aol::Matrix22<RealType> &/*mat*/ ) const {
    cerr << "explicit part";
  }


  void gammaFirstDerivative ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  void gammaFirstDerivative ( const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    RealType denom = ( 1. + aol::Sqr ( _epsilon ) ) * ( aol::Sqr ( z[0] ) + aol::Sqr ( z[1] ) );
    v[0] = ( ( 1. + aol::Sqr ( _epsilon ) ) * z[0] - z[1] ) / sqrt ( denom - 2. * z[0] * z[1] )
           + ( ( 1. + aol::Sqr ( _epsilon ) ) * z[0] + z[1] ) / sqrt ( denom + 2. * z[0] * z[1] );

    v[1] = ( ( 1. + aol::Sqr ( _epsilon ) ) * z[1] - z[0] ) / sqrt ( denom - 2. * z[0] * z[1] )
           + ( ( 1. + aol::Sqr ( _epsilon ) ) * z[1] + z[0] ) / sqrt ( denom + 2. * z[0] * z[1] );
    //cerr << "gamma_z( " << z << ") = " << v << endl;
  }

  RealType gammaNorm ( const aol::Vec2<RealType> &z ) const {
    return sqrt ( aol::Sqr ( z[0] - z[1] ) + aol::Sqr ( _epsilon ) * z.normSqr() ) +
           sqrt ( aol::Sqr ( z[0] + z[1] ) + aol::Sqr ( _epsilon ) * z.normSqr() );

  }
};


template <typename RealType>
class Isotropy: public Anisotropy3dGraphInterface<RealType, Isotropy<RealType> > {
public:
  using Anisotropy3dGraphInterface<RealType, Isotropy<RealType> >::gammaNorm;
  using Anisotropy3dGraphInterface<RealType, Isotropy<RealType> >::gammaFirstDerivative;
  using Anisotropy3dGraphInterface<RealType, Isotropy<RealType> >::implicitPart;

private:
  RealType _epsilon;

public:
  Isotropy ( RealType eps ) : Anisotropy3dGraphInterface<RealType, Isotropy<RealType> > ( ), _epsilon ( eps ) { }

  void implicitPart ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z, aol::Matrix22<RealType> &mat ) const {
    implicitPart ( z, mat );
  }

  void explicitPart ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z, aol::Matrix22<RealType> &mat ) const {
    explicitPart ( z, mat );
  }

  void implicitPart ( const aol::Vec2<RealType> &z, aol::Matrix22<RealType> &mat ) const {
    mat[0][0] = mat[1][1] = 1.;
    mat[1][0] = mat[0][1] = 0.;
    mat /= sqrt ( aol::Sqr ( z[0] ) + aol::Sqr ( z[1] ) + aol::Sqr ( _epsilon ) ) ;
  }

  void explicitPart ( const aol::Vec2<RealType> &z, aol::Matrix22<RealType> &mat ) const {
    for ( int i = 0; i < 2; i++ ) for ( int j = 0; j < 2; j++ ) mat[i][j] = z[i] * z[j];
    mat /= pow ( sqrt ( aol::Sqr ( z[0] ) + aol::Sqr ( z[1] ) + aol::Sqr ( _epsilon ) ), 3. ) ;
  }


  void gammaFirstDerivative ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  void gammaFirstDerivative ( const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    v = z;
    v /= gammaNorm ( z );
  }

  RealType gammaNorm ( const aol::Vec2<RealType> &z ) const {
    return sqrt ( aol::Sqr ( z[0] ) + aol::Sqr ( z[1] ) + aol::Sqr ( _epsilon ) ) ;
  }
};


template <typename RealType>
class Isotropy3d : public Anisotropy3dGraphInterface<RealType, Isotropy3d<RealType> >{
public:
  using Anisotropy3dGraphInterface<RealType, Isotropy3d<RealType> >::gammaNorm;
  using Anisotropy3dGraphInterface<RealType, Isotropy3d<RealType> >::gammaFirstDerivative;
  using Anisotropy3dGraphInterface<RealType, Isotropy3d<RealType> >::implicitPart;
  using Anisotropy3dGraphInterface<RealType, Isotropy3d<RealType> >::explicitPart;

private:
  RealType _epsilon;
public:
  Isotropy3d ( RealType eps ) :
    Anisotropy3dGraphInterface<RealType, Isotropy3d<RealType> > ( ), _epsilon ( eps ) { }

  void implicitPart ( const aol::Vec3<RealType> &, const aol::Vec3<RealType> &z, aol::Matrix33<RealType> &mat ) const {
    implicitPart ( z, mat );
  }

  void explicitPart ( const aol::Vec3<RealType> &, const aol::Vec3<RealType> &z, aol::Matrix33<RealType> &mat ) const {
    explicitPart ( z, mat );
  }

  void implicitPart ( const aol::Vec3<RealType> &z, aol::Matrix33<RealType> &mat ) const {
    mat.setIdentity();
    mat /= sqrt ( aol::Sqr ( z[0] ) + aol::Sqr ( z[1] ) + aol::Sqr ( z[2] ) + aol::Sqr ( _epsilon ) ) ;
  }

  void explicitPart ( const aol::Vec3<RealType> &z, aol::Matrix33<RealType> &mat ) const {
    for ( int i = 0; i < 3; i++ ) for ( int j = 0; j < 3; j++ ) mat[i][j] = z[i] * z[j];
    mat /= pow ( sqrt ( aol::Sqr ( z[0] ) + aol::Sqr ( z[1] ) + aol::Sqr ( z[2] ) + aol::Sqr ( _epsilon ) ), 3. ) ;
  }


  void gammaFirstDerivative ( const aol::Vec3<RealType> &, const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  void gammaFirstDerivative ( const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    v = z;
    v /= gammaNorm ( z );
  }

  RealType gammaNorm ( const aol::Vec3<RealType> &z ) const {
    return sqrt ( aol::Sqr ( z[0] ) + aol::Sqr ( z[1] ) + aol::Sqr ( z[2] ) + aol::Sqr ( _epsilon ) ) ;
  }
};


template <typename RealType>
class EllipseAnisotropy {
private:
  RealType _a, _b, _epsilon;
public:
  EllipseAnisotropy ( RealType a, RealType b, RealType eps ) : _a ( a ), _b ( b ), _epsilon ( eps ) { }


  void implicitPart ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z, aol::Matrix22<RealType> &mat ) const {
    implicitPart ( z, mat );
  }

  void explicitPart ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z, aol::Matrix22<RealType> &mat ) const {
    explicitPart ( z, mat );
  }

  void gammaFirstDerivative ( const aol::Vec2<RealType> &, const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }


  void implicitPart ( const aol::Vec2<RealType> &z, aol::Matrix22<RealType> &mat ) const {
    mat[0][0] = aol::Sqr ( _a );
    mat[1][1] = aol::Sqr ( _b );
    mat[0][1] = mat[1][0] = 0.;
    RealType normSqr = aol::Sqr ( z[0] * _a ) + aol::Sqr ( z[1] * _b ) + _epsilon * _epsilon;
    mat /= sqrt ( normSqr );
  }

  void explicitPart ( const aol::Vec2<RealType> &z, aol::Matrix22<RealType> &mat ) const {
    aol::Vec2<RealType> v;
    gammaFirstDerivative ( z, v );
    mat[0][0] = v[0] * v[0];
    mat[1][1] = v[1] * v[1];
    mat[0][1] =  mat[1][0] = v[1] * v[0];
    RealType normSqr = aol::Sqr ( z[0] * _a ) + aol::Sqr ( z[1] * _b ) + _epsilon * _epsilon;
    mat /= sqrt ( normSqr );
  }

  void gammaFirstDerivative ( const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    v[0] = aol::Sqr ( _a ) * z[0];
    v[1] = aol::Sqr ( _b ) * z[1];
    RealType normSqr = aol::Sqr ( z[0] * _a ) + aol::Sqr ( z[1] * _b ) + _epsilon * _epsilon;
    v /= sqrt ( normSqr );
  }

  RealType gammaNorm ( const aol::Vec2<RealType> &z ) const {
    return sqrt ( aol::Sqr ( z[0]*_a ) + aol::Sqr ( z[1] * _b ) + _epsilon*_epsilon );
  }

};


// ------------------------ ELLIPSOID ----------------------------------------------

template <typename RealType>
class EllipsoidAnisotropy : public Anisotropy3dGraphInterface<RealType, EllipsoidAnisotropy<RealType> > {
public:
  using Anisotropy3dGraphInterface<RealType, EllipsoidAnisotropy<RealType> >::gammaNorm;
  using Anisotropy3dGraphInterface<RealType, EllipsoidAnisotropy<RealType> >::gammaFirstDerivative;

private:
  RealType _a, _b, _c, _epsilon;
public:
  EllipsoidAnisotropy ( RealType a, RealType b, RealType c, RealType eps ) :
    Anisotropy3dGraphInterface<RealType, EllipsoidAnisotropy<RealType> > ( ),
      _a ( a ), _b ( b ), _c ( c ), _epsilon ( eps ) { }

  // evaluating the function itself
  RealType gamma ( const aol::Vec3<RealType> &z ) const {
    return sqrt ( aol::Sqr ( z[0]*_a ) + aol::Sqr ( z[1]*_b )
                  + aol::Sqr ( z[2]*_c ) + _epsilon*_epsilon );
  }

  void gammaFirstDerivative ( const aol::Vec3<RealType> &, const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  // evaluating the first derivative
  void gammaFirstDerivative ( const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    v[0] = aol::Sqr ( _a ) * z[0];
    v[1] = aol::Sqr ( _b ) * z[1];
    v[2] = aol::Sqr ( _c ) * z[2];
    RealType normSqr = aol::Sqr ( z[0] * _a ) + aol::Sqr ( z[1] * _b )
                       + aol::Sqr ( z[2] * _c ) + _epsilon * _epsilon;
    v /= sqrt ( normSqr );
  }

  RealType gammaNorm ( const aol::Vec3<RealType> &z ) const {
    return sqrt ( aol::Sqr ( z[0]*_a ) + aol::Sqr ( z[1] * _b )
                  + aol::Sqr ( z[2] * _c ) + _epsilon*_epsilon );
  }

};

// -----------------------------------------------------------------------------------------------


// ------------------------ ROTIERTER ELLIPSOID ----------------------------------------------
// diese Klasse stellt einen Ellipsoid zur Verfuegung, der so rotiert ist,
// dass die x-Achse auf einem gegebenen Vektor v liegt. Die anderen beiden
// Koordinatenachsen sind aufgrund der Symmetrie des Ellipsoids willkuerlich
// gewaehlt.
// -------------------------------------------------------------------------------------------

template <typename RealType>
class RotatedEllipsoidAnisotropy {
private:
  RealType _a, _b, _c, _epsilon;
  aol::Vec3<RealType> _v;        // direction, on which the x-axxis lies
  aol::Matrix33<RealType> _B;    // matrix to perform a basis-transformation
  aol::Matrix33<RealType> _BT;   // the transposed transformation-matrix
  bool _isotrop;                // if v is 0, then calc like isotropic mcm => flag
public:
  RotatedEllipsoidAnisotropy ( RealType a, RealType b, RealType c, RealType eps ) :
      _a ( a ), _b ( b ), _c ( c ), _epsilon ( eps ), _v ( 1., 0., 0. ) {
    calcB();
    _isotrop = false;
  }

  RotatedEllipsoidAnisotropy ( RealType a, RealType b, RealType c, RealType eps, aol::Vec3<RealType> v ) :
      _a ( a ), _b ( b ), _c ( c ), _epsilon ( eps ),  _v ( v ) {
    calcB();
    _isotrop = false;
  }


  // calc the matrix B for transforming the coordinate-system for the x-axxis
  // lying on v
  void calcB() {
    // the first row is the direction itself
    _B.setRow ( 0, _v );
    // the other directions must be orthogonal to _v, but then, they can be
    // chosen arbitrary, first choose y orthogonal to _v
    // TODO: arranging by size to avoid taking the smallest values (=> stability)
    aol::Vec3<RealType> y ( 0, 0, 0 );
    if ( _v[1] != 0 ) {
      if ( _v[2] != 0 ) {
        y[1] = -_v[2];
        y[2] = _v[1];
        y /= y.norm();
      } else y[2] = 1;
    } else y[1] = 1;

    _B.setRow ( 1, y );

    // now z = vector-product of v and y, both normed => z is normed too
    aol::Vec3<RealType> z ( _v.crossProduct ( y ) );
    _B.setRow ( 2, z );

    // finally save the transposed matrix
    _BT = _B;
    _BT.transpose();
  }


  // set the rotate-direction v
  void setRotateDirection ( const aol::Vec3<RealType> &v ) {
    // only anisotropic, if v is not 0
    // if v is not 0, it is very improbable, that v is equal to _v
    // (that means = v from the timestep before)
    // => it's not worthwile to check this
    if ( v[0] != 0 || v[1] != 0 || v[2] != 0 ) {
      _v = v / v.norm();
      calcB();
      _isotrop = false;
    } else _isotrop = true;
  }

  // set the rotate-direction v
  void setRotateDirection ( RealType v0, RealType v1, RealType v2 ) {
    if ( v0 != 0 || v1 != 0 || v2 != 0 ) {
      _v[0] = v0; _v[1] = v1; _v[2] = v2;
      _v /= _v.norm();
      calcB();
      _isotrop = false;
    } else _isotrop = true;
  }



  // evaluating the function itself
  RealType gamma ( const aol::Vec3<RealType> &arg ) const {
    if ( !_isotrop ) {
      // first apply B on the argument => rotation
      aol::Vec3<RealType> _z;
      _B.mult ( arg, _z );

      return sqrt ( aol::Sqr ( _z[0]*_a ) + aol::Sqr ( _z[1]*_b )
                    + aol::Sqr ( _z[2]*_c ) + _epsilon*_epsilon );
    } else return arg.norm();
//     1.;    // |z|
  }

  void gammaFirstDerivative ( const aol::Vec3<RealType> &, const aol::Vec3<RealType> &z, aol::Vec3<RealType> &v ) const {
    gammaFirstDerivative ( z, v );
  }

  // evaluating the first derivative
  void gammaFirstDerivative ( const aol::Vec3<RealType> &arg, aol::Vec3<RealType> &v ) const {
    if ( !_isotrop ) {
      aol::Vec3<RealType> _z, temp;
      _B.mult ( arg, _z );
      temp[0] = aol::Sqr ( _a ) * _z[0];
      temp[1] = aol::Sqr ( _b ) * _z[1];
      temp[2] = aol::Sqr ( _c ) * _z[2];
      RealType normSqr = aol::Sqr ( _z[0] * _a ) + aol::Sqr ( _z[1] * _b )
                         + aol::Sqr ( _z[2] * _c ) + _epsilon * _epsilon;
      temp /= sqrt ( normSqr );

      // due to the chain rule we have to multiply the rotation-matrix once again
      _BT.mult ( temp, v );
    } else {
      v = arg;
      v /= arg.norm();
    }
  }

  RealType gammaNorm ( const aol::Vec3<RealType> &z ) const {
    aol::Vec3<RealType> _z;
    _B.mult ( z, _z );
    return sqrt ( aol::Sqr ( _z[0]*_a ) + aol::Sqr ( _z[1] * _b )
                  + aol::Sqr ( _z[2] * _c ) + _epsilon*_epsilon );
  }

};

// -----------------------------------------------------------------------------------------------




template <typename RealType>
class MixedAnisotropy {
  LpAnisotropy<RealType>      _a11, _a12;
  EllipseAnisotropy<RealType> _a21, _a22;
public:
  MixedAnisotropy ( RealType p, RealType a, RealType b, RealType epsilon ) :
      _a11 ( 3., p, epsilon ), _a12 ( p, p, epsilon ),
      _a21 ( a, b, epsilon ), _a22 ( 1.0, 1.0, epsilon )  {}

  void implicitPart ( const aol::Vec2<RealType> &pt, const aol::Vec2<RealType> &z, aol::Matrix22<RealType> &mat ) const {
    _a12.implicitPart ( z, mat );
    return;

    if ( pt[0] < 0.5 ) {
      if ( pt[1] < 0.5 ) {
        _a11.implicitPart ( z, mat );
      } else {
        _a12.implicitPart ( z, mat );
      }
    } else {
      if ( pt[1] < 0.5 ) {
        _a21.implicitPart ( z, mat );
      } else {
        _a22.implicitPart ( z, mat );
      }
    }
    return;

    aol::Matrix22<RealType> mtmp, munten, moben;

    if ( pt[0] < 0.45 ) {
      _a11.implicitPart ( z, munten );
    } else if ( pt[0] > 0.55 ) {
      _a12.implicitPart ( z, munten );
    } else {
      RealType lambda = ( pt[0] - 0.45 ) / 0.1;
      _a11.implicitPart ( z, mtmp );
      mtmp *= ( 1. - lambda );
      _a12.implicitPart ( z, munten );
      munten *= lambda;
      munten += mtmp;
    }
    if ( pt[0] < 0.45 ) {
      _a21.implicitPart ( z, moben );
    } else if ( pt[0] > 0.55 ) {
      _a22.implicitPart ( z, moben );
    } else {
      RealType lambda = ( pt[0] - 0.45 ) / 0.1;
      _a21.implicitPart ( z, mtmp );
      mtmp *= ( 1. - lambda );
      _a22.implicitPart ( z, moben );
      moben *= lambda;
      moben += mtmp;
    }
    if ( pt[1] < 0.45 ) {
      mat = munten;
      return;
    }
    if ( pt[1] > 0.55 ) {
      mat = moben;
      return;
    }
    RealType lambda = ( pt[1] - 0.45 ) / 0.1;
    munten *= ( 1. - lambda );
    moben *= lambda;
    mat = munten;
    mat += moben;
  }

  void explicitPart ( const aol::Vec2<RealType> &pt, const aol::Vec2<RealType> &z, aol::Matrix22<RealType> &mat ) const {
    _a12.explicitPart ( z, mat );
    return;

    if ( pt[0] < 0.5 ) {
      if ( pt[1] < 0.5 ) {
        _a11.explicitPart ( z, mat );
      } else {
        _a12.explicitPart ( z, mat );
      }
    } else {
      if ( pt[1] < 0.5 ) {
        _a21.explicitPart ( z, mat );
      } else {
        _a22.explicitPart ( z, mat );
      }
    }
    return;

    aol::Matrix22<RealType> mtmp, munten, moben;

    if ( pt[0] < 0.45 ) {
      _a11.explicitPart ( z, munten );
    } else if ( pt[0] > 0.55 ) {
      _a12.explicitPart ( z, munten );
    } else {
      RealType lambda = ( pt[0] - 0.45 ) / 0.1;
      _a11.explicitPart ( z, mtmp );
      mtmp *= ( 1. - lambda );
      _a12.explicitPart ( z, munten );
      munten *= lambda;
      munten += mtmp;
    }
    if ( pt[0] < 0.45 ) {
      _a21.explicitPart ( z, moben );
    } else if ( pt[0] > 0.55 ) {
      _a22.explicitPart ( z, moben );
    } else {
      RealType lambda = ( pt[0] - 0.45 ) / 0.1;
      _a21.explicitPart ( z, mtmp );
      mtmp *= ( 1. - lambda );
      _a22.explicitPart ( z, moben );
      moben *= lambda;
      moben += mtmp;
    }
    if ( pt[1] < 0.45 ) {
      mat = munten;
      return;
    }
    if ( pt[1] > 0.55 ) {
      mat = moben;
      return;
    }
    RealType lambda = ( pt[1] - 0.45 ) / 0.1;
    munten *= ( 1. - lambda );
    moben *= lambda;
    mat = munten;
    mat += moben;
  }

  void gammaFirstDerivative ( const aol::Vec2<RealType> &pt, const aol::Vec2<RealType> &z, aol::Vec2<RealType> &v ) const {
    _a12.gammaFirstDerivative ( z, v );
    return;

    if ( pt[0] < 0.5 ) {
      if ( pt[1] < 0.5 ) {
        _a11.gammaFirstDerivative ( z, v );
      } else {
        _a12.gammaFirstDerivative ( z, v );
      }
    } else {
      if ( pt[1] < 0.5 ) {
        _a21.gammaFirstDerivative ( z, v );
      } else {
        _a22.gammaFirstDerivative ( z, v );
      }
    }
    return;

    aol::Vec2<RealType> mtmp, munten, moben;

    if ( pt[0] < 0.45 ) {
      _a11.gammaFirstDerivative ( z, munten );
    } else if ( pt[0] > 0.55 ) {
      _a12.gammaFirstDerivative ( z, munten );
    } else {
      RealType lambda = ( pt[0] - 0.45 ) / 0.1;
      _a11.gammaFirstDerivative ( z, mtmp );
      mtmp *= ( 1. - lambda );
      _a12.gammaFirstDerivative ( z, munten );
      munten *= lambda;
      munten += mtmp;
    }
    if ( pt[0] < 0.45 ) {
      _a21.gammaFirstDerivative ( z, moben );
    } else if ( pt[0] > 0.55 ) {
      _a22.gammaFirstDerivative ( z, moben );
    } else {
      RealType lambda = ( pt[0] - 0.45 ) / 0.1;
      _a21.gammaFirstDerivative ( z, mtmp );
      mtmp *= ( 1. - lambda );
      _a22.gammaFirstDerivative ( z, moben );
      moben *= lambda;
      moben += mtmp;
    }
    if ( pt[1] < 0.45 ) {
      v = munten;
      return;
    }
    if ( pt[1] > 0.55 ) {
      v = moben;
      return;
    }
    RealType lambda = ( pt[1] - 0.45 ) / 0.1;
    munten *= ( 1. - lambda );
    moben *= lambda;
    v = munten;
    v += moben;
  }

  RealType gammaNorm ( const aol::Vec2<RealType> &z ) const {
    return _a12.gammaNorm ( z );
  }

};

} // namespace qc

#endif

