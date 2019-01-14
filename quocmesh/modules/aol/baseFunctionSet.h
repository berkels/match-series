#ifndef __BASEFUNCTIONSET_H
#define __BASEFUNCTIONSET_H

#include <aol.h>
#include <quoc.h>
#include <vec.h>
#include <smallMat.h>

namespace aol {

//! provides basis functions, derivatives and cached values for a given
//! quadrature rule and dimension.
/*!
 * \author Droske
 */
template <typename RealType, typename VecType, typename DomVecType, int NumBaseFuncs, class QuadRuleType, typename Imp>
class BaseFunctionSetInterface  {
public:
  BaseFunctionSetInterface( ) {}

  enum { numBaseFuncs = NumBaseFuncs };

  int numQuadPoints( ) const {
    return QuadRuleType::numQuadPoints;
  }

  inline RealType getWeight ( int QuadPoint ) const {
    return _quadRule.getWeight ( QuadPoint );
  }

  inline const DomVecType& getRefCoord ( int QuadPoint ) const {
    return _quadRule.getRefCoord ( QuadPoint );
  }

  //! read the cached value of the basis function with number BaseFuncNum at the given
  //! quadrature point
  inline RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return basisQuadValues[BaseFuncNum][QuadPoint];
  }

  inline const VecType& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return basisQuadGradients[BaseFuncNum][QuadPoint];
  }

  inline RealType evaluate ( int BaseFuncNum, const DomVecType &RefCoord ) const {
    return asImp().evaluate ( BaseFuncNum, RefCoord );
  }

  inline void evaluateGradient ( int BaseFuncNum, const DomVecType& RefCoord, VecType& Gradient ) const {
    asImp().evaluateGradient ( BaseFuncNum, RefCoord, Gradient );
  }

protected:
  Imp &asImp() { return static_cast<Imp&> ( *this ); }

  const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

  void initializeQuadCache( ) {
    for ( int b = 0; b < numBaseFuncs; b++ ) {
      for ( int i = 0; i < QuadRuleType::numQuadPoints; i++ ) {
        basisQuadValues[b][i] = evaluate ( b,  _quadRule.getRefCoord ( i ) );
        evaluateGradient ( b, _quadRule.getRefCoord ( i ), basisQuadGradients[b][i] );
      }
    }
  }

  /**** cache the values of the basis functions at the quadrature points ****/
  RealType         basisQuadValues   [numBaseFuncs][QuadRuleType::numQuadPoints];
  VecType          basisQuadGradients[numBaseFuncs][QuadRuleType::numQuadPoints];
  QuadRuleType _quadRule;
};


template <typename RealType, qc::Dimension Dim, class QuadRuleType>
class BaseFunctionSetMultiLin {
};

//! The basefunctionset for bilinear elements in 1d.
template <typename RealType, class QuadRuleType>
class BaseFunctionSetMultiLin<RealType, qc::QC_1D, QuadRuleType> :
      public BaseFunctionSetInterface<RealType, Vec<1, RealType>, Vec<1, RealType>, 2, QuadRuleType, BaseFunctionSetMultiLin<RealType, qc::QC_1D, QuadRuleType> > {

  static RealType _dx_b1   ( const Vec<1, RealType> &/*RefCoord*/ ) { return - 1.; }
  static RealType _dx_b2   ( const Vec<1, RealType> &/*RefCoord*/ ) { return 1.; }

  static RealType _b1   ( const Vec<1, RealType> &RefCoord ) { return ( 1. - RefCoord[0] ); }
  static RealType _b2   ( const Vec<1, RealType> &RefCoord ) { return RefCoord[0]; }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const Vec<1, RealType> &RefCoord );
  BASIS_FUNC_TYPE _deriv_basis[1][2];
  BASIS_FUNC_TYPE _basis[2];
  const RealType _h;

public:
  BaseFunctionSetMultiLin ( const RealType H ) : _h ( H ) {
    _basis[0] = _b1;
    _basis[1] = _b2;

    _deriv_basis[0][0] = _dx_b1;
    _deriv_basis[0][1] = _dx_b2;

    this->initializeQuadCache( );
  }

  enum { numBaseFuncs = 2 };

  void evaluateGradient ( int BaseFuncNum, const Vec<1, RealType> &RefCoord, Vec<1, RealType> &Gradient ) const {
    Gradient[0] = _deriv_basis[0][BaseFuncNum] ( RefCoord ) / _h;
  }

  inline const Vec<1, RealType>& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec<1, RealType>, Vec<1, RealType>, 2, QuadRuleType, BaseFunctionSetMultiLin<RealType, qc::QC_1D, QuadRuleType> >::evaluateGradient ( BaseFuncNum, QuadPoint );
  }

  RealType evaluate ( int BaseFuncNum, const Vec<1, RealType> &RefCoord ) const {
    return _basis[BaseFuncNum] ( RefCoord );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec<1, RealType>, Vec<1, RealType>, 2, QuadRuleType, BaseFunctionSetMultiLin<RealType, qc::QC_1D, QuadRuleType> >::evaluate ( BaseFuncNum, QuadPoint );
  }
protected:

};

//! The basefunctionset for bilinear elements in 2d.
template <typename RealType, class QuadRuleType>
class BaseFunctionSetMultiLin<RealType, qc::QC_2D, QuadRuleType> :
      public BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 4, QuadRuleType, BaseFunctionSetMultiLin<RealType, qc::QC_2D, QuadRuleType> > {

  static RealType _dx_b1   ( const Vec2<RealType> &RefCoord ) { return RefCoord[1] - aol::ZOTrait<RealType>::one; }
  static RealType _dx_b2   ( const Vec2<RealType> &RefCoord ) { return aol::ZOTrait<RealType>::one - RefCoord[1]; }
  static RealType _dx_b3   ( const Vec2<RealType> &RefCoord ) { return -RefCoord[1]; }
  static RealType _dx_b4   ( const Vec2<RealType> &RefCoord ) { return RefCoord[1]; }

  static RealType _dy_b1   ( const Vec2<RealType> &RefCoord ) { return RefCoord[0] - aol::ZOTrait<RealType>::one; }
  static RealType _dy_b2   ( const Vec2<RealType> &RefCoord ) { return -RefCoord[0]; }
  static RealType _dy_b3   ( const Vec2<RealType> &RefCoord ) { return aol::ZOTrait<RealType>::one - RefCoord[0]; }
  static RealType _dy_b4   ( const Vec2<RealType> &RefCoord ) { return RefCoord[0]; }

  static RealType _b1   ( const Vec2<RealType> &RefCoord ) { return ( aol::ZOTrait<RealType>::one - RefCoord[0] ) * ( aol::ZOTrait<RealType>::one - RefCoord[1] ); }
  static RealType _b2   ( const Vec2<RealType> &RefCoord ) { return RefCoord[0]* ( aol::ZOTrait<RealType>::one - RefCoord[1] ); }
  static RealType _b3   ( const Vec2<RealType> &RefCoord ) { return ( aol::ZOTrait<RealType>::one - RefCoord[0] ) *RefCoord[1]; }
  static RealType _b4   ( const Vec2<RealType> &RefCoord ) { return RefCoord[0]*RefCoord[1]; }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const Vec2<RealType> &RefCoord );
  BASIS_FUNC_TYPE _deriv_basis[2][4];
  BASIS_FUNC_TYPE _basis[4];
  const RealType _h;

public:
  BaseFunctionSetMultiLin ( const RealType H ) : _h ( H ) {
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;

    _deriv_basis[0][0] = _dx_b1;
    _deriv_basis[0][1] = _dx_b2;
    _deriv_basis[0][2] = _dx_b3;
    _deriv_basis[0][3] = _dx_b4;

    _deriv_basis[1][0] = _dy_b1;
    _deriv_basis[1][1] = _dy_b2;
    _deriv_basis[1][2] = _dy_b3;
    _deriv_basis[1][3] = _dy_b4;

    this->initializeQuadCache( );
  }

  enum { numBaseFuncs = 4 };

  void evaluateGradient ( int BaseFuncNum, const Vec2<RealType> &RefCoord, Vec2<RealType> &Gradient ) const {
    Gradient[0] = _deriv_basis[0][BaseFuncNum] ( RefCoord ) / _h;
    Gradient[1] = _deriv_basis[1][BaseFuncNum] ( RefCoord ) / _h;
  }

  inline const Vec2<RealType>& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 4, QuadRuleType, BaseFunctionSetMultiLin<RealType, qc::QC_2D, QuadRuleType> >::evaluateGradient ( BaseFuncNum, QuadPoint );
  }

  RealType evaluate ( int BaseFuncNum, const Vec2<RealType> &RefCoord ) const {
    return _basis[BaseFuncNum] ( RefCoord );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 4, QuadRuleType, BaseFunctionSetMultiLin<RealType, qc::QC_2D, QuadRuleType> >::evaluate ( BaseFuncNum, QuadPoint );
  }
protected:

};


//! The basefunctionset for bilinear elements in 3d.
template <typename RealType, class QuadRuleType>
class BaseFunctionSetMultiLin<RealType, qc::QC_3D, QuadRuleType> :
      public BaseFunctionSetInterface<RealType, Vec3<RealType>, Vec3<RealType>, 8, QuadRuleType, BaseFunctionSetMultiLin<RealType, qc::QC_3D, QuadRuleType> > {

  static inline RealType _b1_1d ( RealType x ) { return 1. - x; }
  static inline RealType _b2_1d ( RealType x ) { return x; }

  static inline RealType _d_b1_1d ( RealType ) { return -1.; }
  static inline RealType _d_b2_1d ( RealType ) { return 1.; }

  static RealType _b1   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b1_1d ( RefCoord[2] ); }
  static RealType _b2   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b1_1d ( RefCoord[2] ); }
  static RealType _b3   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b1_1d ( RefCoord[2] ); }
  static RealType _b4   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b1_1d ( RefCoord[2] ); }
  static RealType _b5   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b2_1d ( RefCoord[2] ); }
  static RealType _b6   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b2_1d ( RefCoord[2] ); }
  static RealType _b7   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b2_1d ( RefCoord[2] ); }
  static RealType _b8   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b2_1d ( RefCoord[2] ); }

  static RealType _dx_b1 ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) *_b1_1d ( RefCoord[1] ) *_b1_1d ( RefCoord[2] ); }
  static RealType _dx_b2 ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) *_b1_1d ( RefCoord[1] ) *_b1_1d ( RefCoord[2] ); }
  static RealType _dx_b3 ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) *_b2_1d ( RefCoord[1] ) *_b1_1d ( RefCoord[2] ); }
  static RealType _dx_b4 ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) *_b2_1d ( RefCoord[1] ) *_b1_1d ( RefCoord[2] ); }
  static RealType _dx_b5 ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) *_b1_1d ( RefCoord[1] ) *_b2_1d ( RefCoord[2] ); }
  static RealType _dx_b6 ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) *_b1_1d ( RefCoord[1] ) *_b2_1d ( RefCoord[2] ); }
  static RealType _dx_b7 ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) *_b2_1d ( RefCoord[1] ) *_b2_1d ( RefCoord[2] ); }
  static RealType _dx_b8 ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) *_b2_1d ( RefCoord[1] ) *_b2_1d ( RefCoord[2] ); }

  static RealType _dy_b1 ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) *_d_b1_1d ( RefCoord[1] ) *_b1_1d ( RefCoord[2] ); }
  static RealType _dy_b2 ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) *_d_b1_1d ( RefCoord[1] ) *_b1_1d ( RefCoord[2] ); }
  static RealType _dy_b3 ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) *_d_b2_1d ( RefCoord[1] ) *_b1_1d ( RefCoord[2] ); }
  static RealType _dy_b4 ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) *_d_b2_1d ( RefCoord[1] ) *_b1_1d ( RefCoord[2] ); }
  static RealType _dy_b5 ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) *_d_b1_1d ( RefCoord[1] ) *_b2_1d ( RefCoord[2] ); }
  static RealType _dy_b6 ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) *_d_b1_1d ( RefCoord[1] ) *_b2_1d ( RefCoord[2] ); }
  static RealType _dy_b7 ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) *_d_b2_1d ( RefCoord[1] ) *_b2_1d ( RefCoord[2] ); }
  static RealType _dy_b8 ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) *_d_b2_1d ( RefCoord[1] ) *_b2_1d ( RefCoord[2] ); }

  static RealType _dz_b1 ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) *_b1_1d ( RefCoord[1] ) *_d_b1_1d ( RefCoord[2] ); }
  static RealType _dz_b2 ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) *_b1_1d ( RefCoord[1] ) *_d_b1_1d ( RefCoord[2] ); }
  static RealType _dz_b3 ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) *_b2_1d ( RefCoord[1] ) *_d_b1_1d ( RefCoord[2] ); }
  static RealType _dz_b4 ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) *_b2_1d ( RefCoord[1] ) *_d_b1_1d ( RefCoord[2] ); }
  static RealType _dz_b5 ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) *_b1_1d ( RefCoord[1] ) *_d_b2_1d ( RefCoord[2] ); }
  static RealType _dz_b6 ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) *_b1_1d ( RefCoord[1] ) *_d_b2_1d ( RefCoord[2] ); }
  static RealType _dz_b7 ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) *_b2_1d ( RefCoord[1] ) *_d_b2_1d ( RefCoord[2] ); }
  static RealType _dz_b8 ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) *_b2_1d ( RefCoord[1] ) *_d_b2_1d ( RefCoord[2] ); }

#if 0
  static RealType _dx_b1   ( const Vec3<RealType> &RefCoord ) { return ( RefCoord[1] - 1. ) * ( 1. - RefCoord[2] ); }
  static RealType _dx_b2   ( const Vec3<RealType> &RefCoord ) { return ( 1. - RefCoord[1] ) * ( 1. - RefCoord[2] ); }
  static RealType _dx_b3   ( const Vec3<RealType> &RefCoord ) { return -RefCoord[1]* ( 1. - RefCoord[2] ); }
  static RealType _dx_b4   ( const Vec3<RealType> &RefCoord ) { return RefCoord[1]* ( 1. - RefCoord[2] ); }
  static RealType _dx_b5   ( const Vec3<RealType> &RefCoord ) { return ( RefCoord[1] - 1. ) *RefCoord[2]; }
  static RealType _dx_b6   ( const Vec3<RealType> &RefCoord ) { return ( 1. - RefCoord[1] ) *RefCoord[2]; }
  static RealType _dx_b7   ( const Vec3<RealType> &RefCoord ) { return -RefCoord[1]*RefCoord[2]; }
  static RealType _dx_b8   ( const Vec3<RealType> &RefCoord ) { return RefCoord[1]*RefCoord[2]; }

  static RealType _dy_b1   ( const Vec3<RealType> &RefCoord ) { return ( RefCoord[0] - 1. ) * ( 1. - RefCoord[2] ); }
  static RealType _dy_b2   ( const Vec3<RealType> &RefCoord ) { return -RefCoord[0]* ( 1. - RefCoord[2] ); }
  static RealType _dy_b3   ( const Vec3<RealType> &RefCoord ) { return ( 1. - RefCoord[0] ) * ( 1. - RefCoord[2] ); }
  static RealType _dy_b4   ( const Vec3<RealType> &RefCoord ) { return RefCoord[0]* ( 1. - RefCoord[2] ); }
  static RealType _dy_b5   ( const Vec3<RealType> &RefCoord ) { return ( RefCoord[0] - 1. ) *RefCoord[2]; }
  static RealType _dy_b6   ( const Vec3<RealType> &RefCoord ) { return -RefCoord[0]*RefCoord[2]; }
  static RealType _dy_b7   ( const Vec3<RealType> &RefCoord ) { return ( 1. - RefCoord[0] ) *RefCoord[2]; }
  static RealType _dy_b8   ( const Vec3<RealType> &RefCoord ) { return RefCoord[0]*RefCoord[2]; }

  static RealType _dz_b1   ( const Vec3<RealType> &RefCoord ) { return ( RefCoord[0] - 1. ) * ( 1. - RefCoord[1] ); }
  static RealType _dz_b2   ( const Vec3<RealType> &RefCoord ) { return -RefCoord[0]* ( 1. - RefCoord[1] ); }
  static RealType _dz_b3   ( const Vec3<RealType> &RefCoord ) { return ( RefCoord[0] - 1. ) *RefCoord[1]; }
  static RealType _dz_b4   ( const Vec3<RealType> &RefCoord ) { return -RefCoord[1]*RefCoord[0]; }
  static RealType _dz_b5   ( const Vec3<RealType> &RefCoord ) { return ( 1. - RefCoord[0] ) * ( 1. - RefCoord[1] ); }
  static RealType _dz_b6   ( const Vec3<RealType> &RefCoord ) { return RefCoord[0]* ( 1. - RefCoord[1] ); }
  static RealType _dz_b7   ( const Vec3<RealType> &RefCoord ) { return ( 1. - RefCoord[0] ) *RefCoord[1]; }
  static RealType _dz_b8   ( const Vec3<RealType> &RefCoord ) { return RefCoord[0]*RefCoord[1]; }

  static RealType _b1   ( const Vec3<RealType> &RefCoord ) { return ( 1. - RefCoord[0] ) * ( 1 - RefCoord[1] ) * ( 1. - RefCoord[2] ); }
  static RealType _b2   ( const Vec3<RealType> &RefCoord ) { return RefCoord[0]* ( 1. - RefCoord[1] ) * ( 1. - RefCoord[2] ); }
  static RealType _b3   ( const Vec3<RealType> &RefCoord ) { return ( 1. - RefCoord[0] ) *RefCoord[1]* ( 1. - RefCoord[2] ); }
  static RealType _b4   ( const Vec3<RealType> &RefCoord ) { return RefCoord[0]*RefCoord[1]* ( 1 - RefCoord[2] ); }

  static RealType _b5   ( const Vec3<RealType> &RefCoord ) { return ( 1. - RefCoord[0] ) * ( 1 - RefCoord[1] ) *RefCoord[2]; }
  static RealType _b6   ( const Vec3<RealType> &RefCoord ) { return RefCoord[0]* ( 1. - RefCoord[1] ) *RefCoord[2]; }
  static RealType _b7   ( const Vec3<RealType> &RefCoord ) { return ( 1. - RefCoord[0] ) *RefCoord[1]*RefCoord[2]; }
  static RealType _b8   ( const Vec3<RealType> &RefCoord ) { return RefCoord[0]*RefCoord[1]*RefCoord[2]; }
#endif

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const Vec3<RealType> &RefCoord );
  BASIS_FUNC_TYPE _deriv_basis[3][8];
  BASIS_FUNC_TYPE _basis[8];

  const RealType _h;

public:
  BaseFunctionSetMultiLin( const RealType H ) : _h ( H ) {
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;
    _basis[4] = _b5;
    _basis[5] = _b6;
    _basis[6] = _b7;
    _basis[7] = _b8;

    _deriv_basis[0][0] = _dx_b1;
    _deriv_basis[0][1] = _dx_b2;
    _deriv_basis[0][2] = _dx_b3;
    _deriv_basis[0][3] = _dx_b4;
    _deriv_basis[0][4] = _dx_b5;
    _deriv_basis[0][5] = _dx_b6;
    _deriv_basis[0][6] = _dx_b7;
    _deriv_basis[0][7] = _dx_b8;

    _deriv_basis[1][0] = _dy_b1;
    _deriv_basis[1][1] = _dy_b2;
    _deriv_basis[1][2] = _dy_b3;
    _deriv_basis[1][3] = _dy_b4;
    _deriv_basis[1][4] = _dy_b5;
    _deriv_basis[1][5] = _dy_b6;
    _deriv_basis[1][6] = _dy_b7;
    _deriv_basis[1][7] = _dy_b8;

    _deriv_basis[2][0] = _dz_b1;
    _deriv_basis[2][1] = _dz_b2;
    _deriv_basis[2][2] = _dz_b3;
    _deriv_basis[2][3] = _dz_b4;
    _deriv_basis[2][4] = _dz_b5;
    _deriv_basis[2][5] = _dz_b6;
    _deriv_basis[2][6] = _dz_b7;
    _deriv_basis[2][7] = _dz_b8;

    this->initializeQuadCache( );
  }

  enum { numBaseFuncs = 8 };

  void evaluateGradient ( int BaseFuncNum, const Vec3<RealType> &RefCoord, Vec3<RealType> &Gradient ) const {
    Gradient[0] = _deriv_basis[0][BaseFuncNum] ( RefCoord ) / _h;
    Gradient[1] = _deriv_basis[1][BaseFuncNum] ( RefCoord ) / _h;
    Gradient[2] = _deriv_basis[2][BaseFuncNum] ( RefCoord ) / _h;
  }

  inline const Vec3<RealType>& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec3<RealType>, Vec3<RealType>, 8, QuadRuleType, BaseFunctionSetMultiLin<RealType, qc::QC_3D, QuadRuleType> >::evaluateGradient ( BaseFuncNum, QuadPoint );
  }

  RealType evaluate ( int BaseFuncNum, const Vec3<RealType> &RefCoord ) const {
    return _basis[BaseFuncNum] ( RefCoord );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec3<RealType>, Vec3<RealType>, 8, QuadRuleType, BaseFunctionSetMultiLin<RealType, qc::QC_3D, QuadRuleType> >::evaluate ( BaseFuncNum, QuadPoint );
  }
protected:

};


template <typename RealType, qc::Dimension Dim, class QuadRuleType>
class BaseFunctionSetMultiQuad : public BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 9, QuadRuleType, BaseFunctionSetMultiQuad<RealType, qc::QC_2D, QuadRuleType> > {
};

/**
 *  The basefunctionset for biquadratic elements in 2d.
 */
template <typename RealType, class QuadRuleType>
class BaseFunctionSetMultiQuad<RealType, qc::QC_2D, QuadRuleType> :
      public BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 9, QuadRuleType, BaseFunctionSetMultiQuad<RealType, qc::QC_2D, QuadRuleType> > {



  static inline RealType _b1_1d ( RealType x ) { return 2.*x*x - 3.*x + 1.; }
  static inline RealType _b2_1d ( RealType x ) { return -4.*x*x + 4.*x; }
  static inline RealType _b3_1d ( RealType x ) { return 2.*x*x - x;  }

  static inline RealType _d_b1_1d ( RealType x ) { return 4.*x - 3.; }
  static inline RealType _d_b2_1d ( RealType x ) { return -8.*x + 4.; }
  static inline RealType _d_b3_1d ( RealType x ) { return 4.*x - 1.;  }

  // will be converted below, static const RealType is only valid for C++11 using constexpr keyword
  static const int _d2_b1_1d =  4;
  static const int _d2_b2_1d = -8;
  static const int _d2_b3_1d =  4;


  static RealType _b1   ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ); }
  static RealType _b2   ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ); }
  static RealType _b3   ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ); }

  static RealType _b4   ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ); }
  static RealType _b5   ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ); }
  static RealType _b6   ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ); }

  static RealType _b7   ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ); }
  static RealType _b8   ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ); }
  static RealType _b9   ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ); }


  static RealType _dx_b1   ( const Vec2<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ); }
  static RealType _dx_b2   ( const Vec2<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ); }
  static RealType _dx_b3   ( const Vec2<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ); }
  static RealType _dx_b4   ( const Vec2<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ); }
  static RealType _dx_b5   ( const Vec2<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ); }
  static RealType _dx_b6   ( const Vec2<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ); }
  static RealType _dx_b7   ( const Vec2<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ); }
  static RealType _dx_b8   ( const Vec2<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ); }
  static RealType _dx_b9   ( const Vec2<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ); }

  static RealType _dy_b1   ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ); }
  static RealType _dy_b2   ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ); }
  static RealType _dy_b3   ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ); }
  static RealType _dy_b4   ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ); }
  static RealType _dy_b5   ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ); }
  static RealType _dy_b6   ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ); }
  static RealType _dy_b7   ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ); }
  static RealType _dy_b8   ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ); }
  static RealType _dy_b9   ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ); }


  static RealType _dxx_b1   ( const Vec2<RealType> &RefCoord ) { return _d2_b1_1d * _b1_1d ( RefCoord[1] ); }
  static RealType _dxx_b2   ( const Vec2<RealType> &RefCoord ) { return _d2_b2_1d * _b1_1d ( RefCoord[1] ); }
  static RealType _dxx_b3   ( const Vec2<RealType> &RefCoord ) { return _d2_b3_1d * _b1_1d ( RefCoord[1] ); }
  static RealType _dxx_b4   ( const Vec2<RealType> &RefCoord ) { return _d2_b1_1d * _b2_1d ( RefCoord[1] ); }
  static RealType _dxx_b5   ( const Vec2<RealType> &RefCoord ) { return _d2_b2_1d * _b2_1d ( RefCoord[1] ); }
  static RealType _dxx_b6   ( const Vec2<RealType> &RefCoord ) { return _d2_b3_1d * _b2_1d ( RefCoord[1] ); }
  static RealType _dxx_b7   ( const Vec2<RealType> &RefCoord ) { return _d2_b1_1d * _b3_1d ( RefCoord[1] ); }
  static RealType _dxx_b8   ( const Vec2<RealType> &RefCoord ) { return _d2_b2_1d * _b3_1d ( RefCoord[1] ); }
  static RealType _dxx_b9   ( const Vec2<RealType> &RefCoord ) { return _d2_b3_1d * _b3_1d ( RefCoord[1] ); }

  static RealType _dyy_b1   ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d2_b1_1d; }
  static RealType _dyy_b2   ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d2_b1_1d; }
  static RealType _dyy_b3   ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d2_b1_1d; }
  static RealType _dyy_b4   ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d2_b2_1d; }
  static RealType _dyy_b5   ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d2_b2_1d; }
  static RealType _dyy_b6   ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d2_b2_1d; }
  static RealType _dyy_b7   ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d2_b3_1d; }
  static RealType _dyy_b8   ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d2_b3_1d; }
  static RealType _dyy_b9   ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d2_b3_1d; }

  static RealType _dxy_b1   ( const Vec2<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ); }
  static RealType _dxy_b2   ( const Vec2<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ); }
  static RealType _dxy_b3   ( const Vec2<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ); }
  static RealType _dxy_b4   ( const Vec2<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ); }
  static RealType _dxy_b5   ( const Vec2<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ); }
  static RealType _dxy_b6   ( const Vec2<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ); }
  static RealType _dxy_b7   ( const Vec2<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ); }
  static RealType _dxy_b8   ( const Vec2<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ); }
  static RealType _dxy_b9   ( const Vec2<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ); }




  typedef RealType ( *BASIS_FUNC_TYPE ) ( const Vec2<RealType> &RefCoord );
  BASIS_FUNC_TYPE _deriv_basis[2][9];
  BASIS_FUNC_TYPE _basis[9];
  BASIS_FUNC_TYPE _secondDeriv_basis[3][9];

  const RealType _h;

public:
  BaseFunctionSetMultiQuad( const RealType H ) : _h ( H ) {
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;
    _basis[4] = _b5;
    _basis[5] = _b6;
    _basis[6] = _b7;
    _basis[7] = _b8;
    _basis[8] = _b9;

    _deriv_basis[0][0] = _dx_b1;
    _deriv_basis[0][1] = _dx_b2;
    _deriv_basis[0][2] = _dx_b3;
    _deriv_basis[0][3] = _dx_b4;
    _deriv_basis[0][4] = _dx_b5;
    _deriv_basis[0][5] = _dx_b6;
    _deriv_basis[0][6] = _dx_b7;
    _deriv_basis[0][7] = _dx_b8;
    _deriv_basis[0][8] = _dx_b9;

    _deriv_basis[1][0] = _dy_b1;
    _deriv_basis[1][1] = _dy_b2;
    _deriv_basis[1][2] = _dy_b3;
    _deriv_basis[1][3] = _dy_b4;
    _deriv_basis[1][4] = _dy_b5;
    _deriv_basis[1][5] = _dy_b6;
    _deriv_basis[1][6] = _dy_b7;
    _deriv_basis[1][7] = _dy_b8;
    _deriv_basis[1][8] = _dy_b9;

    _secondDeriv_basis[0][0] = _dxx_b1;
    _secondDeriv_basis[0][1] = _dxx_b2;
    _secondDeriv_basis[0][2] = _dxx_b3;
    _secondDeriv_basis[0][3] = _dxx_b4;
    _secondDeriv_basis[0][4] = _dxx_b5;
    _secondDeriv_basis[0][5] = _dxx_b6;
    _secondDeriv_basis[0][6] = _dxx_b7;
    _secondDeriv_basis[0][7] = _dxx_b8;
    _secondDeriv_basis[0][8] = _dxx_b9;

    _secondDeriv_basis[1][0] = _dyy_b1;
    _secondDeriv_basis[1][1] = _dyy_b2;
    _secondDeriv_basis[1][2] = _dyy_b3;
    _secondDeriv_basis[1][3] = _dyy_b4;
    _secondDeriv_basis[1][4] = _dyy_b5;
    _secondDeriv_basis[1][5] = _dyy_b6;
    _secondDeriv_basis[1][6] = _dyy_b7;
    _secondDeriv_basis[1][7] = _dyy_b8;
    _secondDeriv_basis[1][8] = _dyy_b9;

    _secondDeriv_basis[2][0] = _dxy_b1;
    _secondDeriv_basis[2][1] = _dxy_b2;
    _secondDeriv_basis[2][2] = _dxy_b3;
    _secondDeriv_basis[2][3] = _dxy_b4;
    _secondDeriv_basis[2][4] = _dxy_b5;
    _secondDeriv_basis[2][5] = _dxy_b6;
    _secondDeriv_basis[2][6] = _dxy_b7;
    _secondDeriv_basis[2][7] = _dxy_b8;
    _secondDeriv_basis[2][8] = _dxy_b9;

    this->initializeQuadCache( );
    // BaseFunctionSetInterface caches only values and gradients (second derivatives are zero in the linear case)
    initializeExtendedQuadCache( );
  }

  enum { numBaseFuncs = 9 };

  void evaluateGradient ( int BaseFuncNum, const Vec2<RealType> &RefCoord, Vec2<RealType> &Gradient ) const {
    Gradient[0] = _deriv_basis[0][BaseFuncNum] ( RefCoord ) / _h;
    Gradient[1] = _deriv_basis[1][BaseFuncNum] ( RefCoord ) / _h;
  }

  inline const Vec2<RealType>& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 9, QuadRuleType, BaseFunctionSetMultiQuad<RealType, qc::QC_2D, QuadRuleType> >::evaluateGradient ( BaseFuncNum, QuadPoint );
  }

  RealType evaluate ( int BaseFuncNum, const Vec2<RealType> &RefCoord ) const {
    return _basis[BaseFuncNum] ( RefCoord );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 9, QuadRuleType, BaseFunctionSetMultiQuad<RealType, qc::QC_2D, QuadRuleType> >::evaluate ( BaseFuncNum, QuadPoint );
  }

  void evaluateHessian ( int BaseFuncNum, const Vec2<RealType> &RefCoord, Matrix22<RealType> &Hessian ) const {
    RealType hsqr = Sqr(_h);
    Hessian[0][0] = _secondDeriv_basis[0][BaseFuncNum] ( RefCoord ) / hsqr;
    Hessian[1][1] = _secondDeriv_basis[1][BaseFuncNum] ( RefCoord ) / hsqr;
    Hessian[0][1] = _secondDeriv_basis[2][BaseFuncNum] ( RefCoord ) / hsqr;
    Hessian[1][0] = _secondDeriv_basis[2][BaseFuncNum] ( RefCoord ) / hsqr;
  }

  inline const Matrix22<RealType>& evaluateHessian ( int BaseFuncNum, int QuadPoint ) const {
    return basisQuadHessians[BaseFuncNum][QuadPoint];
  }

protected:
  void initializeExtendedQuadCache( ) {
    for ( int b = 0; b < numBaseFuncs; b++ ) {
      for ( int i = 0; i < QuadRuleType::numQuadPoints; i++ ) {
        evaluateHessian ( b, this->_quadRule.getRefCoord ( i ), basisQuadHessians[b][i] );
      }
    }
  }

  /**** cache extension for Hessians ****/
  Matrix22<RealType>  basisQuadHessians[numBaseFuncs][QuadRuleType::numQuadPoints];
};


template <typename RealType, qc::Dimension Dim, class QuadRuleType>
class BaseFunctionSetMultiQuart : public BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 9, QuadRuleType, BaseFunctionSetMultiQuart<RealType, qc::QC_2D, QuadRuleType> > {
};



/**
 *  The basefunctionset for biquadratic elements in 3d.
 */
template <typename RealType, class QuadRuleType>
class BaseFunctionSetMultiQuad<RealType, qc::QC_3D, QuadRuleType> :
      public BaseFunctionSetInterface<RealType, Vec3<RealType>, Vec3<RealType>, 27, QuadRuleType, BaseFunctionSetMultiQuad<RealType, qc::QC_3D, QuadRuleType> > {

  const static int numDof = 27;

  static inline RealType _b1_1d ( RealType x ) { return 2.*x*x - 3.*x + 1.; }
  static inline RealType _b2_1d ( RealType x ) { return -4.*x*x + 4.*x; }
  static inline RealType _b3_1d ( RealType x ) { return 2.*x*x - x;  }

  static inline RealType _d_b1_1d ( RealType x ) { return 4.*x - 3.; }
  static inline RealType _d_b2_1d ( RealType x ) { return -8.*x + 4.; }
  static inline RealType _d_b3_1d ( RealType x ) { return 4.*x - 1.;  }

  // will be converted below, static const RealType is only valid for C++11 using constexpr keyword
  static const int _d2_b1_1d =  4;
  static const int _d2_b2_1d = -8;
  static const int _d2_b3_1d =  4;


  static RealType _b1   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b1_1d ( RefCoord[2] ); }
  static RealType _b2   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b1_1d ( RefCoord[2] ); }
  static RealType _b3   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b1_1d ( RefCoord[2] ); }
  static RealType _b4   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b1_1d ( RefCoord[2] ); }
  static RealType _b5   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b1_1d ( RefCoord[2] ); }
  static RealType _b6   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b1_1d ( RefCoord[2] ); }
  static RealType _b7   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b1_1d ( RefCoord[2] ); }
  static RealType _b8   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b1_1d ( RefCoord[2] ); }
  static RealType _b9   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b1_1d ( RefCoord[2] ); }

  static RealType _b10   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b2_1d ( RefCoord[2] ); }
  static RealType _b11   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b2_1d ( RefCoord[2] ); }
  static RealType _b12   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b2_1d ( RefCoord[2] ); }
  static RealType _b13   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b2_1d ( RefCoord[2] ); }
  static RealType _b14   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b2_1d ( RefCoord[2] ); }
  static RealType _b15   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b2_1d ( RefCoord[2] ); }
  static RealType _b16   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b2_1d ( RefCoord[2] ); }
  static RealType _b17   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b2_1d ( RefCoord[2] ); }
  static RealType _b18   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b2_1d ( RefCoord[2] ); }
  
  static RealType _b19   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b3_1d ( RefCoord[2] ); }
  static RealType _b20   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b3_1d ( RefCoord[2] ); }
  static RealType _b21   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b3_1d ( RefCoord[2] ); }
  static RealType _b22   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b3_1d ( RefCoord[2] ); }
  static RealType _b23   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b3_1d ( RefCoord[2] ); }
  static RealType _b24   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b3_1d ( RefCoord[2] ); }
  static RealType _b25   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b3_1d ( RefCoord[2] ); }
  static RealType _b26   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b3_1d ( RefCoord[2] ); }
  static RealType _b27   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b3_1d ( RefCoord[2] ); }
  
  
  
  static RealType _dx_b1   ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dx_b2   ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dx_b3   ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dx_b4   ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dx_b5   ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dx_b6   ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dx_b7   ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dx_b8   ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dx_b9   ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  
  static RealType _dx_b10  ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dx_b11  ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dx_b12  ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dx_b13  ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dx_b14  ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dx_b15  ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dx_b16  ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dx_b17  ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dx_b18  ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  
  static RealType _dx_b19   ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dx_b20   ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dx_b21   ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dx_b22   ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dx_b23   ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dx_b24   ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dx_b25   ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dx_b26   ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dx_b27   ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  
  
  
  static RealType _dy_b1   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dy_b2   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dy_b3   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dy_b4   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dy_b5   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dy_b6   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dy_b7   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dy_b8   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dy_b9   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  
  static RealType _dy_b10  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dy_b11  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dy_b12  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dy_b13  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dy_b14  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dy_b15  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dy_b16  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dy_b17  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dy_b18  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  
  static RealType _dy_b19  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dy_b20  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dy_b21  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dy_b22  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dy_b23  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dy_b24  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dy_b25  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dy_b26  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dy_b27  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }


  
  static RealType _dz_b1   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dz_b2   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dz_b3   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dz_b4   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dz_b5   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dz_b6   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dz_b7   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dz_b8   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dz_b9   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  
  static RealType _dz_b10  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dz_b11  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dz_b12  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dz_b13  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dz_b14  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dz_b15  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dz_b16  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dz_b17  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dz_b18  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  
  static RealType _dz_b19  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dz_b20  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dz_b21  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dz_b22  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dz_b23  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dz_b24  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dz_b25  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dz_b26  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dz_b27  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  
  
  
  static RealType _dxx_b1   ( const Vec3<RealType> &RefCoord ) { return _d2_b1_1d * _b1_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dxx_b2   ( const Vec3<RealType> &RefCoord ) { return _d2_b2_1d * _b1_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dxx_b3   ( const Vec3<RealType> &RefCoord ) { return _d2_b3_1d * _b1_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dxx_b4   ( const Vec3<RealType> &RefCoord ) { return _d2_b1_1d * _b2_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dxx_b5   ( const Vec3<RealType> &RefCoord ) { return _d2_b2_1d * _b2_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dxx_b6   ( const Vec3<RealType> &RefCoord ) { return _d2_b3_1d * _b2_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dxx_b7   ( const Vec3<RealType> &RefCoord ) { return _d2_b1_1d * _b3_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dxx_b8   ( const Vec3<RealType> &RefCoord ) { return _d2_b2_1d * _b3_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dxx_b9   ( const Vec3<RealType> &RefCoord ) { return _d2_b3_1d * _b3_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  
  static RealType _dxx_b10  ( const Vec3<RealType> &RefCoord ) { return _d2_b1_1d * _b1_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dxx_b11  ( const Vec3<RealType> &RefCoord ) { return _d2_b2_1d * _b1_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dxx_b12  ( const Vec3<RealType> &RefCoord ) { return _d2_b3_1d * _b1_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dxx_b13  ( const Vec3<RealType> &RefCoord ) { return _d2_b1_1d * _b2_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dxx_b14  ( const Vec3<RealType> &RefCoord ) { return _d2_b2_1d * _b2_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dxx_b15  ( const Vec3<RealType> &RefCoord ) { return _d2_b3_1d * _b2_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dxx_b16  ( const Vec3<RealType> &RefCoord ) { return _d2_b1_1d * _b3_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dxx_b17  ( const Vec3<RealType> &RefCoord ) { return _d2_b2_1d * _b3_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dxx_b18  ( const Vec3<RealType> &RefCoord ) { return _d2_b3_1d * _b3_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  
  static RealType _dxx_b19  ( const Vec3<RealType> &RefCoord ) { return _d2_b1_1d * _b1_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dxx_b20  ( const Vec3<RealType> &RefCoord ) { return _d2_b2_1d * _b1_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dxx_b21  ( const Vec3<RealType> &RefCoord ) { return _d2_b3_1d * _b1_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dxx_b22  ( const Vec3<RealType> &RefCoord ) { return _d2_b1_1d * _b2_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dxx_b23  ( const Vec3<RealType> &RefCoord ) { return _d2_b2_1d * _b2_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dxx_b24  ( const Vec3<RealType> &RefCoord ) { return _d2_b3_1d * _b2_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dxx_b25  ( const Vec3<RealType> &RefCoord ) { return _d2_b1_1d * _b3_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dxx_b26  ( const Vec3<RealType> &RefCoord ) { return _d2_b2_1d * _b3_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dxx_b27  ( const Vec3<RealType> &RefCoord ) { return _d2_b3_1d * _b3_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  
  
  
  static RealType _dyy_b1   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d2_b1_1d * _b1_1d( RefCoord[2] ); }
  static RealType _dyy_b2   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d2_b1_1d * _b1_1d( RefCoord[2] ); }
  static RealType _dyy_b3   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d2_b1_1d * _b1_1d( RefCoord[2] ); }
  static RealType _dyy_b4   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d2_b2_1d * _b1_1d( RefCoord[2] ); }
  static RealType _dyy_b5   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d2_b2_1d * _b1_1d( RefCoord[2] ); }
  static RealType _dyy_b6   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d2_b2_1d * _b1_1d( RefCoord[2] ); }
  static RealType _dyy_b7   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d2_b3_1d * _b1_1d( RefCoord[2] ); }
  static RealType _dyy_b8   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d2_b3_1d * _b1_1d( RefCoord[2] ); }
  static RealType _dyy_b9   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d2_b3_1d * _b1_1d( RefCoord[2] ); }
  
  static RealType _dyy_b10  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d2_b1_1d * _b2_1d( RefCoord[2] ); }
  static RealType _dyy_b11  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d2_b1_1d * _b2_1d( RefCoord[2] ); }
  static RealType _dyy_b12  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d2_b1_1d * _b2_1d( RefCoord[2] ); }
  static RealType _dyy_b13  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d2_b2_1d * _b2_1d( RefCoord[2] ); }
  static RealType _dyy_b14  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d2_b2_1d * _b2_1d( RefCoord[2] ); }
  static RealType _dyy_b15  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d2_b2_1d * _b2_1d( RefCoord[2] ); }
  static RealType _dyy_b16  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d2_b3_1d * _b2_1d( RefCoord[2] ); }
  static RealType _dyy_b17  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d2_b3_1d * _b2_1d( RefCoord[2] ); }
  static RealType _dyy_b18  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d2_b3_1d * _b2_1d( RefCoord[2] ); }
  
  static RealType _dyy_b19  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d2_b1_1d * _b3_1d( RefCoord[2] ); }
  static RealType _dyy_b20  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d2_b1_1d * _b3_1d( RefCoord[2] ); }
  static RealType _dyy_b21  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d2_b1_1d * _b3_1d( RefCoord[2] ); }
  static RealType _dyy_b22  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d2_b2_1d * _b3_1d( RefCoord[2] ); }
  static RealType _dyy_b23  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d2_b2_1d * _b3_1d( RefCoord[2] ); }
  static RealType _dyy_b24  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d2_b2_1d * _b3_1d( RefCoord[2] ); }
  static RealType _dyy_b25  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d2_b3_1d * _b3_1d( RefCoord[2] ); }
  static RealType _dyy_b26  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d2_b3_1d * _b3_1d( RefCoord[2] ); }
  static RealType _dyy_b27  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d2_b3_1d * _b3_1d( RefCoord[2] ); }
  
  

  static RealType _dzz_b1   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d2_b1_1d; }
  static RealType _dzz_b2   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d2_b1_1d; }
  static RealType _dzz_b3   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d2_b1_1d; }
  static RealType _dzz_b4   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d2_b1_1d; }
  static RealType _dzz_b5   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d2_b1_1d; }
  static RealType _dzz_b6   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d2_b1_1d; }
  static RealType _dzz_b7   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d2_b1_1d; }
  static RealType _dzz_b8   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d2_b1_1d; }
  static RealType _dzz_b9   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d2_b1_1d; }
  
  static RealType _dzz_b10  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d2_b2_1d; }
  static RealType _dzz_b11  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d2_b2_1d; }
  static RealType _dzz_b12  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d2_b2_1d; }
  static RealType _dzz_b13  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d2_b2_1d; }
  static RealType _dzz_b14  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d2_b2_1d; }
  static RealType _dzz_b15  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d2_b2_1d; }
  static RealType _dzz_b16  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d2_b2_1d; }
  static RealType _dzz_b17  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d2_b2_1d; }
  static RealType _dzz_b18  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d2_b2_1d; }
  
  static RealType _dzz_b19  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d2_b3_1d; }
  static RealType _dzz_b20  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d2_b3_1d; }
  static RealType _dzz_b21  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d2_b3_1d; }
  static RealType _dzz_b22  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d2_b3_1d; }
  static RealType _dzz_b23  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d2_b3_1d; }
  static RealType _dzz_b24  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d2_b3_1d; }
  static RealType _dzz_b25  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d2_b3_1d; }
  static RealType _dzz_b26  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d2_b3_1d; }
  static RealType _dzz_b27  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d2_b3_1d; }
  
  
  
  static RealType _dxy_b1   ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dxy_b2   ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dxy_b3   ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dxy_b4   ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dxy_b5   ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dxy_b6   ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dxy_b7   ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dxy_b8   ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }
  static RealType _dxy_b9   ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b1_1d( RefCoord[2] ); }

  static RealType _dxy_b10  ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dxy_b11  ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dxy_b12  ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dxy_b13  ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dxy_b14  ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dxy_b15  ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dxy_b16  ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dxy_b17  ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  static RealType _dxy_b18  ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b2_1d( RefCoord[2] ); }
  
  static RealType _dxy_b19  ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dxy_b20  ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dxy_b21  ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dxy_b22  ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dxy_b23  ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dxy_b24  ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dxy_b25  ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dxy_b26  ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }
  static RealType _dxy_b27  ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _b3_1d( RefCoord[2] ); }

 
 
  static RealType _dxz_b1   ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dxz_b2   ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dxz_b3   ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dxz_b4   ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dxz_b5   ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dxz_b6   ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dxz_b7   ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dxz_b8   ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dxz_b9   ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  
  static RealType _dxz_b10  ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dxz_b11  ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dxz_b12  ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dxz_b13  ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dxz_b14  ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dxz_b15  ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dxz_b16  ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dxz_b17  ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dxz_b18  ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  
  static RealType _dxz_b19  ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dxz_b20  ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dxz_b21  ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dxz_b22  ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dxz_b23  ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dxz_b24  ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dxz_b25  ( const Vec3<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dxz_b26  ( const Vec3<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dxz_b27  ( const Vec3<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  
  
  
  static RealType _dyz_b1   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dyz_b2   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dyz_b3   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dyz_b4   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dyz_b5   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dyz_b6   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dyz_b7   ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dyz_b8   ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  static RealType _dyz_b9   ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _d_b1_1d( RefCoord[2] ); }
  
  static RealType _dyz_b10  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dyz_b11  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dyz_b12  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dyz_b13  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dyz_b14  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dyz_b15  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dyz_b16  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dyz_b17  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  static RealType _dyz_b18  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _d_b2_1d( RefCoord[2] ); }
  
  static RealType _dyz_b19  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dyz_b20  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dyz_b21  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dyz_b22  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dyz_b23  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dyz_b24  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dyz_b25  ( const Vec3<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dyz_b26  ( const Vec3<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  static RealType _dyz_b27  ( const Vec3<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ) * _d_b3_1d( RefCoord[2] ); }
  
  
  typedef RealType ( *BASIS_FUNC_TYPE ) ( const Vec3<RealType> &RefCoord );
  BASIS_FUNC_TYPE _deriv_basis[3][27];
  BASIS_FUNC_TYPE _basis[27];
  BASIS_FUNC_TYPE _secondDeriv_basis[6][27];

  const RealType _h;
  
public:
  const static int numDofs = 27;
  
  BaseFunctionSetMultiQuad( const RealType H ) : _h(H) {
    _basis[0] = _b1; _basis[9]  = _b10; _basis[18] = _b19;
    _basis[1] = _b2; _basis[10] = _b11; _basis[19] = _b20;
    _basis[2] = _b3; _basis[11] = _b12; _basis[20] = _b21;
    _basis[3] = _b4; _basis[12] = _b13; _basis[21] = _b22;
    _basis[4] = _b5; _basis[13] = _b14; _basis[22] = _b23;
    _basis[5] = _b6; _basis[14] = _b15; _basis[23] = _b24;
    _basis[6] = _b7; _basis[15] = _b16; _basis[24] = _b25;
    _basis[7] = _b8; _basis[16] = _b17; _basis[25] = _b26;
    _basis[8] = _b9; _basis[17] = _b18; _basis[26] = _b27;

    _deriv_basis[0][0] = _dx_b1; _deriv_basis[0][9]  = _dx_b10; _deriv_basis[0][18] = _dx_b19;
    _deriv_basis[0][1] = _dx_b2; _deriv_basis[0][10] = _dx_b11; _deriv_basis[0][19] = _dx_b20;
    _deriv_basis[0][2] = _dx_b3; _deriv_basis[0][11] = _dx_b12; _deriv_basis[0][20] = _dx_b21;
    _deriv_basis[0][3] = _dx_b4; _deriv_basis[0][12] = _dx_b13; _deriv_basis[0][21] = _dx_b22;
    _deriv_basis[0][4] = _dx_b5; _deriv_basis[0][13] = _dx_b14; _deriv_basis[0][22] = _dx_b23;
    _deriv_basis[0][5] = _dx_b6; _deriv_basis[0][14] = _dx_b15; _deriv_basis[0][23] = _dx_b24;
    _deriv_basis[0][6] = _dx_b7; _deriv_basis[0][15] = _dx_b16; _deriv_basis[0][24] = _dx_b25;
    _deriv_basis[0][7] = _dx_b8; _deriv_basis[0][16] = _dx_b17; _deriv_basis[0][25] = _dx_b26;
    _deriv_basis[0][8] = _dx_b9; _deriv_basis[0][17] = _dx_b18; _deriv_basis[0][26] = _dx_b27;

    _deriv_basis[1][0] = _dy_b1; _deriv_basis[1][9]  = _dy_b10; _deriv_basis[1][18] = _dy_b19;
    _deriv_basis[1][1] = _dy_b2; _deriv_basis[1][10] = _dy_b11; _deriv_basis[1][19] = _dy_b20;
    _deriv_basis[1][2] = _dy_b3; _deriv_basis[1][11] = _dy_b12; _deriv_basis[1][20] = _dy_b21;
    _deriv_basis[1][3] = _dy_b4; _deriv_basis[1][12] = _dy_b13; _deriv_basis[1][21] = _dy_b22;
    _deriv_basis[1][4] = _dy_b5; _deriv_basis[1][13] = _dy_b14; _deriv_basis[1][22] = _dy_b23;
    _deriv_basis[1][5] = _dy_b6; _deriv_basis[1][14] = _dy_b15; _deriv_basis[1][23] = _dy_b24;
    _deriv_basis[1][6] = _dy_b7; _deriv_basis[1][15] = _dy_b16; _deriv_basis[1][24] = _dy_b25;
    _deriv_basis[1][7] = _dy_b8; _deriv_basis[1][16] = _dy_b17; _deriv_basis[1][25] = _dy_b26;
    _deriv_basis[1][8] = _dy_b9; _deriv_basis[1][17] = _dy_b18; _deriv_basis[1][26] = _dy_b27;
    
    _deriv_basis[2][0] = _dz_b1; _deriv_basis[2][9]  = _dz_b10; _deriv_basis[2][18] = _dz_b19;
    _deriv_basis[2][1] = _dz_b2; _deriv_basis[2][10] = _dz_b11; _deriv_basis[2][19] = _dz_b20;
    _deriv_basis[2][2] = _dz_b3; _deriv_basis[2][11] = _dz_b12; _deriv_basis[2][20] = _dz_b21;
    _deriv_basis[2][3] = _dz_b4; _deriv_basis[2][12] = _dz_b13; _deriv_basis[2][21] = _dz_b22;
    _deriv_basis[2][4] = _dz_b5; _deriv_basis[2][13] = _dz_b14; _deriv_basis[2][22] = _dz_b23;
    _deriv_basis[2][5] = _dz_b6; _deriv_basis[2][14] = _dz_b15; _deriv_basis[2][23] = _dz_b24;
    _deriv_basis[2][6] = _dz_b7; _deriv_basis[2][15] = _dz_b16; _deriv_basis[2][24] = _dz_b25;
    _deriv_basis[2][7] = _dz_b8; _deriv_basis[2][16] = _dz_b17; _deriv_basis[2][25] = _dz_b26;
    _deriv_basis[2][8] = _dz_b9; _deriv_basis[2][17] = _dz_b18; _deriv_basis[2][26] = _dz_b27;
    
    _secondDeriv_basis[0][0] = _dxx_b1; _secondDeriv_basis[0][9]  = _dxx_b10; _secondDeriv_basis[0][18] = _dxx_b19;
    _secondDeriv_basis[0][1] = _dxx_b2; _secondDeriv_basis[0][10] = _dxx_b11; _secondDeriv_basis[0][19] = _dxx_b20;
    _secondDeriv_basis[0][2] = _dxx_b3; _secondDeriv_basis[0][11] = _dxx_b12; _secondDeriv_basis[0][20] = _dxx_b21;
    _secondDeriv_basis[0][3] = _dxx_b4; _secondDeriv_basis[0][12] = _dxx_b13; _secondDeriv_basis[0][21] = _dxx_b22;
    _secondDeriv_basis[0][4] = _dxx_b5; _secondDeriv_basis[0][13] = _dxx_b14; _secondDeriv_basis[0][22] = _dxx_b23;
    _secondDeriv_basis[0][5] = _dxx_b6; _secondDeriv_basis[0][14] = _dxx_b15; _secondDeriv_basis[0][23] = _dxx_b24;
    _secondDeriv_basis[0][6] = _dxx_b7; _secondDeriv_basis[0][15] = _dxx_b16; _secondDeriv_basis[0][24] = _dxx_b25;
    _secondDeriv_basis[0][7] = _dxx_b8; _secondDeriv_basis[0][16] = _dxx_b17; _secondDeriv_basis[0][25] = _dxx_b26;
    _secondDeriv_basis[0][8] = _dxx_b9; _secondDeriv_basis[0][17] = _dxx_b18; _secondDeriv_basis[0][26] = _dxx_b27;

    _secondDeriv_basis[1][0] = _dyy_b1; _secondDeriv_basis[1][9]  = _dyy_b10; _secondDeriv_basis[1][18] = _dyy_b19;
    _secondDeriv_basis[1][1] = _dyy_b2; _secondDeriv_basis[1][10] = _dyy_b11; _secondDeriv_basis[1][19] = _dyy_b20;
    _secondDeriv_basis[1][2] = _dyy_b3; _secondDeriv_basis[1][11] = _dyy_b12; _secondDeriv_basis[1][20] = _dyy_b21;
    _secondDeriv_basis[1][3] = _dyy_b4; _secondDeriv_basis[1][12] = _dyy_b13; _secondDeriv_basis[1][21] = _dyy_b22;
    _secondDeriv_basis[1][4] = _dyy_b5; _secondDeriv_basis[1][13] = _dyy_b14; _secondDeriv_basis[1][22] = _dyy_b23;
    _secondDeriv_basis[1][5] = _dyy_b6; _secondDeriv_basis[1][14] = _dyy_b15; _secondDeriv_basis[1][23] = _dyy_b24;
    _secondDeriv_basis[1][6] = _dyy_b7; _secondDeriv_basis[1][15] = _dyy_b16; _secondDeriv_basis[1][24] = _dyy_b25;
    _secondDeriv_basis[1][7] = _dyy_b8; _secondDeriv_basis[1][16] = _dyy_b17; _secondDeriv_basis[1][25] = _dyy_b26;
    _secondDeriv_basis[1][8] = _dyy_b9; _secondDeriv_basis[1][17] = _dyy_b18; _secondDeriv_basis[1][26] = _dyy_b27;

    _secondDeriv_basis[2][0] = _dzz_b1; _secondDeriv_basis[2][9]  = _dzz_b10; _secondDeriv_basis[2][18] = _dzz_b19;
    _secondDeriv_basis[2][1] = _dzz_b2; _secondDeriv_basis[2][10] = _dzz_b11; _secondDeriv_basis[2][19] = _dzz_b20;
    _secondDeriv_basis[2][2] = _dzz_b3; _secondDeriv_basis[2][11] = _dzz_b12; _secondDeriv_basis[2][20] = _dzz_b21;
    _secondDeriv_basis[2][3] = _dzz_b4; _secondDeriv_basis[2][12] = _dzz_b13; _secondDeriv_basis[2][21] = _dzz_b22;
    _secondDeriv_basis[2][4] = _dzz_b5; _secondDeriv_basis[2][13] = _dzz_b14; _secondDeriv_basis[2][22] = _dzz_b23;
    _secondDeriv_basis[2][5] = _dzz_b6; _secondDeriv_basis[2][14] = _dzz_b15; _secondDeriv_basis[2][23] = _dzz_b24;
    _secondDeriv_basis[2][6] = _dzz_b7; _secondDeriv_basis[2][15] = _dzz_b16; _secondDeriv_basis[2][24] = _dzz_b25;
    _secondDeriv_basis[2][7] = _dzz_b8; _secondDeriv_basis[2][16] = _dzz_b17; _secondDeriv_basis[2][25] = _dzz_b26;
    _secondDeriv_basis[2][8] = _dzz_b9; _secondDeriv_basis[2][17] = _dzz_b18; _secondDeriv_basis[2][26] = _dzz_b27;
    
    _secondDeriv_basis[3][0] = _dxy_b1; _secondDeriv_basis[3][9]  = _dxy_b10; _secondDeriv_basis[3][18] = _dxy_b19;
    _secondDeriv_basis[3][1] = _dxy_b2; _secondDeriv_basis[3][10] = _dxy_b11; _secondDeriv_basis[3][19] = _dxy_b20;
    _secondDeriv_basis[3][2] = _dxy_b3; _secondDeriv_basis[3][11] = _dxy_b12; _secondDeriv_basis[3][20] = _dxy_b21;
    _secondDeriv_basis[3][3] = _dxy_b4; _secondDeriv_basis[3][12] = _dxy_b13; _secondDeriv_basis[3][21] = _dxy_b22;
    _secondDeriv_basis[3][4] = _dxy_b5; _secondDeriv_basis[3][13] = _dxy_b14; _secondDeriv_basis[3][22] = _dxy_b23;
    _secondDeriv_basis[3][5] = _dxy_b6; _secondDeriv_basis[3][14] = _dxy_b15; _secondDeriv_basis[3][23] = _dxy_b24;
    _secondDeriv_basis[3][6] = _dxy_b7; _secondDeriv_basis[3][15] = _dxy_b16; _secondDeriv_basis[3][24] = _dxy_b25;
    _secondDeriv_basis[3][7] = _dxy_b8; _secondDeriv_basis[3][16] = _dxy_b17; _secondDeriv_basis[3][25] = _dxy_b26;
    _secondDeriv_basis[3][8] = _dxy_b9; _secondDeriv_basis[3][17] = _dxy_b18; _secondDeriv_basis[3][26] = _dxy_b27;
    
    _secondDeriv_basis[4][0] = _dxz_b1; _secondDeriv_basis[4][9]  = _dxz_b10; _secondDeriv_basis[4][18] = _dxz_b19;
    _secondDeriv_basis[4][1] = _dxz_b2; _secondDeriv_basis[4][10] = _dxz_b11; _secondDeriv_basis[4][19] = _dxz_b20;
    _secondDeriv_basis[4][2] = _dxz_b3; _secondDeriv_basis[4][11] = _dxz_b12; _secondDeriv_basis[4][20] = _dxz_b21;
    _secondDeriv_basis[4][3] = _dxz_b4; _secondDeriv_basis[4][12] = _dxz_b13; _secondDeriv_basis[4][21] = _dxz_b22;
    _secondDeriv_basis[4][4] = _dxz_b5; _secondDeriv_basis[4][13] = _dxz_b14; _secondDeriv_basis[4][22] = _dxz_b23;
    _secondDeriv_basis[4][5] = _dxz_b6; _secondDeriv_basis[4][14] = _dxz_b15; _secondDeriv_basis[4][23] = _dxz_b24;
    _secondDeriv_basis[4][6] = _dxz_b7; _secondDeriv_basis[4][15] = _dxz_b16; _secondDeriv_basis[4][24] = _dxz_b25;
    _secondDeriv_basis[4][7] = _dxz_b8; _secondDeriv_basis[4][16] = _dxz_b17; _secondDeriv_basis[4][25] = _dxz_b26;
    _secondDeriv_basis[4][8] = _dxz_b9; _secondDeriv_basis[4][17] = _dxz_b18; _secondDeriv_basis[4][26] = _dxz_b27;
    
    _secondDeriv_basis[5][0] = _dyz_b1; _secondDeriv_basis[5][9]  = _dyz_b10; _secondDeriv_basis[5][18] = _dyz_b19;
    _secondDeriv_basis[5][1] = _dyz_b2; _secondDeriv_basis[5][10] = _dyz_b11; _secondDeriv_basis[5][19] = _dyz_b20;
    _secondDeriv_basis[5][2] = _dyz_b3; _secondDeriv_basis[5][11] = _dyz_b12; _secondDeriv_basis[5][20] = _dyz_b21;
    _secondDeriv_basis[5][3] = _dyz_b4; _secondDeriv_basis[5][12] = _dyz_b13; _secondDeriv_basis[5][21] = _dyz_b22;
    _secondDeriv_basis[5][4] = _dyz_b5; _secondDeriv_basis[5][13] = _dyz_b14; _secondDeriv_basis[5][22] = _dyz_b23;
    _secondDeriv_basis[5][5] = _dyz_b6; _secondDeriv_basis[5][14] = _dyz_b15; _secondDeriv_basis[5][23] = _dyz_b24;
    _secondDeriv_basis[5][6] = _dyz_b7; _secondDeriv_basis[5][15] = _dyz_b16; _secondDeriv_basis[5][24] = _dyz_b25;
    _secondDeriv_basis[5][7] = _dyz_b8; _secondDeriv_basis[5][16] = _dyz_b17; _secondDeriv_basis[5][25] = _dyz_b26;
    _secondDeriv_basis[5][8] = _dyz_b9; _secondDeriv_basis[5][17] = _dyz_b18; _secondDeriv_basis[5][26] = _dyz_b27;

    this->initializeQuadCache( );
    // BaseFunctionSetInterface caches only values and gradients (second derivatives are zero in the linear case)
    initializeExtendedQuadCache( );
        
  }

  enum { numBaseFuncs = 27 };

  RealType evaluate ( int BaseFuncNum, const Vec3<RealType> &RefCoord ) const {
    return _basis[BaseFuncNum] ( RefCoord );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec3<RealType>, Vec3<RealType>, 27, QuadRuleType, BaseFunctionSetMultiQuad<RealType, qc::QC_3D, QuadRuleType> >::evaluate ( BaseFuncNum, QuadPoint );
  }

  void evaluateGradient ( int BaseFuncNum, const Vec3<RealType> &RefCoord, Vec3<RealType> &Gradient ) const {
    Gradient[0] = _deriv_basis[0][BaseFuncNum] ( RefCoord ) / _h;
    Gradient[1] = _deriv_basis[1][BaseFuncNum] ( RefCoord ) / _h;
    Gradient[2] = _deriv_basis[2][BaseFuncNum] ( RefCoord ) / _h;
  }
  
  inline const Vec3<RealType>& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec3<RealType>, Vec3<RealType>, 27, QuadRuleType, BaseFunctionSetMultiQuad<RealType, qc::QC_3D, QuadRuleType> >::evaluateGradient ( BaseFuncNum, QuadPoint );
  }


void evaluateHessian ( int BaseFuncNum, const Vec3<RealType> &RefCoord, Mat<3,3,RealType> &Hessian ) const {
    RealType h2 = aol::Sqr(_h);
    
    Hessian[0][0] = _secondDeriv_basis[0][BaseFuncNum] ( RefCoord ) / h2;
    Hessian[0][1] = _secondDeriv_basis[3][BaseFuncNum] ( RefCoord ) / h2;
    Hessian[0][2] = _secondDeriv_basis[4][BaseFuncNum] ( RefCoord ) / h2;
    
    Hessian[1][0] = _secondDeriv_basis[3][BaseFuncNum] ( RefCoord ) / h2;
    Hessian[1][1] = _secondDeriv_basis[1][BaseFuncNum] ( RefCoord ) / h2;
    Hessian[1][2] = _secondDeriv_basis[5][BaseFuncNum] ( RefCoord ) / h2;
    
    Hessian[2][0] = _secondDeriv_basis[4][BaseFuncNum] ( RefCoord ) / h2;
    Hessian[2][1] = _secondDeriv_basis[5][BaseFuncNum] ( RefCoord ) / h2;
    Hessian[2][2] = _secondDeriv_basis[2][BaseFuncNum] ( RefCoord ) / h2;
   
  }

  inline const Mat<3,3,RealType>& evaluateHessian ( int BaseFuncNum, int QuadPoint ) const {
    return basisQuadHessians[BaseFuncNum][QuadPoint];
  }

protected:
  void initializeExtendedQuadCache( ) {
    for ( int b = 0; b < numBaseFuncs; b++ ) {
      for ( int i = 0; i < QuadRuleType::numQuadPoints; i++ ) {
        evaluateHessian ( b, this->_quadRule.getRefCoord ( i ), basisQuadHessians[b][i] );
      }
    }
  }
  /**** cache extension for Hessians ****/
  Matrix33<RealType>  basisQuadHessians[numBaseFuncs][QuadRuleType::numQuadPoints];
   
};





/**
 *  The basefunctionset for biquartic elements in 2d.
 */
template <typename RealType, class QuadRuleType>
class BaseFunctionSetMultiQuart<RealType, qc::QC_2D, QuadRuleType> :
public BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 25, QuadRuleType, BaseFunctionSetMultiQuart<RealType, qc::QC_2D, QuadRuleType> > {

  static inline RealType _b1_1d ( RealType x ) { return 1./3. * ( 3. + x * ( -25. + x * (   70. + x * (  -80. +  32. * x ) ) ) ); }
  static inline RealType _b2_1d ( RealType x ) { return 1./3. * (      x * (  48. + x * ( -208. + x * (  288. - 128. * x ) ) ) ); }
  static inline RealType _b3_1d ( RealType x ) { return 1./3. * (      x * ( -36. + x * (  228. + x * ( -384. + 192. * x ) ) ) ); }
  static inline RealType _b4_1d ( RealType x ) { return 1./3. * (      x * (  16. + x * ( -112. + x * (  224. - 128. * x ) ) ) ); }
  static inline RealType _b5_1d ( RealType x ) { return 1./3. * (      x * (  -3. + x * (   22. + x * (  -48. +  32. * x ) ) ) ); }

  static inline RealType _d_b1_1d ( RealType x ) { return 1./3. * ( -25. + x * (  140. + x * (  -240. + 128. * x ) ) ); }
  static inline RealType _d_b2_1d ( RealType x ) { return 1./3. * (  48. + x * ( -416. + x * (   864. - 512. * x ) ) ); }
  static inline RealType _d_b3_1d ( RealType x ) { return 1./3. * ( -36. + x * (  456. + x * ( -1152. + 768. * x ) ) ); }
  static inline RealType _d_b4_1d ( RealType x ) { return 1./3. * (  16. + x * ( -224. + x * (   672. - 512. * x ) ) ); }
  static inline RealType _d_b5_1d ( RealType x ) { return 1./3. * (  -3. + x * (   44. + x * (  -144. + 128. * x ) ) ); }


  static RealType _b1    ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ); }
  static RealType _b2    ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ); }
  static RealType _b3    ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ); }
  static RealType _b4    ( const Vec2<RealType> &RefCoord ) { return _b4_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ); }
  static RealType _b5    ( const Vec2<RealType> &RefCoord ) { return _b5_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ); }

  static RealType _b6    ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ); }
  static RealType _b7    ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ); }
  static RealType _b8    ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ); }
  static RealType _b9    ( const Vec2<RealType> &RefCoord ) { return _b4_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ); }
  static RealType _b10   ( const Vec2<RealType> &RefCoord ) { return _b5_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ); }

  static RealType _b11   ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ); }
  static RealType _b12   ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ); }
  static RealType _b13   ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ); }
  static RealType _b14   ( const Vec2<RealType> &RefCoord ) { return _b4_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ); }
  static RealType _b15   ( const Vec2<RealType> &RefCoord ) { return _b5_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ); }

  static RealType _b16   ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b4_1d ( RefCoord[1] ); }
  static RealType _b17   ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b4_1d ( RefCoord[1] ); }
  static RealType _b18   ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b4_1d ( RefCoord[1] ); }
  static RealType _b19   ( const Vec2<RealType> &RefCoord ) { return _b4_1d ( RefCoord[0] ) * _b4_1d ( RefCoord[1] ); }
  static RealType _b20   ( const Vec2<RealType> &RefCoord ) { return _b5_1d ( RefCoord[0] ) * _b4_1d ( RefCoord[1] ); }

  static RealType _b21   ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _b5_1d ( RefCoord[1] ); }
  static RealType _b22   ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _b5_1d ( RefCoord[1] ); }
  static RealType _b23   ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _b5_1d ( RefCoord[1] ); }
  static RealType _b24   ( const Vec2<RealType> &RefCoord ) { return _b4_1d ( RefCoord[0] ) * _b5_1d ( RefCoord[1] ); }
  static RealType _b25   ( const Vec2<RealType> &RefCoord ) { return _b5_1d ( RefCoord[0] ) * _b5_1d ( RefCoord[1] ); }


  static RealType _dx_b1    ( const Vec2<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ); }
  static RealType _dx_b2    ( const Vec2<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ); }
  static RealType _dx_b3    ( const Vec2<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ); }
  static RealType _dx_b4    ( const Vec2<RealType> &RefCoord ) { return _d_b4_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ); }
  static RealType _dx_b5    ( const Vec2<RealType> &RefCoord ) { return _d_b5_1d ( RefCoord[0] ) * _b1_1d ( RefCoord[1] ); }

  static RealType _dx_b6    ( const Vec2<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ); }
  static RealType _dx_b7    ( const Vec2<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ); }
  static RealType _dx_b8    ( const Vec2<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ); }
  static RealType _dx_b9    ( const Vec2<RealType> &RefCoord ) { return _d_b4_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ); }
  static RealType _dx_b10   ( const Vec2<RealType> &RefCoord ) { return _d_b5_1d ( RefCoord[0] ) * _b2_1d ( RefCoord[1] ); }

  static RealType _dx_b11   ( const Vec2<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ); }
  static RealType _dx_b12   ( const Vec2<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ); }
  static RealType _dx_b13   ( const Vec2<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ); }
  static RealType _dx_b14   ( const Vec2<RealType> &RefCoord ) { return _d_b4_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ); }
  static RealType _dx_b15   ( const Vec2<RealType> &RefCoord ) { return _d_b5_1d ( RefCoord[0] ) * _b3_1d ( RefCoord[1] ); }

  static RealType _dx_b16   ( const Vec2<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b4_1d ( RefCoord[1] ); }
  static RealType _dx_b17   ( const Vec2<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b4_1d ( RefCoord[1] ); }
  static RealType _dx_b18   ( const Vec2<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b4_1d ( RefCoord[1] ); }
  static RealType _dx_b19   ( const Vec2<RealType> &RefCoord ) { return _d_b4_1d ( RefCoord[0] ) * _b4_1d ( RefCoord[1] ); }
  static RealType _dx_b20   ( const Vec2<RealType> &RefCoord ) { return _d_b5_1d ( RefCoord[0] ) * _b4_1d ( RefCoord[1] ); }

  static RealType _dx_b21   ( const Vec2<RealType> &RefCoord ) { return _d_b1_1d ( RefCoord[0] ) * _b5_1d ( RefCoord[1] ); }
  static RealType _dx_b22   ( const Vec2<RealType> &RefCoord ) { return _d_b2_1d ( RefCoord[0] ) * _b5_1d ( RefCoord[1] ); }
  static RealType _dx_b23   ( const Vec2<RealType> &RefCoord ) { return _d_b3_1d ( RefCoord[0] ) * _b5_1d ( RefCoord[1] ); }
  static RealType _dx_b24   ( const Vec2<RealType> &RefCoord ) { return _d_b4_1d ( RefCoord[0] ) * _b5_1d ( RefCoord[1] ); }
  static RealType _dx_b25   ( const Vec2<RealType> &RefCoord ) { return _d_b5_1d ( RefCoord[0] ) * _b5_1d ( RefCoord[1] ); }

  static RealType _dy_b1    ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ); }
  static RealType _dy_b2    ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ); }
  static RealType _dy_b3    ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ); }
  static RealType _dy_b4    ( const Vec2<RealType> &RefCoord ) { return _b4_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ); }
  static RealType _dy_b5    ( const Vec2<RealType> &RefCoord ) { return _b5_1d ( RefCoord[0] ) * _d_b1_1d ( RefCoord[1] ); }

  static RealType _dy_b6    ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ); }
  static RealType _dy_b7    ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ); }
  static RealType _dy_b8    ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ); }
  static RealType _dy_b9    ( const Vec2<RealType> &RefCoord ) { return _b4_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ); }
  static RealType _dy_b10   ( const Vec2<RealType> &RefCoord ) { return _b5_1d ( RefCoord[0] ) * _d_b2_1d ( RefCoord[1] ); }

  static RealType _dy_b11   ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ); }
  static RealType _dy_b12   ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ); }
  static RealType _dy_b13   ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ); }
  static RealType _dy_b14   ( const Vec2<RealType> &RefCoord ) { return _b4_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ); }
  static RealType _dy_b15   ( const Vec2<RealType> &RefCoord ) { return _b5_1d ( RefCoord[0] ) * _d_b3_1d ( RefCoord[1] ); }

  static RealType _dy_b16   ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b4_1d ( RefCoord[1] ); }
  static RealType _dy_b17   ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b4_1d ( RefCoord[1] ); }
  static RealType _dy_b18   ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b4_1d ( RefCoord[1] ); }
  static RealType _dy_b19   ( const Vec2<RealType> &RefCoord ) { return _b4_1d ( RefCoord[0] ) * _d_b4_1d ( RefCoord[1] ); }
  static RealType _dy_b20   ( const Vec2<RealType> &RefCoord ) { return _b5_1d ( RefCoord[0] ) * _d_b4_1d ( RefCoord[1] ); }

  static RealType _dy_b21   ( const Vec2<RealType> &RefCoord ) { return _b1_1d ( RefCoord[0] ) * _d_b5_1d ( RefCoord[1] ); }
  static RealType _dy_b22   ( const Vec2<RealType> &RefCoord ) { return _b2_1d ( RefCoord[0] ) * _d_b5_1d ( RefCoord[1] ); }
  static RealType _dy_b23   ( const Vec2<RealType> &RefCoord ) { return _b3_1d ( RefCoord[0] ) * _d_b5_1d ( RefCoord[1] ); }
  static RealType _dy_b24   ( const Vec2<RealType> &RefCoord ) { return _b4_1d ( RefCoord[0] ) * _d_b5_1d ( RefCoord[1] ); }
  static RealType _dy_b25   ( const Vec2<RealType> &RefCoord ) { return _b5_1d ( RefCoord[0] ) * _d_b5_1d ( RefCoord[1] ); }


  typedef RealType ( *BASIS_FUNC_TYPE ) ( const Vec2<RealType> &RefCoord );
  BASIS_FUNC_TYPE _deriv_basis[2][25];
  BASIS_FUNC_TYPE _basis[25];

  const RealType _h;

public:
  BaseFunctionSetMultiQuart( const RealType H ) : _h ( H ) {
    _basis[ 0] = _b1;
    _basis[ 1] = _b2;
    _basis[ 2] = _b3;
    _basis[ 3] = _b4;
    _basis[ 4] = _b5;
    _basis[ 5] = _b6;
    _basis[ 6] = _b7;
    _basis[ 7] = _b8;
    _basis[ 8] = _b9;
    _basis[ 9] = _b10;
    _basis[10] = _b11;
    _basis[11] = _b12;
    _basis[12] = _b13;
    _basis[13] = _b14;
    _basis[14] = _b15;
    _basis[15] = _b16;
    _basis[16] = _b17;
    _basis[17] = _b18;
    _basis[18] = _b19;
    _basis[19] = _b20;
    _basis[20] = _b21;
    _basis[21] = _b22;
    _basis[22] = _b23;
    _basis[23] = _b24;
    _basis[24] = _b25;

    _deriv_basis[0][ 0] = _dx_b1;
    _deriv_basis[0][ 1] = _dx_b2;
    _deriv_basis[0][ 2] = _dx_b3;
    _deriv_basis[0][ 3] = _dx_b4;
    _deriv_basis[0][ 4] = _dx_b5;
    _deriv_basis[0][ 5] = _dx_b6;
    _deriv_basis[0][ 6] = _dx_b7;
    _deriv_basis[0][ 7] = _dx_b8;
    _deriv_basis[0][ 8] = _dx_b9;
    _deriv_basis[0][ 9] = _dx_b10;
    _deriv_basis[0][10] = _dx_b11;
    _deriv_basis[0][11] = _dx_b12;
    _deriv_basis[0][12] = _dx_b13;
    _deriv_basis[0][13] = _dx_b14;
    _deriv_basis[0][14] = _dx_b15;
    _deriv_basis[0][15] = _dx_b16;
    _deriv_basis[0][16] = _dx_b17;
    _deriv_basis[0][17] = _dx_b18;
    _deriv_basis[0][18] = _dx_b19;
    _deriv_basis[0][19] = _dx_b20;
    _deriv_basis[0][20] = _dx_b21;
    _deriv_basis[0][21] = _dx_b22;
    _deriv_basis[0][22] = _dx_b23;
    _deriv_basis[0][23] = _dx_b24;
    _deriv_basis[0][24] = _dx_b25;

    _deriv_basis[1][ 0] = _dy_b1;
    _deriv_basis[1][ 1] = _dy_b2;
    _deriv_basis[1][ 2] = _dy_b3;
    _deriv_basis[1][ 3] = _dy_b4;
    _deriv_basis[1][ 4] = _dy_b5;
    _deriv_basis[1][ 5] = _dy_b6;
    _deriv_basis[1][ 6] = _dy_b7;
    _deriv_basis[1][ 7] = _dy_b8;
    _deriv_basis[1][ 8] = _dy_b9;
    _deriv_basis[1][ 9] = _dy_b10;
    _deriv_basis[1][10] = _dy_b11;
    _deriv_basis[1][11] = _dy_b12;
    _deriv_basis[1][12] = _dy_b13;
    _deriv_basis[1][13] = _dy_b14;
    _deriv_basis[1][14] = _dy_b15;
    _deriv_basis[1][15] = _dy_b16;
    _deriv_basis[1][16] = _dy_b17;
    _deriv_basis[1][17] = _dy_b18;
    _deriv_basis[1][18] = _dy_b19;
    _deriv_basis[1][19] = _dy_b20;
    _deriv_basis[1][20] = _dy_b21;
    _deriv_basis[1][21] = _dy_b22;
    _deriv_basis[1][22] = _dy_b23;
    _deriv_basis[1][23] = _dy_b24;
    _deriv_basis[1][24] = _dy_b25;

    this->initializeQuadCache( );
  }

  enum { numBaseFuncs = 25 };

  void evaluateGradient ( int BaseFuncNum, const Vec2<RealType> &RefCoord, Vec2<RealType> &Gradient ) const {
    Gradient[0] = _deriv_basis[0][BaseFuncNum] ( RefCoord ) / _h;
    Gradient[1] = _deriv_basis[1][BaseFuncNum] ( RefCoord ) / _h;
  }

  inline const Vec2<RealType>& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 25, QuadRuleType, BaseFunctionSetMultiQuart<RealType, qc::QC_2D, QuadRuleType> >::evaluateGradient ( BaseFuncNum, QuadPoint );
  }

  RealType evaluate ( int BaseFuncNum, const Vec2<RealType> &RefCoord ) const {
    return _basis[BaseFuncNum] ( RefCoord );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 25, QuadRuleType, BaseFunctionSetMultiQuart<RealType, qc::QC_2D, QuadRuleType> >::evaluate ( BaseFuncNum, QuadPoint );
  }

protected:
};


template <typename RealType, qc::Dimension Dim, class QuadRuleType>
class BaseFunctionSetMultiLinBubble :
      public BaseFunctionSetInterface<RealType, Vec3<RealType>, Vec3<RealType>, 5, QuadRuleType, BaseFunctionSetMultiLin<RealType, qc::QC_2D, QuadRuleType> > {
};

//! The basefunctionset for bilinear elements in 2d.
template <typename RealType, class QuadRuleType>
class BaseFunctionSetMultiLinBubble<RealType, qc::QC_2D, QuadRuleType> :
      public BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 5, QuadRuleType, BaseFunctionSetMultiLinBubble<RealType, qc::QC_2D, QuadRuleType> > {

  static RealType _dx_b1   ( const Vec2<RealType> &RefCoord ) { return RefCoord[1] - 1.; }
  static RealType _dx_b2   ( const Vec2<RealType> &RefCoord ) { return 1. - RefCoord[1]; }
  static RealType _dx_b3   ( const Vec2<RealType> &RefCoord ) { return -RefCoord[1]; }
  static RealType _dx_b4   ( const Vec2<RealType> &RefCoord ) { return RefCoord[1]; }
  static RealType _dx_b5   ( const Vec2<RealType> &RefCoord ) { return ( RefCoord[0]*2. - 1. ) * ( Sqr ( RefCoord[1] ) - RefCoord[1] ); }

  static RealType _dy_b1   ( const Vec2<RealType> &RefCoord ) { return RefCoord[0] - 1.; }
  static RealType _dy_b2   ( const Vec2<RealType> &RefCoord ) { return -RefCoord[0]; }
  static RealType _dy_b3   ( const Vec2<RealType> &RefCoord ) { return 1. - RefCoord[0]; }
  static RealType _dy_b4   ( const Vec2<RealType> &RefCoord ) { return RefCoord[0]; }
  static RealType _dy_b5   ( const Vec2<RealType> &RefCoord ) { return ( RefCoord[1]*2. - 1. ) * ( Sqr ( RefCoord[0] ) - RefCoord[0] ); }

  static RealType _b1   ( const Vec2<RealType> &RefCoord ) { return ( 1. - RefCoord[0] ) * ( 1 - RefCoord[1] ); }
  static RealType _b2   ( const Vec2<RealType> &RefCoord ) { return RefCoord[0]* ( 1. - RefCoord[1] ); }
  static RealType _b3   ( const Vec2<RealType> &RefCoord ) { return ( 1. - RefCoord[0] ) *RefCoord[1]; }
  static RealType _b4   ( const Vec2<RealType> &RefCoord ) { return RefCoord[0]*RefCoord[1]; }
  static RealType _b5   ( const Vec2<RealType> &RefCoord ) { return RefCoord[0]* ( 1. - RefCoord[0] ) *RefCoord[1]* ( 1. - RefCoord[1] ); }


  typedef RealType ( *BASIS_FUNC_TYPE ) ( const Vec2<RealType> &RefCoord );
  BASIS_FUNC_TYPE _deriv_basis[2][5];
  BASIS_FUNC_TYPE _basis[5];

public:
  BaseFunctionSetMultiLinBubble( ) {
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;
    _basis[4] = _b5;


    _deriv_basis[0][0] = _dx_b1;
    _deriv_basis[0][1] = _dx_b2;
    _deriv_basis[0][2] = _dx_b3;
    _deriv_basis[0][3] = _dx_b4;
    _deriv_basis[0][4] = _dx_b5;

    _deriv_basis[1][0] = _dy_b1;
    _deriv_basis[1][1] = _dy_b2;
    _deriv_basis[1][2] = _dy_b3;
    _deriv_basis[1][3] = _dy_b4;
    _deriv_basis[1][4] = _dy_b5;

    this->initializeQuadCache( );
  }

  enum { numBaseFuncs = 4 };

  void evaluateGradient ( int BaseFuncNum, const Vec2<RealType> &RefCoord, Vec2<RealType> &Gradient ) const {
    Gradient[0] = _deriv_basis[0][BaseFuncNum] ( RefCoord );
    Gradient[1] = _deriv_basis[1][BaseFuncNum] ( RefCoord );
  }

  RealType evaluate ( int BaseFuncNum, const Vec2<RealType> &RefCoord ) const {
    return _basis[BaseFuncNum] ( RefCoord );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 5, QuadRuleType, BaseFunctionSetMultiLinBubble<RealType, qc::QC_2D, QuadRuleType> >::evaluate ( BaseFuncNum, QuadPoint );
  }

  inline const Vec2<RealType>& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 5, QuadRuleType, BaseFunctionSetMultiLinBubble<RealType, qc::QC_2D, QuadRuleType> >::evaluateGradient ( BaseFuncNum, QuadPoint );
  }

protected:

};

template <typename _RealType, int NumQuadPoints>
class Quadrature1D {
public:
  enum { numQuadPoints = NumQuadPoints };

  Quadrature1D() {}

  inline const Vec<1, _RealType> &getRefCoord ( int QuadPoint ) const {
    return _points [QuadPoint];
  }

protected:
  Vec<1, _RealType> _points [ numQuadPoints ];
};

template <typename _RealType>
class TrapeziumQuadrature : public Quadrature1D<_RealType,2> {
public:
  TrapeziumQuadrature( ) {
    this->_points[0][0] = 0.0;
    this->_points[1][0] = 1.0;
  }

  inline _RealType getWeight ( int /*QuadPoint*/ ) const {
    return .5;
  }
};

template <typename _RealType>
class SimpsonQuadrature : public Quadrature1D<_RealType,3> {
public:
  SimpsonQuadrature( ) {
    this->_points[0][0] = 0.0;
    this->_points[1][0] = 0.5;
    this->_points[2][0] = 1.0;
  }

  inline _RealType getWeight ( int QuadPoint ) const {
    static const _RealType _weights [Quadrature1D<_RealType,3>::numQuadPoints] = { .16666666666666666, .6666666666666666, .1666666666666666 };
    return _weights [QuadPoint];
  }
};

template <typename _RealType, typename Quadrature1DType>
class Quadrature2D {
public:
  typedef _RealType RealType;
  static const qc::Dimension Dim = qc::QC_2D;

  enum { numQuadPoints = Quadrature1DType::numQuadPoints * Quadrature1DType::numQuadPoints };

  Quadrature2D( ) {

    Quadrature1DType quad1D;

    for ( int k = 0, i = 0; i < Quadrature1DType::numQuadPoints; i++ ) {
      for ( int j = 0; j < Quadrature1DType::numQuadPoints; j++, k++ ) {
        _points[k][0] = quad1D.getRefCoord ( i ).get ( 0 );
        _points[k][1] = quad1D.getRefCoord ( j ).get ( 0 );
        _weights[k] = quad1D.getWeight ( i ) * quad1D.getWeight ( j );
      }
    }
  }

  inline const Vec2<_RealType> &getRefCoord ( int QuadPoint ) const {
    return _points[ QuadPoint ];
  }

  inline _RealType getWeight ( int QuadPoint ) const {
    return _weights[ QuadPoint ];
  }

protected:
  Vec2<_RealType> _points [ numQuadPoints ];
  _RealType       _weights[ numQuadPoints ];
};

template <typename _RealType, typename Quadrature1DType>
class Quadrature3D {
public:
  typedef _RealType RealType;
  static const qc::Dimension Dim = qc::QC_3D;

  enum { numQuadPoints = Quadrature1DType::numQuadPoints * Quadrature1DType::numQuadPoints * Quadrature1DType::numQuadPoints };

  Quadrature3D( ) {

    Quadrature1DType quad1D;

    for ( int k = 0, i = 0; i < Quadrature1DType::numQuadPoints; i++ ) {
      for ( int j = 0; j < Quadrature1DType::numQuadPoints; j++ ) {
        for ( int l = 0; l < Quadrature1DType::numQuadPoints; l++, k++ ) {
          _points[k][0] = quad1D.getRefCoord ( i ).get ( 0 );
          _points[k][1] = quad1D.getRefCoord ( j ).get ( 0 );
          _points[k][2] = quad1D.getRefCoord ( l ).get ( 0 );

          _weights[k] = quad1D.getWeight ( i ) * quad1D.getWeight ( j ) * quad1D.getWeight ( l );
        }
      }
    }
  }

  inline const Vec3<_RealType> &getRefCoord ( int QuadPoint ) const {
    return _points[ QuadPoint ];
  }

  inline _RealType getWeight ( int QuadPoint ) const {
    return _weights[ QuadPoint ];
  }

protected:
  Vec3<_RealType> _points [ numQuadPoints ];
  _RealType       _weights[ numQuadPoints ];
};


// 1D Gauss points and weights computed with Maple worksheet from http://oregonstate.edu/~peterseb/mth351/351index.html

template <typename _RealType, qc::Dimension _Dim, int _Order>
class GaussQuadrature {};

template <typename _RealType>
class GaussQuadrature<_RealType, qc::QC_1D, 1> {
public:
  enum { numQuadPoints = 1 };

  GaussQuadrature( ) :
  _point ( Scalar<_RealType> (0.5) ) {}


  inline const Vec<1, _RealType> &getRefCoord ( int /*QuadPoint*/ ) const {
    return _point;
  }

  inline _RealType getWeight ( int /*QuadPoint*/ ) const {
    return 1;
  }

protected:
  Vec<1, _RealType> _point;
};

template <typename _RealType>
class GaussQuadrature<_RealType, qc::QC_1D, 3> {
public:
  enum { numQuadPoints = 2 };

  GaussQuadrature( ) {
    _points[0][0] = 0.2113248654051871;
    _points[1][0] = 0.7886751345948129;
  }

  inline const Vec<1, _RealType> &getRefCoord ( int QuadPoint ) const {
    return _points [QuadPoint];
  }

  inline _RealType getWeight ( int /*QuadPoint*/ ) const {
    return 0.5;
  }

protected:
  Vec<1, _RealType> _points [ numQuadPoints ];
};

template <typename _RealType>
class GaussQuadrature<_RealType, qc::QC_1D, 5> {
public:
  enum { numQuadPoints = 3 };

  GaussQuadrature( ) {
    _points[0][0] = 0.1127016653792583;
    _points[1][0] = 0.5;
    _points[2][0] = 0.8872983346207417;
  }

  inline const Vec<1, _RealType> &getRefCoord ( int QuadPoint ) const {
    return _points [QuadPoint];
  }

  inline _RealType getWeight ( int QuadPoint ) const {
    static const _RealType _weights [numQuadPoints] = { .2777777777777778, .4444444444444444, .2777777777777778 };
    return _weights [QuadPoint];
  }

protected:
  Vec<1, _RealType> _points [ numQuadPoints ];
};

template <typename _RealType>
class GaussQuadrature<_RealType, qc::QC_1D, 7> {
public:
  enum { numQuadPoints = 4 };

  GaussQuadrature( ) {
    _points[0][0] = 0.6943184420297371e-1;
    _points[1][0] = 0.3300094782075719;
    _points[2][0] = 0.6699905217924281;
    _points[3][0] = 0.9305681557970263;
  }

  inline const Vec<1, _RealType> &getRefCoord ( int QuadPoint ) const {
    return _points [QuadPoint];
  }

  inline _RealType getWeight ( int QuadPoint ) const {
    static const _RealType _weights [numQuadPoints] = { .1739274225687269, .3260725774312731, .3260725774312731, .1739274225687269 };
    return _weights [QuadPoint];
  }

protected:
  Vec<1, _RealType> _points [ numQuadPoints ];
};

template <typename _RealType>
class GaussQuadrature<_RealType, qc::QC_1D, 9> {
public:
  enum { numQuadPoints = 5 };

 GaussQuadrature( ) {
    _points[0][0] = 0.4691007703066800e-1;
    _points[1][0] = 0.2307653449471585;
    _points[2][0] = 0.5;
    _points[3][0] = 0.7692346550528415;
    _points[4][0] = 0.9530899229693320;
  }

  inline const Vec<1, _RealType> &getRefCoord ( int QuadPoint ) const {
    return _points [QuadPoint];
  }

  inline _RealType getWeight ( int QuadPoint ) const {
    static const _RealType _weights [numQuadPoints] = { .1184634425280945, .2393143352496832, .2844444444444444, .2393143352496832, .1184634425280945 };
    return _weights [QuadPoint];
  }

protected:
  Vec<1, _RealType> _points [ numQuadPoints ];
};

template <typename _RealType>
class GaussQuadrature<_RealType, qc::QC_1D, 11> {
public:
  enum { numQuadPoints = 6 };

  GaussQuadrature( ) {
    _points[0][0] = 0.3376524289842399e-1;
    _points[1][0] = 0.1693953067668677;
    _points[2][0] = 0.3806904069584015;
    _points[3][0] = 0.6193095930415985;
    _points[4][0] = 0.8306046932331323;
    _points[5][0] = 0.9662347571015760;
  }

  inline const Vec<1, _RealType> &getRefCoord ( int QuadPoint ) const {
    return _points [QuadPoint];
  }

  inline _RealType getWeight ( int QuadPoint ) const {
    static const _RealType _weights [numQuadPoints] = {
                                                       .8566224618958517e-1, .1803807865240693, .2339569672863455,
                                                       .2339569672863455, .1803807865240693, .8566224618958517e-1
                                                     };
    return _weights [QuadPoint];
  }

protected:
  Vec<1, _RealType> _points [ numQuadPoints ];
};

template <typename _RealType>
class GaussQuadrature<_RealType, qc::QC_1D, 13> {
public:
  enum { numQuadPoints = 7 };

  GaussQuadrature( ) {
    _points[0][0] = 0.2544604382862074e-1;
    _points[1][0] = 0.1292344072003028;
    _points[2][0] = 0.2970774243113014;
    _points[3][0] = 0.5;
    _points[4][0] = 0.7029225756886986;
    _points[5][0] = 0.8707655927996972;
    _points[6][0] = 0.9745539561713793;
  }

  inline const Vec<1, _RealType> &getRefCoord ( int QuadPoint ) const {
    return _points [QuadPoint];
  }

  inline _RealType getWeight ( int QuadPoint ) const {
    static const _RealType _weights [numQuadPoints] = {
                                                       .6474248308443485e-1, .1398526957446383, .1909150252525595, .2089795918367347,
                                                       .1909150252525595, .1398526957446383, .6474248308443485e-1
                                                     };
    return _weights [QuadPoint];
  }

protected:
  Vec<1, _RealType> _points [ numQuadPoints ];
};

template <typename _RealType>
class GaussQuadrature<_RealType, qc::QC_1D, 15> {
public:
  enum { numQuadPoints = 8 };

  GaussQuadrature( ) {
    _points[0][0] = 0.1985507175123188e-1;
    _points[1][0] = 0.1016667612931866;
    _points[2][0] = 0.2372337950418355;
    _points[3][0] = 0.4082826787521751;
    _points[4][0] = 0.5917173212478249;
    _points[5][0] = 0.7627662049581645;
    _points[6][0] = 0.8983332387068134;
    _points[7][0] = 0.9801449282487681;
  }

  inline const Vec<1, _RealType> &getRefCoord ( int QuadPoint ) const {
    return _points [QuadPoint];
  }

  inline _RealType getWeight ( int QuadPoint ) const {
    static const _RealType _weights [numQuadPoints] = {
                                                       .5061426814518813e-1, .1111905172266872, .1568533229389436, .1813418916891810,
                                                       .1813418916891810, .1568533229389436, .1111905172266872, .5061426814518813e-1
                                                     };
    return _weights [QuadPoint];
  }

protected:
  Vec<1, _RealType> _points [ numQuadPoints ];
};

template <typename _RealType>
class GaussQuadrature<_RealType, qc::QC_1D, 17> {
public:
  enum { numQuadPoints = 9 };

  GaussQuadrature( ) {
    _points[0][0] = 0.1591988024618696e-1;
    _points[1][0] = 0.8198444633668210e-1;
    _points[2][0] = 0.1933142836497048;
    _points[3][0] = 0.3378732882980955;
    _points[4][0] = 0.5;
    _points[5][0] = 0.6621267117019045;
    _points[6][0] = 0.8066857163502952;
    _points[7][0] = 0.9180155536633179;
    _points[8][0] = 0.9840801197538130;
  }

  inline const Vec<1, _RealType> &getRefCoord ( int QuadPoint ) const {
    return _points [QuadPoint];
  }

  inline _RealType getWeight ( int QuadPoint ) const {
    static const _RealType _weights [numQuadPoints] = {
                                                       .4063719418078721e-1, .9032408034742870e-1, .1303053482014677,
                                                       .1561735385200014, .1651196775006299, .1561735385200014,
                                                       .1303053482014677, .9032408034742870e-1, .4063719418078721e-1
                                                     };
    return _weights [QuadPoint];
  }

protected:
  Vec<1, _RealType> _points [ numQuadPoints ];
};

template <typename _RealType>
class GaussQuadrature<_RealType, qc::QC_1D, 19> {
public:
  enum { numQuadPoints = 10 };

  GaussQuadrature( ) {
    _points[0][0] = 0.1304673574141414e-1;
    _points[1][0] = 0.6746831665550774e-1;
    _points[2][0] = 0.1602952158504878;
    _points[3][0] = 0.2833023029353764;
    _points[4][0] = 0.4255628305091844;
    _points[5][0] = 0.5744371694908156;
    _points[6][0] = 0.7166976970646236;
    _points[7][0] = 0.8397047841495122;
    _points[8][0] = 0.9325316833444923;
    _points[9][0] = 0.9869532642585859;
  }

  inline const Vec<1, _RealType> &getRefCoord ( int QuadPoint ) const {
    return _points [QuadPoint];
  }

  inline _RealType getWeight ( int QuadPoint ) const {
    static const _RealType _weights [numQuadPoints] = {
                                                       .3333567215434407e-1, .7472567457529030e-1, .1095431812579910, .1346333596549982,
                                                       .1477621123573764, .1477621123573764, .1346333596549982,
                                                       .1095431812579910, .7472567457529030e-1, .3333567215434407e-1
                                                     };
    return _weights [QuadPoint];
  }

protected:
  Vec<1, _RealType> _points [ numQuadPoints ];
};

template <typename _RealType>
class GaussQuadrature<_RealType, qc::QC_1D, 21> {
public:
  enum { numQuadPoints = 11 };

  GaussQuadrature( ) {
    _points[0][0] = 0.1088567092697150e-1;
    _points[1][0] = 0.5646870011595235e-1;
    _points[2][0] = 0.1349239972129753;
    _points[3][0] = 0.2404519353965941;
    _points[4][0] = 0.3652284220238275;
    _points[5][0] = 0.5;
    _points[6][0] = 0.6347715779761725;
    _points[7][0] = 0.7595480646034059;
    _points[8][0] = 0.8650760027870247;
    _points[9][0] = 0.9435312998840476;
    _points[10][0] = 0.9891143290730285;
  }


  inline const Vec<1, _RealType> &getRefCoord ( int QuadPoint ) const {
    return _points [QuadPoint];
  }

  inline _RealType getWeight ( int QuadPoint ) const {
    static const _RealType _weights [numQuadPoints] = {
                                                       .2783428355808683e-1, .6279018473245231e-1, .9314510546386713e-1, .1165968822959952,
                                                       .1314022722551233, .1364625433889503, .1314022722551233, .1165968822959952,
                                                       .9314510546386713e-1, .6279018473245231e-1, .2783428355808683e-1
                                                     };
    return _weights [QuadPoint];
  }

protected:
  Vec<1, _RealType> _points [ numQuadPoints ];
};

template <typename _RealType>
class GaussQuadrature<_RealType, qc::QC_1D, 23> {
public:
  enum { numQuadPoints = 12 };

  GaussQuadrature( ) {
    _points[0][0] = 0.9219682876640375e-2;
    _points[1][0] = 0.4794137181476257e-1;
    _points[2][0] = 0.1150486629028477;
    _points[3][0] = 0.2063410228566913;
    _points[4][0] = 0.3160842505009099;
    _points[5][0] = 0.4373832957442655;
    _points[6][0] = 0.5626167042557345;
    _points[7][0] = 0.6839157494990901;
    _points[8][0] = 0.7936589771433087;
    _points[9][0] = 0.8849513370971523;
    _points[10][0] = 0.9520586281852374;
    _points[11][0] = 0.9907803171233596;
  }

  inline const Vec<1, _RealType> &getRefCoord ( int QuadPoint ) const {
    return _points [QuadPoint];
  }

  inline _RealType getWeight ( int QuadPoint ) const {
    static const _RealType _weights [numQuadPoints] = {
                                                       .2358766819325591e-1, .5346966299765922e-1, .8003916427167311e-1, .1015837133615330,
                                                       .1167462682691774, .1245735229067014, .1245735229067014, .1167462682691774,
                                                       .1015837133615330, .8003916427167311e-1, .5346966299765922e-1, .2358766819325591e-1
                                                     };
    return _weights [QuadPoint];
  }

protected:
  Vec<1, _RealType> _points [ numQuadPoints ];
};

template <typename _RealType, int _Order>
class GaussQuadrature1D {
public:
  typedef _RealType RealType;
  static const qc::Dimension Dim = qc::QC_1D;
  static const int Order = _Order;

  enum { numQuadPoints =
           GaussQuadrature<_RealType, qc::QC_1D, _Order>::numQuadPoints *
         GaussQuadrature<_RealType, qc::QC_1D, _Order>::numQuadPoints };

  GaussQuadrature1D( ) {

    GaussQuadrature<_RealType, qc::QC_1D, _Order> quad1D;

    for ( int k = 0, i = 0; i < GaussQuadrature<_RealType, qc::QC_1D, _Order>::numQuadPoints; i++ ) {
      for ( int j = 0; j < GaussQuadrature<_RealType, qc::QC_1D, _Order>::numQuadPoints; j++, k++ ) {
        _points[k][0] = quad1D.getRefCoord ( i )[0];
        _weights[k] = quad1D.getWeight ( i ) * quad1D.getWeight ( j );
      }
    }
  }

  inline const Vec<1, _RealType> &getRefCoord ( int QuadPoint ) const {
    return _points[ QuadPoint ];
  }

  inline _RealType getWeight ( int QuadPoint ) const {
    return _weights[ QuadPoint ];
  }

protected:
  Vec<1, _RealType> _points [ numQuadPoints ];
  _RealType         _weights[ numQuadPoints ];
};

template <typename _RealType, int _Order>
class GaussQuadrature<_RealType, qc::QC_2D, _Order> : public Quadrature2D<_RealType, GaussQuadrature<_RealType, qc::QC_1D, _Order> > {
public:
  static const int Order = _Order;
  GaussQuadrature( ) : Quadrature2D<_RealType, GaussQuadrature<_RealType, qc::QC_1D, _Order> >() { }
};

template <typename _RealType, int _Order>
class GaussQuadrature<_RealType, qc::QC_3D, _Order> : public Quadrature3D<_RealType, GaussQuadrature<_RealType, qc::QC_1D, _Order> > {
public:
  static const int Order = _Order;
  GaussQuadrature( ) : Quadrature3D<_RealType, GaussQuadrature<_RealType, qc::QC_1D, _Order> >() { }
};


/*!
* \author Horn
*/
template <typename RealType, qc::Dimension Dim, int order >
class TriangleIntegration {};


//! integration over a triangle: area * value at the center of gravity

template <typename RealType>
class TriangleIntegration <RealType, qc::QC_2D, 1> {

public:
  enum { numQuadPoints = 2 };

  TriangleIntegration ( ) {
     _points[0][0] = 1./3.;
     _points[0][1] = 1./3.;
     _points[1][0] = 2./3.;
     _points[1][1] = 2./3.;

     _weights[0] = 0.5;
     _weights[1] = 0.5;
  }

  inline const Vec2<RealType>& getRefCoord ( int QuadPoint ) const {
    return _points[ QuadPoint ];
  }

  inline RealType getWeight ( int QuadPoint ) const {
    return _weights[ QuadPoint ];
  }

protected:
  Vec2<RealType> _points [ numQuadPoints ];
  RealType       _weights[ numQuadPoints ];
};


//! integration over a triangle: area * (sum of weighted values at side middle points)

template <typename RealType>
class TriangleIntegration <RealType, qc::QC_2D, 2> {

public:
  enum { numQuadPoints = 5 };

  TriangleIntegration ( ) {
     _points[0][0] = 0.5;
     _points[0][1] = 0.;
     _points[1][0] = 0.5;
     _points[1][1] = 0.5;
     _points[2][0] = 0.;
     _points[2][1] = 0.5;
     _points[3][0] = 1.;
     _points[3][1] = 0.5;
     _points[4][0] = 0.5;
     _points[4][1] = 1.;

     _weights[0] = 1./6.;
     _weights[1] = 1./3.;
     _weights[2] = 1./6.;
     _weights[3] = 1./6.;
     _weights[4] = 1./6.;
  }

  inline const Vec2<RealType>& getRefCoord ( int QuadPoint ) const {
    return _points[ QuadPoint ];
  }

  inline RealType getWeight ( int QuadPoint ) const {
    return _weights[ QuadPoint ];
  }

protected:
  Vec2<RealType> _points [ numQuadPoints ];
  RealType       _weights[ numQuadPoints ];
};


//! integration over a tetrahedron: test
template <typename RealType>
class TriangleIntegration <RealType, qc::QC_3D, 0> {

public:
  enum { numQuadPoints = 1 };

  TriangleIntegration ( ) {
     _points[0][0] = 0.5;
     _points[0][1] = 0.5;
     _points[0][2] = 0.5;

     _weights[0] = 1.;
  }

  inline const Vec3<RealType>& getRefCoord ( int QuadPoint ) const {
    return _points[ QuadPoint ];
  }

  inline RealType getWeight ( int QuadPoint ) const {
    return _weights[ QuadPoint ];
  }

protected:
  Vec3<RealType> _points [ numQuadPoints ];
  RealType       _weights[ numQuadPoints ];
};

//! center of mass integration over the six tetrahedra which are decomposition of a unit cube
template <typename RealType>
class TriangleIntegration <RealType, qc::QC_3D, 1> {

public:
  enum { numQuadPoints = 6 };

  TriangleIntegration ( ) {
     _points[0][0] = 0.25;
     _points[0][1] = 0.25;
     _points[0][2] = 0.25;

     _points[1][0] = 0.5;
     _points[1][1] = 0.5;
     _points[1][2] = 0.25;

     _points[2][0] = 0.25;
     _points[2][1] = 0.75;
     _points[2][2] = 0.5;

     _points[3][0] = 0.75;
     _points[3][1] = 0.25;
     _points[3][2] = 0.5;

     _points[4][0] = 0.5;
     _points[4][1] = 0.5;
     _points[4][2] = 0.75;

     _points[5][0] = 0.75;
     _points[5][1] = 0.75;
     _points[5][2] = 0.75;

     _weights[0] = 1./6.;
     _weights[1] = 1./6.;
     _weights[2] = 1./6.;
     _weights[3] = 1./6.;
     _weights[4] = 1./6.;
     _weights[5] = 1./6.;
  }

  inline const Vec3<RealType>& getRefCoord ( int QuadPoint ) const {
    return _points[ QuadPoint ];
  }

  inline RealType getWeight ( int QuadPoint ) const {
    return _weights[ QuadPoint ];
  }

protected:
  Vec3<RealType> _points [ numQuadPoints ];
  RealType       _weights[ numQuadPoints ];
};




template <typename RealType, qc::Dimension Dim, class QuadRuleType>
class BaseFunctionSetLinear : public BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 4, QuadRuleType, BaseFunctionSetLinear<RealType, qc::QC_2D, QuadRuleType> > {
};

//! provides linear base functions for a regular triangulation (cross type: bisection of the rectangle from left bottom to right top)
/*!
* \author Horn
*/

template <typename RealType, class QuadRuleType>
class BaseFunctionSetLinear<RealType, qc::QC_2D, QuadRuleType> :
      public BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 4, QuadRuleType,
                                      BaseFunctionSetLinear<RealType, qc::QC_2D, QuadRuleType> > {

  static RealType _dx_b1   ( const Vec2<RealType> &RefCoord ) { return inBottomTriangle(RefCoord) ?  -1. : 0.; }
  static RealType _dx_b2   ( const Vec2<RealType> &RefCoord ) { return inBottomTriangle(RefCoord) ?   1. : 0.; }
  static RealType _dx_b3   ( const Vec2<RealType> &RefCoord ) { return inBottomTriangle(RefCoord) ?   0. : -1.; }
  static RealType _dx_b4   ( const Vec2<RealType> &RefCoord ) { return inBottomTriangle(RefCoord) ?   0. : 1.;}

  static RealType _dy_b1   ( const Vec2<RealType> &RefCoord ) { return inBottomTriangle(RefCoord) ?  -1. : 0.; }
  static RealType _dy_b2   ( const Vec2<RealType> &RefCoord ) { return inBottomTriangle(RefCoord) ?   0. : -1.;}
  static RealType _dy_b3   ( const Vec2<RealType> &RefCoord ) { return inBottomTriangle(RefCoord) ?   1. : 0.; }
  static RealType _dy_b4   ( const Vec2<RealType> &RefCoord ) { return inBottomTriangle(RefCoord) ?   0. : 1.; }


  static RealType _b1   ( const Vec2<RealType> &RefCoord ) { return inBottomTriangle(RefCoord) ?   1. - RefCoord[0] - RefCoord[1] : 0.; }
  static RealType _b2   ( const Vec2<RealType> &RefCoord ) { return inBottomTriangle(RefCoord) ?   RefCoord[0] : 1. - RefCoord[1]; }
  static RealType _b3   ( const Vec2<RealType> &RefCoord ) { return inBottomTriangle(RefCoord) ?   RefCoord[1] : 1. - RefCoord[0]; }
  static RealType _b4   ( const Vec2<RealType> &RefCoord ) { return inBottomTriangle(RefCoord) ?   0.: RefCoord[0] + RefCoord[1] - 1; }

  static bool inBottomTriangle ( const Vec2<RealType> &RefCoord ) {
    return ((RefCoord[0] +  RefCoord[1]) < 1.);
  }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const Vec2<RealType> &RefCoord );
  BASIS_FUNC_TYPE _deriv_basis[2][4];
  BASIS_FUNC_TYPE _basis[4];
  const RealType _h;

public:
  BaseFunctionSetLinear ( const RealType H ) : _h ( H ) {
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;

    _deriv_basis[0][0] = _dx_b1;
    _deriv_basis[0][1] = _dx_b2;
    _deriv_basis[0][2] = _dx_b3;
    _deriv_basis[0][3] = _dx_b4;

    _deriv_basis[1][0] = _dy_b1;
    _deriv_basis[1][1] = _dy_b2;
    _deriv_basis[1][2] = _dy_b3;
    _deriv_basis[1][3] = _dy_b4;

    this->initializeQuadCache( );
  }

  enum { numBaseFuncs = 4 };

  void evaluateGradient ( int BaseFuncNum, const Vec2<RealType> &RefCoord, Vec2<RealType> &Gradient ) const {
    Gradient[0] = _deriv_basis[0][BaseFuncNum] ( RefCoord ) / _h;
    Gradient[1] = _deriv_basis[1][BaseFuncNum] ( RefCoord ) / _h;
  }

  inline const Vec2<RealType>& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 4, QuadRuleType, BaseFunctionSetLinear<RealType, qc::QC_2D, QuadRuleType> >::evaluateGradient ( BaseFuncNum, QuadPoint );
  }

  RealType evaluate ( int BaseFuncNum, const Vec2<RealType> &RefCoord ) const {
    return _basis[BaseFuncNum] ( RefCoord );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec2<RealType>, Vec2<RealType>, 4, QuadRuleType,
                                    BaseFunctionSetLinear<RealType, qc::QC_2D, QuadRuleType> >::evaluate ( BaseFuncNum, QuadPoint );
  }
protected:

};

//! The basefunctionset for linear elements in 3d
template <typename RealType, class QuadRuleType>
class BaseFunctionSetLinear<RealType, qc::QC_3D, QuadRuleType> :
      public BaseFunctionSetInterface<RealType, Vec3<RealType>, Vec3<RealType>, 8, QuadRuleType,
                                      BaseFunctionSetLinear<RealType, qc::QC_3D, QuadRuleType> > {

  static RealType _dx_b1   ( const Vec3<RealType> &RefCoord ) {
    return inWhichSimplex(RefCoord) == 1 ? -1. : 0. ;
  }

  static RealType _dx_b2   ( const Vec3<RealType> &RefCoord ) {
    return inWhichSimplex(RefCoord) == 1 ? 1. : 0. ;
  }

  static RealType _dx_b3   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 2 :
        return -1.;
        break;
      case 3 :
        return -1.;
        break;
      default:
        return 0.;
      }
  }

  static RealType _dx_b4   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 2 :
        return 1.;
        break;
      case 3 :
        return 1.;
        break;
      default:
        return 0.;
      }
  }

  static RealType _dx_b5   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 4 :
        return -1.;
        break;
      case 5 :
        return -1.;
        break;
      default:
        return 0.;
      }
  }
  static RealType _dx_b6   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 4 :
        return 1.;
        break;
      case 5 :
        return 1.;
        break;
      default:
        return 0.;
      }
  }
  static RealType _dx_b7   ( const Vec3<RealType> &RefCoord ) {
    return inWhichSimplex(RefCoord) == 6 ? -1. : 0. ;
  }

  static RealType _dx_b8   ( const Vec3<RealType> &RefCoord ) {
    return inWhichSimplex(RefCoord) == 6 ? 1. : 0. ;
  }

  static RealType _dy_b1   ( const Vec3<RealType> &RefCoord ) {
    return inWhichSimplex(RefCoord) == 1 ? -1. : 0. ;
  }

  static RealType _dy_b2   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 2 :
        return -1.;
        break;
      case 4 :
        return -1.;
        break;
      default:
        return 0.;
    }
  }
  static RealType _dy_b3   ( const Vec3<RealType> &RefCoord ) {
    return inWhichSimplex(RefCoord) == 1 ? 1. : 0. ;
  }
  static RealType _dy_b4   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 2 :
        return 1.;
        break;
      case 4 :
        return 1.;
        break;
      default:
        return 0.;
      }
  }

  static RealType _dy_b5   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 3 :
        return -1.;
        break;
      case 5 :
        return -1.;
        break;
      default:
        return 0.;
      }
  }

  static RealType _dy_b6   ( const Vec3<RealType> &RefCoord ) {
    return inWhichSimplex(RefCoord) == 6 ? -1. : 0. ;
  }

  static RealType _dy_b7   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 3 :
        return 1.;
        break;
      case 5 :
        return 1.;
        break;
      default:
        return 0.;
      }
  }

  static RealType _dy_b8   ( const Vec3<RealType> &RefCoord ) {
    return inWhichSimplex(RefCoord) == 6 ? 1. : 0. ;
  }

  static RealType _dz_b1   ( const Vec3<RealType> &RefCoord ) {
    return inWhichSimplex(RefCoord) == 1 ? -1. : 0. ;
  }

  static RealType _dz_b2   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 2 :
        return -1.;
        break;
      case 4 :
        return -1.;
        break;
      default:
        return 0.;
    }
  }
  static RealType _dz_b3   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 2 :
        return -1.;
        break;
      case 3 :
        return -1.;
        break;
      default:
        return 0.;
      }
  }
  static RealType _dz_b4   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 2 :
        return 1.;
        break;
      case 5 :
        return -1.;
        break;
      case 6 :
        return -1.;
        break;
      default:
        return 0.;
      }
  }
  static RealType _dz_b5   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 1 :
        return 1.;
        break;
      case 2 :
        return 1.;
        break;
      case 5 :
        return -1.;
        break;
      default:
        return 0.;
      }
  }
  static RealType _dz_b6   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 4 :
        return 1.;
        break;
      case 5 :
        return 1.;
        break;
      default:
        return 0.;
      }
  }

  static RealType _dz_b7   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 3 :
        return 1.;
        break;
      case 5 :
        return 1.;
        break;
      default:
        return 0.;
      }
  }

  static RealType _dz_b8   ( const Vec3<RealType> &RefCoord ) {
    return inWhichSimplex(RefCoord) == 6 ? 1. : 0. ;
  }

 /* base functions */
  static RealType _b1   ( const Vec3<RealType> &RefCoord ) {
    return inWhichSimplex(RefCoord) == 1 ? 1. - RefCoord[0] - RefCoord[1] - RefCoord[2] : 0. ;
  }

  static RealType _b2   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 1 :
        return RefCoord[0];
        break;
      case 2 :
        return 1. - RefCoord[1] - RefCoord[2];
        break;
      case 4 :
        return 1. - RefCoord[1] - RefCoord[2];
        break;
      default:
        return 0.;
    }
  }
  static RealType _b3   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 1 :
        return RefCoord[1];
        break;
      case 2 :
        return 1. - RefCoord[0] - RefCoord[2];
        break;
      case 3 :
        return 1. - RefCoord[0] - RefCoord[2];
        break;
      default:
        return 0.;
      }
  }
  static RealType _b4   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 2 :
        return RefCoord[0] + RefCoord[1] + RefCoord[2] - 1.;
        break;
      case 3 :
        return RefCoord[0];
        break;
      case 4 :
        return RefCoord[1];
        break;
      case 5 :
        return 1. - RefCoord[2];
        break;
      case 6 :
        return 1. - RefCoord[2];
        break;
      default:
        return 0.;
      }
  }
  static RealType _b5   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 1 :
        return RefCoord[2];
        break;
      case 2 :
        return RefCoord[2];
        break;
      case 3 :
        return 1. - RefCoord[1];
        break;
      case 4 :
        return 1. - RefCoord[0];
        break;
      case 5 :
        return 2. - RefCoord[0] - RefCoord[1] - RefCoord[2];
        break;
      default:
        return 0.;
      }
  }
  static RealType _b6   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 4 :
        return RefCoord[0] + RefCoord[2] - 1.;
        break;
      case 5 :
        return RefCoord[0] + RefCoord[2] - 1.;
        break;
      case 6 :
        return 1. - RefCoord[1];
        break;
      default:
        return 0.;
      }
  }
  static RealType _b7   ( const Vec3<RealType> &RefCoord ) {
    int simplexNum = inWhichSimplex(RefCoord);
    switch(simplexNum){
      case 3 :
        return RefCoord[1] + RefCoord[2] - 1.;
        break;
      case 5 :
        return RefCoord[1] + RefCoord[2] - 1.;
        break;
      case 6 :
        return 1. - RefCoord[0];
        break;
      default:
        return 0.;
      }
   }
  static RealType _b8   ( const Vec3<RealType> &RefCoord ) {
    return inWhichSimplex(RefCoord) == 6 ? RefCoord[0] + RefCoord[1] + RefCoord[2] - 2. : 0. ;
  }


  static int inWhichSimplex ( const Vec3<RealType> &RefCoord ) {
    RealType barCoord[3];
    int simplexNum = -1;
    barCoord[0] = RefCoord[0];
    barCoord[1] = RefCoord[1];
    barCoord[2] = RefCoord[2];
    if (barCoord[0] >= 0. && barCoord[1] >= 0. && barCoord[2] >= 0. && barCoord[0] + barCoord[1] + barCoord[2] <= 1.){
      return 1;
      }

    barCoord[0] = -RefCoord[0] + 1. - RefCoord[2];
    barCoord[1] = RefCoord[0] - 1. + RefCoord[1] + RefCoord[2];
    barCoord[2] = RefCoord[2];
    if (barCoord[0] >= 0. && barCoord[1] >= 0. && barCoord[2] >= 0. && barCoord[0] + barCoord[1] + barCoord[2] <= 1.){
      return 2;
  }
    barCoord[0] = RefCoord[0];
    barCoord[1] = -RefCoord[1] + 1.;
    barCoord[2] = RefCoord[1] - 1. + RefCoord[2];
    if (barCoord[0] >= 0. && barCoord[1] >= 0. && barCoord[2] >= 0. && barCoord[0] + barCoord[1] + barCoord[2] <= 1.){
      return 3;
    }
    barCoord[0] = RefCoord[1];
    barCoord[1] = -RefCoord[0] + 1.;
    barCoord[2] = RefCoord[0] - 1. + RefCoord[2];
    if (barCoord[0] >= 0. && barCoord[1] >= 0. && barCoord[2] >= 0. && barCoord[0] + barCoord[1] + barCoord[2] <= 1.){
      return 4;
    }
    barCoord[0] = -RefCoord[0] - RefCoord[1] + 2. - RefCoord[2];
    barCoord[1] = RefCoord[0] - 1. + RefCoord[2];
    barCoord[2] = RefCoord[1] - 1. + RefCoord[2];
    if (barCoord[0] >= 0. && barCoord[1] >= 0. && barCoord[2] >= 0. && barCoord[0] + barCoord[1] + barCoord[2] <= 1.){
      return 5;
    }
    barCoord[0] = -RefCoord[1] + 1.;
    barCoord[1] = -RefCoord[0] + 1.;
    barCoord[2] = RefCoord[0] + RefCoord[1] - 2. + RefCoord[2];
    if (barCoord[0] >= 0. && barCoord[1] >= 0. && barCoord[2] >= 0. && barCoord[0] + barCoord[1] + barCoord[2] <= 1.){
      return 6;
    }

    if(true)
      throw aol::Exception( "The RefCoord is outside of all given tetrahedra!", __FILE__, __LINE__);
    return simplexNum;
  }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const Vec3<RealType> &RefCoord );
  BASIS_FUNC_TYPE _deriv_basis[3][8];
  BASIS_FUNC_TYPE _basis[8];
  const RealType _h;

public:
  BaseFunctionSetLinear ( const RealType H ) : _h ( H ) {
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;
    _basis[4] = _b5;
    _basis[5] = _b6;
    _basis[6] = _b7;
    _basis[7] = _b8;

    _deriv_basis[0][0] = _dx_b1;
    _deriv_basis[0][1] = _dx_b2;
    _deriv_basis[0][2] = _dx_b3;
    _deriv_basis[0][3] = _dx_b4;
    _deriv_basis[0][4] = _dx_b5;
    _deriv_basis[0][5] = _dx_b6;
    _deriv_basis[0][6] = _dx_b7;
    _deriv_basis[0][7] = _dx_b8;

    _deriv_basis[1][0] = _dy_b1;
    _deriv_basis[1][1] = _dy_b2;
    _deriv_basis[1][2] = _dy_b3;
    _deriv_basis[1][3] = _dy_b4;
    _deriv_basis[1][4] = _dy_b5;
    _deriv_basis[1][5] = _dy_b6;
    _deriv_basis[1][6] = _dy_b7;
    _deriv_basis[1][7] = _dy_b8;

    _deriv_basis[2][0] = _dz_b1;
    _deriv_basis[2][1] = _dz_b2;
    _deriv_basis[2][2] = _dz_b3;
    _deriv_basis[2][3] = _dz_b4;
    _deriv_basis[2][4] = _dz_b5;
    _deriv_basis[2][5] = _dz_b6;
    _deriv_basis[2][6] = _dz_b7;
    _deriv_basis[2][7] = _dz_b8;

    this->initializeQuadCache( );
  }


  enum { numBaseFuncs = 8 };

  void evaluateGradient ( int BaseFuncNum, const Vec3<RealType> &RefCoord, Vec3<RealType> &Gradient ) const {
    Gradient[0] = _deriv_basis[0][BaseFuncNum] ( RefCoord ) / _h;
    Gradient[1] = _deriv_basis[1][BaseFuncNum] ( RefCoord ) / _h;
    Gradient[2] = _deriv_basis[2][BaseFuncNum] ( RefCoord ) / _h;
  }

  inline const Vec3<RealType>& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec3<RealType>, Vec3<RealType>, 8, QuadRuleType, BaseFunctionSetLinear<RealType, qc::QC_3D, QuadRuleType> >::evaluateGradient ( BaseFuncNum, QuadPoint );
  }

  RealType evaluate ( int BaseFuncNum, const Vec3<RealType> &RefCoord ) const {
    return _basis[BaseFuncNum] ( RefCoord );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return BaseFunctionSetInterface<RealType, Vec3<RealType>, Vec3<RealType>, 8, QuadRuleType,
                                    BaseFunctionSetLinear<RealType, qc::QC_3D, QuadRuleType> >::evaluate ( BaseFuncNum, QuadPoint );
  }
protected:


};


}

#endif
