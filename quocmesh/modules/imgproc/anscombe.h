#ifndef __ANSCOMBE_H
#define __ANSCOMBE_H

#include <aol.h>
#include <op.h>


namespace im {


enum ANSCOMBE_INVERSE_TYPE {
  ANSCOMBE_INV_ASYMPTOTICALLY_UNBIASED,
  ANSCOMBE_INV_EXACT_UNBIASED
};


template<typename RealType>
class VSTOp : public aol::Op<aol::Vector<RealType> > {
public:
  VSTOp ( ) { }
  
  virtual RealType getMinTransformed ( ) const = 0;
};


template<typename RealType>
class AnscombeOp : public VSTOp<RealType> {
public:
  RealType getMinTransformed ( ) const {
    return 2 * sqrt ( 0.0 + 3.0 / 8.0 );
  }
};


template<typename _RealType>
class AnscombeForward : public AnscombeOp<_RealType> {
  typedef _RealType RealType;
  
  const static int anscombeLastIdx;
  const static double anscombeEfz[];
  const static double anscombeEz[];
public:
  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    const RealType a = 3.0 / 8.0;
    for ( int k = 0; k < Arg.size ( ) ; ++k )
      Dest[k] = 2 * sqrt ( Arg[k] + a );
  }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};

template<typename _RealType>
class AnscombeInverse : public AnscombeOp<_RealType> {
  typedef _RealType RealType;
public:
  const static int anscombeLastIdx;
  const static double anscombeEfz[];
  const static double anscombeEz[];
public:
  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    const RealType a = 3.0 / 8.0, b = 1.0 / 8.0;
    for ( int k = 0; k < Arg.size ( ) ; ++k ) {
      if ( Arg[k] < 2 * sqrt ( a ) )
        Dest[k] = 0;
      else if ( Arg[k] > anscombeEfz[anscombeLastIdx] )
        Dest[k] = aol::Sqr<RealType> ( 0.5 * Arg[k] ) - b;
      else
        Dest[k] = interpInverse ( Arg[k] );
    }
  }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void applyInverseOfExactAlgebraicInverse ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) {
    const RealType b = 1.0 / 8.0;
    for ( int k = 0; k < Arg.size ( ) ; ++k ) {
      const RealType y = 2 * sqrt ( Arg[k] + b );
      if ( y  > anscombeEfz[anscombeLastIdx] )
        Dest[k] = y;
      else
        Dest[k] = interp ( Arg[k] );
    }
  }
  
private:
  RealType interpInverse ( const RealType Val ) const {
    // Return extrapolated value if D beyond range
    if ( Val > anscombeEfz[anscombeLastIdx] )
      return ( Val - anscombeEfz[anscombeLastIdx-1] ) * ( anscombeEz[anscombeLastIdx] - anscombeEz[anscombeLastIdx-1] )
           / ( anscombeEfz[anscombeLastIdx] - anscombeEfz[anscombeLastIdx-1] ) + anscombeEz[anscombeLastIdx-1];
    
    // Perform binary search for k s.t. D in [Efz[k], Efz[k+1])
    int a = 0, b = anscombeLastIdx, k = ( a + b ) / 2;
    while ( Val < anscombeEfz[k] || Val > anscombeEfz[k+1] ) {
      if ( Val < anscombeEfz[k] ) {
        b = k-1;
        k = ( b + a ) / 2;
      } else {
        a = k+1;
        k = ( b + a ) / 2;
      }
    }
    
    // Return interpolated value on [Efz[k], Efz[k+1]]
    return ( Val - anscombeEfz[k-1] ) * ( anscombeEz[k] - anscombeEz[k-1] ) / ( anscombeEfz[k] - anscombeEfz[k-1] ) + anscombeEz[k-1];
  }
  
  RealType interp ( const RealType Val ) const {
    // Return extrapolated value if D beyond range
    if ( Val > anscombeEz[anscombeLastIdx] )
      return ( Val - anscombeEz[anscombeLastIdx-1] ) * ( anscombeEfz[anscombeLastIdx] - anscombeEfz[anscombeLastIdx-1] )
      / ( anscombeEz[anscombeLastIdx] - anscombeEz[anscombeLastIdx-1] ) + anscombeEfz[anscombeLastIdx-1];
    
    // Perform binary search for k s.t. D in [Efz[k], Efz[k+1])
    int a = 0, b = anscombeLastIdx, k = ( a + b ) / 2;
    while ( Val < anscombeEz[k] || Val > anscombeEz[k+1] ) {
      if ( Val < anscombeEz[k] ) {
        b = k-1;
        k = ( b + a ) / 2;
      } else {
        a = k+1;
        k = ( b + a ) / 2;
      }
    }
    
    // Return interpolated value on [Efz[k], Efz[k+1]]
    return ( Val - anscombeEz[k-1] ) * ( anscombeEfz[k] - anscombeEfz[k-1] ) / ( anscombeEz[k] - anscombeEz[k-1] ) + anscombeEfz[k-1];
  }
};


template<typename RealType>
class GeneralizedAnscombeOp : public VSTOp<RealType> {
public:
  RealType getMinTransformed ( ) const {
    return 0.0;
  }
  
  static void clampParameters ( RealType &A, RealType &B, const RealType MinInput, const bool Quiet = false ) {
    const RealType minB = -A * ( MinInput + 3.0 / 8.0 * A );
    if ( B < minB ) {
      if ( !Quiet )
        std::cerr << "GeneralizedAnscombeOp: clamped b = " << B << " to min value b_min = " << minB << std::endl;
      B = minB;
    }
  }
};


template<typename _RealType>
class GeneralizedAnscombeForward : public GeneralizedAnscombeOp<_RealType> {
  typedef _RealType RealType;
protected:
  const RealType _alpha, _threshold, _twoAlphaInv, _c;
public:
  GeneralizedAnscombeForward ( const RealType Alpha, const RealType SigmaSqr, const RealType Mu = 0 )
    : _alpha ( Alpha ),
      _threshold ( -3.0 / 8.0 * Alpha - ( SigmaSqr - Alpha * Mu ) / Alpha ),
      _twoAlphaInv ( 2.0 / Alpha ), _c ( 3.0 / 8.0 * aol::Sqr<RealType> ( Alpha ) + SigmaSqr - Alpha * Mu ) { }
  
  inline RealType getTransformedValue ( const RealType Arg ) const {
    if ( Arg > _threshold ) return _twoAlphaInv * sqrt ( _alpha * Arg + _c );
    else return 0.0;
  }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    for ( int k = 0; k < Arg.size ( ) ; ++k )
      Dest[k] = getTransformedValue ( Arg[k] );
  }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};


template<typename _RealType>
class GeneralizedAnscombeInverse : public GeneralizedAnscombeOp<_RealType> {
  typedef _RealType RealType;
protected:
  const RealType _alpha, _c, _c1, _c2, _c3, _c4, _c5;
public:
  GeneralizedAnscombeInverse ( const RealType Alpha, const RealType SigmaSqr, const RealType Mu = 0 )
    : _alpha ( Alpha ), _c ( ( SigmaSqr - Alpha * Mu ) / Alpha ),
      _c1 ( 0.25 ), _c2 ( 0.25 * sqrt ( 1.5 ) ), _c3 ( -11.0 / 8.0 ), _c4 ( 5.0 / 8.0 * sqrt ( 1.5 ) ), _c5 ( -1.0 / 8.0 ) { }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    for ( int k = 0; k < Arg.size ( ) ; ++k )
      Dest[k] = _alpha * ( _c1 * aol::Sqr<RealType> ( Arg[k] ) + _c2 / Arg[k] + _c3 / aol::Sqr<RealType> ( Arg[k] ) + _c4 / ( aol::Sqr<RealType> ( Arg[k] ) * Arg[k] ) + _c5 ) - _c;
  }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};

  
} // namespace im
  
  
#endif
