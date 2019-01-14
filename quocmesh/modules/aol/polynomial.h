#ifndef __POLYNOMIAL_H
#define __POLYNOMIAL_H

#include <matrix.h>

namespace aol {

/* remark:
 * did not care yet about numerical problems for higher degrees.
 * use only for small degrees!
 * later: BSpline class.
 */

/***
 * A class for polynomials.
 * \author Droske
 */
template <class RealType>
class Polynomial {
public:
  explicit Polynomial ( int Degree ) : _deg ( Degree ) {
    _coeff.reallocate ( _deg + 1 );
  }

  ~Polynomial() {
  }

  RealType operator() ( const RealType v ) const {
    RealType r = ZOTrait<RealType>::zero;
    // for ( int i = 0; i < _deg + 1; i++ ) {
    //  r += _coeff->get( i ) * pow( v, i );
    //}
    r = _coeff.get ( _deg );
    for ( int i = _deg - 1; i >= 0; i-- ) {
      r = r * v + _coeff.get ( i );
    }
    return r;
  }

  void setCoeff ( const int I, const RealType V ) {
    _coeff [ I ] = V ;
  }

  RealType getCoeff ( const int I ) const {
    return _coeff.get ( I );
  }

  int getDegree() const {
    return _deg;
  }

  void makeFirstDerivative ( const Polynomial<RealType> &Other ) {
    if ( _deg != Other.getDegree() - 1 ) {
      throw ( Exception ( "wrong degrees. Polynomial::makeFirstDerivative" ) );
    }
    for ( int i = 0; i < _deg + 1; i++ ) {
      _coeff [ i ] = static_cast<RealType> ( ( i + 1 ) ) * Other.getCoeff ( i + 1 );
    }
  }

  RealType firstDerivative ( const RealType v ) const {
    RealType r = ZOTrait<RealType>::zero;
    for ( int i = 1; i < _deg + 1; i++ ) {
      r += _coeff.get ( i ) * static_cast<RealType> ( i ) * pow ( v, i - 1 );
    }
  }

  void L2Projection ( const Vector<RealType> &AbsPts, const Vector<RealType> &Values ) {
    RealType s;
    int n = AbsPts.size();
    if ( n != Values.size() ) {
      throw Exception ( "length of AbsPts and Values. ", "Polynomial::L2Projection" );
    }
    if ( n < _deg + 1 ) {
      throw Exception ( "underdetermined system.", "Polynomial::L2Projection" );
    }
    Vector<RealType> rhs ( _deg + 1 );
    rhs.setZero();
    FullMatrix<RealType> mat ( _deg + 1, _deg + 1 ), inv ( _deg + 1, _deg + 1 );
    for ( int i = 0; i < _deg + 1; i++ ) {
      for ( int j = 0; j < _deg + 1; j++ ) {
        s = ZOTrait<RealType>::zero;
        for ( int k = 0; k < n; k++ ) {
          s += pow ( AbsPts.get ( k ), j + i );
        }
        mat.set ( i, j, s );
      }
      s = ZOTrait<RealType>::zero;
      for ( int k = 0; k < n; k++ ) {
        s += pow ( AbsPts.get ( k ), i ) * Values.get ( k );
      }
      rhs[ i ] = s;
    }
    inv.makeInverse ( mat );
    inv.mult ( rhs, _coeff );
  }

  /***
   * Attention: Since there's (yet) no QR the method is numerically unstable for high degree.
   */
  void interpolate ( const Vector<RealType> &AbsPts, const Vector<RealType> &Values ) {
    if ( AbsPts.size() != _deg + 1 || Values.size() != _deg + 1 ) {
      throw Exception ( "AbsPts or Values has wrong size.. has to be degree + 1", "Polynomial::interpolate" );
    }
    FullMatrix<RealType> mat ( _deg + 1, _deg + 1 ), inv ( _deg + 1, _deg + 1 );
    for ( int i = 0; i < _deg + 1; i++ ) {
      for ( int j = 0; j < _deg + 1; j++ ) {
        mat.set ( i, j, pow ( AbsPts.get ( i ), j ) );
      }
    }
    inv.makeInverse ( mat );
    inv.mult ( Values, _coeff );
  }

  void interpolateWithBoundaryDerivatives ( const Vector<RealType> &AbsPts, const Vector<RealType> &Values, const RealType LeftDer, const RealType RightDer ) {
    if ( static_cast<int> ( AbsPts.size() ) != _deg - 1  || static_cast<int> ( Values.size() ) != _deg - 1 ) {
      throw Exception ( "AbsPts or Values has wrong size.. has to be degree - 1", "Polynomial::interpolateWithBoundaryDerivatives" );
    }
    FullMatrix<RealType> mat ( _deg + 1, _deg + 1 ), inv ( _deg + 1, _deg + 1 )/*, id ( _deg + 1, _deg + 1 )*/;
    Vector<RealType> rhs ( _deg + 1 );

    for ( int i = 0; i < _deg - 1; i++ ) {
      for ( int j = 0; j < _deg + 1; j++ ) {
        mat.set ( i, j, pow ( AbsPts.get ( i ), j ) );
      }
      rhs[ i ] = Values[ i ];
    }

    for ( int j = 1; j < _deg + 1; j++ ) {
      mat.set ( _deg - 1, j, ( static_cast<RealType> ( j ) ) * pow ( AbsPts.get ( 0 ), j - 1 ) );
      mat.set ( _deg, j, ( static_cast<RealType> ( j ) ) * pow ( AbsPts.get ( _deg - 2 ), j - 1 ) );
    }
    mat.set ( _deg - 1, 0, 0. );
    mat.set ( _deg, 0, 0. );

    rhs[ _deg - 1 ] = LeftDer;
    rhs[ _deg ] = RightDer;

    inv.makeInverse ( mat );
    inv.mult ( rhs, _coeff );
  }

  void determineZeros( RealType &X1, RealType &X2 ) const {
    if ( _deg == 2 ) {
      RealType discr = Sqr( _coeff[1] ) - 4 * _coeff[0]*_coeff[2];
      if ( discr < 0 ) {
        throw aol::Exception( "Polynomial<>::determineZeroes: there are no real roots!", __FILE__, __LINE__ );
      } else {
        discr = sqrt( discr );
        X1 = (-_coeff[1] - discr)/(2.f*_coeff[2]);
        X2 = (-_coeff[1] + discr)/(2.f*_coeff[2]);
      }
    } else if ( _deg == 1 ) {
      // X1 = X2 = ( -_coeff[1]/_coeff[0] );
      X1 = X2 = ( -_coeff[0]/_coeff[1]);
    } else {
      throw aol::Exception( "Polynomial<>::determineZeroes: determineZeroes works only for degree 1, 2 so far!", __FILE__, __LINE__ );
    }
  }

  void dump ( ostream &out ) const {
    for ( int i = 0; i < _deg + 1; i++ ) {
      out << _coeff.get ( i ) << "*x**" << i;
      if ( i < _deg && _coeff->get ( i + 1 ) >= 0. ) out << " + ";
    }
  }

protected:
  int _deg;
  Vector<RealType> _coeff;

private:
  // think about proper implementation if needed
  Polynomial();
  Polynomial ( const Polynomial<RealType> &other );
  Polynomial<RealType>& operator= ( const Polynomial<RealType> &other ); // mathematical assignment or assignment only for matching degree?
};

}

#endif
