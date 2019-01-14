#ifndef __ENOCONVECT_H
#define __ENOCONVECT_H

#include <scalarArray.h>

namespace qc {

template <typename RealType, qc::Dimension Dim>
class ENODerivative { };


//! Computes essentially non oscillating derivatives
//! \author droske
template <typename RealType>
class ENODerivative<RealType, qc::QC_2D> {
  const qc::ScalarArray<RealType, qc::QC_2D> &_data;
  const RealType hx, hy;
public:
  ENODerivative ( const qc::ScalarArray<RealType, qc::QC_2D> &data ) : _data ( data ),
      hx ( 1. / static_cast<RealType> ( data.getNumX() - 1 ) ),
      hy ( 1. / static_cast<RealType> ( data.getNumY() - 1 ) ) {}

  inline RealType firstDerivativeXENO ( int centerX, int X, int Y, int order ) {
    RealType d = firstDerX ( X, Y ) / hx;
    if ( order <= 1 ) { return d; }

    RealType d21 = secondDerX ( X, Y ), d22 = secondDerX ( X + 1, Y );
    RealType d2;
    int kstar;
    if ( aol::Abs ( d21 ) <= aol::Abs ( d22 ) ) {
      d2 = d21 / ( 2. * hx );
      kstar = X - 1;
    } else {
      d2 = d22 / ( 2. * hx );
      kstar = X;
    }
    d += d2 * ( 2. * static_cast<RealType> ( centerX - X ) - 1. );
    if ( order <= 2 ) { return d; }

    const RealType d31 = thirdDerX ( kstar, Y ), d32 = thirdDerX ( kstar + 1, Y );
    RealType d3;
    if ( aol::Abs ( d31 ) <= aol::Abs ( d32 ) ) {
      d3 = d31 / ( 6. * hx );
    } else {
      d3 = d32 / ( 6. * hx );
    }
    d += d3 * ( 3. * aol::Sqr ( static_cast<RealType> ( centerX - kstar ) ) - 6. * static_cast<RealType> ( centerX - kstar ) + 2. );
    return d;
  }

  inline RealType firstDerivativeYENO ( int centerY, int X, int Y, int order ) {
    RealType d = firstDerY ( X, Y ) / hy;
    if ( order <= 1 ) { return d; }

    RealType d21 = secondDerY ( X, Y ), d22 = secondDerY ( X, Y + 1 );
    RealType d2;
    int kstar;
    if ( aol::Abs ( d21 ) <= aol::Abs ( d22 ) ) {
      d2 = d21 / ( 2. * hy );
      kstar = Y - 1;
    } else {
      d2 = d22 / ( 2. * hy );
      kstar = Y;
    }
    d += d2 * ( 2. * static_cast<RealType> ( centerY - Y ) - 1. );
    if ( order <= 2 ) { return d; }

    const RealType d31 = thirdDerY ( X, kstar ), d32 = thirdDerY ( X, kstar + 1 );
    RealType d3;
    if ( aol::Abs ( d31 ) <= aol::Abs ( d32 ) ) {
      d3 = d31 / ( 6. * hy );
    } else {
      d3 = d32 / ( 6. * hy );
    }
    d += d3 * ( 3. * aol::Sqr ( static_cast<RealType> ( centerY - kstar ) ) - 6. * static_cast<RealType> ( centerY - kstar ) + 2. );
    return d;
  }

protected:
  RealType thirdDerX ( int X, int Y ) {
    return _data.getClip ( X + 2, Y ) - 3. * _data.getClip ( X + 1, Y ) + 3. * _data.getClip ( X, Y ) - _data.getClip ( X - 1, Y );
  }

  RealType thirdDerY ( int X, int Y ) {
    return _data.getClip ( X, Y + 2 ) - 3. * _data.getClip ( X, Y + 1 ) + 3. * _data.getClip ( X, Y ) - _data.getClip ( X, Y - 1 );
  }

  RealType secondDerX ( int X, int Y ) {
    return _data.getClip ( X + 1, Y ) - 2. * _data.getClip ( X, Y ) + _data.getClip ( X - 1, Y );
  }

  RealType secondDerY ( int X, int Y ) {
    return _data.getClip ( X, Y + 1 ) - 2. * _data.getClip ( X, Y ) + _data.getClip ( X, Y - 1 );
  }

  RealType firstDerX ( int X, int Y ) {
    return _data.getClip ( X + 1, Y ) - _data.getClip ( X, Y );
  }

  RealType firstDerY ( int X, int Y ) {
    return _data.getClip ( X, Y + 1 ) - _data.getClip ( X, Y );
  }


};

/****
 * @brief general interface for level set evolution, provides an higher-order ENO-scheme.
 */
template <typename RealType, qc::Dimension Dim, typename Imp>
class ENOConvectInt {
public:
  ENOConvectInt( ) :
      _data ( NULL ), _update ( NULL ) {
  }

  ~ENOConvectInt( ) {
    if ( _update ) delete _update;
  }

Imp& asImp( ) { return static_cast<Imp&> ( *this ); }
  const Imp& asImp( ) const { return static_cast<const Imp&> ( *this ); }

  void velocity ( const qc::ScalarArray<RealType, qc::QC_2D> &Data, int x, int y, aol::Vec2<RealType> &vec ) const {
    asImp().velocity ( Data, x, y, vec );
  }

  /** member function to specify the array on which the level set image is stored.
   *
   */
  void setData ( qc::ScalarArray<RealType, qc::QC_2D> *Data ) {
    _data = Data;

    hx = 1. / static_cast<RealType> ( _data->getNumX() - 1 );
    hy = 1. / static_cast<RealType> ( _data->getNumY() - 1 );

    delete _update;
    _update = new qc::ScalarArray<RealType, qc::QC_2D> ( Data->getNumX(), Data->getNumY() );
  }

  /** member function to compute a single time step with a third order ENO-scheme.
   * @param Tau the time step size
    * @param order
   */
  void timeStepENO ( RealType Tau, int order = 3 ) {


    ENODerivative<RealType, Dim> der ( *_data );

    const int numX = _data->getNumX();
    const int numY = _data->getNumY();

    for ( int X = 0; X < numX; X++ ) {
      for ( int Y = 0; Y < numY; Y++ ) {
        aol::Vec2<RealType> vel;
        velocity ( *_data, X, Y, vel );

        RealType dxr = der.firstDerivativeXENO ( X, X, Y, order );
        RealType dxl = der.firstDerivativeXENO ( X, X - 1, Y, order );
        RealType dyb = der.firstDerivativeYENO ( Y, X, Y - 1, order );
        RealType dyt = der.firstDerivativeYENO ( Y, X, Y, order );

        RealType update = 0.0;
        if ( vel[0] > 0. ) {
          update += dxl * vel[0];
        } else {
          update += dxr * vel[0];
        }

        //update -= 0.5 * (dxr - dxl) / 2.;

        if ( vel[1] > 0. ) {
          update += dyb * vel[1];
        } else {
          update += dyt * vel[1];
        }

        //update -= 0.5 * (dyt - dyb) / 2.;

        _update->set ( X, Y, -Tau * update ) ;
      }
    }
    *_data += *_update;
  }

protected:


  RealType hx, hy;

  qc::ScalarArray<RealType, qc::QC_2D> *_data;
  qc::ScalarArray<RealType, qc::QC_2D> *_update;
};



template <typename RealType>
class ENOConvectExample : public ENOConvectInt<RealType, qc::QC_2D, ENOConvectExample<RealType> > {
public:
  void velocity ( const qc::ScalarArray<RealType, qc::QC_2D> &/*Data*/, int x, int y, aol::Vec2<RealType> &vec ) const {
    vec[0] = 0.5;
    vec[1] = 0.4;

    RealType X = ( static_cast<RealType> ( x ) * this->hx ) - 0.5;
    RealType Y = ( static_cast<RealType> ( y ) * this->hy ) - 0.5;

    RealType s = sin ( sqrt ( X * X + Y * Y ) * 10 );

    vec[0] = Y;
    vec[1] = -X;

    vec[0] = X;
    vec[1] = Y;
    vec /= sqrt ( vec.normSqr() + 0.01 );
    vec *= s;

    //return asImp().velocity( Data, x, y, vec );
  }
};


template <typename RealType>
class ENOConvectWithVelocityField : public ENOConvectInt<RealType, qc::QC_2D, ENOConvectWithVelocityField<RealType> > {
protected:
  const aol::MultiVector<RealType> &_velocity;
public:
  ENOConvectWithVelocityField ( const aol::MultiVector<RealType> &Velocity ) :
      _velocity ( Velocity ) {}

  void velocity ( const qc::ScalarArray<RealType, qc::QC_2D> &Data, int x, int y, aol::Vec2<RealType> &vec ) const {
    vec[0] = _velocity[0][ y * Data.getNumX() + x ];
    vec[1] = _velocity[1][ y * Data.getNumX() + x ];
  }


};


} // end namespace qc

#endif
