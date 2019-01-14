#ifndef __LEVELSET_H
#define __LEVELSET_H

#include <array.h>
#include <scalarArray.h>

namespace qc {

/** Class for computing viscosity-solution of the problem
 *  \f[ \partial_t \phi + F\|\nabla \phi \| = 0 \f]
 *  \f[  \phi(0,\cdot) = \phi_0 \f]
 * @param Velocity a template argument which should supply the speed function F. Velocity MUST contain a member function RealType veolicity( Array3d<RealType> &data, int X, int Y, int Z )
 */
template <typename RealType, class Imp>
class LevelSet3dInt  {
public:
  LevelSet3dInt() :
      _data ( NULL ), _update ( NULL ) {}

  ~LevelSet3dInt() {
    if ( _update ) delete _update;
  }

  Imp& asImp() {
    return static_cast<Imp&> ( *this );
  }
  const Imp& asImp() const {
    return static_cast<const Imp&> ( *this );
  }

  RealType velocity ( const ScalarArray<RealType, qc::QC_3D> &Data, int x, int y, int z ) const {
    return asImp().velocity ( Data, x, y, z );
  }

  /** member function to specify the array on which the level set image is stored.
   *
   */
  void setData ( ScalarArray<RealType, qc::QC_3D> *Data ) {
    _data = Data;
    if ( _update ) delete _update;
    _update = new ScalarArray<RealType, qc::QC_3D> ( Data->getNumX(), Data->getNumY(), Data->getNumZ() );
  }

  /** member function to compute a single Engquist-Osher update of the level set.
   * @param Tau
   */
  void timeStepEO ( RealType Tau ) {

    const int numX = _data->getNumX();
    const int numY = _data->getNumY();
    const int numZ = _data->getNumZ();

    const RealType hxr = static_cast< RealType > ( numX - 1 );
    const RealType hyr = static_cast< RealType > ( numY - 1 );
    const RealType hzr = static_cast< RealType > ( numZ - 1 );

    int X, Y, Z;

    for ( X = 0; X < numX; X++ ) {
      for ( Y = 0; Y < numY; Y++ ) {
        for ( Z = 0; Z < numZ; Z++ ) {
          RealType v = velocity ( *_data, X, Y, Z );
          RealType nw, cen = _data->get ( X, Y, Z );
          if ( v > 0 ) {
            RealType update = 0.0;
            if ( X > 0 && ( cen > ( nw = _data->get ( X - 1, Y, Z ) ) ) ) {
              update += aol::Sqr ( ( nw - cen ) * hxr );
            }
            if ( X < numX - 1 && ( cen > ( nw = _data->get ( X + 1, Y, Z ) ) ) ) {
              update += aol::Sqr ( ( nw - cen ) * hxr );
            }

            if ( Y > 0 && ( cen > ( nw = _data->get ( X, Y - 1, Z ) ) ) ) {
              update += aol::Sqr ( ( nw - cen ) * hyr );
            }
            if ( Y < numY - 1 && ( cen > ( nw = _data->get ( X, Y + 1, Z ) ) ) ) {
              update += aol::Sqr ( ( nw - cen ) * hyr );
            }

            if ( Z > 0 && ( cen > ( nw = _data->get ( X, Y, Z - 1 ) ) ) ) {
              update += aol::Sqr ( ( nw - cen ) * hzr );
            }
            if ( Z < numZ - 1 && ( cen > ( nw = _data->get ( X, Y, Z + 1 ) ) ) ) {
              update += aol::Sqr ( ( nw - cen ) * hzr );
            }

            _update->set ( X, Y, Z, -Tau * v * sqrt ( update ) );
          }
          if ( v < 0 ) {
            RealType update = 0.0;
            if ( X > 0 && ( cen < ( nw = _data->get ( X - 1, Y, Z ) ) ) ) {
              update += aol::Sqr ( ( nw - cen ) * hxr );
            }
            if ( X < numX - 1 && ( cen < ( nw = _data->get ( X + 1, Y, Z ) ) ) ) {
              update += aol::Sqr ( ( nw - cen ) * hxr );
            }

            if ( Y > 0 && ( cen < ( nw = _data->get ( X, Y - 1, Z ) ) ) ) {
              update += aol::Sqr ( ( nw - cen ) * hyr );
            }
            if ( Y < numY - 1 && ( cen < ( nw = _data->get ( X, Y + 1, Z ) ) ) ) {
              update += aol::Sqr ( ( nw - cen ) * hyr );
            }

            if ( Z > 0 && ( cen < ( nw = _data->get ( X, Y, Z - 1 ) ) ) ) {
              update += aol::Sqr ( ( nw - cen ) * hzr );
            }
            if ( Z < numZ - 1 && ( cen < ( nw = _data->get ( X, Y, Z + 1 ) ) ) ) {
              update += aol::Sqr ( ( nw - cen ) * hzr );
            }

            _update->set ( X, Y, Z, -Tau * v * sqrt ( update ) );
          }
        }
      }
    }
    *_data += *_update;
  }

protected:
  ScalarArray<RealType, qc::QC_3D> *_data;
  ScalarArray<RealType, qc::QC_3D> *_update;
};


/****
 * @brief general interface for level set evolution, provides an higher-order ENO-scheme.
 */
template <typename RealType, class Imp>
class LevelSet2dInt {
public:
  LevelSet2dInt() :
      _data ( NULL ), _update ( NULL ) {}

  ~LevelSet2dInt() {
    if ( _update ) delete _update;
  }

  Imp& asImp() {
    return static_cast<Imp&> ( *this );
  }
  const Imp& asImp() const {
    return static_cast<const Imp&> ( *this );
  }

  RealType velocity ( const ScalarArray<RealType, qc::QC_2D> &Data, int x, int y ) const {
    return asImp().velocity ( Data, x, y );
  }

  /** member function to specify the array on which the level set image is stored.
   *
   */
  void setData ( ScalarArray<RealType, qc::QC_2D> *Data ) {
    _data = Data;

    hx = 1. / static_cast< RealType > ( _data->getNumX() - 1 );
    hy = 1. / static_cast< RealType > ( _data->getNumY() - 1 );

    // cerr << "hx = " << hx << "  hy = " << hy << endl;

    hx2 = aol::Sqr ( hx ), hx3 = aol::Cub ( hx );
    hy2 = aol::Sqr ( hy ), hy3 = aol::Cub ( hy );

    delete _update;
    _update = new ScalarArray<RealType, qc::QC_2D> ( Data->getNumX(), Data->getNumY() );
  }

  /** member function to compute a single Engquist-Osher-Upwinding update of the level set.
   * @param Tau the time step size
   * \author Droske
   */
  void timeStepEO ( RealType Tau ) {

    const int numX = _data->getNumX();
    const int numY = _data->getNumY();

    const RealType hxr = static_cast< RealType > ( numX - 1 );
    const RealType hyr = static_cast< RealType > ( numY - 1 );

    int X, Y;

    for ( X = 0; X < numX; X++ ) {
      for ( Y = 0; Y < numY; Y++ ) {
        RealType v = velocity ( *_data, X, Y );
        RealType nw, cen = _data->get ( X, Y );
        if ( v > 0 ) {
          RealType update = 0.0;
          if ( X > 0 && ( cen > ( nw = _data->get ( X - 1, Y ) ) ) ) {
            update += aol::Sqr ( ( nw - cen ) * hxr );
          }
          if ( X < numX - 1 && ( cen > ( nw = _data->get ( X + 1, Y ) ) ) ) {
            update += aol::Sqr ( ( nw - cen ) * hxr );
          }

          if ( Y > 0 && ( cen > ( nw = _data->get ( X, Y - 1 ) ) ) ) {
            update += aol::Sqr ( ( nw - cen ) * hyr );
          }
          if ( Y < numY - 1 && ( cen > ( nw = _data->get ( X, Y + 1 ) ) ) ) {
            update += aol::Sqr ( ( nw - cen ) * hyr );
          }

          _update->set ( X, Y, -Tau * v * sqrt ( update ) );
        }
        if ( v < 0 ) {
          RealType update = 0.0;
          if ( X > 0 && ( cen < ( nw = _data->get ( X - 1, Y ) ) ) ) {
            update += aol::Sqr ( ( nw - cen ) * hxr );
          }
          if ( X < numX - 1 && ( cen < ( nw = _data->get ( X + 1, Y ) ) ) ) {
            update += aol::Sqr ( ( nw - cen ) * hxr );
          }

          if ( Y > 0 && ( cen < ( nw = _data->get ( X, Y - 1 ) ) ) ) {
            update += aol::Sqr ( ( nw - cen ) * hyr );
          }
          if ( Y < numY - 1 && ( cen < ( nw = _data->get ( X, Y + 1 ) ) ) ) {
            update += aol::Sqr ( ( nw - cen ) * hyr );
          }

          _update->set ( X, Y, -Tau * v * sqrt ( update ) );
        }
      }
    }
    *_data += *_update;
  }

  /** member function to compute a single time step with a third order ENO-scheme.
  * @param Tau the time step size
   * @param order
  * \author Droske
  */
  void timeStepENO ( RealType Tau, int order = 3 ) {

    const int numX = _data->getNumX();
    const int numY = _data->getNumY();

    //const RealType hx = 1./(RealType)(numX - 1);
    //const RealType hy = 1./(RealType)(numY - 1);

    // const RealType hx2 = aol::Sqr( hx ), hx3 = aol::Cub( hx );
    // const RealType hy2 = aol::Sqr( hy ), hy3 = aol::Cub( hy );

    for ( int X = 0; X < numX; X++ ) {
      for ( int Y = 0; Y < numY; Y++ ) {
        RealType v = velocity ( *_data, X, Y );

        RealType dxr = firstDerivativeXENO ( X, X, Y, order );
        RealType dxl = firstDerivativeXENO ( X, X - 1, Y, order );
        RealType dyb = firstDerivativeYENO ( Y, X, Y - 1, order );
        RealType dyt = firstDerivativeYENO ( Y, X, Y, order );

        RealType update = 0.0;
        if ( v > 0. ) {
          update += aol::Max ( aol::Sqr ( aol::Max ( dxl, 0. ) ), aol::Sqr ( aol::Min ( dxr, 0. ) ) );
          update += aol::Max ( aol::Sqr ( aol::Max ( dyb, 0. ) ), aol::Sqr ( aol::Min ( dyt, 0. ) ) );
        } else if ( v < 0. ) {
          update += aol::Max ( aol::Sqr ( aol::Min ( dxl, 0. ) ), aol::Sqr ( aol::Max ( dxr, 0. ) ) );
          update += aol::Max ( aol::Sqr ( aol::Min ( dyb, 0. ) ), aol::Sqr ( aol::Max ( dyt, 0. ) ) );
        }
        _update->set ( X, Y, -Tau * v * sqrt ( update ) );
      }
    }
    // cerr << "min = " << _update->getMinValue() << " max = " << _update->getMaxValue() << endl;

    *_data += *_update;
  }

protected:

  inline RealType firstDerivativeXENO ( int centerX, int X, int Y, int order ) {
    RealType d = firstDerX ( X, Y ) / hx;
    if ( order <= 1 ) {
      return d;
    }

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
    d += d2 * ( 2. * static_cast< RealType > ( centerX - X ) - 1. );
    if ( order == 2 ) {
      return d;
    }

    RealType d31 = thirdDerX ( kstar, Y ), d32 = thirdDerX ( kstar + 1, Y );
    RealType d3;
    if ( aol::Abs ( d31 ) <= aol::Abs ( d32 ) ) {
      d3 = d31 / ( 6. * hx );
    } else {
      d3 = d32 / ( 6. * hx );
    }
    d += d3 * ( 3. * aol::Sqr ( static_cast< RealType > ( centerX - kstar ) ) - 6. * static_cast< RealType > ( centerX - kstar ) + 2. );
    return d;
  }

  inline RealType firstDerivativeYENO ( int centerY, int X, int Y, int order ) {
    RealType d = firstDerY ( X, Y ) / hy;
    if ( order <= 1 ) {
      return d;
    }

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
    d += d2 * ( 2. * static_cast< RealType > ( centerY - Y ) - 1. );
    if ( order == 2 ) {
      return d;
    }

    RealType d31 = thirdDerY ( X, kstar ), d32 = thirdDerX ( X, kstar + 1 );
    RealType d3;
    if ( aol::Abs ( d31 ) <= aol::Abs ( d32 ) ) {
      d3 = d31 / ( 6. * hy );
    } else {
      d3 = d32 / ( 6. * hy );
    }
    d += d3 * ( 3. * aol::Sqr ( static_cast< RealType > ( centerY - kstar ) ) - 6. * static_cast< RealType > ( centerY - kstar ) + 2. );
    return d;
  }

  RealType thirdDerX ( int X, int Y ) {
    return _data->getClip ( X + 2, Y ) - 3. * _data->getClip ( X + 1, Y ) + 3. * _data->getClip ( X, Y ) - _data->getClip ( X - 1, Y );
  }

  RealType thirdDerY ( int X, int Y ) {
    return _data->getClip ( X, Y + 2 ) - 3. * _data->getClip ( X, Y + 1 ) + 3. * _data->getClip ( X, Y ) - _data->getClip ( X, Y - 1 );
  }

  RealType secondDerX ( int X, int Y ) {
    return _data->getClip ( X + 1, Y ) - 2. * _data->getClip ( X, Y ) + _data->getClip ( X - 1, Y );
  }

  RealType secondDerY ( int X, int Y ) {
    return _data->getClip ( X, Y + 1 ) - 2. * _data->getClip ( X, Y ) + _data->getClip ( X, Y - 1 );
  }

  RealType firstDerX ( int X, int Y ) {
    return _data->getClip ( X + 1, Y ) - _data->getClip ( X, Y );
  }

  RealType firstDerY ( int X, int Y ) {
    return _data->getClip ( X, Y + 1 ) - _data->getClip ( X, Y );
  }

  RealType hx, hx2, hx3, hy, hy2, hy3;

  ScalarArray<RealType, qc::QC_2D> *_data;
  ScalarArray<RealType, qc::QC_2D> *_update;
};

/** Class for computing viscosity-solution of the problem
 *  \f[ \partial_t \phi + F\|\nabla \phi \| = 0 \f]
 *  \f[  \phi(0,\cdot) = \phi_0 \f]
 * @param Velocity a template argument which should supply the speed function F. Velocity MUST contain a member function RealType veolicity( ScalarArray<RealType, qc::QC_3D> &data, int X, int Y, int Z )
 */
template <typename RealType, class Velocity>
class LevelSet2d : public Velocity, public qc::LevelSet2dInt<RealType, qc::LevelSet2d<RealType, Velocity> > {
public:
  LevelSet2d() :
      LevelSet2dInt<RealType, qc::LevelSet2d<RealType, Velocity> > () {}

  RealType velocity ( const ScalarArray<RealType, qc::QC_2D> &Data, int x, int y ) const {
    return Velocity::velocity ( Data, x, y );
  }

protected:

};

template <typename RealType>
class LevelSet2dWithSpeedField : public qc::LevelSet2dInt<RealType, qc::LevelSet2dWithSpeedField<RealType> > {
public:
  LevelSet2dWithSpeedField ( const qc::ScalarArray<RealType, qc::QC_2D> &Velocity )
      : _velocity ( Velocity ) {}

  RealType velocity ( const ScalarArray<RealType, qc::QC_2D> &/*Data*/, int x, int y ) const {
    return _velocity.get ( x, y );
  }

  const ScalarArray<RealType, qc::QC_2D> &_velocity;
};


template <typename RealType>
class LevelSet3dWithSpeedField : public qc::LevelSet3dInt<RealType, qc::LevelSet3dWithSpeedField<RealType> > {
public:
  LevelSet3dWithSpeedField ( const qc::ScalarArray<RealType, qc::QC_3D> &velocity )
      : _velocity ( velocity ) {}

  RealType velocity ( const ScalarArray<RealType, qc::QC_3D> &/*Data*/, int x, int y, int z ) const {
    return _velocity.get ( x, y, z );
  }

  const ScalarArray<RealType, qc::QC_3D> &_velocity;
};





}

#endif
