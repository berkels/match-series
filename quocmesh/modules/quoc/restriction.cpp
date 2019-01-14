#include <restriction.h>

template <typename RealType, qc::RestrictType RestrType>
void qc::RestrictOp<RealType, RestrType>::assemble_std_mg_restrict_matrix ( aol::SparseMatrix<RealType> &Mat ) const {

  Mat.setZero();

  // const is complete nonsense here. unfortunately, apply and apply add are const, thus require constness here.
  // could allow more general matrix consisting of rows.

  const int
    n_fi = _fineGrid.getWidth(),
    n_co = _coarseGrid.getWidth();

  switch ( _fineGrid.getDimOfWorld() ) {
  case qc::QC_2D: {
    for ( int i = 0; i < n_fi; ++i ) {
      for ( int j = 0; j < n_fi; ++j ) {
        // integer divisions used on purpose!

        // points in fine grid are also points in coarse grid
        if ( ! ( i % 2 ) && ! ( j % 2 ) ) {
          Mat.add ( n_co * ( ( i / 2 ) )     + ( ( j / 2 )     ) ,   n_fi * i + j,  1.00 );
        }

        // interpolation / restriction in 2nd coordinate
        if ( ! ( i % 2 ) && ( j % 2 ) ) {
          Mat.add ( n_co * ( ( i / 2 ) )     + ( ( j / 2 )     ) ,   n_fi * i + j,  0.50 );
          Mat.add ( n_co * ( ( i / 2 ) )     + ( ( j / 2 ) + 1 ) ,   n_fi * i + j,  0.50 );
        }

        // interpolation / restriction in 1st coordinate
        if ( ( i % 2 ) && ! ( j % 2 ) ) {
          Mat.add ( n_co * ( ( i / 2 ) )     + ( ( j / 2 )     ) ,   n_fi * i + j,  0.50 );
          Mat.add ( n_co * ( ( i / 2 ) + 1 ) + ( ( j / 2 )     ) ,   n_fi * i + j,  0.50 );
        }

        // interpolation / restriction in 1st and 2nd coordinate
        if ( ( i % 2 ) && ( j % 2 ) ) {
          Mat.add ( n_co * ( ( i / 2 ) )     + ( ( j / 2 )     ) ,   n_fi * i + j,  0.25 );
          Mat.add ( n_co * ( ( i / 2 ) )     + ( ( j / 2 ) + 1 ) ,   n_fi * i + j,  0.25 );
          Mat.add ( n_co * ( ( i / 2 ) + 1 ) + ( ( j / 2 )     ) ,   n_fi * i + j,  0.25 );
          Mat.add ( n_co * ( ( i / 2 ) + 1 ) + ( ( j / 2 ) + 1 ) ,   n_fi * i + j,  0.25 );
        }
      }
    }
  } break;

  case qc::QC_3D: {
    const int
      n_fi_2 = n_fi * n_fi,
      n_co_2 = n_co * n_co;

    for ( int i = 0; i < n_fi; ++i ) {
      for ( int j = 0; j < n_fi; ++j ) {
        for ( int k = 0; k < n_fi; ++k ) {
          // integer divisions used on purpose!

          // points in fine grid are also points in coarse grid
          if ( ! ( i % 2 ) && ! ( j % 2 ) && ! ( k % 2 ) ) {
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      1.00 );
          }

          // interpolation in 3rd coordinate
          if ( ! ( i % 2 ) && ! ( j % 2 ) && ( k % 2 ) ) {
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.50 );
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) + 1 ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.50 );
          }

          // interpolation in 2nd coordinate
          if ( ! ( i % 2 ) && ( j % 2 ) && ! ( k % 2 ) ) {
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.50 );
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.50 );
          }

          // interpolation in 1st coordinate
          if ( ( i % 2 ) && ! ( j % 2 ) && ! ( k % 2 ) ) {
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.50 );
            Mat.add ( n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.50 );
          }

          // interpolation in 2nd and 3rd coordinate
          if ( ! ( i % 2 ) && ( j % 2 ) && ( k % 2 ) ) {
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.25 );
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) + 1 ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.25 );
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.25 );
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) + 1 ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.25 );
          }

          // interpolation in 1st and 3rd coordinate
          if ( ( i % 2 ) && ! ( j % 2 ) && ( k % 2 ) ) {
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.25 );
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) + 1 ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.25 );
            Mat.add ( n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.25 );
            Mat.add ( n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) + 1 ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.25 );
          }

          // interpolation in 1st and 2nd coordinate
          if ( ( i % 2 ) && ( j % 2 ) && ! ( k % 2 ) ) {
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.25 );
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.25 );
            Mat.add ( n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.25 );
            Mat.add ( n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.25 );
          }

          // interpolation in 1st, 2nd and 3rd coordinate
          if ( ( i % 2 ) && ( j % 2 ) && ( k % 2 ) ) {
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.125 );
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) + 1 ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.125 );
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.125 );
            Mat.add ( n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) + 1 ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.125 );
            Mat.add ( n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.125 );
            Mat.add ( n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) + 1 ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.125 );
            Mat.add ( n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.125 );
            Mat.add ( n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) + 1 ),
                      n_fi_2 *     i         + n_fi *     j         +     k        ,
                      0.125 );
          }
        }
      }
    }

  } break;

  default:
    throw aol::Exception ( "qc::RestrictOp: Dimension must be 2 or 3", __FILE__, __LINE__ );
  }


}


template <typename RealType, qc::RestrictType RestrType>
void qc::RestrictOp<RealType, RestrType>::assemble_std_quoc_restrict_matrix ( aol::SparseMatrix<RealType> &Mat ) const {

  Mat.setZero();

  assemble_std_mg_restrict_matrix ( Mat );

  // then rescale rows:
  for ( int i = 0; i < Mat.getNumRows(); ++i ) {
    Mat.scaleRow ( i, ( aol::NumberTrait<RealType>::one / ( Mat.rowSum ( i ) ) ) );
  }

}


template <typename RealType, qc::RestrictType RestrType>
void qc::RestrictOp<RealType, RestrType>::assemble_throw_away_restrict_matrix ( aol::SparseMatrix<RealType> &Mat ) const {

  Mat.setZero();

  if ( _fineGrid.getGridDepth() != ( _coarseGrid.getGridDepth() + 1 ) ) {
    cerr << "Incompatible grid depths! Fine level = " << _fineGrid.getGridDepth() << ", coarse level = " << _coarseGrid.getGridDepth() << endl;
    throw aol::Exception ( "qc::RestrictOp::assemble_throw_away_restrict_matrix: Incompatible grid levels", __FILE__, __LINE__ );
  }

  if ( _fineGrid.getDimOfWorld() != _coarseGrid.getDimOfWorld() ) {
    throw aol::Exception ( "qc::RestrictOp::assemble_throw_away_restrict_matrix: Incompatible dimensions of coarse grid", __FILE__, __LINE__ );
  }

  if ( Mat.getNumRows() != _coarseGrid.getNumberOfNodes() || Mat.getNumCols() != _fineGrid.getNumberOfNodes() ) {
    throw aol::Exception ( "qc::RestrictOp::assemble_throw_away_restrict_matrix: Incompatible matrix", __FILE__, __LINE__ );
  }

  const int
    n_fi = _fineGrid.getWidth(),
    n_co = _coarseGrid.getWidth();

  switch ( _fineGrid.getDimOfWorld() ) {
  case qc::QC_2D: {

    for ( int j = 0; j < n_co; ++j ) {
      for ( int i = 0; i < n_co; ++i ) {
        Mat.set ( n_co * i + j,    n_fi * 2 * i + 2 * j,    1 );
      }
    }
  } break;

  case qc::QC_3D: {
    for ( int k = 0; k < n_co; ++k ) {
      for ( int j = 0; j < n_co; ++j ) {
        for ( int i = 0; i < n_co; ++i ) {
          Mat.set ( n_co * n_co * i + n_co * j + k,    n_fi * n_fi * 2 * i + n_fi * 2 * j + 2 * k,    1 );
        }
      }
    }
  } break;

  default:
    throw aol::Exception ( "qc::RestrictOp: Dimension must be 2 or 3", __FILE__, __LINE__ );
  }


}

template <typename RealType, qc::RestrictType RestrType>
void qc::RestrictOp<RealType, RestrType>::std_mg_restrict ( const aol::Vector<RealType> &fineVector,
                                                            aol::Vector<RealType> &coarseVector       ) const {

  switch ( _fineGrid.getDimOfWorld() ) {
  case qc::QC_2D: {

    const int coarseWidth = _coarseGrid.getWidth( );
    const int fineWidth = _fineGrid.getWidth( );

    if ( coarseWidth != ( fineWidth >> 1 ) + 1 ) {
      throw aol::Exception ( "qc::RestrictOp::std_mg_restrict: coarseWidth != ( fineWidth >> 1 ) + 1", __FILE__, __LINE__ );
    }

    const Array<RealType> fineArray ( fineVector, fineWidth, fineWidth, 1, aol::FLAT_COPY ); // 2D

    Array<RealType> coarseArray ( coarseVector, coarseWidth, coarseWidth );

    coarseArray.setZero();

    for ( int X = 0; X < coarseWidth; ++X ) {
      for ( int Y = 0; Y < coarseWidth; ++Y ) {
        const aol::Vec3<short> fc ( X << 1, Y << 1, 0 );
        coarseArray.set ( X, Y, 0, fineArray.get ( fc ) );
      }
    }

    for ( int X = 0; X < coarseWidth - 1; ++X ) {
      for ( int Y = 0; Y < coarseWidth - 1; ++Y ) {
        const aol::Vec3<short> fc ( X << 1, Y << 1, 0 );
        const RealType val = fineArray.get ( fc.x() + 1, fc.y() + 1, 0 ) / 4;

        coarseArray.add ( X,         Y, 0, val );
        coarseArray.add ( X + 1,     Y, 0, val );
        coarseArray.add ( X,     Y + 1, 0, val );
        coarseArray.add ( X + 1, Y + 1, 0, val );
      }
    }

    for ( int X = 0; X < coarseWidth - 1; ++X ) {
      for ( int Y = 0; Y < coarseWidth; ++Y ) {
        const aol::Vec3<short> fc ( X << 1, Y << 1, 0 );
        const RealType val = fineArray.get ( fc.x() + 1, fc.y(), 0 ) / 2;

        coarseArray.add ( X, Y, val );
        coarseArray.add ( X + 1, Y, val );
      }
    }

    for ( int X = 0; X < coarseWidth; ++X ) {
      for ( int Y = 0; Y < coarseWidth - 1; ++Y ) {
        const aol::Vec3<short> fc ( X << 1, Y << 1, 0 );
        const RealType val = fineArray.get ( fc.x(), fc.y() + 1 ) / 2;

        coarseArray.add ( X, Y + 1, val );
        coarseArray.add ( X, Y, val );
      }
    }
  } break;

  case qc::QC_3D: {
    coarseVector.setZero();
    const int nco = _coarseGrid.getWidth( );
    const int nfi = _fineGrid.getWidth( );

    Array<RealType> coarseArray ( coarseVector, nco, nco, nco, aol::FLAT_COPY ); // true = reference.
    const Array<RealType> fineArray ( fineVector, nfi, nfi, nfi, aol::FLAT_COPY );

    // coarse grid points that also exist in the fine grid.
    for ( int i = 0; i < nco; ++i ) {
      for ( int j = 0; j < nco; ++j ) {
        for ( int k = 0; k < nco; ++k ) {
          coarseArray.add ( i, j, k, fineArray.get ( 2*i, 2*j, 2*k ) );
        }
      }
    }

    // interpolate on lines: x, y, z
    for ( int r = 0; r < nco - 1; ++r ) { // interpolating index
      for ( int s = 0; s < nco  ; ++s ) { // indices for lines
        for ( int t = 0; t < nco  ; ++t ) { // dto
          const RealType
            value1 = fineArray.get ( 2 * s     ,  2 * t     , 2 * r + 1 ) / 2,
            value2 = fineArray.get ( 2 * s     ,  2 * r + 1 , 2 * t     ) / 2,
            value3 = fineArray.get ( 2 * r + 1 ,  2 * s     , 2 * t     ) / 2;
          coarseArray.add (  s    ,   t   ,  r    ,  value1 );
          coarseArray.add (  s    ,   t   , r + 1 ,  value1 );
          coarseArray.add (  s    ,   r   ,  t    ,  value2 );
          coarseArray.add (  s    , r + 1 ,  t    ,  value2 );
          coarseArray.add (  r    ,   s   ,  t    ,  value3 );
          coarseArray.add ( r + 1 ,   s   ,  t    ,  value3 );
        }
      }
    }

    // interpolate on planes: xy, xz, yz
    for ( int r = 0; r < nco - 1; ++r ) { // interpolating index
      for ( int s = 0; s < nco - 1; ++s ) { // dto
        for ( int t = 0; t < nco  ; ++t ) { // indices for lines
          const RealType
            value1 = fineArray.get ( 2 * r + 1 , 2 * s + 1 , 2 * t     ) / 4,
            value2 = fineArray.get ( 2 * r + 1 , 2 * t     , 2 * s + 1 ) / 4,
            value3 = fineArray.get ( 2 * t     , 2 * r + 1 , 2 * s + 1 ) / 4;


          coarseArray.add (  r    ,  s    ,  t    ,  value1 );
          coarseArray.add (  r    , s + 1 ,  t    ,  value1 );
          coarseArray.add ( r + 1 , s     ,  t    ,  value1 );
          coarseArray.add ( r + 1 , s + 1 ,  t    ,  value1 );

          coarseArray.add (  r    ,  t    ,  s    ,  value2 );
          coarseArray.add (  r    ,  t    , s + 1 ,  value2 );
          coarseArray.add ( r + 1 ,  t    ,  s    ,  value2 );
          coarseArray.add ( r + 1 ,  t    , s + 1 ,  value2 );

          coarseArray.add (  t    ,  r    ,  s    ,  value3 );
          coarseArray.add (  t    ,  r    , s + 1 ,  value3 );
          coarseArray.add (  t    , r + 1 ,  s    ,  value3 );
          coarseArray.add (  t    , r + 1 , s + 1 ,  value3 );
        }
      }
    }

    // interpolate in cubes: xyz
    for ( int r = 0; r < nco - 1; ++r ) { // interpolating index
      for ( int s = 0; s < nco - 1; ++s ) { // dto
        for ( int t = 0; t < nco - 1; ++t ) { // dto
          const RealType value = fineArray.get ( 2 * r + 1, 2 * s + 1, 2 * t + 1 ) / 8;

          coarseArray.add (  r    ,  s    ,  t    ,  value );
          coarseArray.add (  r    ,  s    , t + 1 ,  value );
          coarseArray.add (  r    , s + 1 ,  t    ,  value );
          coarseArray.add (  r    , s + 1 , t + 1 ,  value );
          coarseArray.add ( r + 1 ,  s    ,  t    ,  value );
          coarseArray.add ( r + 1 ,  s    , t + 1 ,  value );
          coarseArray.add ( r + 1 , s + 1 ,  t    ,  value );
          coarseArray.add ( r + 1 , s + 1 , t + 1 ,  value );
        }
      }
    }
  } break;

  default:
    throw aol::Exception ( "qc::RestrictOp: Dimension must be 2 or 3", __FILE__, __LINE__ );
  }
}

template <typename RealType, qc::RestrictType RestrType>
void qc::RestrictOp<RealType, RestrType>::std_quoc_restrict ( const aol::Vector<RealType> &fineVector,
                                                              aol::Vector<RealType> &coarseVector       ) const {


  switch ( _fineGrid.getDimOfWorld() ) {
  case qc::QC_2D: {

    const int coarseWidth = _coarseGrid.getWidth( );
    const int fineWidth = _fineGrid.getWidth( );
    int X, Y;

    if ( coarseWidth != ( fineWidth >> 1 ) + 1 ) {
      throw aol::Exception( "coarseWidth != ( fineWidth >> 1 ) + 1",
          "qcStandardMultigrid<REAL>::restrict" );
    }

    const qc::Array<RealType> fineArray( fineVector, fineWidth, fineWidth, 1, aol::FLAT_COPY ); // 2D

    qc::Array<RealType> coarseArray( coarseVector, coarseWidth, coarseWidth );

    coarseArray.setZero();

    for ( X = 0; X < coarseWidth; ++X ) {
      for ( Y = 0; Y < coarseWidth; ++Y ) {
        const aol::Vec3<short> fc( X << 1, Y << 1, 0 );
        coarseArray.set( X, Y, 0, fineArray.get( fc ) / 4 );
      }
    }

    for ( X = 0; X < coarseWidth - 1; ++X ) {
      for ( Y = 0; Y < coarseWidth - 1; ++Y ) {
        const aol::Vec3<short> fc( X << 1, Y << 1, 0 );
        const RealType val = fineArray.get( fc.x() + 1, fc.y() + 1, 0 ) / 16 ;

        coarseArray.add( X, Y, 0, val );
        coarseArray.add( X+1, Y,0, val );
        coarseArray.add( X, Y+1,0, val );
        coarseArray.add( X+1, Y+1,0, val );
      }
    }

    for ( X = 0; X < coarseWidth - 1; ++X ) {
      for ( Y = 0; Y < coarseWidth; ++Y ) {
        const aol::Vec3<short> fc( X << 1, Y << 1, 0 );
        const RealType val = fineArray.get( fc.x() + 1, fc.y(), 0 ) / 8;

        coarseArray.add( X, Y, val );
        coarseArray.add( X+1, Y, val );
      }
    }

    for ( X = 0; X < coarseWidth; ++X ) {
      for ( Y = 0; Y < coarseWidth - 1; ++Y ) {
        const aol::Vec3<short> fc( X << 1, Y << 1, 0 );
        const RealType val = fineArray.get( fc.x(), fc.y() + 1 ) / 8;

        coarseArray.add( X, Y+1, val );
        coarseArray.add( X, Y, val );
      }
    }

    for ( X = 1; X < coarseWidth - 1; ++X ) {
      coarseArray.set(               X,               0, 4 * ( coarseArray.get(               X,               0 ) / 3 ) );
      coarseArray.set(               X, coarseWidth - 1, 4 * ( coarseArray.get(               X, coarseWidth - 1 ) / 3 ) );
      coarseArray.set(               0,               X, 4 * ( coarseArray.get(               0,               X ) / 3 ) );
      coarseArray.set( coarseWidth - 1,               X, 4 * ( coarseArray.get( coarseWidth - 1,               X ) / 3 ) );
    }

    coarseArray.set(               0,               0, 16 * ( coarseArray.get(               0,               0 ) / 9 ) );
    coarseArray.set( coarseWidth - 1,               0, 16 * ( coarseArray.get( coarseWidth - 1,               0 ) / 9 ) );
    coarseArray.set(               0, coarseWidth - 1, 16 * ( coarseArray.get(               0, coarseWidth - 1 ) / 9 ) );
    coarseArray.set( coarseWidth - 1, coarseWidth - 1, 16 * ( coarseArray.get( coarseWidth - 1, coarseWidth - 1 ) / 9 ) );

  } break;

  case qc::QC_3D: {
    coarseVector.setZero();

    const int nco = _coarseGrid.getWidth( );
    const int nfi = _fineGrid.getWidth( );

    const qc::Array<RealType> fineArray( fineVector, nfi, nfi, nfi, aol::FLAT_COPY );
    qc::Array<RealType> coarseArray( coarseVector, nco, nco, nco, aol::FLAT_COPY );

    // factor 1/8 in all positions, later only boundary nodes need to be rescaled.

    // coarse grid points that also exist in the fine grid.
    for( int i = 0; i < nco; ++i ){
      for( int j = 0; j < nco; ++j ){
        for( int k = 0; k < nco; ++k ){
          coarseArray.add( i, j, k, fineArray.get( 2*i, 2*j, 2*k ) / 8 );
        }
      }
    }

    // interpolate on lines: x, y, z
    for( int r = 0; r < nco-1; ++r ){ // interpolating index
      for( int s = 0; s < nco  ; ++s ){ // indices for lines
        for( int t = 0; t < nco  ; ++t ){ // dto
          const RealType
            value1 = fineArray.get(  2*s ,   2*t  , 2*r+1 ) / 16,
            value2 = fineArray.get(  2*s ,  2*r+1 ,  2*t  ) / 16,
            value3 = fineArray.get( 2*r+1,   2*s  ,  2*t  ) / 16;
          coarseArray.add(  s ,  t ,  r ,  value1 );
          coarseArray.add(  s ,  t , r+1,  value1 );
          coarseArray.add(  s ,  r ,  t ,  value2 );
          coarseArray.add(  s , r+1,  t ,  value2 );
          coarseArray.add(  r ,  s ,  t ,  value3 );
          coarseArray.add( r+1,  s ,  t ,  value3 );
        }
      }
    }

    // interpolate on planes: xy, xz, yz
    for( int r = 0; r < nco-1; ++r ){ // interpolating index
      for( int s = 0; s < nco-1; ++s ){ // dto
        for( int t = 0; t < nco  ; ++t ){ // indices for lines
          const RealType
            value1 =  fineArray.get( 2*r+1, 2*s+1,  2*t  ) / 32,
            value2 =  fineArray.get( 2*r+1,  2*t , 2*s+1 ) / 32,
            value3 =  fineArray.get(  2*t , 2*r+1, 2*s+1 ) / 32;


          coarseArray.add(  r ,  s ,  t ,  value1 );
          coarseArray.add(  r , s+1,  t ,  value1 );
          coarseArray.add( r+1,  s ,  t ,  value1 );
          coarseArray.add( r+1, s+1,  t ,  value1 );

          coarseArray.add(  r ,  t ,  s ,  value2 );
          coarseArray.add(  r ,  t , s+1,  value2 );
          coarseArray.add( r+1,  t ,  s ,  value2 );
          coarseArray.add( r+1,  t , s+1,  value2 );

          coarseArray.add(  t ,  r ,  s ,  value3 );
          coarseArray.add(  t ,  r , s+1,  value3 );
          coarseArray.add(  t , r+1,  s ,  value3 );
          coarseArray.add(  t , r+1, s+1,  value3 );
        }
      }
    }

    // interpolate in cubes: xyz
    for( int r = 0; r < nco-1; ++r ){ // interpolating index
      for( int s = 0; s < nco-1; ++s ){ // dto
        for( int t = 0; t < nco-1; ++t ){ // dto
          const RealType value = fineArray.get( 2*r + 1, 2*s + 1, 2*t + 1 ) / 64;

          coarseArray.add(  r ,  s ,  t ,  value );
          coarseArray.add(  r ,  s , t+1,  value );
          coarseArray.add(  r , s+1,  t ,  value );
          coarseArray.add(  r , s+1, t+1,  value );
          coarseArray.add( r+1,  s ,  t ,  value );
          coarseArray.add( r+1,  s , t+1,  value );
          coarseArray.add( r+1, s+1,  t ,  value );
          coarseArray.add( r+1, s+1, t+1,  value );
        }
      }
    }

    // rescale:

    //  8 vertices
    coarseArray.set(   0  ,   0  ,   0  ,  64 * ( coarseArray.get(   0  ,   0  ,   0   ) / 27 ) );
    coarseArray.set(   0  ,   0  , nco-1,  64 * ( coarseArray.get(   0  ,   0  , nco-1 ) / 27 ) );
    coarseArray.set(   0  , nco-1,   0  ,  64 * ( coarseArray.get(   0  , nco-1,   0   ) / 27 ) );
    coarseArray.set(   0  , nco-1, nco-1,  64 * ( coarseArray.get(   0  , nco-1, nco-1 ) / 27 ) );
    coarseArray.set( nco-1,   0  ,   0  ,  64 * ( coarseArray.get( nco-1,   0  ,   0   ) / 27 ) );
    coarseArray.set( nco-1,   0  , nco-1,  64 * ( coarseArray.get( nco-1,   0  , nco-1 ) / 27 ) );
    coarseArray.set( nco-1, nco-1,   0  ,  64 * ( coarseArray.get( nco-1, nco-1,   0   ) / 27 ) );
    coarseArray.set( nco-1, nco-1, nco-1,  64 * ( coarseArray.get( nco-1, nco-1, nco-1 ) / 27 ) );

    // 12 lines
    for( int r = 1; r < nco-1; ++r ){
      coarseArray.set(   0  ,   0  ,   r  ,  64 * ( coarseArray.get(   0  ,   0  ,   r   ) / 36 ) );
      coarseArray.set(   0  , nco-1,   r  ,  64 * ( coarseArray.get(   0  , nco-1,   r   ) / 36 ) );
      coarseArray.set( nco-1,   0  ,   r  ,  64 * ( coarseArray.get( nco-1,   0  ,   r   ) / 36 ) );
      coarseArray.set( nco-1, nco-1,   r  ,  64 * ( coarseArray.get( nco-1, nco-1,   r   ) / 36 ) );

      coarseArray.set(   0  ,   r  ,   0  ,  64 * ( coarseArray.get(   0  ,   r  ,   0   ) / 36 ) );
      coarseArray.set(   0  ,   r  , nco-1,  64 * ( coarseArray.get(   0  ,   r  , nco-1 ) / 36 ) );
      coarseArray.set( nco-1,   r  ,   0  ,  64 * ( coarseArray.get( nco-1,   r  ,   0   ) / 36 ) );
      coarseArray.set( nco-1,   r  , nco-1,  64 * ( coarseArray.get( nco-1,   r  , nco-1 ) / 36 ) );

      coarseArray.set(   r  ,   0  ,   0  ,  64 * ( coarseArray.get(   r  ,   0  ,   0   ) / 36 ) );
      coarseArray.set(   r  ,   0  , nco-1,  64 * ( coarseArray.get(   r  ,   0  , nco-1 ) / 36 ) );
      coarseArray.set(   r  , nco-1,   0  ,  64 * ( coarseArray.get(   r  , nco-1,   0   ) / 36 ) );
      coarseArray.set(   r  , nco-1, nco-1,  64 * ( coarseArray.get(   r  , nco-1, nco-1 ) / 36 ) );
    }

    //  6 faces
    for( int r = 1; r < nco-1; ++r ){
      for( int s = 1; s < nco-1; ++s ){
        coarseArray.set(   0  ,   r  ,   s  , 64 * ( coarseArray.get(   0  ,   r  ,   s   ) / 48 ) );
        coarseArray.set( nco-1,   r  ,   s  , 64 * ( coarseArray.get( nco-1,   r  ,   s   ) / 48 ) );
        coarseArray.set(   r  ,   0  ,   s  , 64 * ( coarseArray.get(   r  ,   0  ,   s   ) / 48 ) );
        coarseArray.set(   r  , nco-1,   s  , 64 * ( coarseArray.get(   r  , nco-1,   s   ) / 48 ) );
        coarseArray.set(   r  ,   s  ,   0  , 64 * ( coarseArray.get(   r  ,   s  ,   0   ) / 48 ) );
        coarseArray.set(   r  ,   s  , nco-1, 64 * ( coarseArray.get(   r  ,   s  , nco-1 ) / 48 ) );
      }
    }

    // all inner points: nothing left to do.
  } break;

  default:
    throw aol::Exception ( "qc::RestrictOp: Dimension must be 2 or 3", __FILE__, __LINE__ );
  }
}


template <typename RealType, qc::RestrictType RestrType>
void qc::RestrictOp<RealType, RestrType>::throw_away_restrict ( const aol::Vector<RealType> &fineVector,
                                                                aol::Vector<RealType> &coarseVector       ) const {
  const int
    n_fi = _fineGrid.getWidth(),
    n_co = _coarseGrid.getWidth();

  switch ( _fineGrid.getDimOfWorld() ) {
  case qc::QC_2D: {

    for ( int j = 0; j < n_co; ++j ) {
      for ( int i = 0; i < n_co; ++i ) {
        coarseVector.set ( i*n_co + j,   fineVector.get ( 2*i*n_fi +  2*j ) );
      }
    }
  } break;

  case qc::QC_3D: {
    for ( int k = 0; k < n_co; ++k ) {
      for ( int j = 0; j < n_co; ++j ) {
        for ( int i = 0; i < n_co; ++i ) {
          coarseVector.set ( i*n_co*n_co + j*n_co + k,   fineVector.get ( 2*i*n_fi*n_fi +  2*j*n_fi + 2*k ) );
        }
      }
    }
  } break;

  default:
    throw aol::Exception ( "qc::RestrictOp: Dimension must be 2 or 3", __FILE__, __LINE__ );
  }
}



template <typename RealType, qc::RestrictType RestrType>
void qc::RestrictOp<RealType, RestrType>::periodic_restrict ( const aol::Vector<RealType> &fineVector,
                                                              aol::Vector<RealType> &coarseVector       ) const {

  switch ( _fineGrid.getDimOfWorld() ) {
  case qc::QC_2D: {

    const int coarseWidth = _coarseGrid.getWidth( );
    const int fineWidth = _fineGrid.getWidth( );
    int X, Y;

    if ( coarseWidth != ( fineWidth >> 1 ) + 1 ) {
      throw aol::Exception( "coarseWidth != ( fineWidth >> 1 ) + 1",
          "qcStandardMultigrid<REAL>::restrict" );
    }

    const qc::Array<RealType> fineArray( fineVector, fineWidth, fineWidth, 1, aol::FLAT_COPY ); // 2D

    qc::Array<RealType> coarseArray( coarseVector, coarseWidth, coarseWidth );

    coarseArray.setZero();

    //inner nodes
    for ( X = 1; X < coarseWidth - 1; ++X ) {
      for ( Y = 1; Y < coarseWidth - 1; ++Y ) {
        const aol::Vec3<short> fc( X << 1, Y << 1, 0 );
        coarseArray.set( X, Y, 0, fineArray.get( fc ) / 4 );
      }
    }
    
    //simple boundary nodes
    for ( X = 1; X < coarseWidth - 1; ++X ) {
        const aol::Vec3<short> fc1( X << 1, 0, 0 );
        coarseArray.set(               X,               0, 0, fineArray.get( fc1 ) / 8 );
        const aol::Vec3<short> fc2( X << 1, (coarseWidth - 1) << 1, 0 );
        coarseArray.set(               X, coarseWidth - 1, 0, fineArray.get( fc2 ) / 8 );
        const aol::Vec3<short> fc3( 0, X << 1, 0 );
        coarseArray.set(               0,               X, 0, fineArray.get( fc3 ) / 8 );
        const aol::Vec3<short> fc4( (coarseWidth - 1) << 1, X << 1, 0 );
        coarseArray.set( coarseWidth - 1,               X, 0, fineArray.get( fc4 ) / 8 );
    }
    
    //4 boundary nodes
    const aol::Vec3<short> fc00( 0, 0, 0 );
    coarseArray.set(               0,               0, 0, fineArray.get( fc00 ) / 16 );
    const aol::Vec3<short> fcN0( (coarseWidth - 1) << 1 , 0, 0 );
    coarseArray.set( coarseWidth - 1,               0, 0, fineArray.get( fcN0 ) / 16 );
    const aol::Vec3<short> fc0N( 0, (coarseWidth - 1) << 1, 0 );
    coarseArray.set(               0, coarseWidth - 1, 0, fineArray.get( fc0N ) / 16 );
    const aol::Vec3<short> fcNN( (coarseWidth - 1) << 1, (coarseWidth - 1) << 1 , 0 );
    coarseArray.set( coarseWidth - 1, coarseWidth - 1, 0, fineArray.get( fcNN ) / 16 );
    
    
    //diagonal neighbours (for all types)
    for ( X = 0; X < coarseWidth - 1; ++X ) {
      for ( Y = 0; Y < coarseWidth - 1; ++Y ) {
        const aol::Vec3<short> fc( X << 1, Y << 1, 0 );
        const RealType val = fineArray.get( fc.x() + 1, fc.y() + 1, 0 ) / 16 ;
        coarseArray.add( X, Y, 0, val );
        coarseArray.add( X+1, Y,0, val );
        coarseArray.add( X, Y+1,0, val );
        coarseArray.add( X+1, Y+1,0, val );
      }
    }

    //right neighbours (for all types)
    for ( X = 0; X < coarseWidth - 1; ++X ) {
      
      for ( Y = 1; Y < coarseWidth - 1; ++Y ) {
        const aol::Vec3<short> fc( X << 1, Y << 1, 0 );
        const RealType val = fineArray.get( fc.x() + 1, fc.y(), 0 ) / 8;
        coarseArray.add( X, Y, val );
        coarseArray.add( X+1, Y, val );
      }
      
      const aol::Vec3<short> fcY0( X << 1, 0, 0 );
      const RealType valY0 = fineArray.get( fcY0.x() + 1, fcY0.y(), 0 ) / 16;
      coarseArray.add( X, 0, valY0 );
      coarseArray.add( X+1, Y, valY0 );
      
      const aol::Vec3<short> fcYN( X << 1, (coarseWidth - 1) << 1, 0 );
      const RealType valYN = fineArray.get( fcYN.x() + 1, fcYN.y(), 0 ) / 16;
      coarseArray.add( X, 0, valYN );
      coarseArray.add( X+1, coarseWidth-1, valYN );
      
    }

    //up and down neighbours (for all types)
    for ( Y = 0; Y < coarseWidth - 1; ++Y ) {
      
      for ( X = 1; X < coarseWidth - 1; ++X ) {
        const aol::Vec3<short> fc( X << 1, Y << 1, 0 );
        const RealType val = fineArray.get( fc.x(), fc.y() + 1 ) / 8;
        coarseArray.add( X, Y+1, val );
        coarseArray.add( X, Y, val );
      }
      
      const aol::Vec3<short> fcX0( 0, Y << 1, 0 );
      const RealType valX0 = fineArray.get( fcX0.x(), fcX0.y() + 1 ) / 16;
      coarseArray.add( 0, Y+1, valX0 );
      coarseArray.add( 0, Y, valX0 );
      
      const aol::Vec3<short> fcXN( (coarseWidth - 1) << 1, Y << 1, 0 );
      const RealType valXN = fineArray.get( fcXN.x(), fcXN.y() + 1 ) / 16;
      coarseArray.add( coarseWidth - 1, Y+1, valXN );
      coarseArray.add( coarseWidth - 1, Y, valXN );
      
    }
    
    //identify peridodic nodes
    for ( X = 1; X < coarseWidth - 1; ++X ) {
      coarseArray.add( X, 0, coarseArray.get( X, coarseWidth - 1 ) );
      coarseArray.set( X, coarseWidth - 1, coarseArray.get( X, 0 ) );
      coarseArray.add( 0, X, coarseArray.get( coarseWidth - 1, X ) );
      coarseArray.set( coarseWidth - 1, X, coarseArray.get( 0, X ) );
    }

    coarseArray.add( 0, 0, coarseArray.get( 0, coarseWidth - 1 ) );
    coarseArray.add( 0, 0, coarseArray.get( coarseWidth - 1, 0 ) );
    coarseArray.add( 0, 0, coarseArray.get( coarseWidth - 1, coarseWidth - 1 ) );
    coarseArray.set( coarseWidth - 1, 0, coarseArray.get( 0, 0 ) );
    coarseArray.set( 0, coarseWidth - 1, coarseArray.get( 0, 0 ) );
    coarseArray.set( coarseWidth - 1, coarseWidth - 1, coarseArray.get( 0, 0 ) );

  } break;

  default:
    throw aol::Exception ( "qc::RestrictOp: Dimension must be 2", __FILE__, __LINE__ );
  }
}




template class qc::RestrictOp< float,       qc::STD_MG_RESTRICT >;
template class qc::RestrictOp< double,      qc::STD_MG_RESTRICT >;
template class qc::RestrictOp< long double, qc::STD_MG_RESTRICT >;

template class qc::RestrictOp< float,       qc::STD_QUOC_RESTRICT >;
template class qc::RestrictOp< double,      qc::STD_QUOC_RESTRICT >;
template class qc::RestrictOp< long double, qc::STD_QUOC_RESTRICT >;

template class qc::RestrictOp< float,       qc::THROW_AWAY_RESTRICT >;
template class qc::RestrictOp< double,      qc::THROW_AWAY_RESTRICT >;
template class qc::RestrictOp< long double, qc::THROW_AWAY_RESTRICT >;

template class qc::RestrictOp< float,       qc::PERIODIC_RESTRICT >;
template class qc::RestrictOp< double,      qc::PERIODIC_RESTRICT >;
template class qc::RestrictOp< long double, qc::PERIODIC_RESTRICT >;

