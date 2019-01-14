#include <prolongation.h>


template <typename RealType>
void qc::ProlongOp<RealType>::assemble_std_prolong_matrix ( aol::SparseMatrix<RealType> &Mat ) const {

  if ( _fineGrid.getGridDepth() != ( _coarseGrid.getGridDepth() + 1 ) ) {
    cerr << "Incompatible grid depths! Fine level = " << _fineGrid.getGridDepth() << ", coarse level = " << _coarseGrid.getGridDepth() << endl;
    throw aol::Exception ( "qc::ProlongOp::assemble_std_prolong_matrix(): Incompatible grid levels", __FILE__, __LINE__ );
  }

  if ( _fineGrid.getDimOfWorld() != _coarseGrid.getDimOfWorld() ) {
    throw aol::Exception ( "qc::ProlongOp::assemble_std_prolong_matrix(): Incompatible dimensions of coarse grid", __FILE__, __LINE__ );
  }

  if ( Mat.getNumRows() != _fineGrid.getNumberOfNodes() || Mat.getNumCols() != _coarseGrid.getNumberOfNodes() ) {
    throw aol::Exception ( "qc::ProlongOp::assemble_std_prolong_matrix: Incompatible matrix", __FILE__, __LINE__ );
  }

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
          Mat.add ( n_fi * i + j,                           n_co * ( ( i / 2 ) ) + ( ( j / 2 ) ),           1.0  );
        }

        // interpolation / restriction in 2nd coordinate
        if ( ! ( i % 2 ) && ( j % 2 ) ) {
          Mat.add ( n_fi * i + j,                           n_co * ( ( i / 2 ) ) + ( ( j / 2 ) ),           0.5  );
          Mat.add ( n_fi * i + j,                           n_co * ( ( i / 2 ) ) + ( ( j / 2 ) + 1 ),       0.5  );
        }

        // interpolation / restriction in 1st coordinate
        if ( ( i % 2 ) && ! ( j % 2 ) ) {
          Mat.add ( n_fi * i + j,                           n_co * ( ( i / 2 ) ) + ( ( j / 2 ) ),           0.5  );
          Mat.add ( n_fi * i + j,                           n_co * ( ( i / 2 ) + 1 ) + ( ( j / 2 ) ),       0.5  );
        }

        // interpolation / restriction in 1st and 2nd coordinate
        if ( ( i % 2 ) && ( j % 2 ) ) {
          Mat.add ( n_fi * i + j,                           n_co * ( ( i / 2 ) ) + ( ( j / 2 ) ),           0.25 );
          Mat.add ( n_fi * i + j,                           n_co * ( ( i / 2 ) ) + ( ( j / 2 ) + 1 ),       0.25 );
          Mat.add ( n_fi * i + j,                           n_co * ( ( i / 2 ) + 1 ) + ( ( j / 2 ) ),       0.25 );
          Mat.add ( n_fi * i + j,                           n_co * ( ( i / 2 ) + 1 ) + ( ( j / 2 ) + 1 ),   0.25 );
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
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      1.00 );
          }

          // interpolation in 3rd coordinate
          if ( ! ( i % 2 ) && ! ( j % 2 ) && ( k % 2 ) ) {
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      0.50 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) + 1 ),
                      0.50 );
          }

          // interpolation in 2nd coordinate
          if ( ! ( i % 2 ) && ( j % 2 ) && ! ( k % 2 ) ) {
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      0.50 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) ),
                      0.50 );
          }

          // interpolation in 1st coordinate
          if ( ( i % 2 ) && ! ( j % 2 ) && ! ( k % 2 ) ) {
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      0.50 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      0.50 );
          }

          // interpolation in 2nd and 3rd coordinate
          if ( ! ( i % 2 ) && ( j % 2 ) && ( k % 2 ) ) {
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      0.25 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) + 1 ),
                      0.25 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) ),
                      0.25 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) + 1 ),
                      0.25 );
          }

          // interpolation in 1st and 3rd coordinate
          if ( ( i % 2 ) && ! ( j % 2 ) && ( k % 2 ) ) {
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      0.25 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) + 1 ),
                      0.25 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      0.25 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) + 1 ),
                      0.25 );
          }

          // interpolation in 1st and 2nd coordinate
          if ( ( i % 2 ) && ( j % 2 ) && ! ( k % 2 ) ) {
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      0.25 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) ),
                      0.25 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      0.25 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) ),
                      0.25 );
          }

          // interpolation in 1st, 2nd and 3rd coordinate
          if ( ( i % 2 ) && ( j % 2 ) && ( k % 2 ) ) {
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      0.125 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) + 1 ),
                      0.125 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) ),
                      0.125 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) + 1 ),
                      0.125 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) ),
                      0.125 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) ) + ( ( k / 2 ) + 1 ),
                      0.125 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) ),
                      0.125 );
            Mat.add ( n_fi_2 *     i         + n_fi *     j         +     k        ,
                      n_co_2 * ( ( i / 2 ) + 1 ) + n_co * ( ( j / 2 ) + 1 ) + ( ( k / 2 ) + 1 ),
                      0.125 );
          }

        }
      }
    }

  } break;

  default:
    throw aol::Exception ( "qc::ProlongOp::assembleMatrix: Dimension must be 2 or 3", __FILE__, __LINE__ );
  }

}


template <typename RealType>
void qc::ProlongOp<RealType>::std_mg_prolongate ( const GridDefinition &FineGrid,
                                                  aol::Vector<RealType> &FineVector,
                                                  const GridDefinition &CoarseGrid,
                                                  const aol::Vector<RealType> &CoarseVector ) {

  int coarseCounter = 0, fineCounter = 0, tmpCounter = 0;
  int x, y, i;

  switch ( FineGrid.getDimOfWorld() ) {
  case qc::QC_2D: {
    const int coarseWidth = CoarseGrid.getWidth() - 1;
    const int fineWidth = FineGrid.getWidth();

    // Interpolate even lines
    for ( y = 0; y <= coarseWidth; ++y ) {
      RealType u = CoarseVector.get ( coarseCounter++ );
      for ( x = 0; x < coarseWidth; ++x ) {
        RealType w = CoarseVector.get ( coarseCounter++ );

        FineVector[fineCounter++] = u;
        FineVector[fineCounter++] = ( u + w ) / 2;
        u = w;
      }
      FineVector[fineCounter++] = u;
      fineCounter += fineWidth;
    }

    // Interpolate odd lines
    fineCounter = fineWidth;
    tmpCounter = 0;
    coarseCounter = fineWidth << 1;
    for ( y = 0; y < coarseWidth; ++y ) {
      for ( x = 0; x < fineWidth; ++x ) {
        FineVector[fineCounter++] = ( FineVector[tmpCounter++] + FineVector[coarseCounter++] ) / 2;
      }
      fineCounter   += fineWidth;
      tmpCounter    += fineWidth;
      coarseCounter += fineWidth;
    }

  } break;

  case QC_3D: {
    const int coarseWidth = CoarseGrid.getWidth( );
    const int fineWidth = FineGrid.getWidth( );
    int X, Y, Z;

    Array<RealType> fineArray ( FineVector, fineWidth, fineWidth, fineWidth, aol::FLAT_COPY );
    const Array<RealType> coarseArray ( CoarseVector, coarseWidth, coarseWidth, coarseWidth, aol::FLAT_COPY );

    /* run through all the nodes */
    for ( X = 0; X < coarseWidth; ++X ) {
      for ( Y = 0; Y < coarseWidth; ++Y ) {
        for ( Z = 0; Z < coarseWidth; ++Z ) {
          const aol::Vec3<int> fc ( X << 1, Y << 1, Z << 1 );
          const RealType u = coarseArray.get ( X, Y, Z );
          fineArray.set ( fc[0], fc[1], fc[2], u );
        }
      }
    }

    const aol::Vec3<int> co = coarseArray.getOffset( );

    const int _a3d_off1[ 8 ] = {
                                 coarseArray.getIndexOffset ( 0, 0, 0 ),
                                 coarseArray.getIndexOffset ( 1, 0, 0 ),
                                 coarseArray.getIndexOffset ( 0, 1, 0 ),
                                 coarseArray.getIndexOffset ( 0, 0, 1 ),
                                 coarseArray.getIndexOffset ( 1, 1, 0 ),
                                 coarseArray.getIndexOffset ( 1, 0, 1 ),
                                 coarseArray.getIndexOffset ( 0, 1, 1 ),
                                 coarseArray.getIndexOffset ( 1, 1, 1 )
                               } ;


    /* run through all nodes at the center of elements */
    for ( X = 0; X < coarseWidth - 1; ++X ) {
      for ( Y = 0; Y < coarseWidth - 1; ++Y ) {
        for ( Z = 0; Z < coarseWidth - 1; ++Z ) {
          const aol::Vec3<int> fc ( X << 1, Y << 1, Z << 1 );
          const int in = coarseArray.index ( X, Y, Z );
          RealType v = coarseArray.get ( in );
          for ( i = 1; i < 8; ++i ) {
            v += coarseArray.get ( in + _a3d_off1[i] );
          }
          fineArray.set ( fc.x() + 1, fc.y() + 1, fc.z() + 1, v / 8 );
        }
      }
    }

    /* interpolate all nodes at the center of YZ-aligned faces */
    for ( X = 0; X < coarseWidth; ++X ) {
      for ( Y = 0; Y < coarseWidth - 1; ++Y ) {
        for ( Z = 0; Z < coarseWidth - 1; ++Z ) {
          const aol::Vec3<int> fc ( X << 1, Y << 1, Z << 1 );
          const int in = coarseArray.index ( X, Y, Z );
          fineArray.set ( fc.x(), fc.y() + 1, fc.z() + 1,
                          ( coarseArray.get ( in ) +
                            coarseArray.get ( in + co.y() ) +
                            coarseArray.get ( in + co.z() ) +
                            coarseArray.get ( in + co.y() + co.z() ) ) / 4 );
        }
      }
    }

    /* interpolate all nodes at the center of XZ-aligned faces */
    for ( X = 0; X < coarseWidth - 1; ++X ) {
      for ( Y = 0; Y < coarseWidth; ++Y ) {
        for ( Z = 0; Z < coarseWidth - 1; ++Z ) {
          const aol::Vec3<int> fc ( X << 1, Y << 1, Z << 1 );
          const int in = coarseArray.index ( X, Y, Z );
          fineArray.set ( fc.x() + 1, fc.y(), fc.z() + 1,
                          ( coarseArray.get ( in ) +
                            coarseArray.get ( in + co.x() ) +
                            coarseArray.get ( in + co.z() ) +
                            coarseArray.get ( in + co.x() + co.z() ) ) / 4 );
        }
      }
    }

    /* interpolate all nodes at the center of XY-aligned faces */
    for ( X = 0; X < coarseWidth - 1; ++X ) {
      for ( Y = 0; Y < coarseWidth - 1; ++Y ) {
        for ( Z = 0; Z < coarseWidth; ++Z ) {
          const aol::Vec3<int> fc ( X << 1, Y << 1, Z << 1 );
          const int in = coarseArray.index ( X, Y, Z );
          fineArray.set ( fc.x() + 1, fc.y() + 1, fc.z(),
                          ( coarseArray.get ( in ) +
                            coarseArray.get ( in + co.y() ) +
                            coarseArray.get ( in + co.x() ) +
                            coarseArray.get ( in + co.y() + co.x() ) ) / 4 );
        }
      }
    }

    /* interpolate all nodes at the center of x-aligned edges */
    for ( X = 0; X < coarseWidth - 1; ++X ) {
      for ( Y = 0; Y < coarseWidth; ++Y ) {
        for ( Z = 0; Z < coarseWidth; ++Z ) {
          const aol::Vec3<int> fc ( X << 1, Y << 1, Z << 1 );
          const int in = coarseArray.index ( X, Y, Z );
          fineArray.set ( fc.x() + 1, fc.y(), fc.z(),
                          ( coarseArray.get ( in ) + coarseArray.get ( in + co.x() ) ) / 2 );
        }
      }
    }

    /* interpolate all nodes at the center of y-aligned edges */
    for ( X = 0; X < coarseWidth; ++X ) {
      for ( Y = 0; Y < coarseWidth - 1; ++Y ) {
        for ( Z = 0; Z < coarseWidth; ++Z ) {
          const aol::Vec3<int> fc ( X << 1, Y << 1, Z << 1 );
          const int in = coarseArray.index ( X, Y, Z );
          fineArray.set ( fc.x(), fc.y() + 1, fc.z(),
                          ( coarseArray.get ( in ) + coarseArray.get ( in + co.y() ) ) / 2 );
        }
      }
    }

    /* interpolate all nodes at the center of z-aligned edges */
    for ( X = 0; X < coarseWidth; ++X ) {
      for ( Y = 0; Y < coarseWidth; ++Y ) {
        for ( Z = 0; Z < coarseWidth - 1; ++Z ) {
          const aol::Vec3<int> fc ( X << 1, Y << 1, Z << 1 );
          const int in = coarseArray.index ( X, Y, Z );
          fineArray.set ( fc.x(), fc.y(), fc.z() + 1,
                          ( coarseArray.get ( in ) + coarseArray.get ( in + co.z() ) ) / 2 );
        }
      }
    }

  } break;

  default:
    throw aol::Exception ( "qc::ProlongOp: Dimension must be 2 or 3", __FILE__, __LINE__ );
  }
}

template class qc::ProlongOp<float>;
template class qc::ProlongOp<double>;
template class qc::ProlongOp<long double>;
