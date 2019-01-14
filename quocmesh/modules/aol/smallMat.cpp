#include <smallMat.h>
#include <polynomial.h>

namespace aol {

template <> const Matrix22<float> ZTrait<Matrix22<float> >::zero = Matrix22<float> ( 0, 0, 0, 0 );
template <> const Matrix22<float> ZOTrait<Matrix22<float> >::one = Matrix22<float> ( 1, 0, 0, 1 );
template <> const Matrix22<double> ZTrait<Matrix22<double> >::zero = Matrix22<double> ( 0, 0, 0, 0 );
template <> const Matrix22<double> ZOTrait<Matrix22<double> >::one = Matrix22<double> ( 1, 0, 0, 1 );
template <> const Matrix22<long double> ZTrait<Matrix22<long double> >::zero = Matrix22<long double> ( 0, 0, 0, 0 );
template <> const Matrix22<long double> ZOTrait<Matrix22<long double> >::one = Matrix22<long double> ( 1, 0, 0, 1 );

template <> const Matrix33<float> ZTrait<Matrix33<float> >::zero = Matrix33<float> ( 0, 0, 0, 0, 0, 0, 0, 0, 0 );
template <> const Matrix33<float> ZOTrait<Matrix33<float> >::one = Matrix33<float> ( 1, 0, 0, 0, 1, 0, 0, 0, 1 );
template <> const Matrix33<double> ZTrait<Matrix33<double> >::zero = Matrix33<double> ( 0, 0, 0, 0, 0, 0, 0, 0, 0 );
template <> const Matrix33<double> ZOTrait<Matrix33<double> >::one = Matrix33<double> ( 1, 0, 0, 0, 1, 0, 0, 0, 1 );
template <> const Matrix33<long double> ZTrait<Matrix33<long double> >::zero = Matrix33<long double> ( 0, 0, 0, 0, 0, 0, 0, 0, 0 );
template <> const Matrix33<long double> ZOTrait<Matrix33<long double> >::one = Matrix33<long double> ( 1, 0, 0, 0, 1, 0, 0, 0, 1 );

template <typename _DataType>
void Matrix22<_DataType>::eigenValues ( Vec2<_DataType> &Eig ) const {
  aol::Polynomial<_DataType> p ( 2 );
  p.setCoeff ( 0, this->_row[0][0]*this->_row[1][1] - this->_row[0][1]*this->_row[1][0] );
  p.setCoeff ( 1, -this->_row[0][0] - this->_row[1][1] );
  p.setCoeff ( 2, 1. );
  try {
      p.determineZeros ( Eig[0], Eig[1] );
  } catch ( ... ) {
    throw aol::Exception( "Matrix22<DataType>::eigenValues: no real eigenvalues!", __FILE__, __LINE__ );
  }
  if ( Eig[0] < Eig[1] ) {
    _DataType v = Eig[0];
    Eig[0] = Eig[1];
    Eig[1] = v;
  }
}


namespace {

// Jacobi iteration for the solution of eigenvectors/eigenvalues of a nxn
// real symmetric matrix. Square nxn matrix a; size of matrix in n;
// output eigenvalues in w; and output eigenvectors in v. Resulting
// eigenvalues/vectors are sorted in decreasing order; eigenvectors are
// normalized.

template<class T>
inline void INLINE_ROTATE ( T **a, int &i, int &j, int &k, int &l, T &g, T &h, T &s, T &tau ) {
  g = a[i][j];
  h = a[k][l];
  a[i][j] = g - s * ( h + g * tau );
  a[k][l] = h + s * ( g - h * tau );
}


template<class T>
int jacobi_n ( T **a, int n, T *w, T **v ) {

  const int INLINE_MAX_ROTATIONS = 50;

  int i, j, k, iq, ip, numPos;
  T tresh, theta, tau, t, sm, s, h, g, c, tmp;
  T bspace[4], zspace[4];
  T *b = bspace;
  T *z = zspace;

  // only allocate memory if the matrix is large
  if ( n > 4 ) {  b = new T[n];  z = new T[n];  }

  // initialize
  for ( ip = 0; ip < n; ++ip ) {
    for ( iq = 0; iq < n; ++iq )
      v[ip][iq] = 0.0;
    v[ip][ip] = 1.0;
  }
  for ( ip = 0; ip < n; ++ip ) {
    b[ip] = w[ip] = a[ip][ip]; z[ip] = 0.0;
  }


  for ( i = 0; i < INLINE_MAX_ROTATIONS; ++i ) { // begin rotation sequence
    sm = 0.0;
    for ( ip = 0; ip < n - 1; ++ip )
      for ( iq = ip + 1; iq < n; ++iq )
        sm += fabs ( a[ip][iq] );
    if ( sm == 0.0 ) break;

    if ( i < 3 )  tresh = sm / ( 5 * n * n ); // first 3 sweeps
    else tresh = 0.0;

    for ( ip = 0; ip < n - 1; ++ip ) {
      for ( iq = ip + 1; iq < n; ++iq ) {
        g = 100 * fabs ( a[ip][iq] );

        // after 4 sweeps
        if ( i > 3 && ( fabs ( w[ip] ) + g ) == fabs ( w[ip] ) && ( fabs ( w[iq] ) + g ) == fabs ( w[iq] ) )
          a[ip][iq] = 0.0;
        else if ( fabs ( a[ip][iq] ) > tresh ) {
          h = w[iq] - w[ip];
          if ( ( fabs ( h ) + g ) == fabs ( h ) )
            t = ( a[ip][iq] ) / h;
          else {
            theta = h / ( 2 * a[ip][iq] );
            t = aol::ZOTrait<T>::one / ( fabs ( theta ) + sqrt ( aol::ZOTrait<T>::one + theta * theta ) );
            if ( theta < 0.0 ) t = -t;
          }
          c = aol::ZOTrait<T>::one / sqrt ( 1 + t * t );
          s = t * c;
          tau = s / ( aol::ZOTrait<T>::one + c );
          h = t * a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          w[ip] -= h;
          w[iq] += h;
          a[ip][iq] = 0.0;

          // ip already shifted left by 1 unit
          for ( j = 0;j <= ip - 1; ++j ) {
            INLINE_ROTATE ( a, j, ip, j, iq, g, h, s, tau );
          }
          // ip and iq already shifted left by 1 unit
          for ( j = ip + 1;j <= iq - 1; ++j ) {
            INLINE_ROTATE ( a, ip, j, j, iq, g, h, s, tau );
          }
          // iq already shifted left by 1 unit
          for ( j = iq + 1; j < n; ++j ) {
            INLINE_ROTATE ( a, ip, j, iq, j, g, h, s, tau );
          }
          for ( j = 0; j < n; ++j ) {
            INLINE_ROTATE ( v, j, ip, j, iq, g, h, s, tau );
          }
        }
      }
    }

    for ( ip = 0; ip < n; ++ip ) {
      b[ip] += z[ip];
      w[ip] = b[ip];
      z[ip] = 0.0;
    }
  }

  //// this should be NEVER called
  if ( i >= INLINE_MAX_ROTATIONS ) {
    cout << "Jacobi: Error extracting eigenfunctions" << endl;

    if ( n > 4 ) {   //release if we allocated anything
      delete [] b; delete [] z;
    }

    return 0;
  }

  // sort eigenfunctions                 these changes do not affect accuracy
  for ( j = 0; j < n - 1; ++j ) {        // boundary incorrect
    k = j; tmp = w[k];
    for ( i = j + 1; i < n; ++i ) {        // boundary incorrect, shifted already
      if ( w[i] >= tmp ) {                // why exchage if same?
        k = i; tmp = w[k]; }
    }
    if ( k != j ) {
      w[k] = w[j]; w[j] = tmp;
      for ( i = 0; i < n; ++i ) { tmp = v[i][j]; v[i][j] = v[i][k]; v[i][k] = tmp; }
    }
  }

  // ensure eigenvector consistency (i.e., Jacobi can compute vectors that
  // are negative of one another (.707,.707,0) and (-.707,-.707,0). This can
  // reek havoc in hyperstreamline/other stuff. We will select the most
  // positive eigenvector.
  for ( j = 0; j < n; ++j ) {
    for ( numPos = 0, i = 0; i < n; ++i )
      if ( v[i][j] >= 0.0 ) numPos++;
    if ( numPos < ceil ( double ( n ) / double ( 2.0 ) ) )
      for ( i = 0; i < n; ++i ) v[i][j] *= -1.0;
  }

  if ( n > 4 ) {   //release if we allocated anything
    delete [] b; delete [] z;
  }
  return 1;

}

}
  template <class _DataType>
  void Matrix33Symm<_DataType>::eigenVectors ( Vec3<_DataType> &eigenVals, Matrix33<_DataType> &eigenVecs ) {
    _DataType *mat[3], *vecs[3],
              vals[3],
              mat0[3], mat1[3], mat2[3],
              vecs0[3], vecs1[3], vecs2[3];
              mat[0]  = mat0;
              mat[1]  = mat1;
              mat[2]    = mat2;
              vecs[0] = vecs0;
              vecs[1] = vecs1;
              vecs[2] = vecs2;

    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        mat[i][j] = get ( i, j );
      }
    }

    jacobi_n <_DataType> ( mat, 3, vals, vecs );

    for ( int i = 0; i < 3; ++i ) {
      eigenVals[i] = vals[i];
      for ( int j = 0; j < 3; ++j ) {
        eigenVecs[i][j] = vecs[i][j];
      }
    }
  }

template <typename _DataType> Matrix33<_DataType>::Matrix33 ( const Matrix33Symm<_DataType>& mat ) {
  fill ( mat[0][0], mat [0][1], mat [0][2],
         mat[1][0], mat [1][1], mat [1][2],
         mat[2][0], mat [2][1], mat [2][2] );
}

template class Matrix22<float>;
template class Matrix22<double>;
template class Matrix22<long double>;
template class Matrix33<float>;
template class Matrix33<double>;
template class Matrix33<long double>;
template class Matrix33Symm<float>;
template class Matrix33Symm<double>;
template class Matrix33Symm<long double>;

} // end namespace
