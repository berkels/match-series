#include <matrixInverse.h>

namespace aol {

template <class DataType>
LUInverse<DataType>::LUInverse ( const LUInverse<DataType>& inv )
    : Matrix<DataType> ( inv.getNumCols (), inv.getNumRows () ), _decomposition ( inv._decomposition ), _permut ( inv._permut ) {}

template <class DataType>
LUInverse<DataType>& LUInverse<DataType>::operator = ( const LUInverse<DataType>& inv ) {
  // Beware of self-assignment
  if ( this == &inv ) return *this;

  _decomposition = inv._decomposition;
  _permut = inv._permut;

  return *this;
}

template <class DataType>
LUInverse<DataType>::~LUInverse () {}

template <class DataType>
LUInverse<DataType>::LUInverse ( const Matrix<DataType>& m )
    : Matrix<DataType> ( m.getNumCols (), m.getNumRows () ), _decomposition ( m.getNumCols (), m.getNumRows () ), _permut ( m.getNumRows () ) {
  _decomposition.makeLU ( m, _permut );
}

template <class DataType>
LUInverse<DataType>::LUInverse ( const int N, const aol::FullMatrix<DataType> &Decomp, const aol::PermutationMatrix<DataType> &Permut )
  : Matrix<DataType> ( N, N ), _decomposition ( Decomp ), _permut ( Permut ) {
}

template <class DataType>
void LUInverse<DataType>::apply ( const Vector<DataType>& arg, Vector<DataType>& dest ) const {
  _decomposition.LUSolve ( dest, arg, _permut );
}

template <class DataType>
void LUInverse<DataType>::applyAdd ( const Vector<DataType>& arg, Vector<DataType>& dest ) const {
  aol::Vector<DataType> temp ( dest );
  apply ( arg, dest );
  dest += temp;
}

template <class DataType>
DataType LUInverse<DataType>::get ( int row, int col ) const {
    // Use getFM instead
    cerr << "Warning: LUInverse<DataType>::get is slow." << endl;

    Vector<DataType> colv ( this->_numCols );
    Vector<DataType> resv ( this->_numRows );
    colv.setZero ();
    colv [col] = 1.0;
    apply ( colv, resv );
    return resv [row];
  }

template <class DataType>
FullMatrix<DataType> LUInverse<DataType>::getFM () const {
  FullMatrix<DataType> m ( this->_numRows, this->_numCols );
  Vector<DataType> x ( this->_numRows ), b ( this->_numCols );

  b.setZero ();
  for ( int i = 0; i < this->_numCols; ++i ) {
    b [i] = 1;
    if ( i ) b [i - 1] = 0;
    apply ( b, x );
    for ( int j = 0; j < this->_numRows; ++j )
      m.set ( j, i, x [j] );
  }

  return m;
}

template <class DataType>
void LUInverse<DataType>::set ( int, int, DataType ) {
  throw UnimplementedCodeException ( "Cannot write to LU Inverse", __FILE__, __LINE__ );
}

template <class DataType>
void LUInverse<DataType>::add ( int, int, DataType ) {
  throw UnimplementedCodeException ( "Cannot write to LU Inverse", __FILE__, __LINE__ );
}


template <class DataType>
void LUInverse<DataType>::saveToFiles ( const char* const FileNameMask ) const {
  const aol::Vector<DataType> &decompVec = this->getDecompositionMatrixRef().getDataVectorReference();
  const aol::Vector<int> &permutVec = this->getPermutationMatrixRef().getPermutationVectorRef();
  char filename[1024];
  sprintf ( filename, FileNameMask, "decomp" );
  decompVec.saveToFile ( filename );
  sprintf ( filename, FileNameMask, "permut" );
  permutVec.saveToFile ( filename );
}


template <class DataType>
void LUInverse<DataType>::loadFromFiles ( const char* const FileNameMask ) {
  aol::Vector<DataType> decompVec;
  aol::Vector<int> permutVec;

  char filename[1024];
  sprintf ( filename, FileNameMask, "decomp" );
  decompVec.loadFromFile ( filename );
  sprintf ( filename, FileNameMask, "permut" );
  permutVec.loadFromFile ( filename );

  const int n = permutVec.size();

  aol::FullMatrix<DataType> decomp ( decompVec, n, n );
  aol::PermutationMatrix<DataType> permut ( permutVec );

  this->matrixInit ( n, n );
  this->_decomposition.reallocate ( n, n );
  this->_decomposition = decomp;
  this->_permut.reallocate ( n );
  this->_permut = permut;
}


template class LUInverse<float>;
template class LUInverse<double>;
template class LUInverse<long double>;



template <class DataType>
QRInverse<DataType>::QRInverse ()
  : Matrix<DataType> (), _decomposition (), _lengths () {}

//--- QRInverse copy constructor --------------------------------------------

template <class DataType>
QRInverse<DataType>::QRInverse ( const QRInverse<DataType>& inv )
    : Matrix<DataType> ( inv.getNumRows (), inv.getNumCols () ),
    _decomposition ( inv._decomposition ),
    _lengths ( inv._lengths ) {}

//--- QRInverse assignment operator -----------------------------------------

template <class DataType>
QRInverse<DataType>& QRInverse<DataType>::operator = ( const QRInverse<DataType>& inv ) {
  // Beware of self-assignment
  if ( this == &inv ) return *this;

  _decomposition = inv._decomposition;
  _lengths = inv._lengths;

  return *this;
}

//--- QRInverse constructor from matrix -------------------------------------

template <class DataType>
QRInverse<DataType>::QRInverse ( const Matrix<DataType>& m )
    : Matrix<DataType> ( m.getNumCols (), m.getNumRows () ),
    _decomposition ( m ),
    _lengths ( m.getNumCols () ) {
  if ( m.getNumCols () > m.getNumRows () )
    throw OutOfBoundsException ( "Matrix must be square or higher "
                                 "for QR inversion", __FILE__, __LINE__ );

  int r = m.getNumRows (), c = m.getNumCols ();

  for ( int i = 0; i < c; ++i ) {
    // Consider entries [i;r) of column i
    Vector<DataType> v ( r - i );
    _decomposition.getSubColumn ( i, i, v );

    // Compute vector normal to reflection hyperplane, eliminating entries (i;n)
    DataType alpha;
    computeReflection ( v, alpha );

    // Reflect all columns
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for ( int j = i + 1; j < c; ++j ) { // Result in current column is clear

      doReflectionColumn ( _decomposition, j, v, alpha );
    }

    // Store reflection vector
    writeReflection ( _decomposition, _lengths, i, v, alpha );
  }
}

//--- QRInverse::apply() ----------------------------------------------------

template <class DataType>
void QRInverse<DataType>::apply ( const Vector<DataType> &arg, Vector<DataType> &dest ) const {
  int c = _decomposition.getNumCols ();

  Vector<DataType> rhs ( arg ), v;
  DataType alpha;

  // Reflect right hand side like matrix
  for ( int i = 0; i < c; ++i ) {

    // Reflect subblock
    readReflection ( _decomposition, _lengths, i, v, alpha );
    doReflectionVector ( rhs, v, alpha );
  }

  // Solve triangular system, diagonal is in _lengths
  for ( int i = c - 1; i >= 0; --i ) {

    dest [i] = rhs [i];

    for ( int j = i + 1; j < c; ++j )
      dest [i] -= dest [j] * _decomposition.get ( i, j );

    dest [i] /= _lengths [i];
  }
}

//--- QRInverse::applyAdd() -------------------------------------------------

template <class DataType>
void QRInverse<DataType>::applyAdd ( const Vector<DataType> &arg, Vector<DataType> &dest ) const {
  Vector<DataType> temp ( dest );
  apply ( arg, dest );
  dest += temp;
}

//--- QRInverse::get() ------------------------------------------------------

template <class DataType>
DataType QRInverse<DataType>::get ( int row, int col ) const {
    // Use getFM if all entries are needed
    cerr << "Warning: QRInverse<DataType>::get is slow." << endl;

    Vector<DataType> colv ( this->_numCols );
    Vector<DataType> resv ( this->_numRows );
    colv.setZero ();
    colv [col] = 1.0;
    apply ( colv, resv );
    return resv [row];
  }

//--- QRInverse::getFM() ----------------------------------------------------

template <class DataType>
FullMatrix<DataType> QRInverse<DataType>::getFM () const {
  // Apply to columns of unit matrix
  FullMatrix<DataType> m ( this->_numRows, this->_numCols );
  Vector<DataType> x ( this->_numRows ), b ( this->_numCols );

  b.setZero ();
  for ( int i = 0; i < this->_numCols; ++i ) {
    b [i] = 1;
    if ( i ) b [i - 1] = 0;
    apply ( b, x );
    for ( int j = 0; j < this->_numRows; ++j )
      m.set ( j, i, x [j] );
  }

  return m;
}

//--- QRInverse::getR() -----------------------------------------------------

template <class DataType>
FullMatrix<DataType> QRInverse<DataType>::getR () const {
  FullMatrix<DataType> m ( this->_numCols, this->_numRows );
  m.setZero ();
  for ( int i = 0; i < this->_numRows; ++i )
    m.set ( i, i, _lengths [i] );

  for ( int i = 0; i < this->_numCols; ++i )
    for ( int j = i + 1; j < this->_numRows; ++j )
      m.set ( i, j, _decomposition.get ( i, j ) );

  return m;
}

//--- QRInverse::getQ() -----------------------------------------------------

template <class DataType>
FullMatrix<DataType> QRInverse<DataType>::getQ () const {
  FullMatrix<DataType> m ( this->_numCols, this->_numCols );

  // Applies all reflections successively to the unit matrix
  // Must start from end since inverse is wanted
  m.setIdentity ();
  Vector<DataType> v ( this->_numCols );
  DataType alpha;

  // For all reflections
  for ( int i = this->_numCols - 2; i >= 0; --i ) {

    readReflection ( _decomposition, _lengths, i, v, alpha );

    // For all columns of m where the reflection does something
    for ( int j = i; j < this->_numCols; ++j )
      doReflectionColumn ( m, j, v, alpha );
  }
  return m;
}

//--- QRInverse::set() ------------------------------------------------------

template <class DataType>
void QRInverse<DataType>::set ( int, int, DataType ) {
  throw UnimplementedCodeException ( "Cannot write to QR Inverse", __FILE__, __LINE__ );
}

//--- QRInverse::add() ------------------------------------------------------

template <class DataType>
void QRInverse<DataType>::add ( int, int, DataType ) {
  throw UnimplementedCodeException ( "Cannot write to QR Inverse", __FILE__, __LINE__ );
}


template class QRInverse<float>;
template class QRInverse<double>;
template class QRInverse<long double>;



//--- QRInversePivot constructor from matrix --------------------------------

template <class DataType>
QRInversePivot<DataType>::QRInversePivot ( const Matrix<DataType>& m )
    : Matrix<DataType> ( m.getNumCols (), m.getNumRows () )
    , _decomposition ( m )
    , _permMat ( m.getNumCols () )
    , _lengths ( m.getNumCols () )
    , _rank ( m.getNumCols () )
    , _m ( m ) {
  if ( m.getNumCols () > m.getNumRows () )
    throw OutOfBoundsException ( "Matrix must be square or higher "
                                 "for QR inversion", __FILE__, __LINE__ );

  int r = m.getNumRows (), c = m.getNumCols ();

  int i;
  for ( i = 0; i < c; ++i ) {
    // Consider entries [i;r) of column i
    Vector<DataType> v ( r - i );
    pivote ( i );

    _decomposition.getSubColumn ( i, i, v );

    // Compute vector normal to reflection hyperplane, eliminating entries (i;n)
    DataType alpha;
    computeReflection ( v, alpha );

    // Reflect all columns
    for ( int j = i + 1; j < c; ++j ) { // Result in current column is clear

      doReflectionColumn ( _decomposition, j, v, alpha );
    }

    // Store reflection vector
    writeReflection ( _decomposition, _lengths, i, v, alpha );

    if ( testForSingularity ( i ) )
      break;
  }
  if ( i < c - 1 ) {
    cout << "Rank is at most " << i << ", stopping QR." << endl;
    _rank = i;
  }
}

//--- QRInversePivot::pivote() ----------------------------------------------

template <class DataType>
bool QRInversePivot<DataType>::testForSingularity ( int n ) {
  typedef typename RealTrait<DataType>::RealType RealType;

  // compute subcondition number
  RealType sc = _decomposition.get ( 1, 1 ) / _decomposition.get ( n, n );
  bool isSingular = ( fabs ( sc ) > 1E13 );
  if ( isSingular )
    for ( int i = n + 1; i < _decomposition.getNumRows(); ++i )
      for ( int j = n + 1; j < _decomposition.getNumCols(); ++j )
        _decomposition.set ( i, j, ZOTrait<DataType>::zero );
  return isSingular;
}

//--- QRInversePivot::pivote() ----------------------------------------------

template <class DataType>
void QRInversePivot<DataType>::pivote ( int i ) {
  int index = getSubcolWithMaxNorm ( i );
  if ( i != index ) {
    _permMat.makeSwap ( i, index );
    swapCols ( i, index );
  }
}

//--- QRInversePivot::getSubcolWithMaxNorm() --------------------------------

template <class DataType>
int QRInversePivot<DataType>::getSubcolWithMaxNorm ( int startRowCol ) const {
  typename RealTrait<DataType>::RealType maxNorm = getSubcolNorm ( startRowCol, startRowCol );
  int maxNormIndex = startRowCol;
  for ( int i = startRowCol + 1; i < _decomposition.getNumCols(); ++i )
    if ( getSubcolNorm ( startRowCol, i ) > maxNorm ) {
      maxNormIndex = i;
      maxNorm = getSubcolNorm ( startRowCol, i );
    }
  return maxNormIndex;
}

//--- QRInversePivot::getSubcolNorm() ---------------------------------------

template <class DataType>
typename RealTrait<DataType>::RealType QRInversePivot<DataType>::getSubcolNorm
( int startRow, int col ) const {
  int n = _decomposition.getNumRows() - startRow;
  Vector<DataType> v ( n );
  _decomposition.getSubColumn ( startRow, col, v );
  return v.norm();
}

//--- QRInversePivot::swapRows() --------------------------------------------

template <class DataType>
void QRInversePivot<DataType>::swapCols ( int i, int j ) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for ( int row = 0; row < _decomposition.getNumRows(); ++row )
    swap ( _decomposition.ref ( row, i ), _decomposition.ref ( row, j ) );
}

//--- QRInversePivot::apply() -----------------------------------------------

template <class DataType>
void QRInversePivot<DataType>::applyWithoutPostIteration ( const Vector<DataType> &Arg, Vector<DataType> &Dest ) const {
  int c = _decomposition.getNumCols ();

  Vector<DataType> rhs ( Arg ), v, dest ( Dest, STRUCT_COPY );
  DataType alpha;

  // Reflect right hand side like matrix
  for ( int i = 0; i < c; ++i ) {

    // Reflect subblock
    readReflection ( _decomposition, _lengths, i, v, alpha );
    doReflectionVector ( rhs, v, alpha );
  }

  // Solve triangular system, diagonal is in _lengths
  for ( int i = _rank - 1; i >= 0; --i ) {

    dest [i] = rhs [i];

    for ( int j = i + 1; j < c; ++j )
      dest [i] -= dest [j] * _decomposition.get ( i, j );

    dest [i] /= _lengths [i];
  }
  _permMat.applyTransposed ( dest, Dest );
}

template <class DataType>
void QRInversePivot<DataType>::apply ( const Vector<DataType> &Arg, Vector<DataType> &Dest ) const {

  applyWithoutPostIteration ( Arg, Dest );
  Vector<DataType> res ( Arg, STRUCT_COPY );

  // residual computation
  _m.apply ( Dest, res );
  res -= Arg;
  // cout << "Residual^2 before post iteration: " << res.normSqr() << endl;

  Vector<DataType> correction ( Dest, STRUCT_COPY );
  applyWithoutPostIteration ( res, correction );
  Dest -= correction;

  // residual computation again
  _m.apply ( Dest, res );
  res -= Arg;
  // cout << "Residual^2 after post iteration:  " << res.normSqr() << endl;
}

//--- QRInversePivot::applyAdd() --------------------------------------------

template <class DataType>
void QRInversePivot<DataType>::applyAdd ( const Vector<DataType> &Arg, Vector<DataType> &Dest ) const {
  Vector<DataType> temp ( Dest );
  apply ( Arg, Dest );
  Dest += temp;
}

//--- QRInversePivot::get() -------------------------------------------------

template <class DataType>
DataType QRInversePivot<DataType>::get ( int row, int col ) const {
    // Use getFM if all entries are needed
    cerr << "Warning: QRInverse<DataType>::get is slow." << endl;

    Vector<DataType> colv ( this->_numCols );
    Vector<DataType> resv ( this->_numRows );
    colv.setZero ();
    colv [col] = 1.0;
    apply ( colv, resv );
    return resv [row];
  }

//--- QRInversePivot::set() -------------------------------------------------

template <class DataType>
void QRInversePivot<DataType>::set ( int, int, DataType ) {
  throw UnimplementedCodeException ( "Cannot write to QR Inverse", __FILE__, __LINE__ );
}

//--- QRInversePivot::add() -------------------------------------------------

template <class DataType>
void QRInversePivot<DataType>::add ( int, int, DataType ) {
  throw UnimplementedCodeException ( "Cannot write to QR Inverse", __FILE__, __LINE__ );
}


template class QRInversePivot<float>;
template class QRInversePivot<double>;
template class QRInversePivot<long double>;

}
