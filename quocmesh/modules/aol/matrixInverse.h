#ifndef __MATRIXINVERSE_H
#define __MATRIXINVERSE_H

#include <matrix.h>
#include <bandMatrix.h>
#include <progressBar.h>

namespace aol {

/****************************************************************************
 *
 *  \brief LU decomposition of a matrix
 *  \ingroup directSolver
 *  \author Lenz
 *
 */
template <class DataType>
class LUInverse : public Matrix<DataType> {

protected:

  //! Default constructor not allowed
  LUInverse ();

public:

  //! Copy constructor
  explicit LUInverse ( const LUInverse<DataType>& inv );

  //! Assignment operator
  LUInverse<DataType>& operator = ( const LUInverse<DataType>& inv );

  //! Destructor
  virtual ~LUInverse ();

  //! Construct from matrix via LU decomposition
  explicit LUInverse ( const Matrix<DataType>& m );

  //! construct from given size, decomposition and permutation
  LUInverse ( const int N, const aol::FullMatrix<DataType> &Decomp, const aol::PermutationMatrix<DataType> &Permut );

  //! Apply solves the LU system
  void apply ( const Vector<DataType>& arg, Vector<DataType>& dest ) const;

  //! Add the solution to the previous value of the destination vector
  void applyAdd ( const Vector<DataType>& arg, Vector<DataType>& dest ) const;

  //! Compute entry in inverse matrix, slow
  DataType get ( int row, int col ) const;

  //! Compute direct representation of inverse, for output, debugging etc.
  FullMatrix<DataType> getFM () const;

  //! Write operations to single entries not permitted
  void set ( int, int, DataType );

  //! Write operations to single entries not permitted
  void add ( int, int, DataType );

  const aol::FullMatrix<DataType>& getDecompositionMatrixRef() const {
    return ( _decomposition );
  }

  const aol::PermutationMatrix<DataType>& getPermutationMatrixRef() const {
    return ( _permut );
  }

  //! save decomposition and permutation matrix to two files (FileNameMask should contain %s for decomp/permut)
  void saveToFiles ( const char* const FileNameMask ) const;

  //! load decomposition and permutation matrix from files (FileNameMask should contain %s for decomp/permut)
  void loadFromFiles ( const char* const FileNameMask );

private:

  //! Stores the LU decomposition
  FullMatrix<DataType> _decomposition;

  //! Stores the permutation matrix
  PermutationMatrix<DataType> _permut;

  void setZero() {
    throw aol::UnimplementedCodeException ( "aol::LUInverse::setZero not implemented yet.", __FILE__, __LINE__ );
  }
};


/****************************************************************************
 *
 *  \brief QR decomposition of a matrix via Householder's method
 *
 *  \author Lenz
 *
 */

template <class DataType>
class QRInverse : public Matrix<DataType> {

protected:

  //! Default constructor not allowed
  QRInverse ();

public:

  //! Copy constructor
  explicit QRInverse ( const QRInverse<DataType>& inv );

  //! Assignment operator
  QRInverse<DataType>& operator = ( const QRInverse<DataType>& inv );

  //! Construct from matrix via QR decomposition
  explicit QRInverse ( const Matrix<DataType>& m );

  //! Apply solves the QR system
  void apply ( const Vector<DataType>& arg, Vector<DataType>& dest ) const;

  //! Add the solution to the previous value of the destination vector
  void applyAdd ( const Vector<DataType>& arg, Vector<DataType>& dest ) const;

  //! Compute entry in inverse matrix, slow
  DataType get ( int row, int col ) const;

  //! Compute direct representation of inverse, for output, debugging etc.
  FullMatrix<DataType> getFM () const;

  //! Compute direct representation of reflection
  FullMatrix<DataType> getQ () const;

  //! Show upper triangluar part of decomposition
  FullMatrix<DataType> getR () const;

private:
  //! Write operations to single entries not permitted
  void set ( int, int, DataType );
  //! Write operations to single entries not permitted
  void add ( int, int, DataType );

  //! Stores the QR decomposition
  FullMatrix<DataType> _decomposition;

  //! Stores lengths of reflected vectors
  Vector <DataType> _lengths;

  void setZero() {
    throw aol::UnimplementedCodeException ( "aol::QRInverse::setZero() "
                                            "does not make sense for this class.",
                                            __FILE__, __LINE__ );
  }
};

/****************************************************************************
 *
 *  \brief QR decomposition of a matrix via Householder's method
 *            with column pivoting.
 *
 *  For details on the use of pivoting for QR decomposition and the
 *  role of subcondition numbers cf. Deuflhard/Hohmann: "Num. Math. I",
 *  p. 82 sqq.
 *
 *  \author von Deylen,
 *          Created 2008-04-10
 *          Last change 2008-04-10
 *
 */

template <class DataType>
class QRInversePivot : public Matrix<DataType> {

protected:

  //! Default constructor not allowed
  QRInversePivot ();

public:

  //! Construct from matrix via QR decomposition
  QRInversePivot ( const Matrix<DataType>& m );

  //! Apply solves the QR system
  void apply ( const Vector<DataType>& Arg, Vector<DataType>& Dest ) const;

  //! Add the solution to the previous value of the destination vector
  void applyAdd ( const Vector<DataType>& Arg, Vector<DataType>& Dest ) const;

  //! Compute entry in inverse matrix, slow
  DataType get ( int row, int col ) const;

protected:
  void applyWithoutPostIteration ( const Vector<DataType>& Arg, Vector<DataType>& Dest ) const;
  bool testForSingularity ( int n );
  void pivote ( int startRowCol );
  int getSubcolWithMaxNorm ( int startRowCol ) const;
  typename RealTrait<DataType>::RealType getSubcolNorm ( int startRow, int col ) const;
  void swapCols ( int i, int j );

private:
  //! Write operations to single entries not permitted
  void set ( int, int, DataType );
  //! Write operations to single entries not permitted
  void add ( int, int, DataType );

  //! Stores the QR decomposition
  FullMatrix<DataType> _decomposition;
  PermutationMatrix<DataType> _permMat;

  //! Stores lengths of reflected vectors
  Vector <DataType> _lengths;

  int _rank;
  const Matrix<DataType> & _m;

  void setZero() {
    throw aol::UnimplementedCodeException ( "aol::QRInverse::setZero() "
                                            "does not make sense for this class.",
                                            __FILE__, __LINE__ );
  }
};

/** **************************************************************************
 *
 *  \brief one single Givens rotation
 *
 *  \author Lenz
 *
 */

template <class DataType>
class GivensRotation {

public:

  GivensRotation () {}

  GivensRotation ( DataType eliminator, DataType eliminee ) {

    if ( fabs ( eliminee ) > fabs ( eliminator ) ) {

      DataType t = eliminator / eliminee;
      s = 1. / sqrt ( 1. + t * t );
      c = s * t;
    } else {

      DataType t = eliminee / eliminator;
      c = 1. / sqrt ( 1. + t * t );
      s = c * t;
    }
  }

  void rotate ( DataType& a, DataType& b ) const {

    DataType aa = c * a + s * b;
    b = - s * a + c * b;
    a = aa;
  }

private:

  DataType c, s;
};

/****************************************************************************
 *
 *  \brief list of Givens rotations needed for elimination
 *
 *  \author Lenz
 *
 */

template <class DataType>
class GivensRotations : public vector<GivensRotation<DataType> > {

public:
  GivensRotations ( int n ) : vector<GivensRotation<DataType> > ( n - 1 ) {}};

/** **************************************************************************
 *
 *  \brief QR tridiagonalization of a matrix via Givens rotations
 *         Constructed with a symmetric full matrix A, it computes an
 *         orthogonal matrix Q such that Q^T A Q is tridiagonal.
 *
 *  Bringing symmetric A to tridiagonal form via Householder reflections
 *  or Givens rotations is a 4/3n^3 process.
 *
 *  \author Lenz, von Deylen
 *
 */

template <class MatrixType>
class QRGivensTridiag {

protected:

  //! Default constructor not allowed
  QRGivensTridiag ();

public:
  typedef typename MatrixType::DataType DataType;

  //! Copy constructor
  QRGivensTridiag ( const QRGivensTridiag<MatrixType>& inv );

  //! Assignment operator
  QRGivensTridiag<MatrixType>& operator = ( const QRGivensTridiag<MatrixType>& inv );

  //! Construct for matrix
  QRGivensTridiag ( const MatrixType & m );

  //! Compute direct representation of reflection
  //! for input \f$ M \f$ and tridiagonalization \f$ D \f$, this
  //! functions outputs Q such that \f$ Q^T M Q = D \f$.
  void getQ ( MatrixType & Q ) const;

  //! Export tridiagonalization
  void getTridiagFull ( MatrixType & ret ) const;
  void getTridiagBand ( TriBandMatrix<DataType> & ret ) const;

private:

  //! Stores the Q part of the decomposition
  MatrixType _decomposition;

  //! Stores the three diagonals
  //! Indexed as follows: first index: 1 is diagonal, second index: Column number
  //! _diags [3] is only needed temporarily during the Givens decomposition
  Vector<DataType> _diags [4];

  //! Stores the concatenated Givens rotations
  GivensRotations<DataType> _rotations;
};


/*** FUNCTION DECLARATIONS **************************************************
 *
 *
 */

//! Computes reflection paramters
//! \param normal initially holds the vector to be reflected onto a
//!        multiple of the first standard basis vector, upon return
//!        it holds the direction normal to the reflection hyperplane
template <class DataType>
void computeReflection ( Vector<DataType>& normal, DataType& alpha );

//! Reads reflection from matrix
//! @param offset number of entries below diagonal where stored normal starts
template <class MatrixType>
void readReflection ( MatrixType & decomp,
                      const Vector<typename MatrixType::DataType>& lengths,
                      int i, Vector<typename MatrixType::DataType>& normal,
                      typename MatrixType::DataType& alpha, int offset = 0 );

//! Stores reflection in matrix
template <class MatrixType>
void writeReflection ( MatrixType& decomp,
                       Vector<typename MatrixType::DataType>& lengths,
                       int i, const Vector<typename MatrixType::DataType>& normal,
                       typename MatrixType::DataType alpha );

//! Reflect one vector according to paramters
template <class DataType>
void doReflection ( Vector<DataType>& arg,
                    const Vector<DataType>& normal,
                    DataType alpha );

//! Reflect the end part of one vector according to paramters
template <class DataType>
void doReflectionVector ( Vector<DataType>& arg,
                          const Vector<DataType>& normal,
                          DataType alpha );

//! Reflect the end part of one matrix column according to paramters
template <class MatrixType>
void doReflectionColumn ( MatrixType& arg, int c,
                          const Vector<typename MatrixType::DataType>& normal,
                          typename MatrixType::DataType alpha );

//! Reflect the end part of one matrix row according to paramters, mutiplying
//! from the right with the reflection matrix
template <class MatrixType>
void doReflectionRow ( MatrixType & arg, int r,
                       const Vector<typename MatrixType::DataType>& normal,
                       typename MatrixType::DataType alpha );

//! Rotate subdiagonal
template <class DataType>
void computeRotationSubdiagonal ( GivensRotations<DataType>& rotations,
                                  Vector<DataType> diags [], int col );

//! Rotate subdiagonal
template <class DataType>
void doRotationSubdiagonal ( GivensRotations<DataType>& rotations,
                             Vector<DataType> diags [],
                             int col, int offset );

/*************************************************************************
 *
 *    function implementation
 *
 ************************************************************************* */

//--- computeReflection -----------------------------------------------------

template <class DataType>
inline void computeReflection ( Vector<DataType>& normal, DataType& alpha ) {
  alpha = normal.norm () * signum1at0 ( normal [0] );
  normal [0] -= alpha;
}

//--- readReflection --------------------------------------------------------

template <class MatrixType>
void readReflection ( MatrixType & decomp,
                      const Vector<typename MatrixType::DataType>& lengths,
                      int i, Vector<typename MatrixType::DataType>& normal,
                      typename MatrixType::DataType& alpha, int offset ) {
  normal.resize ( decomp.getNumRows () - i - offset );
  decomp.getSubColumn ( i + offset, i, normal );
  alpha = lengths [i];
}

//--- writeReflection -------------------------------------------------------

template <class MatrixType>
void writeReflection ( MatrixType & decomp,
                       Vector<typename MatrixType::DataType>& lengths,
                       int i, Vector<typename MatrixType::DataType>& normal,
                       typename MatrixType::DataType& alpha ) {
  int p = decomp.getNumRows () - normal.size ();
  decomp.setSubColumn ( p, i, normal );
  lengths [i] = alpha;
}

//---------------------------------------------------------------------------

template <class DataType>
inline void doReflection ( Vector<DataType>& arg,
                           const Vector<DataType>& normal,
                           DataType alpha ) {
  const DataType denom = alpha * normal [0];
  if ( !denom ) return; // Nothing to be done

  arg.addMultiple ( normal, normal * arg / denom );
}

//--- doReflectionVector ----------------------------------------------------

template <class DataType>
inline void doReflectionVector ( Vector<DataType>& arg,
                                 const Vector<DataType>& normal,
                                 DataType alpha ) {
  int l = normal.size (), p = arg.size () - l;
  Vector<DataType> x ( l );
  arg.getBlock ( p, x );
  doReflection ( x, normal, alpha );
  arg.setBlock ( p, x );
}

//--- doReflectionColumn ----------------------------------------------------

template <class MatrixType>
void doReflectionColumn ( MatrixType & arg, int c,
                          const Vector<typename MatrixType::DataType>& normal,
                          typename MatrixType::DataType alpha ) {
  int r = arg.getNumRows (),
          l = normal.size (),
              p = r - l;
  Vector<typename MatrixType::DataType> x ( r );
  x.resize ( l ); // Optimize for vector manager
  arg.getSubColumn ( p, c, x );
  doReflection ( x, normal, alpha );
  arg.setSubColumn ( p, c, x );
}

//--- doReflectionRow -------------------------------------------------------

template <class MatrixType>
void doReflectionRow ( MatrixType & arg, int r,
                       const Vector<typename MatrixType::DataType>& normal,
                       typename MatrixType::DataType alpha ) {
  int l = normal.size (),
          p = arg.getNumCols () - l;
  Vector<typename MatrixType::DataType> x ( l );
  arg.getSubRow ( r, p, x );
  doReflection ( x, normal, alpha );
  arg.setSubRow ( r, p, x );
}

//--- computeRotationSubdiagonal --------------------------------------------

template <class DataType>
inline void computeRotationSubdiagonal ( GivensRotations<DataType>& /*rotations*/,
                                         int /*col*/,
                                         DataType dval, DataType sdval,
                                         DataType & s, DataType & c ) {
  double t = dval / sdval;
  s = 1 / sqrt ( 1 + t * t );
  c = s * t;
}

//--- doRotationSubdiagonal -------------------------------------------------

template <class DataType>
inline void doRotationSubdiagonal ( GivensRotations<DataType>& rotations,
                                    Vector<DataType> diags [], int col ) {
  DataType s, c;
  computeRotationSubdiagonal ( rotations, col, diags [1][col], diags [0][col], s, c );

  // Eliminate subdiagonal
  diags [1][col] = c * diags [1][col] + s * diags [0][col];
  diags [0][col] = 0;

  if ( col + 1 == diags [0].size () ) return;

  // Apply to next column
  diags [2][col+1] = c * diags [2][col+1] + s * diags [1][col+1];
  diags [1][col+1] = c * diags [1][col+1] - s * diags [2][col+1];

  if ( col + 2 == diags [0].size () ) return;

  // Apply to second next column
  // diags [3][col+2] is zero before rotation
  diags [3][col+2] = s * diags [2][col+2];
  diags [2][col+2] = c * diags [2][col+2];

}



/*************************************************************************
 *
 *    QRGivensTridiag implementation
 *
 ************************************************************************* */

//--- QRGivensTridiag constructor -------------------------------------------

template <class MatrixType>
QRGivensTridiag<MatrixType>::QRGivensTridiag ()
    : _decomposition (), _rotations () {}

//--- QRGivensTridiag copy constructor --------------------------------------

template <class MatrixType>
QRGivensTridiag<MatrixType>::QRGivensTridiag ( const QRGivensTridiag<MatrixType>& it )
    : _decomposition ( it._decomposition ), _rotations ( it._rotations ) {
  for ( int i = 0; i < 4; ++i )
    _diags [i] = it._diags [i];
}

//--- QRGivensTridiag assignment operator -----------------------------------

template <class MatrixType>
QRGivensTridiag<MatrixType>& QRGivensTridiag<MatrixType>::operator= ( const QRGivensTridiag<MatrixType>& it ) {
  // Beware of self-assignment
  if ( this == &it ) return *this;

  _decomposition = it._decomposition;
  _rotations = it._rotations;

  for ( int i = 0; i < 4; ++i )
    _diags [i] = it._diags [i];

  return *this;
}

//--- QRGivensTridiag constructor with matrix -------------------------------

template <class MatrixType>
QRGivensTridiag<MatrixType>::QRGivensTridiag ( const MatrixType & m )
    : _decomposition ( m ), _rotations ( m.getNumRows () ) {

  // Preparations
  for ( int i = 0; i < 4; ++i ) {
    _diags [i].resize ( m.getNumRows () );
    _diags [i].setZero ();
  }

  // Tridiagonalize

  int n = m.getNumCols ();

  //  ProgressBar<false, true, int, 50> progressBar ( "matrix tridiagonalization " ); // does not seem to work
  ProgressBar<> progressBar ( "matrix tridiagonalization " );
  progressBar.start ( n - 2 );
  for ( int i = 0; i < n - 2; ++i ) {
    progressBar++;

    // Consider entries [i+1;n) of column i
    Vector<DataType> v ( n - i - 1 );
    _decomposition.getSubColumn ( i + 1, i, v );

    // Compute vector normal to reflection hyperplane, eliminating entries (i;n)
    DataType alpha;
    computeReflection ( v, alpha );

    // Reflect all columns
    for ( int j = i + 1; j < n; ++j )
      doReflectionColumn ( _decomposition, j, v, alpha );

    // Reflect all rows
    for ( int j = i; j < n; ++j )
      doReflectionRow ( _decomposition, j, v, alpha );

    // Store reflection vector
    writeReflection ( _decomposition, _diags [0], i, v, alpha );
  }
  progressBar.finish();

  // Last column is not eliminated, copy value for convenience
  _diags [0].set ( n - 2, _decomposition.get ( n - 1, n - 2 ) );

  // Copy diagonals
  for ( int i = 0; i < n; ++i ) {
    _diags [1][i] = _decomposition.get ( i, i );
    _decomposition.set ( i, i, 0. );
  }
  for ( int i = 0; i < n - 1; ++i ) {
    _diags [2][i+1] = _decomposition.get ( i, i + 1 );
    _decomposition.set ( i, i + 1, 0. );
  }
}

//--- QRGivensTridiag::getQ() -----------------------------------------------

template <class MatrixType>
void QRGivensTridiag<MatrixType>::getQ ( MatrixType & Q ) const {
  int n = _decomposition.getNumRows ();

  // Applies all reflections successively to the unit matrix
  // Must start from end since inverse is wanted
  Q.setIdentity ();
  Vector<DataType> v ( n );
  DataType alpha;

  // For all reflections
  for ( int i = n - 3; i >= 0; --i ) {

    readReflection ( _decomposition, _diags [0], i, v, alpha, 1 );

    // For all columns of m where the reflection does something
    for ( int j = i; j < n; ++j )
      doReflectionColumn ( Q, j, v, alpha );
  }
}

//--- QRGivensTridiag::getTridiagFull() -------------------------------------

template <class MatrixType>
void QRGivensTridiag<MatrixType>::getTridiagFull ( MatrixType & ret ) const {
  int n = _decomposition.getNumRows ();

  ret.setZero ();

  for ( int i = 0; i < n; ++i ) {
    if ( i < n - 1 ) ret.set ( i + 1, i, _diags [0][i] );
    ret.set ( i, i, _diags [1][i] );
    if ( i ) ret.set ( i - 1, i, _diags [2][i] );
  }
}

//--- QRGivensTridiag::getTridiagBand() -------------------------------------

template <class MatrixType>
void QRGivensTridiag<MatrixType>::getTridiagBand ( TriBandMatrix<typename MatrixType::DataType> & ret ) const {
  int n = _decomposition.getNumRows ();

  for ( int i = 0; i < n; ++i ) {
    if ( i < n - 1 ) ret.set ( i + 1, i, _diags [0][i] );
    ret.set ( i, i, _diags [1][i] );
    if ( i ) ret.set ( i - 1, i, _diags [2][i] );
  }
}

//---------------------------------------------------------------------------


} // end of namespace aol.

#endif
