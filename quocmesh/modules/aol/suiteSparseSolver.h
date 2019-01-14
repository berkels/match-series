#ifndef __SUITESPARSESOLVER_H
#define __SUITESPARSESOLVER_H

#include <solver.h>

#include <sparseMatrices.h>

#ifdef USE_EXTERNAL_SUITESPARSE

#include <cholmod.h>
#include <umfpack.h>
#include <SuiteSparseQR.hpp>

namespace aol {

/**
* \brief a direct solver for matrices which are symmetric and positive definite
* \note to compile this class the external cholmod is necessary
*
* \author Wirth
* \ingroup directSolver
*/
template <typename RealType, typename MatrixType>
class CholeskyBlockInverseOp :
  public aol::InverseOp<aol::MultiVector<RealType> > {

private:
  mutable cholmod_common _settings;
  cholmod_sparse* _matrix;
  cholmod_factor* _lowerDiagonalFactor;

  int _maxNumEntriesPerBlockRow;  //!< The max. number of non-zero entries per row in a block.
  bool _matrixPositiveDefinite;

  inline void factorize( cholmod_triplet* triplet, const int numNonZeros )  {
    cholmod_free_sparse( &_matrix, &_settings );
    _matrix = cholmod_triplet_to_sparse( triplet, numNonZeros, &_settings );
#ifdef VERBOSE
    _settings.print = 5;
    cholmod_check_triplet( triplet, &_settings );
    cerr << endl << "Check status " << _settings.status << endl;
    cholmod_print_triplet( triplet, "", &_settings );
    cerr << "Print status " << _settings.status << endl;
    cholmod_print_sparse( _matrix, "", &_settings );
#endif

    if ( _lowerDiagonalFactor == NULL )
      _lowerDiagonalFactor = cholmod_analyze( _matrix, &_settings );
    cholmod_factorize( _matrix, _lowerDiagonalFactor, &_settings );

    
    switch (_settings.status)
    {
      case 0:
        break;
      case 1:
        cerr << "CholeskyBlockInverseOp::factorize: Warning: CHOLMOD_NOT_POSDEF (1)" << endl;
        break;
      case 2:
        cerr << "CholeskyBlockInverseOp::factorize: Warning: CHOLMOD_DSMALL (2)" << endl;
        break;
      case -1:
        cerr << "CholeskyBlockInverseOp::factorize: Failure: CHOLMOD_NOT_INSTALLED (-1)" << endl;
        break;
      case -2:
        cerr << "CholeskyBlockInverseOp::factorize: Failure: CHOLMOD_OUT_OF_MEMORY (-2)" << endl;
        break;
      case -3:
        cerr << "CholeskyBlockInverseOp::factorize: Failure: CHOLMOD_TOO_LARGE (-3)" << endl;
        break;
      case -4:
        cerr << "CholeskyBlockInverseOp::factorize: Failure: CHOLMOD_INVALID (-4)" << endl;
        break;
      case -5:
        cerr << "CholeskyBlockInverseOp::factorize: Failure: CHOLMOD_GPU_PROBLEM (-5)" << endl;
        break;
      default:
        cerr << "CholeskyBlockInverseOp::factorize: undefined status value (" << _settings.status << ")" << endl;
    }
    _matrixPositiveDefinite = ( _settings.status == CHOLMOD_OK );
    if(!_matrixPositiveDefinite)
      throw aol::Exception("CholeskyBlockInverseOp::factorize: Matrix is not positive definite!", __FILE__, __LINE__ );
  }

public:
  /*!
   * \brief Constructs the Op, initializes cholmod
   *
   * _maxNumEntriesPerBlockRow is only used in setMatrix and should be set in that call
   * the argument here is for compatibility only
   */
  CholeskyBlockInverseOp( const int MaxNumEntriesPerBlockRow = 0 ) :
    aol::InverseOp<aol::MultiVector<RealType> > ( ),
    _matrix( NULL ),
    _lowerDiagonalFactor( NULL ),
    _maxNumEntriesPerBlockRow( MaxNumEntriesPerBlockRow ),
    _matrixPositiveDefinite ( false ) {
    cholmod_start ( &_settings );
    if ( aol::RealTrait<RealType>::ALIAS == aol::FLOAT )
      _settings.dtype = CHOLMOD_DOUBLE; // CHOLMOD_SINGLE not yet supported
    else
      _settings.dtype = CHOLMOD_DOUBLE;
  }

  virtual ~CholeskyBlockInverseOp() {
    cholmod_free_factor( &_lowerDiagonalFactor, &_settings );
    cholmod_free_sparse( &_matrix, &_settings );
    cholmod_finish ( &_settings );
  }

  /*!
   * \brief Constructs the Op, assumes d by d blocks (d=dimension).
   * \param[in] MaxNumEntriesPerBlockRow The max. number of non-zero entries in a row of a block (taken over all blocks).
   *
   * The constructor takes the number of non-zero entries per row in a block, which can be set to \f$3^d\f$,
   * (\f$d = \f$ dimension) if the standard stencil in a uniform FE grid and standard FE operators are used.
   */
  inline void setMatrix ( const aol::BlockOpBase<RealType,MatrixType> &BlockMat, const int MaxNumEntriesPerBlockRow = 0  ) {

    if ( _maxNumEntriesPerBlockRow == 0 )  {
      if ( MaxNumEntriesPerBlockRow == 0 )
        throw aol::Exception( "You have to set MaxNumEntriesPerBlockRow in setMatrix call!", __FILE__, __LINE__ );
      else _maxNumEntriesPerBlockRow = MaxNumEntriesPerBlockRow;
    }

    // check whether all blocks are allocated
    for ( int i = 0; i < BlockMat.getNumRows(); i++ )
      for ( int j = 0; j < BlockMat.getNumCols(); j++ )
        if( !BlockMat.getPointer(i, j) )
          throw aol::Exception( "CholeskyBlockInverseOp::setMatrix: all blocks must be allocated!" );

    // compute number of rows (and columns)
    int n = 0;
    for ( int i = 0; i < BlockMat.getNumRows(); i++ )
      n += BlockMat.getReference( i, 0 ).getNumRows();

    // only save the entries below the diagonal (if all entries were stored, the matrix would be treated unsymmetric, yielding bad results)
    cholmod_triplet* triplet = cholmod_allocate_triplet( n, n, n * ( ( _maxNumEntriesPerBlockRow * BlockMat.getNumCols() ) / 2 + 1 ), 1, CHOLMOD_REAL, &_settings );
    int numberOfEntries = 0;
    for ( int row = 0, totalRow = 0; row < BlockMat.getNumRows(); row++ ) {
      for ( int localRow = 0; localRow < BlockMat.getReference( row, 0 ).getNumRows(); localRow++, totalRow++ ) {
        int colOffset = 0;
        for ( int col = 0; col <= row; col++ ) {
          vector<typename aol::Row<RealType>::RowEntry > matRow;
          BlockMat.getReference( row, col ).makeRowEntries( matRow, localRow );
          for ( typename vector<typename aol::Row<RealType>::RowEntry >::iterator it = matRow.begin(); it != matRow.end(); ++it )
            if ( col < row || it->col <= localRow ) {
              reinterpret_cast<int*>( triplet->i )[numberOfEntries] = totalRow;
              reinterpret_cast<int*>( triplet->j )[numberOfEntries] = colOffset + it->col;
              reinterpret_cast<double*>( triplet->x )[numberOfEntries] = it->value;
              numberOfEntries++;
            }
          colOffset += BlockMat.getReference( row, col ).getNumCols();
        }
      }
    }
    triplet->nnz = numberOfEntries;
    factorize( triplet, n * _maxNumEntriesPerBlockRow * BlockMat.getNumCols() );
    cholmod_free_triplet( &triplet, &_settings );
  }

  inline void setMatrix ( const aol::TripletMatrix<RealType> & tmat )  {
    cholmod_triplet ctriplet;
    ctriplet.nrow  = tmat.getNumRows();
    ctriplet.ncol  = tmat.getNumCols();
    ctriplet.nzmax = tmat.getValueReference().size();
    ctriplet.nnz   = tmat.getValueReference().size(); // better way to do this?
    ctriplet.i     = tmat.getRowIndexReference().getData();
    ctriplet.j     = tmat.getColIndexReference().getData();
    ctriplet.x     = tmat.getValueReference().getData();
    ctriplet.z     = NULL;
    ctriplet.stype = 1;  // only upper triangular part is stored
    ctriplet.itype = CHOLMOD_INT;
    ctriplet.xtype = CHOLMOD_REAL;
    ctriplet.dtype = CHOLMOD_DOUBLE;
    factorize( &ctriplet, ctriplet.nzmax ); // number of nonzeros is much smaller!
  }

  /*!
   * \brief setMatrix for matrices in CSCMatrix format
   * \author Toelkes
   * \param[in] mat The system matrix
   *
   * The upper triangle of mat will be copied and used for the Cholesky decomposition.
   * mat is assumed to be symmetric positive definite and to have either all entries
   * of the system matrix or just the upper triangle. In the first case, all entries
   * below the diagonal will be ignored.
   */
  inline void setMatrix ( const aol::CSCMatrix<RealType> &mat ) {
    const aol::Vector<int> &colPtr = mat.getColumnPointerReference ();
    const aol::Vector<int> &rowIdx = mat.getRowIndexReference ();
    const aol::Vector<RealType> &val = mat.getValueReference ();

    // Count entries in upper triangle of system matrix
    // (system matrix is assumed to be symmetric pos. def.
    int upperTriangEntries = 0;
    for ( int j = 0; j < mat.getNumCols (); ++j ) {
      for ( int idx = colPtr[j]; idx < colPtr[j + 1]; ++idx ) {
        if ( rowIdx[idx] <= j )
          ++upperTriangEntries;
      }
    }

    // Allocoate cholmod_sparse matrix
    _matrix = cholmod_allocate_sparse ( mat.getNumRows (), mat.getNumCols (),
        upperTriangEntries, 1, 1, 1, CHOLMOD_REAL, &_settings );

    // Copy upper triangle to _matrix
    static_cast<int*> ( _matrix->p )[0] = 0;
    int c = 0;
    for ( int j = 0; j < mat.getNumCols (); ++j ) {
      // Instead of destIdx, c could be used when the columns are ordered and no gaps
      // exists between them (which should both be true for CSCMatrix).
      // Handle the more general case here anyways (although almost nothing will work
      // if a non-standard format is used).
      int destIdx = static_cast<int*> ( _matrix->p )[j];
      for ( int idx = colPtr[j]; idx < colPtr[j + 1]; ++idx ) {
        if ( rowIdx[idx] > j )
          continue;

        static_cast<int*> ( _matrix->i )[destIdx] = rowIdx[idx];
        static_cast<RealType*> ( _matrix->x )[destIdx] = val[idx];

        ++destIdx;
        ++c;
      }

      // Set p to correct values (less values per column than in full matrix)
      static_cast<int*> ( _matrix->p )[j + 1] = c;
    }

#ifdef VERBOSE
    cerr << "Wrote " << c << " entries"  << endl;
    cerr << "Print status " << _settings.status << endl;
    _settings.print = 5;
    cholmod_print_sparse( _matrix, "", &_settings );
#endif

    // Factorize
    if ( _lowerDiagonalFactor == NULL )
      _lowerDiagonalFactor = cholmod_analyze ( _matrix, &_settings );
    else {
      // This should not happen (reset _lowerDiagonalFactor when resetting matrix)
      throw aol::Exception( "_lowerDiagonalFactor already exists", __FILE__, __LINE__ );
    }
    cholmod_factorize ( _matrix, _lowerDiagonalFactor, &_settings );
    _matrixPositiveDefinite = ( _settings.status == CHOLMOD_OK );

    if ( !_matrixPositiveDefinite ) {
      // It does not make sense to try using a cholesky decomposition
      // when the matrix is not positive definite.
      throw aol::Exception( "Matrix not positive definite", __FILE__, __LINE__ );
    }
  }

  // exports lower diagonal matrix in Matrix Market format
  void exportMatrix( const char* Filename ) const {
    FILE* matrix = fopen( Filename, "w+" );
    if ( matrix == NULL )
      cerr<<"error occured while opening the file"<<endl;
    cerr<<cholmod_write_sparse( matrix, _matrix, NULL, "", &_settings );
    fclose( matrix );
  }

  virtual void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::MultiVector<RealType> result ( Dest, aol::STRUCT_COPY );
    apply ( Arg, result );
    Dest += result;
  }

  virtual void apply ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    if ( !_matrixPositiveDefinite )
      Dest = Arg;
    else {
      cholmod_dense *rhs = cholmod_allocate_dense( Arg.getTotalSize(), 1, Arg.getTotalSize(), CHOLMOD_REAL, &_settings );
      for ( int i = 0, index = 0; i < Arg.numComponents(); i++ )
        for ( int j = 0; j < Arg[i].size(); j++, index++ )
          reinterpret_cast<double*>( rhs->x )[index] = Arg[i][j];

      cholmod_dense *x = cholmod_solve( CHOLMOD_A, _lowerDiagonalFactor, rhs, &_settings );
      for ( int i = 0, index = 0; i < Dest.numComponents(); i++ )
        for ( int j = 0; j < Dest[i].size(); j++, index++ )
          Dest[i][j] = reinterpret_cast<double*>( x->x )[index];

      cholmod_free_dense( &x, &_settings );
      cholmod_free_dense( &rhs, &_settings );
    }
  }
};

/**
* \brief A direct solver for block matrices \f$ B \f$ without zero eigenvalue. The matrix neither needs to be symmetric nor positive definite.
*        When calling apply(r,x) the normal equations \f$ B^TB x = B^T r \f$ are solved for \f$ x \f$ by means of aol::CholeskyBlockInverseOp<>.
* \note To compile this class the external cholmod is necessary.
* \note Each time the matrix changes one has to create a new instance of this solver.
*
* \author Teusner
* \ingroup directSolver
*/
template <typename RealType, typename MatrixType, int NumRowBlocks, int NumColBlocks = NumRowBlocks >
class CholeskyBiBlockInverseOp :
  public aol::InverseOp<aol::MultiVector<RealType> > {
protected:
  mutable aol::SparseBlockMatrix<MatrixType> *_pBlockMatTransposed;
  mutable aol::SparseBlockMatrix<aol::SparseMatrix<RealType> > *_pMatSqr;
  mutable aol::CholeskyBlockInverseOp<RealType, aol::SparseMatrix<RealType> > *_pCholeskySolver;

public:
  CholeskyBiBlockInverseOp( const aol::SparseBlockMatrix<MatrixType> &BlockMat ) :
    aol::InverseOp<aol::MultiVector<RealType> > ( ),
    _pBlockMatTransposed ( NULL ),
    _pMatSqr( new aol::SparseBlockMatrix<aol::SparseMatrix<RealType> > (NumColBlocks, NumColBlocks) ),
    _pCholeskySolver ( NULL ) {

      // note: allocating B^TB matrix is not necessary as this is done automatically in addMatrixProduct()!

      _pBlockMatTransposed = new aol::SparseBlockMatrix<MatrixType> ( NumColBlocks, NumRowBlocks );
      BlockMat.transposeTo(*_pBlockMatTransposed);
          
      _pMatSqr->addMatrixProduct(*_pBlockMatTransposed, BlockMat);
      
      _pCholeskySolver = new aol::CholeskyBlockInverseOp<RealType, aol::SparseMatrix<RealType> > ( _pMatSqr->maxNumNonZeroesPerRow() );
      _pCholeskySolver->setMatrix(*_pMatSqr);
  }

  virtual ~CholeskyBiBlockInverseOp() {
    if (_pBlockMatTransposed) delete _pBlockMatTransposed;
    if (_pMatSqr) delete _pMatSqr;
    if (_pCholeskySolver) delete _pCholeskySolver;
  }

  virtual void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::MultiVector<RealType> result ( Dest, aol::STRUCT_COPY );
    apply ( Arg, result );
    Dest += result;
  }

  virtual void apply ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::MultiVector<RealType> NewArg (Arg, aol::STRUCT_COPY);
    _pBlockMatTransposed->apply(Arg, NewArg);

    _pCholeskySolver->apply(NewArg, Dest);
  }
};
/**
* \brief A direct LU-factorization solver for square sparse unsymmetric matrices, using UMFPACK
*
* (see <a href="https://www.cise.ufl.edu/research/sparse/umfpack/UMFPACK/Doc/QuickStart.pdf">quickstart documentation</a>)
* (a <a href="https://www.cise.ufl.edu/research/sparse/umfpack/UMFPACK/Doc/UserGuide.pdf">more detailed documentation</a>)
*
* \author Wirth, modified by Perl, switched from umfpack_di_* (int 32-bit, only 2GB memory) to umfpack_dl_* (long int 64 bit)
* \ingroup directSolver
*/
template <typename RealType, typename MatrixType>
class UMFPACKBlockInverseOp :
  public aol::InverseOp<aol::MultiVector<RealType> > {

private:
// Only recent versions of SuiteSparse define SuiteSparse_long, but "long int" is not
// what SuiteSparse is using as long data type on all platforms. If the types differ
// we have to resort to SuiteSparse_long.
#if defined(__MINGW32__)
  typedef SuiteSparse_long SuiteSparseIntType;
#else
  typedef long int SuiteSparseIntType;
#endif

  // the system matrix in umfpack format (number of rows and columns, indices and values of compressed column format)
  SuiteSparseIntType _n;
  aol::Vector<SuiteSparseIntType> _Ap, _Ai;
  aol::Vector<RealType> _Ax;
  // umfpack control objects (the LU factors of the matrix, control and info variables)
  void *_umfpackNumeric;
  double _umfpackControl[UMFPACK_CONTROL];
  mutable double _umfpackInfo[UMFPACK_INFO];
  bool _matrixSet;

public:
  UMFPACKBlockInverseOp( const aol::BlockOpBase<RealType,MatrixType> &BlockMat ) :
    aol::InverseOp<aol::MultiVector<RealType> > (),
    _umfpackNumeric( NULL ),
    _matrixSet( false ){
    // set desired parameters (currently only default parameters implemented)
    umfpack_dl_defaults( _umfpackControl );
#ifdef DEBUG
    _umfpackControl [UMFPACK_PRL] = 6 ;
#endif
    umfpack_dl_report_control ( _umfpackControl ) ;
    // check for double RealType (UMFPACK only supports double)
    if ( aol::RealTrait<RealType>::ALIAS != aol::DOUBLE )
      throw aol::Exception( "UMFPACKBlockInverseOp currently only implemented for double." );
    setMatrix( BlockMat );
  }

  UMFPACKBlockInverseOp() :
    aol::InverseOp<aol::MultiVector<RealType> > (),
    _umfpackNumeric( NULL ),
    _matrixSet( false ) {
    // set desired parameters (currently only default parameters implemented)
    umfpack_dl_defaults( _umfpackControl );
#ifdef DEBUG
    _umfpackControl [UMFPACK_PRL] = 6 ;
#endif
    umfpack_dl_report_control ( _umfpackControl ) ;    
    // check for double RealType (UMFPACK only supports double)
    if ( aol::RealTrait<RealType>::ALIAS != aol::DOUBLE )
      throw aol::Exception( "UMFPACKBlockInverseOp currently only implemented for double." );
  }

  virtual ~UMFPACKBlockInverseOp() {
    if ( _umfpackNumeric != NULL )
      // free numeic factorization
      umfpack_dl_free_numeric( &_umfpackNumeric );
  }

  //! \warning Do not call this method with matrices without entries!
  inline void setMatrix ( const aol::BlockOpBase<RealType,MatrixType> &BlockMat ) {
    // read in all values and put them into the triplet vectors
    aol::Vector<SuiteSparseIntType> rowInd, colInd;
    aol::Vector<RealType> val;

    SuiteSparseIntType totalRow = 0;
    for ( int row = 0; row < BlockMat.getNumRows(); row++ ) {
      for ( int localRow = 0; localRow < BlockMat.getReference( row, 0 ).getNumRows(); localRow++ ) {
        SuiteSparseIntType colOffset = 0;
        for ( int col = 0; col < BlockMat.getNumCols(); col++ ) {
          vector<typename aol::Row<RealType>::RowEntry > matRow;
          BlockMat.getReference( row, col ).makeRowEntries( matRow, localRow );
          for ( typename vector<typename aol::Row<RealType>::RowEntry >::iterator it = matRow.begin(); it != matRow.end(); ++it ) {
            rowInd.pushBack( totalRow );
            // this cast should always be possible!
            SuiteSparseIntType totalCol = colOffset + static_cast<SuiteSparseIntType>(it->col);
            colInd.pushBack( totalCol );
            val.pushBack( it->value );
          }
          colOffset += BlockMat.getReference( row, col ).getNumCols();
        }
        totalRow++;
      }
    }

    // convert triplet format into required format
    _n = static_cast<SuiteSparseIntType>( aol::Max( rowInd.getMaxValue(), colInd.getMaxValue() ) + 1 );
    _Ap.resize( _n + 1 );
    _Ai.resize( static_cast<SuiteSparseIntType>(val.size()) );
    _Ax.resize( static_cast<SuiteSparseIntType>(val.size()) );
    umfpack_dl_triplet_to_col( _n, _n, val.size(), rowInd.getData(), colInd.getData(), val.getData(), _Ap.getData(), _Ai.getData(), _Ax.getData(), NULL );

    // perform symbolic and numeric factorization
    if ( _umfpackNumeric != NULL )
      umfpack_dl_free_numeric( &_umfpackNumeric );
    void *umfpackSymbolic;

    // check factorization status
    int status;

    // symbolic factorization
    status = umfpack_dl_symbolic( _n, _n, _Ap.getData(), _Ai.getData(), _Ax.getData(), &umfpackSymbolic, _umfpackControl, _umfpackInfo );
    if ( status != UMFPACK_OK ){
      umfpack_dl_report_status ( _umfpackControl, status) ;
      cerr << aol::color::red << "umfpack_di_symbolic failed!" << aol::color::reset << endl;
    }
    // numeric factorization
    status = umfpack_dl_numeric( _Ap.getData(), _Ai.getData(), _Ax.getData(), umfpackSymbolic, &_umfpackNumeric, _umfpackControl, _umfpackInfo );
    if ( status != UMFPACK_OK ){
      umfpack_dl_report_status ( _umfpackControl, status) ;
      cerr << aol::color::red << "umfpack_di_numeric failed!" << aol::color::reset << endl;
    }

    // free symbolic factorization
    umfpack_dl_free_symbolic( &umfpackSymbolic );
    _matrixSet = true;
  }

  // standard applyAdd
  virtual void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::MultiVector<RealType> result ( Dest, aol::STRUCT_COPY );
    apply ( Arg, result );
    Dest += result;
  }

  virtual void apply ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    
    if( !_matrixSet )
      throw aol::Exception( "UMFPACKBlockInverseOp::apply(): matrix not set." );
    
    // transform rhs into appropriate format
    aol::Vector<RealType> rhs( Arg.getTotalSize () ), sol( Dest.getTotalSize () );
    rhs.copyUnblockedFrom( Arg );

    // solve and copy solution into result vector

    umfpack_dl_solve( UMFPACK_A, _Ap.getData(), _Ai.getData(), _Ax.getData(), sol.getData(), rhs.getData(), _umfpackNumeric, _umfpackControl, _umfpackInfo );
    // if umfpackControl [UMFPACK_PRL] is set, then this we give you al the information you need
    umfpack_dl_report_info (_umfpackControl, _umfpackInfo) ;

    Dest.copySplitFrom( sol );
  }
  
  bool isMatrixSet() const {
    return _matrixSet;
  }
};
/**
* \brief Version of the UMFPACKInverseOp for usage with aol::TripletMatrix.
*
* A direct LU-factorization solver for square sparse unsymmetric matrices, using UMFPACK
* (see quickstart documentation, http://www.cise.ufl.edu/research/sparse/umfpack/current/UMFPACK/Doc/QuickStart.pdf)
*
* Version using aol::TripletMatrix.
*
* \author Toelkes, heavily based on work by Wirth.
* \ingroup directSolver
*/
template <typename RealType>
class TripletMatrixUMFPACKInverseOp :
  public aol::InverseOp<aol::Vector<RealType> > {

private:
  // the system matrix in umfpack format (number of rows and columns, indices and values of compressed column format)
  int _n;
  aol::Vector<int> _Ap, _Ai;
  aol::Vector<RealType> _Ax;
  // umfpack control objects (the LU factors of the matrix, control and info variables)
  void *_umfpackNumeric;
  double _umfpackControl[UMFPACK_CONTROL];
  mutable double _umfpackInfo[UMFPACK_INFO];

  // The constructor should always be given a matrix.
  TripletMatrixUMFPACKInverseOp() :
    aol::InverseOp<aol::MultiVector<RealType> > (),
    _n ( 0 ), _umfpackNumeric( NULL ) {
  }

public:
  //! \brief Standard constructor.
  //! \param[in] mat The system matrix in triplet format.
  TripletMatrixUMFPACKInverseOp ( const aol::TripletMatrix<RealType> &mat ) :
    aol::InverseOp<aol::Vector<RealType> > (),
    _umfpackNumeric( NULL ) {
    // set desired parameters (currently only default parameters implemented)
    umfpack_di_defaults( _umfpackControl );
    // check for double RealType (UMFPACK only supports double)
    if ( aol::RealTrait<RealType>::ALIAS != aol::DOUBLE )
      throw aol::Exception( "UMFPACKBlockInverseOp currently only implemented for double." );

    // convert triplet format into required format
    if ( mat.getNumRows () != mat.getNumCols () )
      throw aol::Exception ( "Matrix is not quadratic", __FILE__, __LINE__, __FUNCTION__ );

    // Use arrays, that are already present in the TripletMatrix.
    _n = mat.getNumRows ();
    _Ap.resize( _n + 1 );
    _Ai.resize( mat.getValueReference ().size () );
    _Ax.resize( mat.getValueReference ().size () );
    if ( umfpack_di_triplet_to_col ( _n, _n, mat.getValueReference ().size (),
        mat.getRowIndexReference ().getData (), mat.getColIndexReference ().getData(), mat.getValueReference ().getData(),
        _Ap.getData(), _Ai.getData(), _Ax.getData(), NULL ) != UMFPACK_OK )
      cerr << aol::color::red << "umfpack_di_triplet_to_col failed!" << aol::color::reset << endl;

    // perform symbolic and numeric factorization
    if ( _umfpackNumeric != NULL )
      umfpack_di_free_numeric( &_umfpackNumeric );
    void *umfpackSymbolic;
    if ( umfpack_di_symbolic( _n, _n, _Ap.getData(), _Ai.getData(), _Ax.getData(), &umfpackSymbolic, _umfpackControl, _umfpackInfo ) != UMFPACK_OK )
      cerr << aol::color::red << "umfpack_di_symbolic failed!" << aol::color::reset << endl;
    if ( umfpack_di_numeric( _Ap.getData(), _Ai.getData(), _Ax.getData(), umfpackSymbolic, &_umfpackNumeric, _umfpackControl, _umfpackInfo )
        != UMFPACK_OK )
      cerr << aol::color::red << "umfpack_di_numeric failed!" << aol::color::reset << endl;
    umfpack_di_free_symbolic( &umfpackSymbolic );
  }

  TripletMatrixUMFPACKInverseOp ( const aol::CSCMatrix<RealType> &mat )
  : aol::InverseOp<aol::Vector<RealType> > (), _Ap ( mat.getColumnPointerReference (), aol::FLAT_COPY ),
    _Ai ( mat.getRowIndexReference (), aol::FLAT_COPY ), _Ax ( mat.getValueReference (), aol::FLAT_COPY ),
    _umfpackNumeric( NULL ) {
    // set desired parameters (currently only default parameters implemented)
    umfpack_di_defaults( _umfpackControl );
    // check for double RealType (UMFPACK only supports double)
    if ( aol::RealTrait<RealType>::ALIAS != aol::DOUBLE )
      throw aol::Exception( "UMFPACKBlockInverseOp currently only implemented for double." );

    // convert triplet format into required format
    if ( mat.getNumRows () != mat.getNumCols () )
      throw aol::Exception ( "Matrix is not quadratic", __FILE__, __LINE__, __FUNCTION__ );

    // Use arrays, that are already present in the TripletMatrix.
    _n = mat.getNumRows ();

    // perform symbolic and numeric factorization
    if ( _umfpackNumeric != NULL )
      umfpack_di_free_numeric( &_umfpackNumeric );
    void *umfpackSymbolic;

    int status = umfpack_di_symbolic( _n, _n, _Ap.getData (), _Ai.getData (), _Ax.getData (),
        &umfpackSymbolic, _umfpackControl, _umfpackInfo );
    if ( status != UMFPACK_OK )
      cerr << aol::color::red << "umfpack_di_symbolic failed! (" << status << ")" << aol::color::reset << endl;
    status = umfpack_di_numeric( _Ap.getData (), _Ai.getData (), _Ax.getData (),
        umfpackSymbolic, &_umfpackNumeric, _umfpackControl, _umfpackInfo );
    if ( status != UMFPACK_OK )
      cerr << aol::color::red << "umfpack_di_numeric failed! (" << status << ")" << aol::color::reset << endl;
    umfpack_di_free_symbolic( &umfpackSymbolic );
  }

  //! \brief Destructor.
  virtual ~TripletMatrixUMFPACKInverseOp() {
    if ( _umfpackNumeric != NULL )
      umfpack_di_free_numeric( &_umfpackNumeric );
  }

  //! \brief Solve the system and add the solution to dest.
  //! \param[in] arg The right hand side.
  //! \param[out] dest The solution is added to this vector.
  virtual void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> result ( Dest, aol::STRUCT_COPY );
    apply ( Arg, result );
    Dest += result;
  }

  //! \brief Solve the system.
  //! \param[in] arg The right hand side.
  //! \param[out] dest The solution.
  virtual void apply ( const aol::Vector<RealType> &arg, aol::Vector<RealType> &dest ) const {
    // solve and copy solution into result vector
    umfpack_di_solve( UMFPACK_A, _Ap.getData(), _Ai.getData(), _Ax.getData(), dest.getData(), arg.getData(), _umfpackNumeric, _umfpackControl, _umfpackInfo );
  }
};

/**
* \brief Version of the UMFPACKBlockInverseOp for usage with aol::CSCMatrix. Takes a unblocked matrix, but MultiVectors as rhs and solution.
*
* A direct LU-factorization solver for square sparse unsymmetric matrices, using UMFPACK
* (see quickstart documentation, http://www.cise.ufl.edu/research/sparse/umfpack/current/UMFPACK/Doc/QuickStart.pdf)
*
* Version using an unblocked aol::CSCMatrix (because otherwise it has to be converted for UMFPACK) and aol::MultiVector.
*
* \author Toelkes, heavily based on work by Wirth.
* \ingroup directSolver
*/
template <typename RealType>
class CSCMatrixUMFPACKBlockInverseOp :
  public aol::InverseOp<aol::MultiVector<RealType> > {

private:
  // the system matrix in umfpack format (number of rows and columns, indices and values of compressed column format)
  int _n;
  const aol::CSCMatrix<RealType> &_cscMat;

  // umfpack control objects (the LU factors of the matrix, control and info variables)
  void *_umfpackNumeric;
  double _umfpackControl[UMFPACK_CONTROL];
  mutable double _umfpackInfo[UMFPACK_INFO];

  // The constructor should always be given a matrix.
  CSCMatrixUMFPACKBlockInverseOp() :
    aol::InverseOp<aol::MultiVector<RealType> > (),
    _n ( 0 ), _umfpackNumeric( NULL ) {
  }

public:
  //! \brief Standard constructor.
  //! \param[in] mat The system matrix in triplet format.
  CSCMatrixUMFPACKBlockInverseOp ( const aol::CSCMatrix<RealType> &mat ) :
    aol::InverseOp<aol::MultiVector<RealType> > (), _cscMat ( mat ),
    _umfpackNumeric( NULL ) {
    // set desired parameters (currently only default parameters implemented)
    umfpack_di_defaults( _umfpackControl );
    // check for double RealType (UMFPACK only supports double)
    if ( aol::RealTrait<RealType>::ALIAS != aol::DOUBLE )
      throw aol::Exception( "UMFPACKBlockInverseOp currently only implemented for double." );

    // convert triplet format into required format
    if ( mat.getNumRows () != mat.getNumCols () )
      throw aol::Exception ( "Matrix is not quadratic", __FILE__, __LINE__, __FUNCTION__ );

    // Use arrays, that are already present in the TripletMatrix.
    _n = mat.getNumRows ();

    // perform symbolic and numeric factorization
    if ( _umfpackNumeric != NULL )
      umfpack_di_free_numeric( &_umfpackNumeric );
    void *umfpackSymbolic;

    int status = umfpack_di_symbolic( _n, _n, _cscMat.getColumnPointerReference ().getData (), _cscMat.getRowIndexReference ().getData (),
        _cscMat.getValueReference ().getData (), &umfpackSymbolic, _umfpackControl, _umfpackInfo );
    if ( status != UMFPACK_OK )
      cerr << aol::color::red << "umfpack_di_symbolic failed! (" << status << ")" << aol::color::reset << endl;

    status = umfpack_di_numeric( _cscMat.getColumnPointerReference ().getData (), _cscMat.getRowIndexReference ().getData (),
        _cscMat.getValueReference ().getData (), umfpackSymbolic, &_umfpackNumeric, _umfpackControl, _umfpackInfo );
    if ( status != UMFPACK_OK )
      cerr << aol::color::red << "umfpack_di_numeric failed! (" << status << ")" << aol::color::reset << endl;

    umfpack_di_free_symbolic( &umfpackSymbolic );
  }

  //! \brief Destructor.
  virtual ~CSCMatrixUMFPACKBlockInverseOp() {
    if ( _umfpackNumeric != NULL )
      umfpack_di_free_numeric( &_umfpackNumeric );
  }

  //! \brief Solve the system and add the solution to dest.
  //! \param[in] Arg The right hand side.
  //! \param[out] Dest The solution is added to this vector.
  virtual void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::MultiVector<RealType> result ( Dest, aol::STRUCT_COPY );
    apply ( Arg, result );
    Dest += result;
  }

  //! \brief Solve the system.
  //! \param[in] arg The right hand side.
  //! \param[out] dest The solution.
  virtual void apply ( const aol::MultiVector<RealType> &arg, aol::MultiVector<RealType> &dest ) const {
    aol::Vector<RealType> rhs( _n ), sol( _n );
    rhs.copyUnblockedFrom ( arg );

    // solve and copy solution into result vector
    int status = umfpack_di_solve( UMFPACK_A, _cscMat.getColumnPointerReference ().getData (), _cscMat.getRowIndexReference ().getData (),
        _cscMat.getValueReference ().getData (), sol.getData(), rhs.getData(), _umfpackNumeric, _umfpackControl, _umfpackInfo );

    if ( status != UMFPACK_OK )
      cerr << aol::color::red << "umfpack_di_solve failed (" << status << ")" << aol::color::reset << endl;

    dest.copySplitFrom ( sol );
  }
};

/**
* \brief a direct QR-solver for square matrices which can be umsymmetric and any kind
* 
* \author Perl, based on the Cholesky code of Wirth. BUT here always long int (indicated by cholmod_l...),
*  i.e. cholmod_i* uses int 32-bit (only 2GB memory) and cholmod_l* uses long int 64 bit
* Details:
*         - http://www.cise.ufl.edu/research/sparse/SPQR/SPQR/Doc/spqr.pdf
*         - http://www.cise.ufl.edu/research/sparse/SPQR/SPQR/Doc/spqr_user_guide.pdf
* 
* \ingroup directSolver
*/
template <typename RealType, typename MatrixType>
class SparseQRBlockInverseOp : public aol::InverseOp<aol::MultiVector<RealType> > {

private:
  // settings
  mutable cholmod_common _settings;
  // matrix in cholmod datatype
  cholmod_sparse* _matrix;

  int _maxNumEntriesPerBlockRow;  //!< The max. number of non-zero entries per row in a block.

  inline void factorize( cholmod_triplet* triplet, const int numNonZeros )  {
    cholmod_l_free_sparse( &_matrix, &_settings );
    _matrix = cholmod_l_triplet_to_sparse( triplet, numNonZeros, &_settings );
#ifdef VERBOSE
    _settings.print = 5;
    cholmod_l_check_triplet( triplet, &_settings );
    cerr << endl << "Check status " << _settings.status << endl;
    cholmod_l_print_triplet( triplet, "", &_settings );
    cerr << "Print status " << _settings.status << endl;
    cholmod_l_print_sparse( _matrix, "", &_settings );
#endif

  }

public:
  /*!
   * \brief Constructs the Op, initializes cholmod
   *
   * _maxNumEntriesPerBlockRow is only used in setMatrix and should be set in that call
   * the argument here is for compatibility only
   */
  SparseQRBlockInverseOp( const int MaxNumEntriesPerBlockRow = 0 ) 
    : aol::InverseOp<aol::MultiVector<RealType> > ( ),
      _matrix( NULL ),
      _maxNumEntriesPerBlockRow( MaxNumEntriesPerBlockRow ){
    // start cholmod with long int
    cholmod_l_start ( &_settings );    
    // CHOLMOD_SINGLE not yet supported
    if ( aol::RealTrait<RealType>::ALIAS == aol::FLOAT )
      _settings.dtype = CHOLMOD_DOUBLE; // CHOLMOD_SINGLE not yet supported
    else
      _settings.dtype = CHOLMOD_DOUBLE;
  }

  SparseQRBlockInverseOp( const aol::BlockOpBase<RealType,MatrixType> &BlockMat, 
                          const int MaxNumEntriesPerBlockRow = 0 )
    : aol::InverseOp<aol::MultiVector<RealType> > ( ),
      _matrix( NULL ),
      _maxNumEntriesPerBlockRow( MaxNumEntriesPerBlockRow ){
      // start cholmod with long int
      cholmod_l_start ( &_settings );    
      // CHOLMOD_SINGLE not yet supported
      if ( aol::RealTrait<RealType>::ALIAS == aol::FLOAT )
        _settings.dtype = CHOLMOD_DOUBLE; // CHOLMOD_SINGLE not yet supported
      else
        _settings.dtype = CHOLMOD_DOUBLE;
      setMatrix( BlockMat, MaxNumEntriesPerBlockRow );
  }

  virtual ~SparseQRBlockInverseOp() {
    // destroy 
    cholmod_l_free_sparse( &_matrix, &_settings );
    cholmod_l_finish ( &_settings );
  }

  /*!
   * \brief Constructs the Op, assumes d by d blocks (d=dimension).
   * \param[in] MaxNumEntriesPerBlockRow The max. number of non-zero entries in a row of a block (taken over all blocks).
   *
   * The constructor takes the number of non-zero entries per row in a block, which can be set to \f$3^d\f$,
   * (\f$d = \f$ dimension) if the standard stencil in a uniform FE grid and standard FE operators are used.
   */
  inline void setMatrix ( const aol::BlockOpBase<RealType,MatrixType> &BlockMat, const int MaxNumEntriesPerBlockRow = 0  ) {

    // MaxNumEntriesPerBlockRow == 0 makes no sense!!!
    if ( _maxNumEntriesPerBlockRow == 0 )  {
      if ( MaxNumEntriesPerBlockRow == 0 )
        throw aol::Exception( "You have to set MaxNumEntriesPerBlockRow in setMatrix call!" );
      else _maxNumEntriesPerBlockRow = MaxNumEntriesPerBlockRow;
    }

    // check whether all blocks are allocated
    for ( int i = 0; i < BlockMat.getNumRows(); i++ )
      for ( int j = 0; j < BlockMat.getNumCols(); j++ )
        if( !BlockMat.getPointer(i, j) )
          throw aol::Exception( "SparseQRBlockInverseOp::setMatrix: all blocks must be allocated!" );

    // compute number of total rows (and columns)
    int n = 0;
    for ( int i = 0; i < BlockMat.getNumRows(); i++ )
      n += BlockMat.getReference( i, 0 ).getNumRows();

    // allocate new matrix
    // Todo: cases
    cholmod_triplet* triplet = cholmod_l_allocate_triplet( n, n, n * _maxNumEntriesPerBlockRow * BlockMat.getNumCols(), 0, CHOLMOD_REAL, &_settings );
    // number of entries
    long int numberOfEntries = 0;
    long int totalRow = 0;
    // iterate over block rows
    for ( int row = 0; row < BlockMat.getNumRows(); row++ ) {
      // iterate over rows in the block
      for ( int localRow = 0; localRow < BlockMat.getReference( row, 0 ).getNumRows(); localRow++ ) {
        // see below for explanation
        long int colOffset = 0;
        // iterate over all columns
        for ( int col = 0; col < BlockMat.getNumCols(); col++ ) {
          // create a new row
          vector< typename aol::Row<RealType>::RowEntry > matRow;
          BlockMat.getReference( row, col ).makeRowEntries( matRow, localRow );
          // iterate over the block row
          for ( typename vector<typename aol::Row<RealType>::RowEntry >::iterator it = matRow.begin(); it != matRow.end(); ++it ){
            reinterpret_cast<long int*>( triplet->i )[numberOfEntries] = totalRow;
            // it->col is the current column, BUT these entries can be unsorted in aol::Row!!!
            long int totalCol = colOffset + static_cast<long int>( it->col);
            reinterpret_cast<long int*>( triplet->j )[numberOfEntries] = totalCol;
            reinterpret_cast<double*>( triplet->x )[numberOfEntries] = it->value;
            // increase
            numberOfEntries++;
          }
          // colOffset = 0, numCols, 2*numCols, ...
          colOffset += BlockMat.getReference( row, col ).getNumCols();
        }
        totalRow++;
      }
    }
    triplet->nnz = numberOfEntries;
    factorize( triplet, n * _maxNumEntriesPerBlockRow * BlockMat.getNumCols() );
    cholmod_l_free_triplet( &triplet, &_settings );
  }

/* TODO: Is it enough to cast to long int????
  inline void setMatrix ( const aol::TripletMatrix<RealType> & tmat )  {
    cholmod_triplet ctriplet;
    ctriplet.nrow  = tmat.getNumRows();
    ctriplet.ncol  = tmat.getNumCols();
    ctriplet.nzmax = tmat.getValueReference().size();
    ctriplet.nnz   = tmat.getValueReference().size(); // better way to do this?
    ctriplet.i     = tmat.getRowIndexReference().getData();
    ctriplet.j     = tmat.getColIndexReference().getData();
    ctriplet.x     = tmat.getValueReference().getData();
    ctriplet.z     = NULL;
    ctriplet.stype = 0;  // store unsymmetric matrix
    ctriplet.itype = CHOLMOD_INT;
    ctriplet.xtype = CHOLMOD_REAL;
    ctriplet.dtype = CHOLMOD_DOUBLE;
    factorize( &ctriplet, ctriplet.nzmax ); // number of nonzeros is much smaller!
  }
*/

  // exports lower diagonal matrix in Matrix Market format
  void exportMatrix( const char* Filename ) const {
    FILE* matrix = fopen( Filename, "w+" );
    if ( matrix == NULL )
      cerr<<"error occured while opening the file"<<endl;
    cerr << cholmod_l_write_sparse( matrix, _matrix, NULL, "", &_settings );
    fclose( matrix );
  }

  virtual void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::MultiVector<RealType> result ( Dest, aol::STRUCT_COPY );
    apply ( Arg, result );
    Dest += result;
  }

  virtual void apply ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {

    // write Arg to rhs
    cholmod_dense *rhs = cholmod_l_allocate_dense( Arg.getTotalSize(), 1, Arg.getTotalSize(), CHOLMOD_REAL, &_settings );
    long int index = 0;
    for ( int i = 0; i < Arg.numComponents(); i++ ){
      for ( int j = 0; j < Arg[i].size(); j++){
        reinterpret_cast<double*>( rhs->x )[index] = Arg[i][j];
        index++;
      }
    }
    
    // solve
    cholmod_dense *x = SuiteSparseQR<RealType> ( _matrix, rhs, &_settings );

    // write solution to Dest
    index = 0;
    for ( int i = 0; i < Dest.numComponents(); i++ ){
      for ( int j = 0; j < Dest[i].size(); j++){
        Dest[i][j] = reinterpret_cast<double*>( x->x )[index];
        index++;
      }
    }

    // free x and rhs
    cholmod_l_free_dense( &x, &_settings );
    cholmod_l_free_dense( &rhs, &_settings );
  }
};

} // end namespace

#else

namespace aol {

template <typename RealType, typename MatrixType>
class CholeskyBlockInverseOp : public aol::InverseOp<aol::MultiVector<RealType> > {

public:
  CholeskyBlockInverseOp( const int /*MaxNumEntriesPerBlockRow*/ = 0 ) { }

  void setMatrix ( const aol::BlockOpBase<RealType,MatrixType> &/*BlockMat*/, const int /*MaxNumEntriesPerBlockRow*/ = 0 ) {
    throw aol::Exception ( "aol::CholeskyBlockInverseOp needs suitesparse! Compile with -DUSE_SUITESPARSE", __FILE__, __LINE__ );
  }

  virtual void applyAdd( const aol::MultiVector<RealType> &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/ ) const {
    throw aol::Exception ( "aol::CholeskyBlockInverseOp needs suitesparse! Compile with -DUSE_SUITESPARSE", __FILE__, __LINE__ );
  }
};

template <typename ConfiguratorType, typename MatrixType, int NumRowBlocks, int NumColBlocks = NumRowBlocks >
class CholeskyBiBlockInverseOp : public aol::InverseOp<aol::MultiVector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;

public:
  CholeskyBiBlockInverseOp( const typename ConfiguratorType::InitType &/*Initializer*/,
                            const aol::SparseBlockMatrix<MatrixType> &/*BlockMat*/ ) { }

  virtual void applyAdd( const aol::MultiVector<RealType> &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/ ) const {
    throw aol::Exception ( "aol::CholeskyBiBlockInverseOp needs suitesparse! Compile with -DUSE_SUITESPARSE", __FILE__, __LINE__ );
  }
};

template <typename RealType, typename MatrixType>
class UMFPACKBlockInverseOp :
  public aol::InverseOp<aol::MultiVector<RealType> > {
public:
  UMFPACKBlockInverseOp( const aol::BlockOpBase<RealType,MatrixType> & ) {}
  UMFPACKBlockInverseOp() {}
  virtual ~UMFPACKBlockInverseOp() {}

  //! \warning Do not call this method with matrices without entries!
  inline void setMatrix ( const aol::BlockOpBase<RealType,MatrixType> & ) {
    throw aol::Exception ( "aol::UMFPACKBlockInverseOp needs suitesparse! Compile with -DUSE_SUITESPARSE", __FILE__, __LINE__ );
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &, aol::MultiVector<RealType> & ) const {
    throw aol::Exception ( "aol::UMFPACKBlockInverseOp needs suitesparse! Compile with -DUSE_SUITESPARSE", __FILE__, __LINE__ );
  }
  bool isMatrixSet() const {
    throw aol::Exception ( "aol::UMFPACKBlockInverseOp needs suitesparse! Compile with -DUSE_SUITESPARSE", __FILE__, __LINE__ );
  }
};

template <typename RealType>
class TripletMatrixUMFPACKInverseOp : public aol::InverseOp<aol::Vector<RealType> > {
public:
  TripletMatrixUMFPACKInverseOp ( const aol::TripletMatrix<RealType> & ) {}
  TripletMatrixUMFPACKInverseOp ( const aol::CSCMatrix<RealType> & ) {}

  virtual void applyAdd( const aol::Vector<RealType> &, aol::Vector<RealType> & ) const {
    throw aol::Exception ( "aol::TripletMatrixUMFPACKInverseOp needs suitesparse! Compile with -DUSE_SUITESPARSE", __FILE__, __LINE__ );
  }
};

template <typename RealType>
class CSCMatrixUMFPACKBlockInverseOp : public aol::InverseOp<aol::MultiVector<RealType> > {
public:
  CSCMatrixUMFPACKBlockInverseOp ( const aol::CSCMatrix<RealType> & ) {}

  virtual void applyAdd( const aol::MultiVector<RealType> &, aol::MultiVector<RealType> & ) const {
    throw aol::Exception ( "aol::CSCMatrixUMFPACKBlockInverseOp needs suitesparse! Compile with -DUSE_SUITESPARSE", __FILE__, __LINE__ );
  }
};

template <typename RealType, typename MatrixType>
class SparseQRBlockInverseOp :
  public aol::InverseOp<aol::MultiVector<RealType> > {
public:
  SparseQRBlockInverseOp( const aol::BlockOpBase<RealType,MatrixType> & , const int MaxNumEntriesPerBlockRow ) {}
  SparseQRBlockInverseOp( const int MaxNumEntriesPerBlockRow ) {}
  virtual ~SparseQRBlockInverseOp() {}

  //! \warning Do not call this method with matrices without entries!
  inline void setMatrix ( const aol::BlockOpBase<RealType,MatrixType> &BlockMat ) {
    throw aol::Exception ( "aol::SparseQRBlockInverseOp needs suitesparse! Compile with -DUSE_SUITESPARSE", __FILE__, __LINE__ );
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &, aol::MultiVector<RealType> & ) const {
    throw aol::Exception ( "aol::SparseQRBlockInverseOp needs suitesparse! Compile with -DUSE_SUITESPARSE", __FILE__, __LINE__ );
  }
};

} // end namespace

#endif

namespace aol {

/**
 * \brief Wrapper that allows to use aol::CholeskyBlockInverseOp on non-block matrices.
 *
 * \author Berkels
 * \ingroup directSolver
 */
template <typename RealType, typename MatrixType>
class CholeskyInverseOp : public aol::InverseOp<aol::Vector<RealType> > {
  aol::BlockOp<RealType, MatrixType> _blockMatrix;
  aol::CholeskyBlockInverseOp<RealType, MatrixType> _choleskyBlockSolver;
  aol::MVecToVecOp<RealType> _solver;
public:
  CholeskyInverseOp ( )
    : _blockMatrix( 1, 1 ),
      _solver ( _choleskyBlockSolver ) { }

  void setMatrix ( MatrixType &Mat, const int MaxNumEntriesPerRow ) {
    _blockMatrix.setReference( 0, 0, Mat );
    _choleskyBlockSolver.setMatrix ( _blockMatrix, MaxNumEntriesPerRow );
  }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    _solver.applyAdd ( Arg, Dest );
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    _solver.apply ( Arg, Dest );
  }
};

} // end namespace

#endif
