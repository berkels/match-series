#ifndef __TRUSTREGIONMETHOD_H
#define __TRUSTREGIONMETHOD_H

#include <aol.h>
#include <ringBuffer.h>
#include <suiteSparseSolver.h>

namespace aol {

/**
 * Class to create aol::MultiVector whose single entries of aol::Vector are stored right one after another.
 * Needed for a simple conversion between quocmesh and cholmod format.
 * Attention: Problematic with USE_SSE due to alignment issues.
 */
template <typename _DataType>
class MultiVectorWithConcatData :
  public aol::MultiVector<_DataType> {

protected:
  aol::Vector<_DataType> _dataVector;
  const aol::Vector<int> _componentSizes;

public:
  typedef _DataType DataType;

  //! Split one long array to MultiVector with components of specified size.
  MultiVectorWithConcatData ( DataType* Data, const aol::Vector<int>& ComponentSizes, aol::CopyFlag CopyFlag = aol::FLAT_COPY ) :
    _dataVector( Data, ComponentSizes.sum(), CopyFlag ),
    _componentSizes( ComponentSizes, aol::DEEP_COPY ) {
    int k = 0;
    this->vecs.reserve ( ComponentSizes.size() );
    for ( int i = 0; i < ComponentSizes.size(); ++i ) {
      try {
        this->vecs.push_back ( typename aol::MultiVector<_DataType>::vec_entry ( new aol::Vector<DataType> ( &_dataVector[k], ComponentSizes[i], aol::FLAT_COPY ) ) );
      }
      catch ( aol::Exception& ex ) {
        ex.consume ();
        throw  Exception ( "aol::MultiVectorWithConcatData<DataType>::MultiVectorWithConcatData ( DataType* Data, const aol::Vector<int>& ComponentSizes, aol::CopyFlag CopyFlag = aol::FLAT_COPY ) With USE_SSE all entries of ComponentSizes must be multiples of 16.\n", __FILE__, __LINE__ );
      }
      k += ComponentSizes[i];
    }
  }

  //! Split one long Vector to MultiVector with shorter components of specified size.
  MultiVectorWithConcatData ( aol::Vector<DataType>& Data, const aol::Vector<int>& ComponentSizes, aol::CopyFlag CopyFlag = aol::FLAT_COPY ) :
    _dataVector( Data, CopyFlag ),
    _componentSizes( ComponentSizes, aol::DEEP_COPY ) {
    int k = 0;
    this->vecs.reserve ( ComponentSizes.size() );
    for ( int i = 0; i < ComponentSizes.size(); ++i ) {
      try {
        this->vecs.push_back ( typename aol::MultiVector<_DataType>::vec_entry ( new aol::Vector<DataType> ( &_dataVector[k], ComponentSizes[i], aol::FLAT_COPY ) ) );
      }
      catch ( aol::Exception& ex ) {
        ex.consume ();
        throw  Exception ( "aol::MultiVectorWithConcatData<DataType>::MultiVectorWithConcatData ( DataType* Data, const aol::Vector<int>& ComponentSizes, aol::CopyFlag CopyFlag = aol::FLAT_COPY ) With USE_SSE all entries of ComponentSizes must be multiples of 16.\n", __FILE__, __LINE__ );
      }
      k += ComponentSizes[i];
    }
  }

  //! copy constructor
  explicit MultiVectorWithConcatData ( const MultiVectorWithConcatData<DataType>& Data, aol::CopyFlag CopyFlag = aol::FLAT_COPY ) :
    aol::MultiVector<DataType>( 0, 0 ),
    _dataVector( Data.getVectorRef(), CopyFlag ),
    _componentSizes( Data._componentSizes, aol::DEEP_COPY ) {
    int k = 0;
    this->vecs.reserve ( _componentSizes.size() );
    for ( int i = 0; i < _componentSizes.size(); ++i ) {
      try {
        this->vecs.push_back ( typename aol::MultiVector<_DataType>::vec_entry ( new aol::Vector<DataType> ( &_dataVector[k], _componentSizes[i], aol::FLAT_COPY ) ) );
      }
      catch ( aol::Exception& ex ) {
        ex.consume ();
        throw  Exception ( "aol::MultiVectorWithConcatData<DataType>::MultiVectorWithConcatData ( DataType* Data, const aol::Vector<int>& ComponentSizes, aol::CopyFlag CopyFlag = aol::FLAT_COPY ) With USE_SSE all entries of ComponentSizes must be multiples of 16.\n", __FILE__, __LINE__ );
      }
      k += _componentSizes[i];
    }
  }

  MultiVectorWithConcatData<DataType> &operator= ( const MultiVectorWithConcatData<DataType> &Other ) {
    if ( _componentSizes != Other._componentSizes )
      throw aol::Exception( "MultiVectorWithConcatData<DataType>::operator= : dimensions don't match.", __FILE__, __LINE__ );
    _dataVector = Other._dataVector;
    return *this;
  }
  
  const DataType& get( const int I ) const {
    return _dataVector[I];
  }

  void set( const int I, const DataType Value ) {
    _dataVector[I] = Value;
  }

  inline int size() const {
    return _dataVector.size();
  }

  inline DataType* getData() const {
    return _dataVector.getData();
  }

  inline const aol::Vector<DataType>& getVectorRef() const {
    return _dataVector;
  }

  inline aol::Vector<DataType>& getVectorRef() {
    return _dataVector;
  }

private:
  void setVectorRef ( int /*Index*/, aol::Vector<_DataType> &/*Vec*/, bool /*deleteFlag*/=false ) {}
  void reallocate ( const aol::Vector<int> &/*size*/ ) {}
  void reallocate ( const int /*NumComponents*/, const int /*SizeOfComponents*/ ) {}
  void reallocate ( const aol::MultiVector<_DataType> &/*other*/ ) {}
  void reallocate ( const qc::GridStructure &/*grid*/ ) {}
};

template<typename VectorType, typename DataType = typename VectorType::DataType>
class ConcatDataTrait {};

template<typename DataType>
class ConcatDataTrait<aol::Vector<DataType>, DataType> {
public:
  typedef aol::Vector<DataType> ConcatDataType;
  typedef int SizeType;
  static void getSize( const aol::Vector<DataType> &Vec, int& size ) {
    size = Vec.size();
  }
};

template<typename DataType>
class ConcatDataTrait<aol::MultiVector<DataType>, DataType> {
public:
  typedef MultiVectorWithConcatData<DataType> ConcatDataType;
  typedef aol::Vector<int> SizeType;
  static void getSize( const aol::MultiVector<DataType> &Vec, aol::Vector<int>& size ) {
    Vec.getSizes (size);
  }
};

template<typename MatrixType>
class MatrixWrapperForMultiVectorWithConcatData : public MatrixType {
public:
  using MatrixType::apply;
  using MatrixType::applyAdd;
  
  MatrixWrapperForMultiVectorWithConcatData<MatrixType> ( const unsigned int NumRows, const unsigned int NumCols ) :
    MatrixType( NumRows, NumCols ) {}
  
  explicit MatrixWrapperForMultiVectorWithConcatData<MatrixType> ( const MatrixWrapperForMultiVectorWithConcatData<MatrixType> &Other, aol::CopyFlag Flag = aol::DEEP_COPY ) :
    MatrixType( Other, Flag ) {}
  
  void apply( const MultiVectorWithConcatData<typename MatrixType::DataType> &Arg, MultiVectorWithConcatData<typename MatrixType::DataType> &Dest ) const {
    MatrixType::apply( Arg.getVectorRef(), Dest.getVectorRef() );
  }

  //! inefficient implementation, but sufficient for aol::SecondDerivativeValidator
  void apply( const MultiVector<typename MatrixType::DataType> &Arg, MultiVector<typename MatrixType::DataType> &Dest ) const {
    static bool performanceWarningPrinted = false;
    if ( ! performanceWarningPrinted ) {
      cerr << aol::color::red << "aol::MatrixWrapperForMultiVectorWithConcatData::apply on aol::MultiVector is not implemented efficiently." << aol::color::reset << endl;
      performanceWarningPrinted = true;
    }

    const aol::Vector<typename MatrixType::DataType> vecArg ( Arg );
    aol::Vector<typename MatrixType::DataType> vecDest ( vecArg, aol::STRUCT_COPY );
    MatrixType::apply( vecArg, vecDest );
    Dest.copySplitFrom ( vecDest );
  }

  void applyAdd( const MultiVectorWithConcatData<typename MatrixType::DataType> &Arg, MultiVectorWithConcatData<typename MatrixType::DataType> &Dest ) const {
    MatrixType::applyAdd( Arg.getVectorRef(), Dest.getVectorRef() );
  }

  //! dummy function to make second derivative validator work, cf. aol::SecondDerivativeValidator<>
  //! \todo re-consider this!
  void transpose() const {
    throw aol::Exception("aol::MatrixWrapperForMultiVectorWithConcatData::transpose: not implemented", __FILE__, __LINE__ );
  }
};

#ifdef USE_EXTERNAL_SUITESPARSE

/**
 * \brief Newton trust region method
 *
 * Implements the algorithms 6.1.1 and 7.3.1-4 in Conn, Gould, Toint: Trust-region methods.
 * Is implemented for _VectorType = aol::Vector or aol::MultiVector;
 * in the latter case, the _SecondDerivativeType should be a block matrix
 * instead of a normal matrix, and the type of each block is given as _SubMatrixType.
 * If an iterative solver instead of a direct one is to be used, also specify the _PrecondType
 * (by default, the preconditioner is set to the _SecondDerivativeType, which makes no sense).
 * \ingroup Optimization
 */
template<typename _RealType, typename _VectorType, typename _SecondDerivativeType, typename _PrecondType = _SecondDerivativeType, typename SubMatrixType = _SecondDerivativeType>
// SubMatrixType = SecondDerivativeType for SecondDerivativeType being normal matrices, but not for block-matrices
class TrustRegionMethod :
  public aol::Op<_VectorType,_VectorType> {

public:
  //typedef typename VectorType::DataType RealType;
  typedef _RealType RealType;
  typedef _VectorType VectorType;
  typedef _VectorType DerivativeType;
  typedef _SecondDerivativeType SecondDerivativeType;
  typedef typename ConcatDataTrait<VectorType>::ConcatDataType ConcatDataType;
  typedef typename ConcatDataTrait<VectorType>::SizeType SizeType;
  typedef enum _TrustRegionStatus {
    METHOD_HAS_NOT_RUN = -1,
    ACCURACY_CONDITION_MET = 0,
    EXCEEDED_MAXITERATIONS = 1,
    EXCEEDED_DELTABOUND = 2
  } TrustRegionStatus;

protected:
  // the energy to be minimized and its first and second derivatives
  const aol::Op<VectorType,aol::Scalar<RealType> > &_e;
  const aol::Op<VectorType,DerivativeType> &_de;
  const aol::Op<VectorType,SecondDerivativeType> &_d2e;

  const int _maxSteps;
  const RealType _maxAccuracy;
  const int _maxLanczosSteps;
  const RealType _eta1, _eta2, _gamma1, _gamma2, _theta, _kappa1, _kappa2;
  bool _iterativeMode, _quiet, _writeSteps;

  SecondDerivativeType * const _pMatrixTemplate;
  mutable cholmod_common _cholmodSetting;
  const int _dtype; // cholmod datatype (float or double); float not yet implemented in cholmod
  typedef double CholmodRealType; // if cholmod implements float, replace "double" by RealType

  // value of the function to be minimized
  mutable aol::Scalar<RealType> _f;
  
  mutable TrustRegionStatus _lastStatus;

public:
  inline int countMatrixEntries( const SubMatrixType &Matrix ) const {
    int result = 0;
    for ( int row = 0; row < Matrix.getNumRows(); row++ ) {
      vector<typename aol::Row<RealType>::RowEntry > matRow;
      Matrix.makeRowEntries( matRow, row );
      result += matRow.size();
    }
    return result;
  }

  inline int countMatrixEntries( const aol::BlockOpBase<RealType,SubMatrixType> &BlockMat ) const {
    int result = 0;
    for ( int row = 0; row < BlockMat.getNumRows(); row++ )
      for ( int col = 0; col < BlockMat.getNumCols(); col++ )
        result += countMatrixEntries( BlockMat.getReference( row, col ) );
    return result;
  }

  inline int getNumRows( const SubMatrixType &Matrix ) const {
    return Matrix.getNumRows();
  }

  inline int getNumRows( const aol::BlockOpBase<RealType,SubMatrixType> &BlockMat ) const {
    int result = 0;
    for ( int row = 0; row < BlockMat.getNumRows(); row++ )
      result += getNumRows( BlockMat.getReference( row, 0 ) );
    return result;
  }

  TrustRegionMethod( const aol::Op<VectorType,aol::Scalar<RealType> > &E,
                     const aol::Op<VectorType,DerivativeType> &DE,
                     const aol::Op<VectorType,SecondDerivativeType> &D2E,
                     SecondDerivativeType * const MatrixTemplate, // will be owned by the class
                     const int MaxSteps = 50,
                     const RealType MaxAccuracy = 1.e-6, // |gradient|~|energy difference|/|step|; with |step|~|gradient|
                                                         // and min|energy difference|~eps_machine we have max accuracy~|gradient|~\sqrt(eps_machine)
                     const int MaxLanczosSteps = 3, // default so small due to the instability of the Lanczos iteration
                     const RealType Eta1 = .1,
                     const RealType Eta2 = .9,
                     const RealType Gamma1 = 1. / 32.,
                     const RealType Gamma2 = 2.,
                     const RealType Theta = .01,
                     const RealType Kappa1 = .1,
                     const RealType Kappa2 = .2,
                     const bool IterativeMode = false,
                     const bool Quiet = false,
                     const bool WriteSteps = false ) :
    _e( E ),
    _de( DE ),
    _d2e( D2E ),
    _maxSteps( MaxSteps ),
    _maxAccuracy( MaxAccuracy ),
    _maxLanczosSteps( MaxLanczosSteps ),
    _eta1( Eta1 ),
    _eta2( Eta2 ),
    _gamma1( Gamma1 ),
    _gamma2( Gamma2 ),
    _theta( Theta ),
    _kappa1( Kappa1 ),
    _kappa2( Kappa2 ),
    _iterativeMode( IterativeMode ),
    _quiet( Quiet ),
    _writeSteps( WriteSteps ),
    _pMatrixTemplate( MatrixTemplate ),
    // float not yet implemented in cholmod; if it is implemented in the future replace first CHOLMOD_DOUBLE by CHOLMOD_SINGLE
    _dtype( RealTrait<RealType>::ALIAS == FLOAT ? CHOLMOD_DOUBLE : CHOLMOD_DOUBLE ),
    _lastStatus( METHOD_HAS_NOT_RUN ) {

    cholmod_start ( &_cholmodSetting );
    _cholmodSetting.final_ll = true; // always compute LL^T factorization instead of LDL^T
    /*_cholmodSetting.supernodal = CHOLMOD_SIMPLICIAL; // always compute simplicial factorization (very slow; better each time copy the factor and transform it)*/
    _cholmodSetting.print = 0; // do not print error messages (e.g. that H is not positive definite)
    /*// to switch off the ordering (strongly decreases sparsity of L and increases computation time, but in the original algorithm, the exact L is needed sometimes)
    _cholmodSetting.nmethods = 1;
    _cholmodSetting.method [0].ordering = CHOLMOD_NATURAL;
    _cholmodSetting.postorder = false;*/
    _cholmodSetting.dtype = _dtype;
  }

  ~TrustRegionMethod() {
    cholmod_finish ( &_cholmodSetting );
    delete _pMatrixTemplate;
  }
  
  void setWriteSteps( const bool WriteSteps = true ) {
    _writeSteps = WriteSteps;
  }

  void setQuietMode( const bool Quiet = true ) {
    _quiet = Quiet;
  }

  void setIterativeMode( const bool IterativeMode = true ) {
    _iterativeMode = IterativeMode;
  }

  //! returns the value of the function to be minimized
  RealType getFunctionValue() {
    return _f.v;
  }
  
  TrustRegionStatus getLastStatus() const {
    return _lastStatus;
  }

  template<typename InputVectorType>
  cholmod_dense* transformVectorIntoCholmodDense( const InputVectorType &Vector ) const {
    // if in the future cholmod also supports float, simply remove whole if-statement by the body of the else-part
    if ( RealTrait<RealType>::ALIAS != DOUBLE ) {
      cholmod_dense* result = cholmod_allocate_dense( Vector.size(), 1, Vector.size(), CHOLMOD_REAL, &_cholmodSetting );
      for ( int row = 0; row < Vector.size(); row++ )
        reinterpret_cast<CholmodRealType*>( result->x )[row] = Vector.get( row );
      return result;
    } else {
      // define all struct fields
      cholmod_dense* result = new cholmod_dense;
      result->nrow = Vector.size();
      result->ncol = 1;
      result->nzmax = Vector.size();
      result->d = Vector.size();
      result->x = Vector.getData();
      result->z = NULL;
      result->xtype = CHOLMOD_REAL;
      result->dtype = _dtype;
      return result;
    }
  }

  cholmod_sparse* transformMatrixIntoCholmodSparse( const SubMatrixType &Matrix ) const {
    cholmod_triplet* triplet = cholmod_allocate_triplet( Matrix.getNumRows(), Matrix.getNumCols(), countMatrixEntries( Matrix ), 1, CHOLMOD_REAL, &_cholmodSetting );
    int numberOfEntries = 0;
    for ( int row = 0; row < Matrix.getNumRows(); row++ ) {
      vector<typename aol::Row<RealType>::RowEntry > matRow;
      Matrix.makeRowEntries( matRow, row );
      for ( typename vector<typename aol::Row<RealType>::RowEntry >::iterator it = matRow.begin(); it != matRow.end(); ++it )
        if ( ( it->col <= row ) && ( it->value != aol::ZTrait<RealType>::zero ) ) { // only nonzero subdiagonal entries
          reinterpret_cast<int*>( triplet->i )[numberOfEntries] = row;
          reinterpret_cast<int*>( triplet->j )[numberOfEntries] = it->col;
          reinterpret_cast<CholmodRealType*>( triplet->x )[numberOfEntries] = it->value;
          numberOfEntries++;
        }
    }
    triplet->nnz = numberOfEntries;
    cholmod_sparse* result = cholmod_triplet_to_sparse( triplet, numberOfEntries, &_cholmodSetting );
    cholmod_free_triplet( &triplet, &_cholmodSetting );
    cholmod_sort( result, &_cholmodSetting );
    return result;
  }

  cholmod_sparse* transformMatrixIntoCholmodSparse( const TripletMatrix<RealType> &Matrix ) const {
    if ( RealTrait<RealType>::ALIAS != DOUBLE )
      throw aol::Exception ( "aol::TrustRegionMethod.transformMatrixIntoCholmodSparse< TripletMatrix<RealType> > not implemented for RealType = float", __FILE__, __LINE__ );
    else {
      cholmod_triplet* triplet = new cholmod_triplet;
      triplet->nrow  = Matrix.getNumRows();
      triplet->ncol  = Matrix.getNumCols();
      triplet->nzmax = Matrix.getValueReference().size();
      triplet->nnz   = Matrix.getValueReference().size(); // better way to do this?
      triplet->i     = Matrix.getRowIndexReference().getData();
      triplet->j     = Matrix.getColIndexReference().getData();
      triplet->x     = Matrix.getValueReference().getData();
      triplet->z     = NULL;
      triplet->stype = 0;  // potentially unsymmetric matrix; (1=upper triangular part is stored; -1=lower triangular part is stored)
      triplet->itype = CHOLMOD_INT;
      triplet->xtype = CHOLMOD_REAL;
      triplet->dtype = CHOLMOD_DOUBLE;
      cholmod_sparse* result = cholmod_triplet_to_sparse( triplet, triplet->nnz, &_cholmodSetting );
      delete triplet;
      cholmod_sort( result, &_cholmodSetting );
      result->stype = 1;  // only use upper triangular part and ignore lower triangular part, since Hessian has to be symmetric anyway
      return result;
    }
  }

  cholmod_sparse* transformMatrixIntoCholmodSparse( const aol::BlockOpBase<RealType,SubMatrixType> &BlockMat ) const {
    // only save the entries below the diagonal (if all entries were stored, the matrix would be treated unsymmetric, yielding bad results)
    const int n = getNumRows( BlockMat );
    cholmod_triplet* triplet = cholmod_allocate_triplet( n, n, countMatrixEntries( BlockMat ), 1, CHOLMOD_REAL, &_cholmodSetting );
    int numberOfEntries = 0;
    for ( int row = 0, totalRow = 0; row < BlockMat.getNumRows(); row++ ) {
      for ( int localRow = 0; localRow < BlockMat.getReference( row, 0 ).getNumRows(); localRow++, totalRow++ ) {
        int colOffset = 0;
        for ( int col = 0; col <= row; col++ ) {
          vector<typename aol::Row<RealType>::RowEntry > matRow;
          BlockMat.getReference( row, col ).makeRowEntries( matRow, localRow );
          for ( typename vector<typename aol::Row<RealType>::RowEntry >::iterator it = matRow.begin(); it != matRow.end(); ++it )
            if ( ( col < row || it->col <= localRow ) && ( it->value != aol::ZTrait<RealType>::zero ) ) { // only nonzero subdiagonal entries
              reinterpret_cast<int*>( triplet->i )[numberOfEntries] = totalRow;
              reinterpret_cast<int*>( triplet->j )[numberOfEntries] = colOffset + it->col;
              reinterpret_cast<CholmodRealType*>( triplet->x )[numberOfEntries] = it->value;
              numberOfEntries++;
            }
          colOffset += BlockMat.getReference( row, col ).getNumCols();
        }
      }
    }
    triplet->nnz = numberOfEntries;
    cholmod_sparse* result = cholmod_triplet_to_sparse( triplet, numberOfEntries, &_cholmodSetting );
    cholmod_free_triplet( &triplet, &_cholmodSetting );
    cholmod_sort( result, &_cholmodSetting );
    return result;
  }

  inline RealType dotProd( const cholmod_dense * const Mat1, const cholmod_dense * const Mat2 ) const {
    aol::Vector<CholmodRealType> mat1( reinterpret_cast<CholmodRealType*>( Mat1->x ), Mat1->ncol * Mat1->nrow, aol::FLAT_COPY );
    aol::Vector<CholmodRealType> mat2( reinterpret_cast<CholmodRealType*>( Mat2->x ), Mat2->ncol * Mat2->nrow, aol::FLAT_COPY );
    return mat1 * mat2;
  }

  inline void addMultiple( const cholmod_dense * const Mat1, const cholmod_dense * const Mat2, const RealType alpha ) const {
    aol::Vector<CholmodRealType> mat1( reinterpret_cast<CholmodRealType*>( Mat1->x ), Mat1->ncol * Mat1->nrow, aol::FLAT_COPY );
    aol::Vector<CholmodRealType> mat2( reinterpret_cast<CholmodRealType*>( Mat2->x ), Mat2->ncol * Mat2->nrow, aol::FLAT_COPY );
    mat1.addMultiple( mat2, alpha );
  }

  inline void multSparseDense( cholmod_sparse * const A, cholmod_dense * const X, cholmod_dense * const Y, const RealType Alpha, const RealType Beta ) const {
    // Y = alpha A X + beta Y
    double alpha[2] = { 0., 0. }, beta[2] = { 0., 0. };
    alpha[0] = Alpha;
    beta[0] = Beta;
    cholmod_sdmult( A, false, alpha, beta, X, Y, &_cholmodSetting ) ;
  }

  void computeGershgorinEstimates( const cholmod_sparse * const Matrix, RealType &MinDiagEntry, RealType &LowerEstimate, RealType &UpperEstimate ) const {
    // assume packed format!
    aol::Vector<RealType> diagEntries( Matrix->ncol );
    aol::Vector<RealType> offDiag1Norm( Matrix->ncol );
    for ( unsigned int j = 0; j < Matrix->ncol; j++ )
      for ( int i = reinterpret_cast<int*>( Matrix->p )[j]; i < reinterpret_cast<int*>( Matrix->p )[j+1]; i++ )
        if ( static_cast<int>( j ) == reinterpret_cast<int*>( Matrix->i )[i] ) // diagonal entry
          diagEntries[j] = reinterpret_cast<CholmodRealType*>( Matrix->x )[i];
        else { // offdiagonal entry
          offDiag1Norm[j] += aol::Abs( reinterpret_cast<CholmodRealType*>( Matrix->x )[i] );
          offDiag1Norm[reinterpret_cast<int*>( Matrix->i )[i]] += aol::Abs( reinterpret_cast<CholmodRealType*>( Matrix->x )[i] );
        }
    if ( Matrix->stype == 0 )
      offDiag1Norm *= .5;
    MinDiagEntry = diagEntries.getMinValue();
    diagEntries += offDiag1Norm;
    UpperEstimate = diagEntries.getMaxValue();
    diagEntries.addMultiple( offDiag1Norm, -2. );
    LowerEstimate = diagEntries.getMinValue();
  }

  RealType getFrobeniusNormSqr( const cholmod_sparse * const Matrix ) const {
    RealType result = 0.;
    // assumes packed format and that only half of the symmetric matrix is stored
    for ( unsigned int col = 0; col < Matrix->ncol; col++ )
      for ( int pos = reinterpret_cast<int*>( Matrix->p )[col]; pos < reinterpret_cast<int*>( Matrix->p )[col+1]; pos++ )
        if ( reinterpret_cast<int*>( Matrix->i )[pos] == static_cast<int>( col ) ) // diagonal entry
          result += aol::Sqr( reinterpret_cast<CholmodRealType*>( Matrix->x )[pos] );
        else // off-diagonal entry
          result += 2 * aol::Sqr( reinterpret_cast<CholmodRealType*>( Matrix->x )[pos] );
    return result;
  }

  cholmod_dense* linpackMethod( cholmod_factor* L ) const {
    // Conn, Gould, Toint p.191
    /*// first transpose the factor so that the rows of L become the columns of lT
    cholmod_factor* aux = cholmod_copy_factor( L, &_cholmodSetting );
    cholmod_sparse* l = cholmod_factor_to_sparse( aux, &_cholmodSetting );
    cholmod_free_factor( &aux, &_cholmodSetting );
    cholmod_sparse* lT = cholmod_transpose( l, 1, &_cholmodSetting );
    cholmod_free_sparse( &l, &_cholmodSetting );*/

    // convert L into row-major order
    // Lp, Li, and Lx have same interpretation as in cholmod, just in row-major order
    aol::Vector<int> entriesPerRow( L->n ), Lp( L->n + 1 );
    // first count nonzeros per row
    if ( !L->is_super )
      for ( unsigned int col = 0; col < L->n; col++ )
        for ( int pos = reinterpret_cast<int*>( L->p )[col]; pos < reinterpret_cast<int*>( L->p )[col] + reinterpret_cast<int*>( L->nz )[col]; pos++ )
          entriesPerRow[reinterpret_cast<int*>( L->i )[pos]]++;
    else
      for ( unsigned int superNode = 0; superNode < L->nsuper; superNode++ )
        for ( int col = reinterpret_cast<int*>( L->super )[superNode], valIndex = reinterpret_cast<int*>( L->px )[superNode]; col < reinterpret_cast<int*>( L->super )[superNode+1]; col++ )
          for ( int rowIndex = reinterpret_cast<int*>( L->pi )[superNode]; rowIndex < reinterpret_cast<int*>( L->pi )[superNode+1]; rowIndex++, valIndex++ )
            if ( reinterpret_cast<CholmodRealType*>( L->x )[valIndex] != 0. )
              entriesPerRow[reinterpret_cast<int*>( L->s )[rowIndex]]++;
    // define Lp
    for ( unsigned int k = 1; k <= L->n; k++ )
      Lp[k] = Lp[k-1] + entriesPerRow[k-1];
    // define Li and Lx
    int nzmax = entriesPerRow.sum();
    aol::Vector<int> Li( nzmax );
    aol::Vector<RealType> Lx( nzmax );
    entriesPerRow.setZero();
    if ( !L->is_super )
      for ( unsigned int col = 0; col < L->n; col++ )
        for ( int pos = reinterpret_cast<int*>( L->p )[col]; pos < reinterpret_cast<int*>( L->p )[col] + reinterpret_cast<int*>( L->nz )[col]; pos++ )  {
          int row = reinterpret_cast<int*>( L->i )[pos];
          Li[Lp[row]+entriesPerRow[row]] = col;
          Lx[Lp[row]+entriesPerRow[row]] = reinterpret_cast<CholmodRealType*>( L->x )[pos];
          entriesPerRow[row]++;
        }
    else
      for ( unsigned int superNode = 0; superNode < L->nsuper; superNode++ )
        for ( int col = reinterpret_cast<int*>( L->super )[superNode], valIndex = reinterpret_cast<int*>( L->px )[superNode]; col < reinterpret_cast<int*>( L->super )[superNode+1]; col++ )
          for ( int rowIndex = reinterpret_cast<int*>( L->pi )[superNode]; rowIndex < reinterpret_cast<int*>( L->pi )[superNode+1]; rowIndex++, valIndex++ )
            if ( reinterpret_cast<CholmodRealType*>( L->x )[valIndex] != 0. ) {
              int row = reinterpret_cast<int*>( L->s )[rowIndex];
              Li[Lp[row]+entriesPerRow[row]] = col;
              Lx[Lp[row]+entriesPerRow[row]] = reinterpret_cast<CholmodRealType*>( L->x )[valIndex];
              entriesPerRow[row]++;
            }

    // forward substitution for w to make L^{-1}v as large as possible
    cholmod_dense* w = cholmod_allocate_dense( L->n, 1, L->n, CHOLMOD_REAL, &_cholmodSetting );
    for ( unsigned int row = 0; row < L->n; row++ ) { // row k of L
      // compute \sum_{i=1}^{k-1} L_{ki}w_i
      RealType sum = 0.;
      for ( int i = Lp[row]; i < Lp[row+1] - 1; i++ ) // assuming the diagonal entry to be i = Lp[row+1]-1
        sum += Lx[i] * reinterpret_cast<CholmodRealType*>( w->x )[Li[i]];
      // set w_k
      if ( sum > 0. )
        reinterpret_cast<CholmodRealType*>( w->x )[row] = - ( 1. + sum ) / Lx[Lp[row+1]-1];
      else
        reinterpret_cast<CholmodRealType*>( w->x )[row] = ( 1. - sum ) / Lx[Lp[row+1]-1];
    }

    // copmute L^{-T}w/||L^{-T}w||
    /*cholmod_free_sparse( &lT, &_cholmodSetting );*/
    cholmod_dense* u = cholmod_solve( CHOLMOD_Lt, L, w, &_cholmodSetting );
    RealType uScaling = 1. / cholmod_norm_dense( u, 2, &_cholmodSetting );
    for ( unsigned int i = 0; i < u->nrow; i++ )
      reinterpret_cast<CholmodRealType*>( u->x )[i] *= uScaling;

    // permute u according to the cholmod-permutation
    for ( unsigned int i = 0; i < L->n; i++ )
      reinterpret_cast<CholmodRealType*>( w->x )[reinterpret_cast<int*>( L->Perm )[i]] = reinterpret_cast<CholmodRealType*>( u->x )[i];
    cholmod_free_dense( &u, &_cholmodSetting );
    return w;
  }

  RealType partialFactorizationEigenvalueBound( cholmod_factor* L, cholmod_sparse* H, RealType Lambda ) const {
    // Conn, Gould, Toint p.191

    // find (permuted) kth diagonal entry of H+\lambda I, for which the Cholesky factorization stopped
    int k = reinterpret_cast<int*>( L->Perm )[L->minor];
    RealType h_kk = 0.;
    int pos = 0;
    if ( H->sorted && H->packed ) {
      pos = reinterpret_cast<int*>( H->p )[k]-1;
      while ( reinterpret_cast<int*>( H->i )[++pos] < k ){}
    } else
      throw aol::Exception( "unexpected matrix type in partial factorization for trust region method", __FILE__, __LINE__ );
    if ( reinterpret_cast<int*>( H->i )[pos] == k )
      h_kk = reinterpret_cast<CholmodRealType*>( H->x )[pos] + Lambda;
    else // kth diagonal entry does not exist in H (i.e. it is 0 and is thus not stored)
      h_kk = Lambda;

    // compute delta
    RealType delta = 0.;
    if ( !L->is_super ) {
      for ( unsigned int col = 0; col < L->minor; col++ )
        for ( int pos = reinterpret_cast<int*>( L->p )[col]; pos < reinterpret_cast<int*>( L->p )[col] + reinterpret_cast<int*>( L->nz )[col]; pos++ )
          if ( reinterpret_cast<int*>( L->i )[pos] == static_cast<int>( L->minor ) )
            delta += aol::Sqr( reinterpret_cast<CholmodRealType*>( L->x )[pos] );
    } else
      for ( unsigned int superNode = 0; superNode < L->nsuper; superNode++ )
        for ( int col = reinterpret_cast<int*>( L->super )[superNode], valIndex = reinterpret_cast<int*>( L->px )[superNode]; col < reinterpret_cast<int*>( L->super )[superNode+1]; col++ )
          for ( int rowIndex = reinterpret_cast<int*>( L->pi )[superNode]; rowIndex < reinterpret_cast<int*>( L->pi )[superNode+1]; rowIndex++, valIndex++ )
            if ( ( reinterpret_cast<int*>( L->s )[rowIndex] == static_cast<int>( L->minor ) ) && ( col < static_cast<int>( L->minor ) ) )
              delta += aol::Sqr( reinterpret_cast<CholmodRealType*>( L->x )[valIndex] );
    delta -= h_kk;

    // replace L->minor'th and following diagonals by 1
    if ( !L->is_super )
      for ( unsigned int col = L->minor; col < L->n; col++ )
        reinterpret_cast<CholmodRealType*>( L->x )[reinterpret_cast<int*>( L->p )[col]] = 1.;
    else // in supernodal case simply replace all columns by 1 instead of just the diagonals, since it is not fully clear how the matrix is stored
      for ( unsigned int superNode = 0; superNode < L->nsuper; superNode++ )
        for ( int col = reinterpret_cast<int*>( L->super )[superNode], valIndex = reinterpret_cast<int*>( L->px )[superNode]; col < reinterpret_cast<int*>( L->super )[superNode+1]; col++ )
          for ( int rowIndex = reinterpret_cast<int*>( L->pi )[superNode]; rowIndex < reinterpret_cast<int*>( L->pi )[superNode+1]; rowIndex++, valIndex++ )
            if ( col >= static_cast<int>( L->minor ) )
              if ( reinterpret_cast<int*>( L->s )[rowIndex] >= col )
                reinterpret_cast<CholmodRealType*>( L->x )[valIndex] = 1.;

    // compute v
    cholmod_dense* rhs = cholmod_zeros( L->n, 1, CHOLMOD_REAL, &_cholmodSetting );
    reinterpret_cast<CholmodRealType*>( rhs->x )[L->minor] = 1.;
    cholmod_dense* v = cholmod_solve( CHOLMOD_Lt, L, rhs, &_cholmodSetting );
    cholmod_free_dense( &rhs, &_cholmodSetting );

    // // test correctness of computation
    // cholmod_dense* aux = cholmod_allocate_dense( L->n, 1, L->n, CHOLMOD_REAL, &_cholmodSetting );
    // for ( unsigned int i = 0; i < L->n; i++ )
    //   reinterpret_cast<CholmodRealType*>( aux->x )[reinterpret_cast<int*>( L->Perm )[i]] = reinterpret_cast<CholmodRealType*>( v->x )[i];
    // for ( unsigned int i = 0; i < L->n; i++ )
    //   reinterpret_cast<CholmodRealType*>( v->x )[i] = reinterpret_cast<CholmodRealType*>( aux->x )[i];
    // multSparseDense( H, v, aux, 1., Lambda );
    // reinterpret_cast<CholmodRealType*>( aux->x )[k] += delta * reinterpret_cast<CholmodRealType*>( v->x )[k];
    // if ( aol::Abs( dotProd( aux, v ) ) > 1.e-15 )
    //   throw aol::Exception( "incorrect computation of eigenvalue estimate for trust region method" );
    // cholmod_free_dense( &aux, &_cholmodSetting );

    // return new lower bound on - negative eigenvalue
    RealType vNorm = cholmod_norm_dense( v, 2, &_cholmodSetting );
    cholmod_free_dense( &v, &_cholmodSetting );
    return Lambda + delta / aol::Sqr( vNorm );
  }

  void solveSubProblem( const ConcatDataType &Gradient, const SecondDerivativeType &Hessian, const RealType Delta, ConcatDataType &Step ) const {
    // transfer the data into cholmod format
    cholmod_dense* g = transformVectorIntoCholmodDense( Gradient );
    cholmod_dense* s = transformVectorIntoCholmodDense( Step );
    cholmod_sparse* H = transformMatrixIntoCholmodSparse( Hessian );
    solveSubProblem( g, H, Delta, s );

    if ( RealTrait<RealType>::ALIAS != DOUBLE ) // cholmod currently only supports double, hence copy s back (is unnecessary else)
      for ( int i = 0; i < Step.size(); i++ )
        Step.set( i, reinterpret_cast<CholmodRealType*>( s->x )[i] );

    // free the dynamic variables (not via cholmod_free_dense, since this would delete the data)
    delete g;
    delete s;
    cholmod_free_sparse( &H, &_cholmodSetting );
  }

  void solveSubProblem( cholmod_dense * const Gradient, cholmod_sparse * const Hessian, const RealType Delta, cholmod_dense * const Step ) const {
    // algorithm 7.3.4-5 in Conn, Gould, Toint
    cholmod_dense * const g = Gradient;
    cholmod_dense * const s = Step;
    cholmod_sparse * const H = Hessian;
    bool positiveDefinite, interiorStep; // the sets N, L, G correspond to {false,-}, {true,false}, {true,true}
    RealType stepNorm; // ||s||_2
    cholmod_sparse* identity = cholmod_speye( H->nrow, H->ncol, CHOLMOD_REAL, &_cholmodSetting );
    identity->stype = H->stype;
    cholmod_factor* choleskyFactor = NULL;

    // initialize lambda and lower and upper bounds for it
    RealType minDiagEntry, lambdaMin, lambdaMax;
    computeGershgorinEstimates( H, minDiagEntry, lambdaMin, lambdaMax );
    RealType gNorm = cholmod_norm_dense( g, 2, &_cholmodSetting );
    RealType H2Norm = sqrt( getFrobeniusNormSqr( H ) );
    RealType HinfNorm = cholmod_norm_sparse( H, 0, &_cholmodSetting );
    RealType lambdaL = aol::Max( aol::ZTrait<RealType>::zero, aol::Max( -minDiagEntry, gNorm / Delta - aol::Min( lambdaMax, aol::Min( H2Norm, HinfNorm ) ) ) );
    RealType lambdaU = aol::Max( aol::ZTrait<RealType>::zero, gNorm / Delta + aol::Min( -lambdaMin, aol::Min( H2Norm, HinfNorm ) ) );
    RealType lambda = aol::Max( aol::ZTrait<RealType>::zero, lambdaL );

#ifdef VERBOSE
    cerr<<endl;
#endif

    do {
#ifdef VERBOSE
      cerr<<"         lambda: "<<lambda;
#endif

      // factorize H+lambda*I = L L^T
      double one[2] = { 1., 0. }, lam[2] = { 0., 0. };
      lam[0] = lambda;
      cholmod_sparse* H_lambda = cholmod_add( H, identity, one, lam, true, false, &_cholmodSetting );
      if ( choleskyFactor == NULL ) // since sparsity structure stays the same!
        choleskyFactor = cholmod_analyze( H_lambda, &_cholmodSetting );
      cholmod_factorize( H_lambda, choleskyFactor, &_cholmodSetting );
      positiveDefinite = ( _cholmodSetting.status == CHOLMOD_OK );
      cholmod_free_sparse( &H_lambda, &_cholmodSetting );
#ifdef VERBOSE
      cerr<<"      pos def: "<<positiveDefinite;
#endif

      if ( positiveDefinite ) { // lambda \in L \cup G
        // solve L L^T s = -g
        cholmod_dense* aux = cholmod_solve( CHOLMOD_A, choleskyFactor, g, &_cholmodSetting );
        for ( unsigned int i = 0; i < aux->nzmax; i++ )
          reinterpret_cast<CholmodRealType*>( s->x )[i] = - reinterpret_cast<CholmodRealType*>( aux->x )[i];
        cholmod_free_dense( &aux, &_cholmodSetting );

        // check whether the minimum lies inside the trust region
        stepNorm = cholmod_norm_dense( s, 2, &_cholmodSetting );
        interiorStep = ( stepNorm < Delta );
        if ( ( interiorStep && ( lambda == 0. ) ) // Hessian positive definite and Newton step inside trust region
          || ( aol::Abs( stepNorm - Delta ) <= _kappa1 * Delta ) ) // optimal step lies at the boundary and has sufficiently been reached
          break;
#ifdef VERBOSE
        cerr<<"      int step: "<<interiorStep;
#endif

        // update upper and lower bounds for lambda
        if ( interiorStep )
          lambdaU = lambda;
        else
          lambdaL = lambda; // not really necessary, since from now on we have monotone convergence of lambda anyway

        // solve -L L^T s' = s and compute -s'\cdot s (instead of Lw=s, which would prohibit reordering in the Cholesky factorization (unless one takes the permutation into account))
        cholmod_dense* sPrimeMinus = cholmod_solve( CHOLMOD_A, choleskyFactor, s, &_cholmodSetting );
        RealType sTimessPrimeMinus = 0.;
        for ( unsigned int i = 0; i < s->nrow * s->ncol; i++ )
          sTimessPrimeMinus += reinterpret_cast<CholmodRealType*>( s->x )[i] * reinterpret_cast<CholmodRealType*>( sPrimeMinus->x )[i];
        cholmod_free_dense( &sPrimeMinus, &_cholmodSetting );

        // Linpack method
        if ( interiorStep ) {
          cholmod_dense* u = linpackMethod( choleskyFactor );
          cholmod_dense* aux = cholmod_copy_dense( u, &_cholmodSetting );
          multSparseDense( H, u, aux, 1., lambda );
          RealType uHu = dotProd( u, aux ), su = dotProd( s, u );
          multSparseDense( H, s, aux, 1., 0. );
          addMultiple( aux, s, lambda );
          RealType sHs = dotProd( s, aux );
          lambdaL = aol::Max( lambdaL, lambda - uHu );
          RealType alphas[2], qs[2];
          alphas[0] = - su + sqrt( aol::Sqr( su ) - aol::Sqr( stepNorm ) + aol::Sqr( Delta ) );
          alphas[1] = -alphas[0] - 2 * su;
          for ( int i = 0; i < 2; i++ ) {
            cholmod_dense* sNew = cholmod_copy_dense( s, &_cholmodSetting );
            addMultiple( sNew, u, alphas[i] );
            multSparseDense( H, sNew, aux, 1., 0. );
            qs[i] = .5 * dotProd( aux, sNew ) + dotProd( sNew, g );
            cholmod_free_dense( &sNew, &_cholmodSetting );
          }
          cholmod_free_dense( &aux, &_cholmodSetting );
          int index = ( qs[0] > qs[1] );
          addMultiple( s, u, alphas[index] );
          cholmod_free_dense( &u, &_cholmodSetting );
          if ( aol::Sqr( alphas[index] ) * uHu <= _kappa2 * ( sHs + lambda * aol::Sqr( Delta ) ) )
            break;
        }

        // lambda = lambda + (||s||-Delta)/Delta*||s||^2/(-s'\cdot s)  (-s'\cdot s replaced ||w||^2, since L is reordered)
        lambda += ( stepNorm - Delta ) / Delta * aol::Sqr( stepNorm ) / sTimessPrimeMinus;
        lambda = aol::Max( lambdaL, lambda ); // in case lambda \in G before

      } else { // lambda \in N
        // update value and lower bound of lambda
        lambdaL = aol::Max( lambdaL, partialFactorizationEigenvalueBound( choleskyFactor, H, lambda ) );
        lambdaL = aol::Max( lambdaL, lambda );
        lambda = aol::Max( sqrt( lambdaL * lambdaU ), lambdaL + _theta * ( lambdaU - lambdaL ) );
        // If H+lambda*I is positive semi-definite (i.e. the model minimizer lies at the trust region boundary),
        // then lambda is exactly -lambda1. If g is orthogonal to the zero-eigenvector ("hard case"), this means that
        // g s + s^T (H + lambda I) s has a whole line of minimizers, going right through the trust region,
        // parallel to the zero-eigenvector. If we simply increase lambda a little, we will have lambda \in G. The
        // minimizing step s will then be found near the middle of that line, and lambda will be identified as too
        // large and will be decreased successively, converging to -lambda1, while s converges to the middle of the line.
        // The termination methods for lambda \in G will then become active, so that we do not have to do anything here.
      }
#ifdef VERBOSE
      cerr<<endl;
#endif

      // termination is done via breaks
    } while ( true );

    cholmod_free_sparse( &identity, &_cholmodSetting );
    cholmod_free_factor( &choleskyFactor, &_cholmodSetting );

#ifdef VERBOSE
    cerr<<endl;
#endif
  }

  void solveSubProblemLanczos( const ConcatDataType &Gradient, const SecondDerivativeType &Hessian, const RealType Delta, ConcatDataType &Step ) const {
    // algorithm 7.5.2 in Conn, Gould, Toint
    VectorType s_k( Step, aol::FLAT_COPY );
    s_k.setZero();

    _PrecondType MInv( Hessian );

    ConcatDataType g_k( Gradient, aol::DEEP_COPY ),
                   v_k( g_k, aol::STRUCT_COPY );
    MInv.apply( g_k, v_k );

    RealType gv = g_k * v_k,
             gamma0 = sqrt( gv ),
             alpha_k = 1.,
             beta_k = 0.,
             s_k_MNormSqr = 0.,
             p_k_MNormSqr = 0.,
             s_k_M_p_k = 0.;

    ConcatDataType p_k( v_k, aol::DEEP_COPY ),
                   H_p_k( v_k, aol::STRUCT_COPY );
    p_k *= -1.;

    aol::MultiVector<RealType> T_k( 2, Gradient.size() ); // first component contains diagonal, second subdiagonal
    aol::Vector<RealType> h_k( Gradient.size() );

#ifdef VERBOSE
    cerr<<endl;
#endif

    bool interior = true;
    int k = 0;
    while ( ( !interior || ( /*gv*/ g_k.normSqr() > aol::Sqr( _maxAccuracy ) ) )
          && ( interior || ( k == 0 ) || ( T_k[1][k-1] * aol::Abs( h_k[k-1] ) > _maxAccuracy ) ) && ( k < _maxLanczosSteps ) ) {
    // gv or g_k.normSqr() is neither monotonic, nor do they decrease fast (in fact, very slowly)
    // we use g_k.normSqr() instead of the proposed gv since this is also the stopping criterion in the outer problem;
    // the subproblem has to reach an accuracy at least as high as the accuracy for which the outer problem terminates!
    // more than k = 3 Lanczos steps make no sense, since the solution of the subproblem might become positive due to the instability of Lanczos

#ifdef VERBOSE
      cerr<<"      step "<<k<<";  interior: "<<interior<<";  gradient-norm: "<<gv<<endl;
#endif

      Hessian.apply( p_k, H_p_k );
      RealType alpha_k_minus_1 = alpha_k;
      alpha_k = gv / ( p_k * H_p_k );

      T_k[0][k] = ( 1. / alpha_k ) + ( beta_k / alpha_k_minus_1 );

      if ( interior ) {
        s_k.addMultiple( p_k, alpha_k );
        // the M-norm of s_k is computed according to p.206 in Conn, Gould, Toint
        s_k_M_p_k = beta_k * ( s_k_M_p_k + alpha_k_minus_1 * p_k_MNormSqr );
        p_k_MNormSqr = gv + p_k_MNormSqr * aol::Sqr( beta_k );
        s_k_MNormSqr += 2 * alpha_k * s_k_M_p_k + aol::Sqr( alpha_k ) * p_k_MNormSqr;
        interior = false;//( alpha_k > 0. ) && ( s_k_MNormSqr <= aol::Sqr( Delta ) );
      }

      if ( !interior ) {
        cholmod_dense* grad = cholmod_zeros( k + 1, 1, CHOLMOD_REAL, &_cholmodSetting );
        reinterpret_cast<CholmodRealType*>( grad->x )[0] = gamma0;

        cholmod_dense* step = transformVectorIntoCholmodDense( h_k );
        step->nrow = k + 1;
        step->nzmax = k + 1;

        cholmod_triplet* triplet = cholmod_allocate_triplet( k + 1, k + 1, 2 * k + 1, 1, CHOLMOD_REAL, &_cholmodSetting );
        for ( int col = 0, index = 0; col < k; col++ ) {
          reinterpret_cast<int*>( triplet->i )[index] = col;
          reinterpret_cast<int*>( triplet->j )[index] = col;
          reinterpret_cast<CholmodRealType*>( triplet->x )[index++] = T_k[0][col];
          reinterpret_cast<int*>( triplet->i )[index] = col + 1;
          reinterpret_cast<int*>( triplet->j )[index] = col;
          reinterpret_cast<CholmodRealType*>( triplet->x )[index++] = T_k[1][col];
        }
        reinterpret_cast<int*>( triplet->i )[2*k] = k;
        reinterpret_cast<int*>( triplet->j )[2*k] = k;
        reinterpret_cast<CholmodRealType*>( triplet->x )[2*k] = T_k[0][k];
        triplet->nnz = 2 * k + 1;
        cholmod_sparse* hessian = cholmod_triplet_to_sparse( triplet, 2 * k + 1, &_cholmodSetting );
        cholmod_free_triplet( &triplet, &_cholmodSetting );

        solveSubProblem( grad, hessian, Delta, step );

        cholmod_free_dense( &grad, &_cholmodSetting );
        cholmod_free_sparse( &hessian, &_cholmodSetting );

        if ( RealTrait<RealType>::ALIAS != DOUBLE ) // cholmod currently only supports double, hence copy step back (is unnecessary else)
          for ( int i = 0; i < h_k.size(); i++ )
            h_k.set( i, reinterpret_cast<CholmodRealType*>( step->x )[i] );
        delete step;
      }

      RealType gvOld = gv;
      g_k.addMultiple( H_p_k, alpha_k );
      MInv.apply( g_k, v_k );
      gv = g_k * v_k;

      beta_k = gv / gvOld;
      p_k *= beta_k;
      p_k -= v_k;

      T_k[1][k] = sqrt( beta_k ) / aol::Abs( alpha_k );

      k++;
    }

    if ( !interior ) {
      // apply algorithm 5.2.4 in Conn, Gould, Toint to recover Q_k
      ConcatDataType t_k( Gradient, aol::DEEP_COPY );
      ConcatDataType w_k( Step, aol::STRUCT_COPY ),
                     w_k_minus_1( Step, aol::STRUCT_COPY ),
                     q_k( Step, aol::STRUCT_COPY );
      s_k.setZero();

      for ( int j = 0; j < k; j++ ) {
        // until it is scaled by gamma_k, q_k is called y_k in Conn, Gould, Toint
        // note that the Lanczos algorithm is not very stable, so the q_k become non-M-orthonormal very quickly
        MInv.apply( t_k, q_k );

        RealType gamma_k = sqrt( t_k * q_k );

        w_k_minus_1 = w_k;

        w_k = t_k;
        w_k /= gamma_k;

        q_k /= gamma_k;

        Hessian.apply( q_k, t_k );
        RealType delta_k = q_k * t_k;

        t_k.addMultiple( w_k, -delta_k );
        t_k.addMultiple( w_k_minus_1, -gamma_k );

        s_k.addMultiple( q_k, h_k[j] );
      }
    }
  }
  
  virtual void writeUpdate( const VectorType &/*x*/, const VectorType &/*gradient*/, const SecondDerivativeType &/*hessian*/ ) const { }

  virtual void apply( const VectorType &Arg, VectorType &Dest ) const {
    // algorithms 6.1.1 and 7.3.1-4 in Conn, Gould, Toint
    // notation: x are the iterates, s are the steps, g and H are the derivative and the Hessian, Delta the trust region radius, f the energy values
    aol::Scalar<RealType> fNew;
    VectorType &x = Dest;
    aol::MultiVector<RealType> arg;
    arg.appendReference( Arg );
    aol::Vector<RealType> sVec( arg.getTotalSize() ), gVec( sVec, aol::STRUCT_COPY );
    SizeType argsize; 
    ConcatDataTrait<VectorType>::getSize( Arg, argsize );
    ConcatDataType s( sVec.getData(), argsize, aol::FLAT_COPY );
    ConcatDataType g( gVec.getData(), argsize, aol::FLAT_COPY );
    ConcatDataType Hs( g, aol::STRUCT_COPY );
    SecondDerivativeType &H = *_pMatrixTemplate;

    // initialization
    x = Arg;
    int k = 0;
    RealType Delta = 1024.;

    // define the (quadratic) model
    _e.apply( x, _f );
    _de.apply( x, g );
    _d2e.apply( x, H );

    // prepare non-monotone Armijo condition (sometimes, in the exact minimum the numerical value is higher than at the last step;
    // hence rho is formed in each step by comparison with the maximum of the last two steps)
    aol::RingBuffer<RealType> fOld( 2 );
    fOld.push_front( _f.v );
    fOld.push_front( _f.v + 1e-12 * aol::Abs( _f.v ) );

    while ( ( k < _maxSteps ) && ( gVec.norm() > _maxAccuracy ) && ( Delta > 1.e-9 ) ) {
      if ( !_quiet )
        cerr << aol::strprintf ( "step %4d: Delta %8.3g, grad-norm %10.5g, func-value %6g   \r", k, static_cast < double > ( Delta ), static_cast < double > ( g.norm() ), static_cast < double > ( _f.v ) );
      // static_cast < double > is to prevent warnings when using long double

      // solve the trust region subproblem
      if ( _iterativeMode )
        solveSubProblemLanczos( g, H, Delta, s );
      else
        solveSubProblem( g, H, Delta, s );

      // compute energy and model decrease
      x += s;
      _e.apply( x, fNew );
      H.apply( s, Hs );
      RealType rho = ( fNew.v - fOld.getMaxValue()/*_f.v*/ ) / ( g * s + .5 * ( s * Hs ) );
      if ( aol::isNaN( rho ) || aol::isInf( rho ) )
        rho = 0;

      // update the iterate and the trust region
      if ( rho < _eta1 ) {  // unsuccessful
        x -= s;
        // decrease trust region radius; careful: in iterative mode, this is radius in preconditioned metric
        Delta *= _gamma1;
      } else {              // successful
        // update the (quadratic) model
        _f = fNew;
        fOld.push_front( _f.v );
        if ( ++k < _maxSteps ) {
          _de.apply( x, g );
          _d2e.apply( x, H );
          if ( _writeSteps )
            writeUpdate( x, g, H );
          if ( ( rho >= _eta2 ) && ( Delta < 1.e7 ) ) // very successful
            Delta *= _gamma2;
        }
      }
    } // while: terminate if desired accuracy reached
    
    if ( gVec.norm() <= _maxAccuracy ) {
      _lastStatus = ACCURACY_CONDITION_MET;
    } else if ( Delta <= 1.e-9 ) {
      _lastStatus = EXCEEDED_DELTABOUND;
    } else {
      _lastStatus = EXCEEDED_MAXITERATIONS;
    }

    if ( !_quiet )
      cerr << aol::strprintf ( "step %4d: Delta %8.3g, grad-norm %10.5g, func-value %6g   \n", k, static_cast < double > ( Delta ), static_cast < double > ( g.norm() ), static_cast < double > ( _f.v ) );
    // static_cast < double > is to prevent warnings when using long double
  }

  virtual void applyAdd( const VectorType &Arg, VectorType &Dest ) const {
    VectorType result ( Dest, aol::STRUCT_COPY );
    apply ( Arg, result );
    Dest += result;
  }

  virtual void applySingle( VectorType &ArgDest ) const {
    VectorType result ( ArgDest, aol::STRUCT_COPY );
    apply ( ArgDest, result );
    ArgDest = result;
  }
};

#else

template<typename _RealType, typename _VectorType, typename _SecondDerivativeType, typename _PrecondType = _SecondDerivativeType, typename SubMatrixType = _SecondDerivativeType>
class TrustRegionMethod : public aol::Op<_VectorType,_VectorType> {
public:
  typedef _RealType RealType;
  typedef _VectorType VectorType;
  typedef _VectorType DerivativeType;
  typedef _SecondDerivativeType SecondDerivativeType;

  TrustRegionMethod( const aol::Op<VectorType,aol::Scalar<RealType> > &/*E*/,
                     const aol::Op<VectorType,DerivativeType> &/*DE*/,
                     const aol::Op<VectorType,SecondDerivativeType> &/*D2E*/,
                     SecondDerivativeType * const /*MatrixTemplate*/,
                     const int /*MaxSteps*/ = 50,
                     const RealType /*MaxAccuracy*/ = 1.e-6,
                     const int /*MaxLanczosSteps*/ = 3, // default so small due to the instability of the Lanczos iteration
                     const RealType /*Eta1*/ = .1,
                     const RealType /*Eta2*/ = .9,
                     const RealType /*Gamma1*/ = 1. / 32.,
                     const RealType /*Gamma2*/ = 2.,
                     const RealType /*Theta*/ = .01,
                     const RealType /*Kappa1*/ = .1,
                     const RealType /*Kappa2*/ = .2,
                     const bool /*IterativeMode*/ = false,
                     const bool /*Quiet*/ = false ) { }

  void setQuietMode( const bool /*Quiet*/ = true ) {
    throw aol::Exception ( "aol::TrustRegionMethod needs suitesparse! Compile with -DUSE_SUITESPARSE", __FILE__, __LINE__ );
  }

  void setIterativeMode( const bool /*IterativeMode*/ = true ) {
    throw aol::Exception ( "aol::TrustRegionMethod needs suitesparse! Compile with -DUSE_SUITESPARSE", __FILE__, __LINE__ );
  }

  virtual void applyAdd( const VectorType &/*Arg*/, VectorType &/*Dest*/ ) const {
    throw aol::Exception ( "aol::TrustRegionMethod needs suitesparse! Compile with -DUSE_SUITESPARSE", __FILE__, __LINE__ );
  }

  virtual void applySingle( VectorType &/*ArgDest*/ ) const {
    throw aol::Exception ( "aol::TrustRegionMethod needs suitesparse! Compile with -DUSE_SUITESPARSE", __FILE__, __LINE__ );
  }
};

#endif

} // namespace aol

#endif
