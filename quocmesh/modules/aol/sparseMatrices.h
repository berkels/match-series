#ifndef __SPARSEMATRICES_H
#define __SPARSEMATRICES_H

#include <rectangularGrid.h>
#include <simplexGrid.h>
#include <rows.h>
#include <matrix.h>
#include <quoc.h>
#include <sparseMatrixRowIterator.h>
#ifdef USE_EXTERNAL_GMM
#include <gmmIncludes.h>
#endif

namespace aol {

template <class T> class Matrix;

/** Basis class for sparse matrix classes that are organized row-wise
 *  Note that pointers to rows are stored and that NULL pointers mean identity rows. It is not clear whether this matrix always behaves as one would expect!
 */
template <typename DataType>
class GenSparseMatrix : public Matrix<DataType> {

  GenSparseMatrix() {}
  // Copying must be done in subclass where Row type is known
  explicit GenSparseMatrix ( const GenSparseMatrix& /*mat*/ ) {}
  // Assignment must be done in subclass where Row type is known
  GenSparseMatrix& operator= ( const GenSparseMatrix& /*mat*/ ) {}

public:
  typedef BitVector MaskType;

  GenSparseMatrix ( int Rows, int Columns )
      : Matrix<DataType> ( Rows, Columns ), rows ( Rows ),
      _diagEntry ( aol::ZOTrait<DataType>::one ) {
    for ( int i = 0; i < this->getNumRows(); ++i )
      rows[i] = NULL;
  }

  template <typename GridType>
  explicit GenSparseMatrix ( const GridType &Grid )
      : Matrix<DataType> ( Grid.getNumberOfNodes(), Grid.getNumberOfNodes() ),
      rows ( Grid.getNumberOfNodes() ),
      _diagEntry ( aol::ZOTrait<DataType>::one ) {
    init( );
  }

  //! resize Matrix, deleting old contents
  void reallocate ( const int Rows, const int Columns ) {
    int i;
    for ( i = 0; i < this->_numRows; ++i ) {
      if ( rows[i] ) {
        delete rows[i];
        rows[i] = NULL;
      }
    }
    this->setZero();

    rows.resize ( Rows );
    for ( i = this->_numRows; i < Rows; ++i ) {
      rows[i] = NULL;
    }
    this->_numRows = Rows;
    this->_numCols = Columns;
  }

  virtual ~GenSparseMatrix() {}

  virtual Row<DataType>* newDefaultRow() const = 0;

  // *** other methods ***
  virtual DataType get ( int I, int J ) const {
#ifdef BOUNDS_CHECK
      if ( I < 0 || J < 0 || I >= static_cast<int> ( this->getNumRows() ) || J >= static_cast<int> ( this->getNumCols() ) ) {
        cerr << I << " " << J << " is out of bounds: " << this->getNumRows() << " " << this->getNumCols() << endl;
        throw aol::Exception ( "aol::GenSparseMatrix::get: Index out of bounds", __FILE__, __LINE__ );
      }
#endif
      if ( rows[I] ) {
        return ( rows[I]->get ( I, J ) );
      } else {
        if ( I == J ) {
          return _diagEntry;
        } else {
          return aol::ZOTrait<DataType>::zero;
        }
      }
    }

  virtual void set ( int I, int J, DataType Value ) {
#ifdef BOUNDS_CHECK
    if ( I < 0 || J < 0 || I >= static_cast<int> ( this->getNumRows() ) || J >= static_cast<int> ( this->getNumCols() ) ) {
      cerr << I << " " << J << " is out of bounds: " << this->getNumRows() << " " << this->getNumCols() << endl;
      throw aol::Exception ( "aol::GenSparseMatrix::set: Index out of bounds", __FILE__, __LINE__ );
    }
#endif
    if ( rows[I] )
      rows[I]->set ( I, J, Value );
    else if ( ( ( I == J ) && Value == _diagEntry ) || Value == aol::NumberTrait<DataType>::zero ) {
      // do nothing
#ifdef VERBOSE
      cerr << "aol::GenSparseMatrix<DataType>::set: setting non-existent entry (" << I << ", " << J << ") to "
      << Value << endl;
#endif
    } else
      throw Exception ( "Row does not exist.", __FILE__, __LINE__ );
  }

  virtual void add ( int I, int J, DataType Value ) {
#ifdef BOUNDS_CHECK
    if ( I < 0 || J < 0 || I >= static_cast<int> ( this->getNumRows() ) || J >= static_cast<int> ( this->getNumCols() ) ) {
      cerr << I << " " << J << " is out of bounds: " << this->getNumRows() << " " << this->getNumCols() << endl;
      throw aol::Exception ( "aol::GenSparseMatrix::add: Index out of bounds", __FILE__, __LINE__ );
    }
#endif
    if ( rows[I] ) {
      rows[I]->add ( I, J, Value );
    } else if ( Value == aol::NumberTrait<DataType>::zero ) {
      // do nothing
#ifdef VERBOSE
      cerr << "aol::GenSparseMatrix<DataType>::add: adding zero to diagonal entry" << endl;
#endif
    } else {
      cerr << rows[I] << " " << I << endl;
      throw Exception ( "Row does not exist.", __FILE__, __LINE__ );
    }
  }

  //! Adds vec1 \f$ \otimes \f$ vec2 to this.
  //! \author Toelkes
  inline virtual void addTensorProduct ( const Vector < DataType > &vec1, const Vector < DataType > &vec2 ) {
    addTensorProductMultiple ( vec1, vec2, 1.0 );
  }

  //! Adds factor * vec1 \f$ \otimes \f$ vec2 to this.
  //! \author Toelkes
  virtual void addTensorProductMultiple ( const Vector < DataType > &vec1, const Vector < DataType > &vec2, const DataType factor ) {
    for ( int i = 0; i < this->getNumRows (); ++i ) {
      for ( int j = 0; j < this->getNumCols (); ++j ) {
        this->add ( i, j, factor * vec1[i] * vec2[j] );
      }
    }
  }

  //! Optimized version of sparse matrix-vector-multiplication.
  /*! warning default behaviour: if there is no row at given index, apply assumes _diagEntry on the diagonal. */
  void applyAdd ( const Vector<DataType> &Arg, Vector<DataType> &Dest ) const {
    if ( this->getNumRows() != Dest.size() || this->getNumCols() != Arg.size() ) {
      string msg = strprintf ( "aol::GenSparseMatrix::applyAdd: Cannot applyAdd %d by %d matrix from vector of size %d to vector of size %d.", this->getNumRows(), this->getNumCols(), Arg.size(), Dest.size() );
      throw ( Exception ( msg, __FILE__, __LINE__ ) );
    }

    for ( int i = 0; i <  this->getNumRows() ; ++i )
      if ( rows[i] )
        Dest[i] += rows[i]->mult ( Arg, i );
      else
        Dest [i] += _diagEntry * Arg[i];
  }

  void applyAddMultiple( const Vector<DataType> &Arg, Vector<DataType> &Dest, const DataType Factor = aol::NumberTrait<DataType>::one ) const {
    for ( int i = 0; i <  this->getNumRows() ; ++i )
      if ( rows[i] )
        Dest[i] += Factor * rows[i]->mult ( Arg, i );
      else
        Dest[i] += Factor * _diagEntry * Arg[i];
  }

  //! Matrix-vector multiplication with masking functionality.
  //! Differently from apply and applyAdd, this function is not re-implemented
  //! in subclasses.
  void applyAddMasked ( const Vector<DataType> &Arg, Vector<DataType> &Dest,
                        const BitVector & Mask, IncludeWriteMode applyMode ) const {

    // check dimensions
    if ( this->getNumRows() != Dest.size() || this->getNumCols() != Arg.size() ) {
      char errmsg [ 1024 ];
      sprintf ( errmsg, "aol::GenSparseMatrix::apply: Cannot apply %d by %d matrix from vector of size %d to vector of size %d.", this->getNumRows(), this->getNumCols(), Arg.size(), Dest.size() );
      throw ( Exception ( errmsg, __FILE__, __LINE__ ) );
    }

    switch ( applyMode ) {
    case INCLUDE_ALL_WRITE_ALL:
      applyAddMasked<BitMaskFunctorTrue, &Row<DataType>::multMaskedFunctorTrue > ( Arg, Dest, Mask );
      break;
    case INCLUDE_BD_WRITE_INT:
      applyAddMasked<BitMaskFunctorIdentity, &Row<DataType>::multMaskedFunctorNegate > ( Arg, Dest, Mask );
      break;
    case INCLUDE_INT_WRITE_ALL:
      applyAddMasked<BitMaskFunctorTrue, &Row<DataType>::multMaskedFunctorIdentity > ( Arg, Dest, Mask );
      break;
    case INCLUDE_ALL_WRITE_INT:
      applyAddMasked<BitMaskFunctorIdentity, &Row<DataType>::multMaskedFunctorTrue > ( Arg, Dest, Mask );
      break;
    case INCLUDE_INT_WRITE_INT:
      applyAddMasked<BitMaskFunctorIdentity, &Row<DataType>::multMaskedFunctorIdentity > ( Arg, Dest, Mask );
      break;
    default:
      throw aol::Exception ( "aol::GenSparseMatrix::applyAddMasked: unknown IncludeWriteMode", __FILE__, __LINE__ );
    }
  }

  template <typename BitMaskFunctorType, DataType ( Row<DataType>::* multMasked ) ( const Vector<DataType> &, int Row, const BitVector & ) >
  void applyAddMasked ( const Vector<DataType> &Arg, Vector<DataType> &Dest,
                        const BitVector & Mask ) const {
    BitMaskFunctorType imFunctor;
    // now multiply only on desired nodes and write the result only to desired nodes
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i <  this->getNumRows() ; ++i )
      if ( imFunctor ( Mask[i] ) ) {
        if ( rows[i] )
          Dest[i] += ( rows[i]->*multMasked ) ( Arg, i, Mask );
        else
          Dest [i] += _diagEntry * Arg[i];
      }
  }

  //! Optimized version of sparse matrix-vector-multiplication.
  /*! warning default behaviour: if there is no row at given index, apply assumes _diagEntry on the diagonal. */
  void apply ( const Vector<DataType> &Arg, Vector<DataType> &Dest ) const {
    if ( this->getNumRows() != Dest.size() || this->getNumCols() != Arg.size() ) {
      string msg = strprintf ( "aol::GenSparseMatrix::apply: Cannot apply %d by %d matrix from vector of size %d to vector of size %d.", this->getNumRows(), this->getNumCols(), Arg.size(), Dest.size() );
      throw ( Exception ( msg, __FILE__, __LINE__ ) );
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i <  this->getNumRows() ; ++i )
      if ( rows[i] )
        Dest[i] = rows[i]->mult ( Arg, i );
      else
        Dest [i] = _diagEntry * Arg[i];
  }


  //! Matrix-vector multiplication with masking functionality.
  //! Differently from apply and applyAdd, this function is not re-implemented
  //! in subclasses.
  //!
  //! Using a INCLUDE_*_WRITE_INT mode, actually a multiple of the identity is
  //! applied for bdry nodes, i. e. values are just copies from arg
  //! into "untouched" nodes.
  void applyMasked ( const Vector<DataType> &Arg, Vector<DataType> &Dest,
                     const BitVector & Mask, IncludeWriteMode applyMode ) const {

    // check dimensions
    if ( this->getNumRows() != Dest.size() || this->getNumCols() != Arg.size() ) {
      char errmsg [ 1024 ];
      sprintf ( errmsg, "aol::GenSparseMatrix::apply: Cannot apply %d by %d matrix from vector of size %d to vector of size %d.", this->getNumRows(), this->getNumCols(), Arg.size(), Dest.size() );
      throw ( Exception ( errmsg, __FILE__, __LINE__ ) );
    }

    switch ( applyMode ) {
    case INCLUDE_ALL_WRITE_ALL:
      applyMasked<BitMaskFunctorTrue, &Row<DataType>::multMaskedFunctorTrue > ( Arg, Dest, Mask );
      break;
    case INCLUDE_BD_WRITE_INT:
      applyMasked<BitMaskFunctorIdentity, &Row<DataType>::multMaskedFunctorNegate> ( Arg, Dest, Mask );
      break;
    case INCLUDE_ALL_WRITE_INT:
      applyMasked<BitMaskFunctorIdentity, &Row<DataType>::multMaskedFunctorTrue > ( Arg, Dest, Mask );
      break;
    case INCLUDE_INT_WRITE_ALL:
      applyMasked<BitMaskFunctorTrue, &Row<DataType>::multMaskedFunctorIdentity > ( Arg, Dest, Mask );
      break;
    case INCLUDE_INT_WRITE_INT:
      applyMasked<BitMaskFunctorIdentity, &Row<DataType>::multMaskedFunctorIdentity > ( Arg, Dest, Mask );
      break;
    default:
      throw aol::Exception ( "aol::GenSparseMatrix::applyMasked: unknown IncludeWriteMode", __FILE__, __LINE__ );
    }
  }

  template<typename BitMaskFunctorType, DataType ( Row<DataType>::* multMasked ) ( const Vector<DataType> &, int Row, const BitVector & ) >
  void applyMasked ( const Vector<DataType> &Arg, Vector<DataType> &Dest,
                     const BitVector & Mask ) const {
    BitMaskFunctorType imFunctor;
    // now multiply only on desired nodes and write the result only to desired nodes
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i <  this->getNumRows() ; ++i )
      if ( imFunctor ( Mask[i] ) ) {
        if ( rows[i] )
          Dest[i] = ( rows[i]->*multMasked ) ( Arg, i, Mask );
        else
          Dest [i] = _diagEntry * Arg[i];
      }

  }


  DataType multRow ( const Vector<DataType> &Arg, int RowNum ) const {
    if ( rows[ RowNum ] ) {
      return rows[ RowNum ]->mult ( Arg, RowNum );
    } else {
      return 0.;
    }
  }


  DataType rowSum ( const int I ) const {
    if ( rows[I] ) {
      return ( rows[I]->sum ( I ) );
#ifdef VERBOSE
    } else {
      cerr << "aol::GenSparseMatrix::rowSum: implicite identity row " << I << ", returning rowSum = 0 anyway." << endl;
#endif
    }
    return 0;
  }

  /** matrix is mapped to factor * matrix, this does not affect implicite diagonal (formerly identity) rows (identified by NULL pointer)
   */
  GenSparseMatrix<DataType>& operator*= ( const DataType factor ) {
    for ( int i = 0; i < this->_numRows; ++i ) {
      if ( rows[i] ) {
        rows[i]->scale ( i, factor );
#ifdef VERBOSE
      } else {
        cerr << "aol::GenSparseMatrix::scale: not scaling implicite identity row " << i << endl;
#endif
      }
    }

    return *this;
  }

  void newRow ( int I, Row<DataType> *NewRow ) {
    if ( rows[I] ) delete rows[I];
    rows[I] = NewRow;
  }

  //! clears all rows but keeps instances
  void setZero() {
    for ( int i = 0; i < this->_numRows; ++i ) {
      if ( rows[i] ) {
        rows[i]->setZero();
      }
    }
  }

  virtual void setRowToZero ( const int I ) {
    if ( rows[I] ) {
      rows[I]->setZero();
    }
  }


  //! deleteRow clears row and deletes instance, makes this row an implicit identity row
  void deleteRow ( int I ) {
    delete ( rows[I] );
    rows[I] = NULL;
  }

  //! scaleRow scales Ith row with factor unless it is implicit identity row
  void scaleRow ( const int RowNum, const DataType Factor ) {
    if ( rows[RowNum] ) {
      rows[RowNum]->scale ( RowNum, Factor );
#ifdef VERBOSE
    } else {
      cerr << "aol::GenSparseMatrix::scale: not scaling implicite identity row " << RowNum << endl;
#endif
    }
  }

  GenSparseMatrix<DataType>& addMultiple ( const GenSparseMatrix<DataType> &Matrix, const DataType factor ) {
    for ( int i = 0; i < this->_numRows; ++i ) {
      if ( !this->rows[ i ] && Matrix.rows[ i ] ) {
        this->rows[ i ] = newDefaultRow();
      }
      if ( this->rows[ i ] && Matrix.rows[ i] ) {
        this->rows[ i ]->addMultiple ( i, *Matrix.rows[ i ], factor );
      }
    }
    return ( *this );
  }

  using aol::Matrix<DataType>::operator+=;
  GenSparseMatrix<DataType>& operator+= ( const GenSparseMatrix<DataType> &Matrix ) {
    addMultiple ( Matrix, aol::NumberTrait<DataType>::one );
    return ( *this );
  }

  using aol::Matrix<DataType>::operator-=;
  GenSparseMatrix<DataType>& operator-= ( const GenSparseMatrix<DataType> &Matrix ) {
    addMultiple ( Matrix, -aol::NumberTrait<DataType>::one );
    return ( *this );
  }

  //! Return vector of row entries. Entries need not be sorted with respect to column index and zeros may be contained.
  void makeRowEntries ( vector<typename Row<DataType>::RowEntry > &vec, const int RowNum ) const {
    if ( rows[RowNum] ) {
      rows[RowNum]->makeRowEntries ( vec, RowNum );
    } else {
      vec.resize ( 1 );
      vec[0].col = RowNum;
      vec[0].value = _diagEntry;
    }
  }

  //! Same as makeRowEntries, only that the entries have to be sorted wrt column index.
  void makeSortedRowEntries ( vector<typename Row<DataType>::RowEntry > &vec, const int RowNum ) const {
    if ( rows[RowNum] ) {
      rows[RowNum]->makeSortedRowEntries ( vec, RowNum );
    } else {
      vec.resize ( 1 );
      vec[0].col = RowNum;
      vec[0].value = _diagEntry;
    }
  }

  void getRow ( int i, Vector<DataType> &v ) const {
    aol::Matrix<DataType>::getRow ( i, v );
  }

  const Row<DataType>& getRow ( int I ) const {
    return *rows[I];
  }


  bool checkForNANsAndINFs() const {
    for ( int i = 0; i < this->_numRows; ++i ) {
      if ( rows[i] && rows[i]->checkForNANsAndINFs() ) {
        return true;
      }
    }
    return false;
  }

  bool isSymmetric ( typename RealTrait<DataType>::RealType tol = 0 ) const {
    for ( int i = 0; i < this->_numRows; ++i )
      if ( rows[i] ) {
        vector< typename Row< DataType >::RowEntry > entries;
        makeRowEntries ( entries, i );
        for ( int j = 0; j < static_cast<int>(entries.size()); ++j )
          if ( Abs (entries[j].value - get ( entries[j].col, i ) ) > tol )
            return false;
      }
    return true;
  }

  /** Approximate comparison
   */
  bool isApproxEqual ( const GenSparseMatrix<DataType> &other, DataType epsilon ) {
    for ( int i = 0; i < this->_numRows; ++i ) {
      // check first if row only exists in one of them
      // (use XOR). Testing " != NULL " is
      // just to obtain a bool. One could equally cast
      // the pointers into booleans.
      if ( ( rows[i] != NULL ) ^ ( other.rows[i] != NULL ) )
        return false;

      // only test equality on rows if they exists
      if ( rows[ i ] && rows[ i ]->isApproxEqual ( i, *other.rows[ i ], epsilon ) == false )
        return false;
      else
        if ( _diagEntry != other._diagEntry )
          return false;
    }
    return true;
  }

  int numNonZeroes() const {
    int num = 0;
    for ( int i = 0; i < this->_numRows; ++i ) {
      if ( rows[i] )
        num += rows[i]->numNonZeroes();
      else
        if ( _diagEntry != aol::ZOTrait<DataType>::zero )
          num++;
    }
    return num;
  }

  int numStoredEntries() const {
    int num = 0;
    for ( int i = 0; i < this->_numRows; ++i ) {
      if ( rows[i] )
        num += rows[i]->numStoredEntries();
      else
        num++;
    }
    return num;
  }

  int numNonZeroRows() const {
    int num = 0;
    for ( int i = 0; i < this->_numRows; ++i )
      if ( rows[i] )
        num++;
      else
        if ( _diagEntry != aol::ZOTrait<DataType>::zero )
          num++;
    return num;
  }

  virtual int numNonZeroes ( int I ) const {
    if ( rows[I] )
      return rows[I]->numNonZeroes();
    else
      if ( _diagEntry != aol::ZOTrait<DataType>::zero )
        return 1;
    return 0;
  }

  virtual int numStoredEntries ( const int I ) const {
    if ( rows[I] )
      return rows[I]->numStoredEntries();
    else
      return 1;
  }

  //! reimplementation of transposition based on makeSortedRowEntries
  virtual void transposeTo ( aol::Matrix<DataType> &other_mat ) const {
    other_mat.setZero();
    for ( int i = 0; i < this->_numRows; ++i ) {
      vector<typename Row<DataType>::RowEntry > vec;
      makeSortedRowEntries ( vec, i );
      for ( typename vector<typename Row<DataType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
        other_mat.set ( it->col, i, it->value );
      }
    }
  }
  
  //
  void transpose () {
    throw UnimplementedCodeException ( "This function is not implemented yet.", __FILE__, __LINE__ );
  }

  //! get diagonal entry
  DataType getDiag ( int i ) const {
    return ( get ( i, i ) );
  }

  DataType getUnsetRowsDiagEntry() const {
    return _diagEntry;
  }
  void setUnsetRowsDiagEntry ( DataType diagEntry ) {
    _diagEntry = diagEntry;
  }

protected:
  vector< Row<DataType>* > rows;

private:
  void init() {
    for ( int i = 0; i < this->getNumRows(); ++i )
      rows[i] = NULL;
  }
  DataType _diagEntry;
};


template < typename DataType > class SparseMatrixRowIterator;

/** \brief a general, unstructured sparse matrix
 *  \ingroup Matrix
 */
template <typename DataType>
class SparseMatrix : public GenSparseMatrix<DataType> {
  friend class SparseMatrixRowIterator < DataType >;

  bool _deleteRows;

public:
  SparseMatrix ( ) : GenSparseMatrix<DataType> ( 0, 0 ), _deleteRows ( true ) {}

  SparseMatrix ( int Rows, int Columns )
      : GenSparseMatrix<DataType> ( Rows, Columns ), _deleteRows ( true ) {
    for ( int i = 0; i < Rows; ++i )
      this->rows[i] = new SparseRow<DataType>;
  }

  template <typename GridType>
  explicit SparseMatrix ( const GridType &Grid )
      : GenSparseMatrix<DataType> ( Grid ), _deleteRows ( true ) {
    init();
  }

  template <qc::Dimension Dim>
  explicit SparseMatrix ( const qc::GridSize<Dim> &GridSize )
      : GenSparseMatrix<DataType> ( GridSize ), _deleteRows ( true ) {
    init();
  }

  //! Copying
  explicit SparseMatrix ( const SparseMatrix& mat, CopyFlag copyFlag = DEEP_COPY )
      : GenSparseMatrix<DataType> ( mat.getNumRows (), mat.getNumCols () ) {
    switch ( copyFlag ) {
    case STRUCT_COPY: {
      init();
      _deleteRows = true;
    }
    break;
    case DEEP_COPY: {
      this->rows.resize ( this->_numRows );
      for ( int i = 0; i < this->_numRows; ++i ) {
        const SparseRow<DataType>* oldrow = static_cast<const SparseRow<DataType>* > ( mat.rows [i] );
        this->rows [i] = new SparseRow<DataType> ( *oldrow );
      }
      _deleteRows = true;
    }
    break;
    case FLAT_COPY: {
      this->rows.resize ( this->_numRows );
      for ( int i = 0; i < this->_numRows; ++i ) {
        this->rows[i] = static_cast<SparseRow<DataType>* > ( mat.rows [i] );
      }
      _deleteRows = false;
    }
    break;
    default: {
      throw UnimplementedCodeException ( "This CopyFlag is not implemented yet.", __FILE__, __LINE__ );
    }
    break;
    }
  }

  Matrix<DataType>* clone ( CopyFlag copyFlag = DEEP_COPY ) const {
    SparseMatrix *clone = new SparseMatrix ( *this, copyFlag );
    return clone;
  }

  //! Assignment of another SparseMatrix of the same size.
  //! If necessary, this matrix has to be resized before assignment.
  SparseMatrix& operator= ( const SparseMatrix& mat ) {
    // Beware of self-assignment
    if ( this == &mat ) return *this;

    if ( mat.getNumRows() != this->getNumRows() || mat.getNumCols() != this->getNumCols() )
      throw Exception ( "SparseMatrix::operator= : dimensions don't match.", __FILE__, __LINE__ ); // and we do not want to change the size

    for ( int i = 0; i < this->_numRows; ++i ) {
      if ( !this->rows[i] && mat.rows[i] ) {
        const SparseRow<DataType>* oldrow = static_cast<const SparseRow<DataType>* > ( mat.rows [i] );
        this->rows [i] = new SparseRow<DataType> ( *oldrow );
      } else {
        * ( dynamic_cast< SparseRow<DataType>* > ( ( this->rows [i] ) ) ) = * ( dynamic_cast< SparseRow<DataType>* > ( ( mat.rows [i] ) ) );
      }
    }

    return *this;
  }

  // change size of matrix, destroying old contents.
  void reallocate ( const int new_rows, const int new_cols ) {
    GenSparseMatrix<DataType>::reallocate ( new_rows, new_cols );
    init();
  }

  //! resize matrix, keeping old contents (as far as possible)
  void resize ( const int new_rows, const int new_cols ) {
    if ( new_rows > this->_numRows ) {
      // increase number of rows

      this->rows.resize ( new_rows );
      for ( int i = this->_numRows; i < new_rows; ++i )
        this->rows[i] = new SparseRow<DataType>;

    } else if ( new_rows < this->_numRows ) {
      // decrease number of rows

      for ( int i = ( this->_numRows ) - 1 ; i >= new_rows; --i )
        delete ( this->rows[i] );
      this->rows.resize ( new_rows );

    } // else number of rows stays the same

    if ( new_cols > this->_numCols ) {

      // rows do not know how large they are so we needn't tell them

    } else if ( new_cols < this->_numCols ) {

      throw aol::UnimplementedCodeException ( "aol::SparseMatrix<DataType>::resize: Decreasing the number of columns not implemented yet.", __FILE__, __LINE__ );
      // rows do not know about their size, but we must not leave invalid entries in them.

    } // else number of columns stays the same.

    // if we got this far, we can set the new size.
    this->_numRows = new_rows;
    this->_numCols = new_cols;

  }

  //! \brief Resizes matrix without checking if old content can be kept.
  void destructiveResize ( const int new_rows, const int new_cols ) {
    if ( new_rows > this->_numRows ) {
      // increase number of rows
      this->rows.resize ( new_rows );
      for ( int i = this->_numRows; i < new_rows; ++i )
        this->rows[i] = new SparseRow<DataType>;
    }
    else if ( new_rows < this->_numRows ) {
      // decrease number of rows
      for ( int i = ( this->_numRows ) - 1 ; i >= new_rows; --i )
        delete ( this->rows[i] );
      this->rows.resize ( new_rows );

    } // else number of rows stays the same

    this->_numRows = new_rows;
    this->_numCols = new_cols;
  }

  //! \brief Destroy (remove) row i.
  void destroyRow ( const int i ) {
    delete this->rows[i];
    this->rows.erase ( this->rows.begin () + i );
    --( this->_numRows );
  }

  //! \brief Insert a new row at position i.
  void insertRow ( const int i ) {
    this->rows.insert ( this->rows.begin () + i, new SparseRow<DataType>() );
    ++( this->_numRows );
  }

  virtual SparseRow<DataType>* newDefaultRow() const {
    return new SparseRow<DataType>();
  }

  virtual ~SparseMatrix() {
    destroy();
  }

  //! Delete entries that are exactly zero if these are stored for some reason
  void eraseZeroEntries () {
    for ( int i = 0; i < this->_numRows; ++i ) {
      dynamic_cast< aol::SparseRow<DataType>* >( this->rows[i] )->eraseZeroEntries();
    }
  }

  //! Adds multiple of one row to another
  void addMultipleRowToRow ( const int from, const int to, const DataType multiple ) {
    QUOC_ASSERT( from != to );

    // iterate over from-row and to-row simultaneously
    typename std::vector < typename aol::SparseRow<DataType>::qcCurMatrixEntry >::iterator it_from;
    typename std::vector < typename aol::SparseRow<DataType>::qcCurMatrixEntry >::iterator it_to;
    for ( it_from =  dynamic_cast< aol::SparseRow<DataType> * >(this->rows[from])->row.begin(),
          it_to =    dynamic_cast< aol::SparseRow<DataType> * >(this->rows[to])->row.begin();
          it_from != dynamic_cast< aol::SparseRow<DataType> * >(this->rows[from])->row.end(); ++it_from )  {
      // advance until column number in to-row is not smaller anymore or no entries are left
      while ( ( it_to != dynamic_cast< aol::SparseRow<DataType> * >(this->rows[to])->row.end() ) && ( it_to->col < it_from->col ) ) ++it_to;
      // no more entries in to-row or no match, insert before
      if ( ( it_to == dynamic_cast< aol::SparseRow<DataType> * >(this->rows[to])->row.end() ) || ( it_to->col != it_from-> col ) )
        // reassign iterator after insertion
        it_to = dynamic_cast< aol::SparseRow<DataType> * >(this->rows[to])->row.insert ( it_to, typename aol::SparseRow<DataType>::qcCurMatrixEntry ( it_from->col, multiple * it_from->value ) );
      // match
      else it_to->value += multiple * it_from->value;
    }
  }

  //! Adds multiple of one column to another
  void addMultipleColToCol ( const int from, const int to, const DataType multiple ) {
    QUOC_ASSERT( from != to );

    const int first  = aol::Min( from, to );
    const int second = aol::Max( from, to );

    // have to do it for all rows
    for ( int i = 0; i < this->getNumRows(); ++i )  {
      typename std::vector < typename aol::SparseRow<DataType>::qcCurMatrixEntry >::iterator it_first  = dynamic_cast< aol::SparseRow<DataType> * >(this->rows[i])->row.begin();
      typename std::vector < typename aol::SparseRow<DataType>::qcCurMatrixEntry >::iterator it_second = dynamic_cast< aol::SparseRow<DataType> * >(this->rows[i])->row.begin();
      // iterate until first entry is passed
      while ( ( it_first != dynamic_cast< aol::SparseRow<DataType> * >(this->rows[i])->row.end() ) && ( it_first->col < first ) ) it_first++;
      // iterate further until second entry is passed
      it_second = it_first;
      while ( ( it_second != dynamic_cast< aol::SparseRow<DataType> * >(this->rows[i])->row.end() ) && ( it_second->col < second ) ) it_second++;

      // it_first shall now correspond to from and it_second to to
      if ( to < from ) std::swap( it_first, it_second );

      // if from-entry was not found there is nothing to do
      if ( ( it_first != dynamic_cast< aol::SparseRow<DataType> * >(this->rows[i])->row.end() ) && ( it_first->col == from ) )  {
        if ( ( it_second == dynamic_cast< aol::SparseRow<DataType> * >(this->rows[i])->row.end() ) || ( it_second->col != to ) )
        // column number of it_second must be larger or end of row vector, insert before
          dynamic_cast< aol::SparseRow<DataType> * >(this->rows[i])->row.insert ( it_second, typename aol::SparseRow<DataType>::qcCurMatrixEntry ( to, multiple * it_first->value ) );
        // match
        else it_second->value += multiple * it_first->value;
      }
    }
  }

  //! Set the whole row and col to zero except for the diagonal entry
  void setRowColToDiagonal ( const int index, const DataType diagEntry = 1. )  {
    // set row to zero
    this->rows[index]->setZero();
    // set col to zero
    for ( int i = 0; i < this->getNumRows(); ++i )  {
      typename std::vector < typename aol::SparseRow<DataType>::qcCurMatrixEntry >::iterator it;
      for ( it = dynamic_cast< aol::SparseRow<DataType> * >(this->rows[i])->row.begin(); it != dynamic_cast< aol::SparseRow<DataType> * >(this->rows[i])->row.end(); ++it )  {
        if ( it->col == index ) dynamic_cast< aol::SparseRow<DataType> * >(this->rows[i])->row.erase( it );
        if ( it->col >= index ) break;
      }
    }
    // set diagonal entry
    dynamic_cast< aol::SparseRow<DataType> * >(this->rows[index])->set( index,diagEntry );
  }

  //! Collapse a pair of rows and cols by adding multiples of one to the other and then setting the former row / col to zero except for the diagonal entry
  void collapseRowCol ( const int from, const int to, const DataType multiple, const DataType diagEntry = 1. )  {
    addMultipleRowToRow( from, to, multiple );
    addMultipleColToCol( from, to, multiple );
    setRowColToDiagonal( from, diagEntry );
  }

  //! Loads a matrix in the Harwell-Boeing format. Requires the external gmm.
  //! \author Berkels
  void loadHarwellBoeing ( const char *FileName ) {
#ifdef USE_EXTERNAL_GMM
    gmm::csc_matrix<DataType> gmmMat;
    gmm::Harwell_Boeing_load ( FileName, gmmMat );

    this->reallocate ( gmmMat.nrows(), gmmMat.ncols() );
    // These loops are not very efficient, so don't load matrices in performance critical code or
    // make this more efficient by only iterating of the non-zero entries of loaded matrix.
    for ( int i = 0; i < this->getNumRows(); ++i )
      for ( int j = 0; j < this->getNumCols(); ++j )
        if ( gmmMat ( i, j ) != 0 )
          this->set( i, j, gmmMat ( i, j ) );
#else
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( FileName );
    throw aol::Exception ( "Reading matrices in the Harwell-Boeing format requires the external gmm.", __FILE__, __LINE__ );
#endif // USE_EXTERNAL_GMM
  }

  // applies transposed matrix to argument
  void applyAddTransposed( const aol::Vector<DataType>& Arg, aol::Vector<DataType>& Dest ) const {
#ifdef BOUNDS_CHECK
    if ( Arg.size() != this->getNumRows() )
      throw aol::Exception ( "applyAddTransposed():: Arg has wrong size", __FILE__, __LINE__, __FUNCTION__ );
    if ( Dest.size() !=  this->getNumCols () )
      throw aol::Exception ( "applyAddTransposed():: Dest has wrong size", __FILE__, __LINE__, __FUNCTION__ );
#endif
    
    for ( int i = 0; i < this->getNumRows(); ++i ){
      vector<typename aol::Row<DataType>::RowEntry > rowEntries;
      this->makeRowEntries ( rowEntries, i );
      for ( typename vector<typename aol::Row<DataType>::RowEntry >::iterator it = rowEntries.begin(); it != rowEntries.end(); ++it )
        Dest[ it->col ] += it->value * Arg[i];
    }
    
  }
  
  
private:
  void init( ) {
    for ( int i = 0; i < this->getNumRows(); ++i )
      this->rows[i] = new SparseRow<DataType>;
  }

  void destroy() {
    if ( this->_deleteRows )
      for ( unsigned int i = 0; i < this->rows.size(); ++i ) {
        delete this->rows[ i ];
      }
  }
};


template <class RealType, class SparseOpType>
class RowEntryOp : public  aol::Op<aol::Vector<RealType> > {
  const SparseOpType &_sparseOp;
  const int _numRows;
public:
  RowEntryOp ( const SparseOpType &SparseOp, int NumRows ) : _sparseOp ( SparseOp ), _numRows ( NumRows ) {}

  int getNumRows() const {
    return _numRows;
  }

  void applyAdd ( const aol::Vector<RealType> &ArgVec, aol::Vector<RealType> &DestVec ) const {
    const RealType * Arg  = ArgVec.getData();
    RealType * Dest = DestVec.getData();
    vector<typename aol::Row<RealType>::RowEntry > rowEntries;
    for ( int i = 0; i < _numRows; ++i ) {
      _sparseOp.makeRowEntries ( rowEntries, i );
      for ( typename vector<typename aol::Row<RealType>::RowEntry >::iterator it = rowEntries.begin(); it != rowEntries.end(); ++it ) {
        Dest[ i ] += Arg[ it->col ] * it->value;
      }
    }
  }
};

/**
 * \brief A sparse matrix in triplet format.
 * \author Toelkes
 * \ingroup Matrix
 *
 * A sparse matrix class which stores its entries in a triplet format. If an entry appears more than once, the values are summed.
 * Useful for assembling, slow for arithmetic operations.
 */
template <typename DataType, template<typename> class VectorType = aol::Vector>
class TripletMatrix : public Matrix<DataType> {
public:
  typedef VectorType<int> IndexVectorType;
  typedef VectorType<DataType> ValueVectorType;
protected:
  VectorType<int> _rowIndex;
  VectorType<int> _colIndex;
  VectorType<DataType> _value;

  template<typename, template <typename>class> friend class TripletMatrixOffset;
  template<typename, template <typename>class> friend class TripletMatrixOffsetUpperTriangle;

  //! \brief Remove entries that have value zero.
  void removeZeroEntries () {
    for ( int i = 0; i < _value.size (); ++i ) {
      if ( _value[i] == 0.0 ) {
        _rowIndex.erase ( i );
        _colIndex.erase ( i );
        _value.erase ( i );
      }
    }
  }

  //! \brief Set a matrix entry to zero by setting _value for all its appearances to 0.0.
  inline void setEntryToZero ( int row, int col ) {
    for ( int k = 0; k < _rowIndex.size (); ++k ) {
      if ( _rowIndex[k] == row && _colIndex[k] == col ) {
        _value[k] = 0.0;
      }
    }
  }

  //! \brief Helper class for sorting indices by matrix position, first by row then by column.
  //! \author Toelkes
  class RowColIndexSorter {
  protected:
    const VectorType<int> &_rowIndex;
    const VectorType<int> &_colIndex;

  public:
    explicit RowColIndexSorter ( const VectorType<int> &rowIndex, const VectorType<int> &colIndex )
    : _rowIndex ( rowIndex ), _colIndex ( colIndex ) {}

    bool operator() ( int i, int j ) {
      return ( ( _rowIndex[i] < _rowIndex[j] ) || ( _rowIndex[i] == _rowIndex[j] && _colIndex[i] < _colIndex[j] ) );
    }
  };

  //! \brief Helper class for sorting indices by matrix position, first by column then by row.
  //! \author Toelkes
  class ColRowIndexSorter {
  protected:
    const VectorType<int> &_rowIndex;
    const VectorType<int> &_colIndex;

  public:
    explicit ColRowIndexSorter ( const VectorType<int> &rowIndex, const VectorType<int> &colIndex )
    : _rowIndex ( rowIndex ), _colIndex ( colIndex ) {}

    bool operator() ( int i, int j ) {
      return ( ( _colIndex[i] < _colIndex[j] ) || ( _colIndex[i] == _colIndex[j] && _rowIndex[i] < _rowIndex[j] ) );
    }
  };
public:
  //! \brief Standard constructor.
  explicit TripletMatrix ( const unsigned int numRows, const unsigned int numCols )
  : Matrix<DataType> ( numRows, numCols ) {}

  //! \brief Copy constructor.
  explicit TripletMatrix ( const TripletMatrix &mat, CopyFlag copyFlag = DEEP_COPY )
  : Matrix<DataType> ( mat.getNumRows (), mat.getNumCols () ) {
    switch ( copyFlag ) {
    case DEEP_COPY:
      _rowIndex = mat._rowIndex;
      _colIndex = mat._colIndex;
      _value = mat._value;
      break;
    case STRUCT_COPY:
      break;
    default:
      throw aol::UnimplementedCodeException ( "Copy flag not implemented", __FILE__, __LINE__ );
      break;
    }
  }

  //! \brief Standard destructor.
  virtual ~TripletMatrix () {}

  //! \brief Add value to entry at position (i, j).
  virtual void add ( int i, int j, DataType value ) {
#ifdef BOUNDS_CHECK
    if ( i > this->getNumRows () )
      throw aol::Exception ( "Row index is out of bounds", __FILE__, __LINE__, __FUNCTION__ );
    if ( j > this->getNumCols () )
      throw aol::Exception ( "Col index is out of bounds", __FILE__, __LINE__, __FUNCTION__ );
#endif

    _rowIndex.pushBack ( i );
    _colIndex.pushBack ( j );
    _value.pushBack ( value );
  }

  //! \brief Remove the i-th row and column.
  void removeRowCol ( unsigned int i ) {
    for ( int k = 0; k < _rowIndex.size (); ++k ) {
      if ( _rowIndex[k] == i || _colIndex[k] == i ) {
        _rowIndex.erase ( k );
        _colIndex.erase ( k );
        _value.erase ( k );
      }
      else {
        if ( _rowIndex[k] > i )
          --( _rowIndex[k] );
        if ( _colIndex[k] > i )
          --( _colIndex[k] );
      }
    }
  }

  //! \brief Set the i-th row to zero
  void setRowToZero ( const int i ) {
    for ( int k = 0; k < _rowIndex.size (); ++k ) {
      if ( _rowIndex[k] == i ) {
        _value[k] = 0.0;
      }
    }
  }

  //! \brief Set the i-th column to zero
  void setColToZero ( const int j ) {
    for ( int k = 0; k < _rowIndex.size (); ++k ) {
      if ( _colIndex[k] == j ) {
        _value[k] = 0.0;
      }
    }
  }

  //! \brief Apply method for compatibility with other matrices.
  virtual void apply ( const aol::Vector<DataType> &Arg, aol::Vector<DataType> &Dest ) const {
    Dest.setZero();
    for ( int k = 0; k < _rowIndex.size(); ++k )
      Dest[_rowIndex[k]] += Arg[_colIndex[k]] * _value[k];
  }

  //! \brief Apply method to be used when only the upper triangular part was stored
  virtual void applyUpperTriangle ( const aol::Vector<DataType> &Arg, aol::Vector<DataType> &Dest ) const {
    Dest.setZero();
    for ( int k = 0; k < _rowIndex.size(); ++k )  {
      Dest[_rowIndex[k]] += Arg[_colIndex[k]] * _value[k];
      if ( _rowIndex[k] != _colIndex[k] )
        Dest[_colIndex[k]] += Arg[_rowIndex[k]] * _value[k];
    }
  }

  //! \brief applyAdd method for compatibility with other matrices.
  //! \warning Not implemented.
  virtual void applyAdd ( const aol::Vector<DataType> &, aol::Vector<DataType> & ) const {
    throw aol::UnimplementedCodeException ( "ApplyAdd not implemented...", __FILE__, __LINE__ );
  }

  //! \brief Set the matrix to zero.
  virtual void setZero () {
    _rowIndex.reallocateClear ( 0 );
    _colIndex.reallocateClear ( 0 );
    _value.reallocateClear ( 0 );
  }

  void reallocate ( const int Rows, const int Columns ) {
    setZero();
    this->_numRows = Rows;
    this->_numCols = Columns;
  }

  //! \brief Get matrix entry (row, col).
  virtual DataType get ( int row, int col ) const {
    DataType ret = static_cast<DataType> ( 0 );

    for ( int i = 0; i < _rowIndex.size (); ++i ) {
      if ( _rowIndex[i] == row && _colIndex[i] == col )
        ret += _value[i];
    }

    return ret;
  }

  //! \brief Set matrix entry (row, col) to value.
  virtual void set ( int row, int col, DataType value ) {
    setEntryToZero ( row, col );
    add ( row, col, value );
  }

  //! \brief Erase value at (i, j).
  void eraseValue ( int row, int col ) {
    int k = 0;
    while ( k < _rowIndex.size () ) {
      if ( _rowIndex[k] == row && _colIndex[k] == col ) {
        _rowIndex.erase ( k );
        _colIndex.erase ( k );
        _value.erase ( k );
      }
      // Increment only, if the kth element has not been removed.
      // Otherwise, the former k+1st element is now the kth element.
      else {
        ++k;
      }
    }
  }

  //! \brief Find duplicate entries and sum these entries into one entry.
  void sumDuplicates () {
    std::vector<unsigned int> index ( _rowIndex.size () );
    for ( unsigned int i = 0; i < index.size (); ++i )
      index[i] = i;

    RowColIndexSorter rcSorter ( _rowIndex, _colIndex );
    std::sort ( index.begin (), index.end (), rcSorter );

    unsigned int i = 0;
    unsigned int j;
    while ( i < index.size () - 1 ) {
      j = i + 1;
      while ( _rowIndex[index[j]] == _rowIndex[index[i]] && _colIndex[index[j]] == _colIndex[index[i]] ) {
        _value[index[i]] += _value[index[j]];
        _value[index[j]] = static_cast<DataType> ( 0 );

        // If the last index has been processed, break.
        if ( ++j >= index.size () )
          break;
      }
      i += ( j - i );
    }

    removeZeroEntries ();
  }

  void getRowColIndexSorting ( std::vector<unsigned int> &index ) const {
    getIndexSorting<RowColIndexSorter> ( index );
  }

  void getColRowIndexSorting ( std::vector<unsigned int> &index ) const {
    getIndexSorting<ColRowIndexSorter> ( index );
  }

  template<typename SorterType>
  void getIndexSorting ( std::vector<unsigned int> &index ) const {
    index.resize ( _rowIndex.size () );
    for ( unsigned int i = 0; i < index.size (); ++i )
      index[i] = i;

    SorterType sorter ( _rowIndex, _colIndex );
    std::sort ( index.begin (), index.end (), sorter );
  }

  //! \brief Get const reference to row index vector.
  const VectorType<int>& getRowIndexReference () const {
    return _rowIndex;
  }

  //! \brief Get reference to row index vector.
  VectorType<int>& getRowIndexReference () {
    return _rowIndex;
  }

  //! \brief Get const reference to column index vector.
  const VectorType<int>& getColIndexReference () const {
    return _colIndex;
  }

  //! \brief Get reference to column index vector.
  VectorType<int>& getColIndexReference () {
    return _colIndex;
  }

  //! \brief Get const reference to row value vector.
  const VectorType<DataType>& getValueReference () const {
    return _value;
  }

  //! \brief Get reference to row value vector.
  VectorType<DataType>& getValueReference () {
    return _value;
  }

  //! \brief Converts the triplet matrix into a SparseMatrix.
  //! \param[out] sparseMatrix Will have the same entries as this, after the conversion (existing entries are truncated).
  void toSparseMatrix ( aol::SparseMatrix<DataType> &sparseMatrix ) {
    sparseMatrix.resize ( this->getNumRows (), this->getNumCols () );
    sparseMatrix.setZero ();

    for ( int i = 0; i < _value.size (); ++i ) {
      sparseMatrix.add ( _rowIndex[i], _colIndex[i], _value[i] );
    }
  }

  //! \brief Add another TripletMatrix.
  virtual TripletMatrix<DataType>& operator+= ( const TripletMatrix<DataType> &Other ) {
    for ( int i = 0; i < Other._value.size(); ++i )
      this->add ( Other._rowIndex[i], Other._colIndex[i], Other._value[i] );
    return *this;
  }
  
  using Matrix<DataType>::operator+=;

  //! \brief Add matrix block at a particular position.
  virtual void addMultipleAtPosition ( int I, int J, const TripletMatrix<DataType> &Other, DataType Factor = 1 ) {
    for ( int i = 0; i < Other._value.size(); ++i )
      this->add ( Other._rowIndex[i]+I, Other._colIndex[i]+J, Factor * Other._value[i] );
  }
  
  using Matrix<DataType>::addMultipleAtPosition;

  //! Adds multiple of one row to another
  void addMultipleRowToRow ( const int from, const int to, const DataType multiple ) {
    QUOC_ASSERT( from != to );
    for ( int k = 0; k < _rowIndex.size (); ++k ) {
      if ( _rowIndex[k] == from ) {
        add( to, _colIndex[k], multiple * _value[k] );
      }
    }
  }

  //! Adds multiple of one column to another
  void addMultipleColToCol ( const int from, const int to, const DataType multiple ) {
    QUOC_ASSERT( from != to );
    for ( int k = 0; k < _rowIndex.size (); ++k ) {
      if ( _colIndex[k] == from ) {
        add( _rowIndex[k], to, multiple * _value[k] );
      }
    }
  }

  //! Set the whole row and col to zero except for the diagonal entry
  //! \warning naive and slow implementation
  void setRowColToDiagonal ( const int index, const DataType diagEntry = 1. )  {
    setRowToZero( index );
    setColToZero( index );
    add( index, index, diagEntry );
  }

  //! Set whole rows and cols given by an index set to zero except for the diagonal entry
  //! In principle it should be faster to traverse the values list only once and query the index set each time
  //! templatized to be usuable with std::unordered_set and previous versions thereof
  template <typename HashSetType>
  void setRowsColsToDiagonal ( const HashSetType & hashset, const DataType diagEntry = 1. )  {
    for ( int k = 0; k < _rowIndex.size (); ++k ) {
      if ( ( hashset.find ( _rowIndex[k] ) != hashset.end() ) || ( hashset.find ( _colIndex[k] ) != hashset.end() ) )
        _value[k] = 0.0;
    }
    for ( typename HashSetType::const_iterator it = hashset.begin(); it != hashset.end(); ++it )
      add ( *it, *it, diagEntry );
  }
  
  /**
   * \brief transpose a triplet matrix to other_mat
   * \author Tatano
   */
  void transposeTo ( aol::TripletMatrix<DataType> &other_mat ) const {
    other_mat._rowIndex.resize(this->_colIndex.size());
    other_mat._colIndex.resize(this->_rowIndex.size());
    other_mat._value.resize(this->_value.size());
    
    other_mat._rowIndex = this->_colIndex;
    other_mat._colIndex = this->_rowIndex;
    other_mat._value = this->_value;
  }

  bool checkForNANsAndINFs() const {
    return _value.checkForNANsAndINFs();
  }
};

/**
 * \brief Wrapper for using the aol::TripletMatrix in a block matrix.
 * \author Toelkes
 *
 * This class takes a TripletMatrix, sizes and offsets and maps matrix accesses
 * to the sub block given by these values.
 */
template <typename DataType, template<typename> class VectorType = aol::Vector>
class TripletMatrixOffset : public Matrix<DataType> {
private:
  //! \brief Standard constructor forbidden
  //! \warning Should never be called, might cause bad behavior!
  TripletMatrixOffset () : _mat ( TripletMatrix<DataType, VectorType>( ) ) {}
  //! \brief Constructor for compatibility with usual matrices
  //! \warning Should never be called, might cause bad behavior!
  TripletMatrixOffset ( const int numRows, const int numCols ) : _mat ( TripletMatrix<DataType, VectorType>( ) ) {}

protected:
  TripletMatrix<DataType, VectorType> &_mat;
  int _rowOffset;
  int _colOffset;
  bool _setWholeMatToZero;  //!< If true, setZero sets _mat to zero instead of only the block (much faster).

public:
  //! \brief Constructor.
  //! \param[in] mat The matrix.
  //! \param[in] numRows The number of rows of the block.
  //! \param[in] numCols The number of columns of the block.
  //! \param[in] rowOffset The row (in the full matrix mat) at which the block starts.
  //! \param[in] colOffset The column in the full matrix mat) at which the block starts.
  explicit TripletMatrixOffset ( TripletMatrix<DataType> &mat, int numRows, int numCols,
      const int rowOffset, const int colOffset )
  : Matrix<DataType> ( numRows, numCols ), _mat ( mat ), _rowOffset ( rowOffset ), _colOffset ( colOffset ), _setWholeMatToZero ( false ) {}

  //! \brief Copy constructor.
  explicit TripletMatrixOffset ( const TripletMatrixOffset &other, CopyFlag copyFlag = DEEP_COPY )
  : Matrix<DataType> ( other.getNumRows (), other.getNumCols () ), _mat ( other._mat ), _setWholeMatToZero ( other._setWholeMatToZero ) {
    switch ( copyFlag ) {
    case DEEP_COPY:
      _rowOffset = other._rowOffset;
      _colOffset = other._colOffset;
      break;
    default:
      throw aol::UnimplementedCodeException ( "Copy flag not implemented", __FILE__, __LINE__ );
      break;
    }
  }

  //! \brief Destructor.
  virtual ~TripletMatrixOffset () {}

  //! \brief Returns entry (row, col) of the block.
  virtual DataType get ( int row, int col ) const {
    return _mat.get ( row + _rowOffset, col + _colOffset );
  }

  //! \brief Get reference to underlying TripletMatrix
  TripletMatrix<DataType> &getTripletMatrixRef () const  {
    return _mat;
  }

  //! \brief Get row offset
  inline int getRowOffset( ) const  {
    return _rowOffset;
  }

  //! \brief Get col offset
  inline int getColOffset( ) const  {
    return _colOffset;
  }

  //! \brief Sets entry (row, col) of the block to value.
  virtual void set ( int row, int col, DataType value ) {
    _mat.set ( row + _rowOffset, col + _colOffset, value );
  }

  //! \brief Adds value to entry (row, col) of the block.
  virtual void add ( int row, int col, DataType value ) {
    _mat.add ( row + _rowOffset, col + _colOffset, value );
  }

  void add ( const TripletMatrixOffset<DataType> &other ) {
    for ( int i = 0; i < other.getTripletMatrixRef ()._rowIndex.size (); ++i ) {
      getTripletMatrixRef ()._rowIndex.pushBack ( other.getTripletMatrixRef ()._rowIndex[i] );
      getTripletMatrixRef ()._colIndex.pushBack ( other.getTripletMatrixRef ()._colIndex[i] );
      getTripletMatrixRef ()._value.pushBack ( other.getTripletMatrixRef ()._value[i] );
    }
  }

  //! \brief Sets the block to zero.
  //! \attention setZero on this class is very slow!
  virtual void setZero () {
    cerr << aol::color::red << "Slow!" << aol::color::reset << endl;
    // Fast and dangerous
    if ( _setWholeMatToZero ) {
      getTripletMatrixRef ().setZero ();
    }
    // slow and safe
    else {
      for ( int i = 0; i < this->getNumRows (); ++i ) {
        for ( int j = 0; j < this->getNumCols (); ++j ) {
          _mat.eraseValue ( i + _rowOffset, j + _colOffset );
        }
      }
    }
  }

  //! \brief Remove the i-th row and column of the block.
  void removeRowCol ( unsigned int i ) {
    if ( _rowOffset != _colOffset )
      throw aol::Exception ( "Only supported if rowOffset == colOffset", __FILE__, __LINE__, __FUNCTION__ );

    _mat.removeRowCol ( i + _rowOffset );
  }

  //! \brief Set the i-th row to zero
  virtual void setRowToZero ( const int i ) {
    for ( int k = 0; k < _mat._rowIndex.size (); ++k ) {
      if ( _mat._rowIndex[k] == i + _rowOffset &&
          _mat._colIndex[k] >= _colOffset && _mat._colIndex[k] < this->getNumCols () + _colOffset ) {
        _mat._value[k] = 0.0;
      }
    }
  }

  //! \brief Set the i-th column to zero
  virtual void setColToZero ( const int j ) {
    for ( int k = 0; k < _mat._colIndex.size (); ++k ) {
      if ( _mat._colIndex[k] == j + _colOffset &&
          _mat._rowIndex[k] >= _rowOffset && _mat._rowIndex[k] < this->getNumRows () + _rowOffset ) {
        _mat._value[k] = 0.0;
      }
    }
  }

  //! \brief Set the i-th row and column to zero
  virtual void setRowColToZero ( const int i ) {
    // Call two seperate methods (this is slow anyway)
    setRowToZero ( i );
    setColToZero ( i );
  }

  //! \brief Set the i-th row and column to zero and the diagonal entry to diag
  void setRowColToDiag ( int i, DataType diag ) {
    // Call two seperate methods (this is slow anyway)
    setRowColToZero ( i );
    set ( i, i, diag );
  }

  //! \brief Apply function for compatibility with other matrices.
  virtual void apply ( const aol::Vector<DataType> & Arg, aol::Vector<DataType> & Dest ) const {
    Dest.setZero();
    applyAdd( Arg, Dest );
  }

  //! \brief applyAdd method for compatibility with other matrices.
  virtual void applyAdd ( const aol::Vector<DataType> & Arg, aol::Vector<DataType> & Dest ) const {
    for ( int k = 0; k < _mat._rowIndex.size(); ++k )
      if (    ( _mat._rowIndex[k] >= _rowOffset ) && ( _mat._rowIndex[k] < _rowOffset + this->getNumRows() )
           && ( _mat._colIndex[k] >= _colOffset ) && ( _mat._colIndex[k] < _colOffset + this->getNumCols() ) )
        Dest[ _mat._rowIndex[k] - _rowOffset ] += Arg[ _mat._colIndex[k] - _colOffset ] * _mat._value[k];
  }

  //! \brief Determine if setZero sets the whole underlying matrix to zero.
  //! \note Setting this to true will make setZero much faster but also make it set not only the handled block but everything to zero.
  //! \warning Use with care!
  void setWholeMatToZero ( bool b ) {
    _setWholeMatToZero = b;
  }
};

/**
 * \brief  Special TripletMatrixOffset wrapper
 * \author geihe
 *
 * Fills only the upper triangular part of the specified block, all other write requests are discarded
 */
template <typename DataType, template<typename> class VectorType = aol::Vector>
class TripletMatrixOffsetUpperTriangle : public TripletMatrixOffset<DataType, VectorType> {
private:
  //! \brief Standard constructor forbidden
  TripletMatrixOffsetUpperTriangle () {}

public:
  //! \brief Constructor.
  //! \param[in] mat The matrix.
  //! \param[in] numRows The number of rows of the block.
  //! \param[in] numCols The number of columns of the block.
  //! \param[in] rowOffset The row (in the full matrix mat) at which the block starts.
  //! \param[in] colOffset The column in the full matrix mat) at which the block starts.
  explicit TripletMatrixOffsetUpperTriangle ( TripletMatrix<DataType> &mat, int numRows, int numCols,
                                              const int rowOffset, const int colOffset )
  : TripletMatrixOffset<DataType, VectorType> ( mat, numRows, numCols, rowOffset, colOffset )
  {}

  //! \brief Destructor.
  virtual ~TripletMatrixOffsetUpperTriangle () {}

  //! \brief Copy constructor.
  explicit TripletMatrixOffsetUpperTriangle ( const TripletMatrixOffsetUpperTriangle &other, CopyFlag copyFlag = DEEP_COPY )
  : TripletMatrixOffset<DataType, VectorType> ( other, copyFlag )  {}

  //! \brief Adds value to entry (row, col) of the block only if within upper triangular part of the block
  virtual void add ( int row, int col, DataType value ) {
    // only the upper triangular part should be filled, filter out the rest
    if ( row <= col ) this->_mat.add ( row + this->_rowOffset, col + this->_colOffset, value );
  }

  //! \brief Apply function for compatibility with other matrices.
  virtual void apply ( const aol::Vector<DataType> & Arg, aol::Vector<DataType> & Dest ) const {
    Dest.setZero();
    applyAdd( Arg, Dest );
  }

  //! \brief applyAdd method for compatibility with other matrices.
  virtual void applyAdd ( const aol::Vector<DataType> & Arg, aol::Vector<DataType> & Dest ) const {
    for ( int k = 0; k < this->_mat._rowIndex.size(); ++k )
      if (    ( this->_mat._rowIndex[k] >= this->_rowOffset ) && ( this->_mat._rowIndex[k] < this->_rowOffset + this->getNumRows() )
           && ( this->_mat._colIndex[k] >= this->_colOffset ) && ( this->_mat._colIndex[k] < this->_colOffset + this->getNumCols() ) )  {
        Dest[ this->_mat._rowIndex[k] - this->_rowOffset ] += Arg[ this->_mat._colIndex[k] - this->_colOffset ] * this->_mat._value[k];
        if ( this->_mat._rowIndex[k] - this->_rowOffset != this->_mat._colIndex[k] - this->_mat._colIndex[k] )
          // we are not on the diagonal of the current block
          Dest[ this->_mat._colIndex[k] - this->_colOffset ] += Arg[ this->_mat._rowIndex[k] - this->_rowOffset ] * this->_mat._value[k];
      }
  }
};

//! \brief Base class for compressed sparse (row, column) matrices
//! \author Toelkes
//! \ingroup Matrix
template <typename DataType, typename IndexType>
class CSMatrix : public Matrix<DataType> {
protected:
  aol::Vector<IndexType> _index;
  aol::Vector<IndexType> _indPointer;
  aol::Vector<DataType>  _value;

public:
  explicit CSMatrix ()
  : Matrix<DataType> ( 0, 0 ), _index ( 0 ), _indPointer ( 0 ), _value ( 0 ) {}

  explicit CSMatrix ( const CSMatrix<DataType, IndexType> &other, CopyFlag copyFlag = DEEP_COPY )
  : Matrix<DataType> ( other.getNumRows(), other.getNumCols() ), _index ( other._index, copyFlag ),
    _indPointer ( other._indPointer, copyFlag ), _value ( other._value, copyFlag ) {
    switch ( copyFlag ) {
    case DEEP_COPY:
    case FLAT_COPY:
      // Nothing more to do.
      break;
    case STRUCT_COPY:
      // Set size to zero, s.t. all methods work as expected.
      this->_index.resize ( 0 );
      this->_indPointer.resize ( 0 );
      this->_value.resize ( 0 );
      break;
    default:
      throw aol::UnimplementedCodeException ( "Copy flag not implemented", __FILE__, __LINE__ );
      break;
    }
  }

  void getSubMatrix ( IndexType minIndex, IndexType maxIndex,
                      IndexType minPointed, IndexType maxPointed,
                      bool pointedIsRow,
                      CSMatrix<DataType, IndexType> &subMatrix ) const {
    // Set sub matrix dimensions
    if ( pointedIsRow ) {
      subMatrix._numRows = maxPointed - minPointed + 1;
      subMatrix._numCols = maxIndex - minIndex + 1;
    }
    else {
      subMatrix._numRows = maxIndex - minIndex + 1;
      subMatrix._numCols = maxPointed - minPointed + 1;
    }

    // Initialize sub matrix data structures
    subMatrix._index.resize ( 0 );
    subMatrix._indPointer.resize ( maxPointed - minPointed + 2 );
    subMatrix._value.resize ( 0 );

    // Copy matrix entries
    for ( IndexType p = minPointed; p <= maxPointed; ++p ) {
      if ( p - minPointed > 0 )
        subMatrix._indPointer[p - minPointed + 1] = subMatrix._indPointer[p - minPointed];

      for ( IndexType k = this->_indPointer[p]; k < this->_indPointer[p + 1]; ++k ) {
        if ( this->_index[k] >= minIndex && this->_index[k] <= maxIndex ) {
          ++( subMatrix._indPointer[p - minPointed + 1] );
          subMatrix._index.pushBack ( this->_index[k] - minIndex );
          subMatrix._value.pushBack ( this->_value[k] );
        }
      }
    }
    // Should not be necessary
    subMatrix._indPointer[subMatrix._indPointer.size () - 1] = subMatrix._value.size ();
  }
};

/**
 * \brief A compressed sparse column matrix.
 * \attention Changes to the sparsity structure of the matrix (via add, set, ...) are slow!
 * \warning If IndexType is set to a type bigger than int (or unsigned int), triplet to csc conversion does no longer work (see CSCMatrix::setFromTriplet).
 * \author Toelkes
 * \ingroup Matrix
 *
 * A compressed sparse column matrix. Fast for matrix-vector multiplication and certain arithmetic operations.
 * Changes to the sparsity structure are slow. Use the TripletMatrix for assembling and then convert it into a CSCMatrix.
 */
template <typename DataType, typename IndexType = int>
class CSCMatrix : public CSMatrix<DataType, IndexType> {
  // In this class, _index lists row indices and _indPointer holds column pointers.
protected:
  void setFromTriplet ( const TripletMatrix<DataType> &tripletMatrix ) {
    const Vector<int> &tripletRow = tripletMatrix.getRowIndexReference ();
    const Vector<int> &tripletCol = tripletMatrix.getColIndexReference ();
    const Vector<DataType> &tripletVal = tripletMatrix.getValueReference ();

    // helper variables, make t long int, because it needs to be signed. If IndexType = unsigned long int, this does not work!
    aol::Vector<int64_t> t ( aol::Max ( this->getNumRows (), this->getNumCols () ) );
    int k;

    // Convert tripletMatrix into a (compressed) row matrix with duplicate entries and unsorted rows:
    aol::Vector<IndexType> rowPointer ( this->getNumRows () + 1 );
    aol::Vector<IndexType> colIndex ( tripletCol.size () );
    aol::Vector<DataType> val ( tripletVal.size () );

    // Count num. of entries per row
    for ( IndexType i = 0; i < tripletVal.size (); ++i ) {
      ++( t[tripletRow[i]] );
    }

    // Set row pointers
    for ( int j = 0; j < this->getNumRows (); ++j ) {
      rowPointer[j + 1] = rowPointer[j] + t[j];
      t[j] = rowPointer[j];
    }

    // Fill matrix
    for ( int i = 0; i < tripletVal.size (); ++i ) {
      k = ( t[tripletRow[i]] )++;
      colIndex[k] = tripletCol[i];
      val[k] = tripletVal[i];
    }

    t.setAll ( -1 );

    IndexType p1;
    IndexType p2;
    IndexType pd;
    IndexType pj;
    aol::Vector<IndexType> rowCount ( this->getNumRows () );

    // Sum up duplicate entries
    for ( int i = 0; i < this->getNumRows (); ++i ) {
      p1 = rowPointer[i];
      p2 = rowPointer[i + 1];
      pd = p1;
      for ( IndexType p = p1; p < p2; ++p ) {
        k = colIndex[p];
        pj = t[k];

        // New entry
        if ( t[k] < p1 ) {
          t[k] = pd;
          colIndex[pd] = k;
          val[pd] = val[p];
          ++pd;
        }
        // Entry already seen.
        else {
          val[pj] += val[p];
        }
      }

      rowCount[i] = pd - p1;
    }

    t.setZero ();

    // Convert into a compressed column form matrix
    // Count num. of entries per column.
    for ( int i = 0; i < this->getNumRows (); ++i ) {
      for ( IndexType p = rowPointer[i]; p < rowPointer[i] + rowCount[i]; ++p ) {
        ++( t[colIndex[p]] );
      }
    }

    IndexType numEntries = 0;
    for ( int i = 0; i < t.size (); ++i )
      numEntries += t[i];

    this->_index.resize ( numEntries );
    this->_value.resize ( numEntries );

    // Set column pointers
    this->_indPointer.resize ( this->getNumCols () + 1);
    this->_indPointer[0] = static_cast<IndexType> ( 0 );
    for ( int j = 0; j < this->getNumCols (); ++j ) {
      this->_indPointer[j + 1] = this->_indPointer[j] + t[j];
      t[j] = this->_indPointer[j];
    }

    // Fill matrix
    for ( int i = 0; i < this->getNumRows (); ++i ) {
      for ( IndexType p = rowPointer[i]; p < rowPointer[i] + rowCount[i]; ++p ) {
        k = ( t[colIndex[p]] )++;
        this->_index[k] = i;
        this->_value[k] = val[p];
      }
    }
  }

  //! \warning Untested!
  //! \note Argument is not const because of SparseMatrixRowIterator
  void setFromSparse ( SparseMatrix<DataType> &sparseMatrix ) {
    // Set index pointer to correct size
    this->_indPointer.resize ( sparseMatrix.getNumCols () + 1);

    // helper variables, make t long int, because it needs to be signed. If IndexType = unsigned long int, this does not work!
    aol::Vector<long int> t ( aol::Max ( this->getNumRows (), this->getNumCols () ) );
    int k;

    // Convert into a compressed column form matrix
    // Count num. of entries per column.
    for ( int i = 0; i < this->getNumRows (); ++i ) {
      for ( aol::SparseMatrixRowIterator<DataType> it ( sparseMatrix, i ); it.notAtEnd (); ++it ) {
        ++( t[it->col] );
      }
    }

    IndexType numEntries = 0;
    for ( int i = 0; i < t.size (); ++i )
      numEntries += t[i];

    this->_index.resize ( numEntries );
    this->_value.resize ( numEntries );

    // Set column pointers
    this->_indPointer.resize ( this->getNumCols () + 1);
    this->_indPointer[0] = static_cast<IndexType> ( 0 );
    for ( int j = 0; j < this->getNumCols (); ++j ) {
      this->_indPointer[j + 1] = this->_indPointer[j] + t[j];
      t[j] = this->_indPointer[j];
    }

    // Fill matrix
    for ( int i = 0; i < this->getNumRows (); ++i ) {
      for ( aol::SparseMatrixRowIterator<DataType> it ( sparseMatrix, i ); it.notAtEnd (); ++it ) {
        k = ( t[it->col] )++;
        this->_index[k] = i;
        this->_value[k] = it->value;
      }
    }
  }

public:
  //! \brief Constructor taking the number of rows and columns.
  explicit CSCMatrix ( IndexType numRows, IndexType numCols ) {
    // Set matrix dimensions
    this->_numRows = numRows;
    this->_numCols = numCols;

    // Set column index vector size
    this->_indPointer.resize ( numCols + 1 );
  }

  //! \brief Standard constructor.
  explicit CSCMatrix ()
  : CSMatrix<DataType, IndexType> () {}

  //! \brief Copy constructor.
  CSCMatrix ( const CSCMatrix &other, CopyFlag copyFlag = DEEP_COPY )
  : CSMatrix<DataType, IndexType> ( other, copyFlag ) {}

  //! \brief Constructor taking dimensions and three vectors for column pointers, row indices and values.
  CSCMatrix ( IndexType numRows, IndexType numCols,
      aol::Vector<IndexType> &columnPointer, aol::Vector<IndexType> &rowIndex, aol::Vector<DataType> &value )
  : CSMatrix<DataType, IndexType> () {
    // Set matrix dimensions
    this->_numRows = numRows;
    this->_numCols = numCols;

    this->_indPointer.resize ( columnPointer.size () );
    this->_indPointer = columnPointer;
    this->_index.resize ( rowIndex.size () );
    this->_index = rowIndex;
    this->_value.resize ( value.size () );
    this->_value = value;
  }

  /** \brief Constructor that converts a TripletMatrix into a CSCMatrix.
   *
   * Converts tripletMatrix into a CSCMatrix. If the TripletMatrix has duplicate
   * entries, the contributions are summed up (so it is not necessary
   * to call sumDuplicates on the TripletMatrix).
   */
  CSCMatrix ( const TripletMatrix<DataType> &tripletMatrix ) {
    // Set matrix dimensions
    this->_numRows = tripletMatrix.getNumRows ();
    this->_numCols = tripletMatrix.getNumCols ();

    // Set matrix entries
    setFromTriplet ( tripletMatrix );
  }

  CSCMatrix ( SparseMatrix<DataType> &sparseMatrix ) {
    // Set matrix dimensions
    this->_numRows = sparseMatrix.getNumRows ();
    this->_numCols = sparseMatrix.getNumCols ();

    // Set matrix entries
    setFromSparse ( sparseMatrix );
  }

  //! \brief Destructor.
  virtual ~CSCMatrix () {}

  //! \brief Set the matrix to zero.
  virtual void setZero () {
    this->_index.reallocateClear ( 0 );
    this->_indPointer.reallocateClear ( this->getNumCols () + 1 );
    this->_value.reallocateClear ( 0 );
  }

  CSCMatrix& operator= ( const TripletMatrix<DataType> &tripletMatrix ) {
    this->_numRows = tripletMatrix.getNumRows ();
    this->_numCols = tripletMatrix.getNumCols ();
    this->setZero ();
    setFromTriplet ( tripletMatrix );

    return *this;
  }

  //! \brief Get matrix entry (row, col).
  virtual DataType get ( int row, int col ) const {
    // Look for an entry in the column. The end of column col is this->_indPointer[col + 1]
    for ( IndexType i = this->_indPointer[col]; i < this->_indPointer[col + 1]; ++i ) {
      if ( this->_index[i] == static_cast<IndexType> ( row ) ) {
        return this->_value[i];
      }
    }

    // If no value has been found, the entry is zero.
    return static_cast<DataType> ( 0 );
  }

  void getSubMatrix ( IndexType minRow, IndexType maxRow,
                      IndexType minCol, IndexType maxCol,
                      CSCMatrix<DataType, IndexType> &subMatrix ) const {
    aol::CSMatrix<DataType, IndexType>::getSubMatrix ( minRow, maxRow, minCol, maxCol, false, subMatrix );
  }

  //! \brief Sets entry (row, col) of the block to value.
  //! \attention Changes to the sparsity structure of the matrix are (very) slow!
  virtual void set ( int row, int col, DataType val ) {
    // Look for an entry in the column. The end of column col is this->_indPointer[col + 1].
    for ( IndexType i = this->_indPointer[col]; i < this->_indPointer[col + 1]; ++i ) {
      // If the entry already exists, set it and exit the method.
      if ( this->_index[i] == static_cast<IndexType> ( row ) ) {
        this->_value[i] = val;
        return;
      }
    }

    // Do not insert a new entry if val is 0.
    if ( val == static_cast<DataType> ( 0 ) )
      return;

    // If no value has been found, insert a new entry.
    // Find position to insert new value
    IndexType pos = this->_indPointer[col];
    while ( pos < this->_indPointer[col + 1] ) {
      if ( this->_index[pos] > row )
        break;

      ++pos;
    }

    // And insert new value at this position
    this->_index.insert ( pos, row );
    this->_value.insert ( pos, val );

    // Following column pointers have to be shifted by 1
    for ( IndexType c = col + 1; c < static_cast<IndexType> ( this->_indPointer.size () ); ++c ) {
      ++( this->_indPointer[c] );
    }
  }

  //! \brief Add value to entry (row, col).
  //! \attention Changes to the sparsity structure of the matrix are (very) slow!
  virtual void add ( int row, int col, DataType val ) {
    // Look for an entry in the column. The end of column col is this->_indPointer[col + 1].
    for ( IndexType i = this->_indPointer[col]; i < this->_indPointer[col + 1]; ++i ) {
      // If the entry already exists, set it and exit the method.
      if ( this->_index[i] == static_cast<IndexType> ( row ) ) {
        this->_value[i] += val;
        return;
      }
    }

    // If no value has been found, insert a new entry.
    // Find position to insert new value
    IndexType pos = this->_indPointer[col];
    while ( pos < this->_indPointer[col + 1] ) {
      if ( this->_index[pos] > row )
        break;

      ++pos;
    }

    // Do not insert a new entry if val is 0.
    if ( val == static_cast<DataType> ( 0 ) )
      return;

    // And insert new value at this position
    this->_index.insert ( pos, row );
    this->_value.insert ( pos, val );

    // Following column pointers have to be shifted by 1
    for ( IndexType c = col + 1; c < static_cast<IndexType> ( this->_indPointer.size () ); ++c ) {
      ++( this->_indPointer[c] );
    }
  }
  
  /**
   * \author Tatano
   */
  aol::CSCMatrix< DataType > & operator*= ( const DataType alpha ){
    this->_value *= alpha;
    return *this;
  }
  
  using aol::Matrix<DataType>::operator+=;

  /**
   * \author Tatano
   */
  aol::CSCMatrix< DataType > & operator+= ( const aol::CSCMatrix< DataType > &other ){
    for ( int col = 0; col < other.getNumCols (); ++col ) {
      for ( IndexType i = other._indPointer[col]; i < other._indPointer[col + 1]; ++i ) {
        this->add ( other._index[i], col, other._value[i] ) ;
      }
    }
    return *this;
  }
  
  /**
   * \brief return number of nonzero entries in the matrix (which is NOT the number of stored entries)
   * \author Tatano
   */
  virtual int numNonZeroes ( int I ) const {
    int nNonZeroes = 0;
    for ( int col = 0; col < this->getNumCols (); ++col ) {
      for ( IndexType k = this->_indPointer[col]; k < this->_indPointer[col + 1]; ++k ) {
        if ( this->_index[k] == I && this->_value[k] != 0 )
          ++nNonZeroes;
      }
    }
    return nNonZeroes;
  }

  //! \brief applyAdd method.
  virtual void applyAdd ( const aol::Vector<DataType> &arg, aol::Vector<DataType> &dest ) const {
    // Traverse the matrix column-wise, as this is faster for the CSCMatrix.
    for ( int col = 0; col < this->getNumCols (); ++col ) {
      for ( IndexType i = this->_indPointer[col]; i < this->_indPointer[col + 1]; ++i ) {
        dest[this->_index[i]] += ( this->_value[i] * arg[col] );
      }
    }
  }

  //! \brief apply method.
  virtual void apply ( const aol::Vector<DataType> &arg, aol::Vector<DataType> &dest ) const {
    dest.setZero ();

    // Traverse the matrix column-wise, as this is faster for the CSCMatrix.
    for ( int col = 0; col < this->getNumCols (); ++col ) {
      for ( IndexType i = this->_indPointer[col]; i < this->_indPointer[col + 1]; ++i ) {
        dest[this->_index[i]] += ( this->_value[i] * arg[col] );
      }
    }
  }

  void applyAdd ( const aol::MultiVector<DataType> &arg, aol::MultiVector<DataType> &dest ) const {
    aol::Vector<DataType> rhs( arg.getTotalSize () ), sol( dest.getTotalSize () );
    rhs.copyUnblockedFrom ( arg );
    sol.copyUnblockedFrom ( dest );

    this->applyAdd ( rhs, sol );

    aol::Vector<int> sizes;
    dest.getSizes ( sizes ); 
    dest.copySplitFrom ( sol, sizes );
  }

  void apply ( const aol::MultiVector<DataType> &arg, aol::MultiVector<DataType> &dest ) const {
    aol::Vector<DataType> rhs( arg.getTotalSize () ), sol( dest.getTotalSize () );
    rhs.copyUnblockedFrom ( arg );

    this->apply ( rhs, sol );

    aol::Vector<int> sizes;
    dest.getSizes ( sizes ); 
    dest.copySplitFrom ( sol, sizes );
  }

  const aol::Vector<IndexType>& getRowIndexReference () const {
    return this->_index;
  }

  const aol::Vector<IndexType>& getColumnPointerReference () const {
    return this->_indPointer;
  }

  const aol::Vector<DataType>& getValueReference () const {
    return this->_value;
  }

  template<typename MatrixType>
  void copyTo ( MatrixType &dest ) const {
    dest.setZero ();

    // Traverse the matrix column-wise, as this is faster for the CSCMatrix.
    for ( int col = 0; col < this->getNumCols (); ++col ) {
      for ( IndexType i = this->_indPointer[col]; i < this->_indPointer[col + 1]; ++i ) {
        dest.add ( this->_index[i], col, this->_value[i] );
      }
    }
  }

  void copyTo ( aol::SparseBlockMatrix<aol::SparseMatrix<DataType> > &sparseBlockMatrix ) {
    sparseBlockMatrix.setZero ();

    for ( int col = 0; col < this->getNumCols (); ++col ) {
      for ( IndexType i = this->_indPointer[col]; i < this->_indPointer[col + 1]; ++i ) {
        int rowSum = sparseBlockMatrix.getReference ( 0, 0 ).getNumRows ();
        int br = 0;
        while ( this->_index[i] >= rowSum ) {
          rowSum += sparseBlockMatrix.getReference ( ++br, 0 ).getNumRows ();
        }

        int colSum = sparseBlockMatrix.getReference ( 0, 0 ).getNumCols ();
        int bc = 0;
        while ( col >= colSum ) {
          colSum += sparseBlockMatrix.getReference ( 0, ++bc ).getNumCols ();
        }

        sparseBlockMatrix.getReference ( br, bc ).add ( this->_index[i] - ( rowSum - sparseBlockMatrix.getReference ( br, 0 ).getNumRows () ),
            col - ( colSum - sparseBlockMatrix.getReference ( 0, bc ).getNumCols () ), this->_value[i] );
      }
    }
  }
  
  /*!
   * \brief Sets the matrix to the matrix product of two other CSCMatrices. Dimensions have to match.
   * \note This method uses a TripletMatrix to store the entries of the result matrix and then
   * converts it to CSCMatrix. This is not optimal in terms of memory consumption and
   * probably computation time.
   *
   * \author Luethen, Toelkes
   */
  void makeProduct ( const CSCMatrix<DataType, IndexType> &A, const CSCMatrix<DataType, IndexType> &B ) {
#ifdef BOUNDS_CHECK
    if ( ( A.getNumRows() != this->getNumRows() ) || ( A.getNumCols() != B.getNumRows() ) || ( B.getNumCols() != this->getNumCols() ) )
      throw aol::Exception( "CSCMatrix::makeProduct: Matrix dimensions do not match!", __FILE__, __LINE__ );
#endif
    setZero ();

    // Use pointers for easier/faster access
    // ap, bp column pointer
    // ai, bi row indices
    // ax, bx values
    IndexType *ap = A._indPointer.getData ();
    IndexType *ai = A._index.getData ();
    DataType  *ax = A._value.getData ();

    IndexType *bp = B._indPointer.getData ();
    IndexType *bi = B._index.getData ();
    DataType  *bx = B._value.getData ();

    // Save the result into a TripletMatrix
    TripletMatrix<DataType> res ( this->getNumRows (), this->getNumCols () );
    
    // Traverse 2nd factor column-wise
    for ( int cb = 0; cb < B.getNumCols (); ++cb ) {
      // For each non-zero element...
      for ( IndexType j = bp[cb]; j < bp[cb + 1]; ++j ) {
        IndexType ca = bi[j];
        // ...loop over the corresponding COLUMN of the first factor...
        for ( IndexType i = ap[ca]; i < ap[ca + 1]; ++i ) {
          // ...and sum up the entries in the result column of the result matrix
          res.add ( ai[i], cb, ax[i] * bx[j] );
        }
      }
    }

    // Finally, convert result back to CSCMatrix
    *this = res;
  }
  
  //! \brief Sets the matrix to the transposed of another CSCMatrix. Dimensions have to match.
  //! \author Luethen
  void transposeFrom ( const CSCMatrix<DataType, IndexType> &A ) {
#ifdef BOUNDS_CHECK
    if ( ( A.getNumRows() != this->getNumCols () ) || ( A.getNumCols () != this->getNumRows () ) )
      throw aol::Exception( "CSCMatrix::transposeFrom: Matrix dimensions do not match!", __FILE__, __LINE__ );
#endif
    setZero ();
    
    const IndexType indexSize = A._index.size ();
    //Go through row index vector of A and construct the column index pointer for the transposed matrix
    for ( IndexType i = 0; i < indexSize; ++i )
      ++( this->_indPointer[A._index[i] + 1] );
    for ( IndexType i = 1; i < this->_indPointer.size (); ++i )
      this->_indPointer[i] += this->_indPointer[i-1];
    
    //row index and value vector have the same size in the transposed matrix
    this->_index.resize ( indexSize );
    this->_value.resize ( indexSize );
    
    //Using the old column indices, go through the old row index vector again.
    //Use the old row indices to find the right column in the transposed matrix, and set the new index to the old column index.
    //offset stores how many new row indices belonging to a new column have been set already
    aol::Vector<IndexType> offset ( this->_indPointer.size () - 1 ); 
    for ( IndexType oldColInd = 0; oldColInd < A._indPointer.size () - 1; ++oldColInd )
      for ( IndexType k = A._indPointer[oldColInd]; k < A._indPointer[oldColInd+1]; ++k ) {
        const IndexType newColInd = A._index[k]; //= oldRowInd
        this->_index[this->_indPointer[newColInd] + offset[newColInd]] = oldColInd;
        this->_value[this->_indPointer[newColInd] + offset[newColInd]] = A._value[k];
        ++offset[newColInd];
      }
  }
};

/**
 * \brief A compressed sparse row matrix.
 * \attention Changes to the sparsity structure of the matrix (via add, set, ...) are slow!
 * \warning If IndexType is set to a type bigger than int (or unsigned int), triplet to csr conversion does no longer work (see CSRMatrix::setFromTriplet).
 * \author Toelkes
 * \ingroup Matrix
 *
 * A compressed sparse row matrix. Fast for matrix-vector multiplication and certain arithmetic operations.
 * Changes to the sparsity structure are slow. Use the TripletMatrix for assembling and then convert it into a CSRMatrix.
 */
template <typename DataType, typename IndexType = int>
class CSRMatrix : public CSMatrix<DataType, IndexType> {
  // In this class, _index lists column indices and _indPointer holds row pointers.
protected:
  void setFromTriplet ( const TripletMatrix<DataType> &tripletMatrix ) {
    const Vector<int> &tripletRow = tripletMatrix.getRowIndexReference ();
    const Vector<int> &tripletCol = tripletMatrix.getColIndexReference ();
    const Vector<DataType> &tripletVal = tripletMatrix.getValueReference ();

    // helper variables, make t long int, because it needs to be signed. If IndexType = unsigned long int, this does not work!
    aol::Vector<int64_t> t ( aol::Max ( this->getNumRows (), this->getNumCols () ) );
    int k;

    // Convert tripletMatrix into a (compressed) column matrix with duplicate entries and unsorted columns:
    aol::Vector<IndexType> rowIndex ( tripletRow.size () );
    aol::Vector<IndexType> colPointer ( this->getNumCols () + 1 );
    aol::Vector<DataType> val ( tripletVal.size () );

    // Count number of entries per column
    for ( IndexType j = 0; j < tripletVal.size (); ++j ) {
      ++( t[tripletCol[j]] );
    }

    // Set column pointers
    for ( int i = 0; i < this->getNumCols (); ++i ) {
      colPointer[i + 1] = colPointer[i] + t[i];
      t[i] = colPointer[i];
    }

    // Fill matrix
    for ( int i = 0; i < tripletVal.size (); ++i ) {
      k = ( t[tripletCol[i]] )++;
      rowIndex[k] = tripletRow[i];
      val[k] = tripletVal[i];
    }

    t.setAll ( -1 );

    IndexType p1;
    IndexType p2;
    IndexType pd;
    IndexType pj;
    aol::Vector<IndexType> colCount ( this->getNumCols () );

    // Sum up duplicate entries
    for ( int j = 0; j < this->getNumCols (); ++j ) {
      p1 = colPointer[j];
      p2 = colPointer[j + 1];
      pd = p1;
      for ( IndexType p = p1; p < p2; ++p ) {
        k = rowIndex[p];
        pj = t[k];

        // New entry
        if ( t[k] < p1 ) {
          t[k] = pd;
          rowIndex[pd] = k;
          val[pd] = val[p];
          ++pd;
        }
        // Entry already seen.
        else {
          val[pj] += val[p];
        }
      }

      colCount[j] = pd - p1;
    }

    t.setZero ();

    // Convert into a compressed row form matrix
    // Count num. of entries per row.
    for ( int j = 0; j < this->getNumCols (); ++j ) {
      for ( IndexType p = colPointer[j]; p < colPointer[j] + colCount[j]; ++p ) {
        ++( t[rowIndex[p]] );
      }
    }

    IndexType numEntries = 0;
    for ( int i = 0; i < t.size (); ++i )
      numEntries += t[i];

    this->_index.resize ( numEntries );
    this->_value.resize ( numEntries );

    // Set column pointers
    this->_indPointer.resize ( this->getNumRows () + 1);
    this->_indPointer[0] = static_cast<IndexType> ( 0 );
    for ( int i = 0; i < this->getNumRows (); ++i ) {
      this->_indPointer[i + 1] = this->_indPointer[i] + t[i];
      t[i] = this->_indPointer[i];
    }

    // Fill matrix
    for ( int j = 0; j < this->getNumCols (); ++j ) {
      for ( IndexType p = colPointer[j]; p < colPointer[j] + colCount[j]; ++p ) {
        k = ( t[rowIndex[p]] )++;
        this->_index[k] = j;
        this->_value[k] = val[p];
      }
    }
  }
  
  void setFromSparse ( const aol::SparseMatrix<DataType> &sparseMatrix ) {
    vector<typename aol::Row<DataType>::RowEntry > rowEntries;

    this->_indPointer.resize( sparseMatrix.getNumRows() + 1 );
    this->_indPointer[0] = static_cast<DataType> ( 0 );

    int num = 0;
    for ( int i = 0; i < sparseMatrix.getNumRows(); ++i ) {
      sparseMatrix.makeRowEntries(rowEntries, i);
      for ( typename vector<typename aol::Row<DataType>::RowEntry >::iterator it = rowEntries.begin(); it != rowEntries.end(); ++it ) {
        if ( it->value != 0.0 ) {
          this->_value.pushBack( it->value );
          this->_index.pushBack( it->col );
          ++num;
        }
      }
      this->_indPointer[i+1] = num;
    }
  }

public:
  //! \brief Constructor taking the number of rows and columns.
  CSRMatrix ( IndexType numRows, IndexType numCols ) {
    // Set matrix dimensions
    this->_numRows = numRows;
    this->_numCols = numCols;

    // Set row index vector size
    this->_indPointer.resize ( numRows + 1 );
  }

  //! \brief Standard constructor.
  CSRMatrix ()
  : CSMatrix<DataType, IndexType> () {}

  //! \brief Copy constructor.
  CSRMatrix ( const CSRMatrix &other, CopyFlag copyFlag = DEEP_COPY )
  : CSMatrix<DataType, IndexType> ( other, copyFlag ) {}

  /** \brief Constructor that converts a TripletMatrix into a CSRMatrix.
   *
   * Converts tripletMatrix into a CSRMatrix. If the TripletMatrix has duplicate
   * entries, the contributions are summed up (so it is not necessary
   * to call sumDuplicates on the TripletMatrix).
   */
  CSRMatrix ( const TripletMatrix<DataType> &tripletMatrix ) {
    // Set matrix dimensions
    this->_numRows = tripletMatrix.getNumRows ();
    this->_numCols = tripletMatrix.getNumCols ();

    // Set matrix entries
    setFromTriplet ( tripletMatrix );
  }
  
  /** \brief Constructor that converts a SparseMatrix into a CSRMatrix
   */
  CSRMatrix ( const SparseMatrix<DataType> &sparseMatrix ) {
    // Set matrix dimensions
    this->_numRows = sparseMatrix.getNumRows ();
    this->_numCols = sparseMatrix.getNumCols ();

    // Set matrix entries
    setFromSparse ( sparseMatrix );
  }

  //! \brief Destructor.
  virtual ~CSRMatrix () {}

  //! \brief Set the matrix to zero.
  virtual void setZero () {
    this->_index.reallocateClear ( 0 );
    this->_indPointer.reallocateClear ( this->getNumRows () + 1 );
    this->_value.reallocateClear ( 0 );
  }

  CSRMatrix& operator= ( const TripletMatrix<DataType> &tripletMatrix ) {
    this->setZero ();
    this->_numRows = tripletMatrix.getNumRows ();
    this->_numCols = tripletMatrix.getNumCols ();
    setFromTriplet ( tripletMatrix );

    return *this;
  }

  //! \brief Get matrix entry (row, col).
  virtual DataType get ( int row, int col ) const {
    // Look for an entry in the row. The end of row row is this->_indPointer[row + 1]
    for ( IndexType j = this->_indPointer[row]; j < this->_indPointer[row + 1]; ++j ) {
      if ( this->_index[j] == static_cast<IndexType> ( col ) ) {
        return this->_value[j];
      }
    }

    // If no value has been found, the entry is zero.
    return static_cast<DataType> ( 0 );
  }

  //! \brief Sets entry (row, col) of the matrix to value.
  //! \attention Changes to the sparsity structure of the matrix are (very) slow!
  virtual void set ( int row, int col, DataType val ) {
    // Look for an entry in the row. The end of row row is this->_indPointer[row + 1]
    for ( IndexType j = this->_indPointer[row]; j < this->_indPointer[row + 1]; ++j ) {
      // If the entry already exists, set it and exit the method.
      if ( this->_index[j] == static_cast<IndexType> ( col ) ) {
        this->_value[j] = val;
        return;
      }
    }

    // If no value has been found, insert a new entry.
    // Find position to insert new value
    IndexType pos = this->_indPointer[row];
    while ( pos < this->_indPointer[row + 1] ) {
      if ( this->_index[pos] > col )
        break;

      ++pos;
    }

    // And insert new value at this position
    this->_index.insert ( pos, col );
    this->_value.insert ( pos, val );

    // Following column pointers have to be shifted by 1
    for ( IndexType c = row + 1; c < static_cast<IndexType> ( this->_indPointer.size () ); ++c ) {
      ++( this->_indPointer[c] );
    }
  }

  //! \brief Add value to entry (row, col).
  //! \attention Changes to the sparsity structure of the matrix are (very) slow!
  virtual void add ( int row, int col, DataType val ) {
    // Look for an entry in the row. The end of row row is this->_indPointer[row + 1]
    for ( IndexType j = this->_indPointer[row]; j < this->_indPointer[row + 1]; ++j ) {
      // If the entry already exists, set it and exit the method.
      if ( this->_index[j] == static_cast<IndexType> ( col ) ) {
        this->_value[j] += val;
        return;
      }
    }

    // If no value has been found, insert a new entry.
    // Find position to insert new value
    IndexType pos = this->_indPointer[row];
    while ( pos < this->_indPointer[row + 1] ) {
      if ( this->_index[pos] > col )
        break;

      ++pos;
    }

    // And insert new value at this position
    this->_index.insert ( pos, col );
    this->_value.insert ( pos, val );

    // Following column pointers have to be shifted by 1
    for ( IndexType c = row + 1; c < static_cast<IndexType> ( this->_indPointer.size () ); ++c ) {
      ++( this->_indPointer[c] );
    }
  }

  //! \brief applyAdd method.
  virtual void applyAdd ( const aol::Vector<DataType> &arg, aol::Vector<DataType> &dest ) const {
    // Traverse the matrix row-wise.
    for ( int row = 0; row < this->getNumRows (); ++row ) {
      for ( IndexType j = this->_indPointer[row]; j < this->_indPointer[row + 1]; ++j ) {
        dest[row] += ( this->_value[j] * arg[this->_index[j]] );
      }
    }
  }

  //! \brief apply method.
  virtual void apply ( const aol::Vector<DataType> &arg, aol::Vector<DataType> &dest ) const {
    // Traverse the matrix row-wise.
    for ( int row = 0; row < this->getNumRows (); ++row ) {
      DataType s = static_cast<DataType> ( 0 );
      for ( IndexType j = this->_indPointer[row]; j < this->_indPointer[row + 1]; ++j ) {
        s += ( this->_value[j] * arg[this->_index[j]] );
      }
      dest[row] = s;
    }
  }

  void applyAdd ( const aol::MultiVector<DataType> &arg, aol::MultiVector<DataType> &dest ) const {
    aol::Vector<DataType> rhs( arg.getTotalSize () ), sol( dest.getTotalSize () );
    rhs.copyUnblockedFrom ( arg );
    sol.copyUnblockedFrom ( dest );

    this->applyAdd ( rhs, sol );

    aol::Vector<int> sizes;
    arg.getSizes ( sizes );
    dest.copySplitFrom ( sol, sizes );
  }

  void apply ( const aol::MultiVector<DataType> &arg, aol::MultiVector<DataType> &dest ) const {
    aol::Vector<DataType> rhs( arg.getTotalSize () ), sol( dest.getTotalSize () );
    rhs.copyUnblockedFrom ( arg );

    this->apply ( rhs, sol );

    aol::Vector<int> sizes;
    arg.getSizes ( sizes );
    dest.copySplitFrom ( sol, sizes );
  }

  const aol::Vector<IndexType>& getRowPointerReference () const {
    return this->_indPointer;
  }
  
  aol::Vector<IndexType>& getRowPointerReference () {
    return this->_indPointer;
  }

  const aol::Vector<IndexType>& getColumnIndexReference () const {
    return this->_index;
  }
  
  aol::Vector<IndexType>& getColumnIndexReference () {
    return this->_index;
  }

  const aol::Vector<DataType>& getValueReference () const {
    return this->_value;
  }
  
  aol::Vector<DataType>& getValueReference () {
    return this->_value;
  }
  
  //! \author Tatano
  void makeRowEntries ( vector<typename Row<DataType>::RowEntry > &vec, const int RowNum ) const {
    vec.resize ( this->_indPointer[RowNum+1] - this->_indPointer[RowNum] );
    for ( IndexType j = this->_indPointer[RowNum]; j < this->_indPointer[RowNum + 1]; j++ ) {
      vec[j-this->_indPointer[RowNum]].col = this->_index[j];
      vec[j-this->_indPointer[RowNum]].value = this->_value[j];
    }
  }
  
  //! \author Tatano
  aol::CSRMatrix< DataType > & operator*= ( const DataType alpha ){
    this->_value *= alpha;
    return *this;
  }
  
  //! \author Tatano
  void getVectorWithNormRows (aol::Vector<DataType> &Dvec) const {
    for ( int row = 0; row < this->getNumRows (); ++row ) {
      vector<typename Row<DataType>::RowEntry > vec;
      this->makeRowEntries(vec, row);
      Dvec[row] = vec.norm();
    }
  }
  
  //! \author Tatano
  void scaleMatrixRowwiseTo(aol::CSRMatrix< DataType > &scalMat) const {
    aol::Vector<DataType> Dvec;
    this->getVectorWithNormRows (Dvec);
    for ( int row = 0; row < this->getNumRows (); ++row ) {
      for ( IndexType j = this->_indPointer[row]; j < this->_indPointer[row + 1]; j++ ) {
        this->_value[j] *= 1./Dvec[row];
      }
    }
  }
  
  //! \author Tatano
  void scaleRow (const int Row, const DataType factor){
    for ( IndexType j = this->_indPointer[Row]; j < this->_indPointer[Row + 1]; j++ ) {
      this->_value[j] *= factor;
    }
  }
  
  //! \author Tatano
  void scaleCol (const int Col, const DataType factor) {
    for ( int row = 0; row < this->getNumRows (); ++row )  {
      if (this->_indPointer[row] == Col) {
        this->_value[Col] *= factor;
      }
    }
  }
  
  //! \author Tatano
  aol::CSRMatrix< DataType > & addMultiple ( const aol::CSRMatrix< DataType > &Matrix, const DataType factor ) {
    cerr << this->_index.size() << endl;
    for ( int row = 0; row < Matrix.getNumRows(); ++row ) {
      for ( IndexType j = Matrix._indPointer[row]; j < Matrix._indPointer[row + 1]; j++ ) {
        DataType value = 0;
        if (this->_indPointer[row] == j)
          value = this->_value[j] + factor * Matrix._value[j];
        else
          value = factor * Matrix._value[j];
        this->set(row, this->_indPointer[row], value);
      }
    }
    return ( *this );
  }
  
  //! \author Tatano
  template<typename MatrixType>
  void copyTo ( MatrixType &dest ) const {
    dest.setZero ();
    for ( int row = 0; row < this->getNumRows (); ++row ) {
      for ( IndexType j = this->_indPointer[row]; j < this->_indPointer[row + 1]; ++j ) {
        dest.add ( row, this->_index[j], this->_value[j] );
      }
    }
  }
};

}

#endif
