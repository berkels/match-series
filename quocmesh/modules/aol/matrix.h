#ifndef __MATRIX_H
#define __MATRIX_H

#include <aol.h>
#include <vec.h>
#include <multiVector.h>
#include <op.h>
#include <rows.h>
#include <qmException.h>
#include <vectorExtensions.h>

// forward declaration
namespace qc {
template <typename DataType, Dimension Dim> class ScalarArray;
}

namespace aol {

/**
 * Tranposes matrix A using makeRowEntries and stores A^T into B.
 *
 * \author Berkels
 */
template <typename MatrixType, typename DataType>
void transposeAToB ( const MatrixType &A, MatrixType &B ) {
  B.setZero();
  vector<typename Row<DataType>::RowEntry > vec;
  for ( int i = 0; i < A.getNumRows(); ++i ) {
    A.makeRowEntries ( vec, i );
    for ( typename vector<typename Row<DataType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
      B.set ( it->col, i, it->value );
    }
  }
}

/**
 * \author Berkels
 */
template <typename MatrixType, typename DataType>
int maxNumNonZeroesPerRow ( const MatrixType &A ) {
  int maxNum = 0;
  for ( int i = 0; i < A.getNumRows(); ++i )
    maxNum = aol::Max ( maxNum, A.numNonZeroes ( i ) );
  return maxNum;
}

/**
 * Writes a matrix into an ostream in the Octave sparse format.
 *
 * \author Berkels, Schwen
 */
template <typename MatrixType, typename DataType>
ostream& printSparseOctave ( const MatrixType &A, ostream& os ) {
  os << "A = sparse(" << A.getNumRows() << "," <<  A.getNumCols() << ");" << endl;
  vector<typename Row<DataType>::RowEntry > vec;
  for ( int i = 0; i < A.getNumRows(); ++i ) {
    A.makeRowEntries ( vec, i );
    for ( typename vector<typename Row<DataType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
      if ( it->value != ZOTrait<DataType>::zero ) {
        if ( isFinite ( it->value ) ) {
          // note that octave's indexing starts from 1
          os << "A (" << i + 1 << ", " << it->col + 1 << ") = " << longScientificFormat ( it->value ) << ";" << endl;
        } else
          cerr << it->value << " found at position (" << i << "," << it->col << ")!\n" << endl;
      }
    }
    os << endl;
  }
  return os;
}


/**
 * Calculates the condition number of a matrix using Octave.
 *
 * \author Berkels, Schwen
 */
template <typename MatrixType, typename DataType>
void computeConditionNumberViaOctave ( const MatrixType &Mat, const bool CalcTime = false ) {
  char octaveDatTempFileName[1024];
  sprintf ( octaveDatTempFileName, "octave.datXXXXXX" );
  ofstream octavedat;
  generateTemporaryFile ( octaveDatTempFileName, octavedat );

  printSparseOctave<MatrixType, DataType> ( Mat, octavedat );

  octavedat << "cond ( A ) " << endl;
  octavedat.close();

  char systemCommand[1024];
  sprintf ( systemCommand, "%soctave -q %s", ( CalcTime ? "time " : "" ), octaveDatTempFileName );
  system ( systemCommand );

  remove ( octaveDatTempFileName );
}

/**
 * Checks whether any of the entries are finite or not. Returns true
 * if any non finite entries (i.e. inf or nan) are found.
 *
 * \author Berkels
 */
template <typename MatrixType, typename DataType>
bool checkForNANsAndINFs ( const MatrixType &A ) {
  vector<typename Row<DataType>::RowEntry > vec;
  for ( int i = 0; i < A.getNumRows(); ++i ) {
    A.makeRowEntries ( vec, i );
    for ( typename vector<typename Row<DataType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
      if ( !isFinite ( it->value ) )
        return true;
    }
  }
  return false;
}


/**
 * Checks whether any of the entries are finite or not. Returns true
 * if any non finite entries (i.e. inf or nan) are found.
 *
 * \note: Calls get on every matrix entry. So only use this for
 * non-sparse matrices.
 *
 * \author Berkels
 */
template <typename MatrixType>
bool checkForNANsAndINFsInFullMatrix ( const MatrixType &A ) {
  for ( int i = 0; i < A.getNumRows(); ++i ) {
    for ( int j = 0; j < A.getNumCols(); ++j ) {
      if ( !aol::isFinite ( A.get ( i, j ) ) )
        return true;
    }
  }
  return false;
}

/** determine Linf difference of two Ops via makeRowEntries
 *  \author Schwen
 */
template< typename OpType1, typename OpType2, typename RealType >
RealType rowwiseOpLinfDifference ( const OpType1 &Op1, const OpType2 &Op2, const int NumRows ) {
  RealType ret = - aol::NumberTrait<RealType>::Inf;

  for ( int i = 0; i < NumRows; ++i ) {
    std::map< int, RealType > rowEntryDifference;
    vector< typename aol::Row<RealType>::RowEntry > rowEntries1, rowEntries2;

    Op1.makeRowEntries ( rowEntries1, i );
    Op2.makeRowEntries ( rowEntries2, i );

    for ( typename vector< typename aol::Row<RealType>::RowEntry >::const_iterator it = rowEntries1.begin(); it != rowEntries1.end(); ++it )
      rowEntryDifference[ it->col ] += it->value;

    for ( typename vector< typename aol::Row<RealType>::RowEntry >::const_iterator it = rowEntries2.begin(); it != rowEntries2.end(); ++it )
      rowEntryDifference[ it->col ] -= it->value;

    for ( typename std::map< int, RealType >::const_iterator it = rowEntryDifference.begin(); it != rowEntryDifference.end(); ++it )
      if ( fabs ( it->second ) > ret )
        ret = fabs ( it->second );
  }

  return ( ret );
}

/** check whether a matrix is symmetric up to a threshold
 *  \author Schwen
 */
template <typename MatrixType, typename DataType>
bool isMatrixSymmetric ( const MatrixType &A, const DataType thres = aol::NumberTrait<DataType>::zero ) {
  vector<typename Row<DataType>::RowEntry > vec;
  for ( int rowN = 0; rowN < A.getNumRows(); ++rowN ) {
    A.makeRowEntries ( vec, rowN );
    for ( typename vector<typename Row<DataType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
      if ( fabs ( it->value - A.get ( it->col, rowN ) ) > thres ) {
        return ( false );
      }
    }
  }

  return ( true );
}


template <class DomType, class RangeType> class Op;
template <typename DataType> class SparseMatrix;

template <typename RealType>
class GenSparseOp : public aol::Op<aol::Vector<RealType> > {
public:
  GenSparseOp( const int Rows, const int Columns ) : _numRows ( Rows ), _numCols ( Columns ) { }
  virtual ~GenSparseOp () {}
  //! Return vector of row entries. This function must be overloaded on derived matrix classes.
  //! Entries need not be sorted with respect to column index and zeros may be contained.
  virtual void makeRowEntries ( vector<typename aol::Row<RealType>::RowEntry > &vec, const int RowNum ) const = 0;

  //! Returns number of rows.
  int getNumRows() const {
    return _numRows;
  }

  //! Returns number of columns.
  int getNumCols() const {
    return _numCols;
  }

protected:
  int _numRows, _numCols;
};

/**
 * Similar to aol::CompositeOp but allows to use non-quadratic ops at the price of only working on
 * aol::GenSparseOps and not on general aol::Ops.
 *
 * \todo Fix code duplication with aol::CompositeOp.
 *
 * \author Berkels
 */
template <class DataType>
class CompositeGenSparseOp : public Op<aol::Vector<DataType> > {
public:
  CompositeGenSparseOp() {}

  virtual ~CompositeGenSparseOp() {}

  void appendReference ( const GenSparseOp<DataType> &Op ) {
    ops.push_back ( &Op );
  }

  void reset() {
    ops.erase ( ops.begin(), ops.end() );
  }

  virtual void applyAdd ( const aol::Vector<DataType> &Arg, aol::Vector<DataType> &Dest ) const {
    aol::Vector<DataType> tmp ( Dest );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  virtual void apply ( const aol::Vector<DataType> &Arg, aol::Vector<DataType> &Dest ) const {
    aol::Vector<DataType> tmp ( Dest );

    Dest = Arg;
    for ( typename vector<const GenSparseOp<DataType>* >::const_iterator it = ops.begin(); it != ops.end(); ++it ) {
      if ( tmp.size() != Dest.size() )
        tmp.reallocate ( Dest );

      tmp = Dest;

      if ( ( *it )->getNumRows() != Dest.size() )
        Dest.reallocate ( ( *it )->getNumRows() );

      ( *it )->apply ( tmp, Dest );
    }
  }
protected:
  vector<const GenSparseOp<DataType>* > ops;
};

template <class _DataType> class Matrix;

template <typename RealType>
class MatrixAbstractBase : public aol::GenSparseOp<RealType> {
public:
  typedef RealType DataType;
  MatrixAbstractBase( const int Rows, const int Columns ) : aol::GenSparseOp<RealType>( Rows, Columns ) {}

  virtual ~MatrixAbstractBase() {}

  virtual void makeRowEntries ( vector<typename aol::Row<RealType>::RowEntry > &/*vec*/, const int /*RowNum*/ ) const{
    throw aol::Exception ( "aol::MatrixAbstractBase::makeRowEntries must be overloaded on derived classes where you want to use it.", __FILE__, __LINE__ );
  }
  //! Sets matrix entries to zeros.
  virtual void setZero ( ) = 0;
  virtual int numNonZeroes ( int I ) const = 0;
  virtual int numStoredEntries ( const int /*I*/ ) const {
    throw ( aol::UnimplementedCodeException ( "aol::MatrixAbstractBase::numStoredEntries() not implemented", __FILE__, __LINE__ ) );
  }

  virtual RealType getDiag ( int I ) const = 0;
};


//! Abstract Matrix class.
template <class _DataType>
class Matrix : public MatrixAbstractBase<_DataType> {
public:
  typedef _DataType DataType;
public:
  Matrix ( const int Rows, const int Columns ) : MatrixAbstractBase<_DataType>( Rows, Columns ) { }

  Matrix() : MatrixAbstractBase<_DataType>( 0, 0 ){ }

  virtual Matrix* clone ( CopyFlag ) const {
    throw aol::Exception ( "aol::Matrix::clone must be overloaded on derived classes.", __FILE__, __LINE__ );
  }

  virtual ~Matrix() { }

  //! Set row to zero.
  virtual void setRowToZero ( const int row );

  //! Set column to zero.
  virtual void setColToZero ( const int col );

  //! Set row and column to zero.
  virtual void setRowColToZero ( const int rowCol );

  //! Returns a matrix entry.
  virtual DataType get ( int i, int j ) const = 0;

  //! Returns diagonal entry. Should be reimplemented for fast access to diagonal entries.
  virtual DataType getDiag ( int i ) const {
    return ( get ( i, i ) );
  }

  //! Sets a matrix entry.
  virtual void set ( int i, int j, DataType value ) = 0;

  //! Set diagonal entry. Should be reimplemented for fast access to diagonal entries.
  virtual void setDiag ( int i, DataType value ) {
    set ( i, i, value );
  }

  //! Adds a value to a matrix element
  virtual void add ( int i, int j, DataType value ) {
    set ( i, j, get ( i, j ) + value );
  }

  //! Writes the matrix, uses get.
  virtual ostream& print ( ostream& out ) const;

  //! Prints entries of matrix (output may contain zero entries)
  virtual ostream& printSparse ( ostream& out = cout, const aol::Format DataFormatter = aol::mixedFormat ) const;

  //! Prints entries of matrix in octave sparse matrix format (output may contain zero entries)
  virtual ostream& printSparseOctave ( ostream& out = cout ) const;

  //! Reads a matrix, uses reallocate & set.
  virtual istream& read ( istream& in );

  //! Writes the matrix header with memory allocation data
  virtual void printHead ( ostream& out ) const;

  //! Reads the matrix Header, reallocating memory
  virtual void readHead ( istream& in );

  //! Writes the j-th Column into the Vector v.
  virtual void getColumn ( int j, Vector<DataType> &v ) const;

  //! Writes the c-th column into the Vector consisting of the concatenation of all components of Column.
  virtual void getColumn ( int c, MultiVector<DataType>& Column ) const;

  //! Writes the i-th row into the Vector v.
  virtual void getRow ( int i, Vector<DataType> &v ) const;

  //! Writes the j-th column from the Vector v.
  virtual void setColumn ( int j, const Vector<DataType> &v );

  //! Writes the c-th column from the Vector consisting of the concatenation of all components of Column.
  virtual void setColumn ( int c, const MultiVector<DataType>& Column );

  //! Adds the Vector v multiplied by Factor to the j-th column.
  virtual void addColumn ( int j, const Vector<DataType> &v, const DataType Factor = ZOTrait<DataType>::one );

  //! Adds the Vector consisting of the concatenation of all components of Column multiplied by Factor to the c-th column.
  virtual void addColumn ( int c, const MultiVector<DataType>& Column, const DataType Factor = ZOTrait<DataType>::one  );

  //! Writes the i-th row from the Vector v.
  virtual void setRow ( int i, const Vector<DataType> &v );

  //! Slow matrix-vector multiplication - redefine on sparse representations!
  void applyAdd ( const Vector<DataType> &Arg, Vector<DataType> &Dest ) const;

  //! returns maximal entry in matrix
  //! \todo reimplement on FullMatrix
  DataType getMaxValue() const ;

  //! Returns the square of the Frobenius norm, i.e. the sum of the squared matrix entries.
  DataType getFrobeniusNormSqr() const;

  /** Creates a non-zero pattern in an image
   *  Zeros have a zeroGray dot whereas non-zeros have a nonZeroGray dot.
   */
  qc::ScalarArray<unsigned char, qc::QC_2D>* createNonZeroPattern ( const int zoomFactor = 1,
                                                                    const int zeroGray = 255,
                                                                    const int nonZeroGray = 0 ) const ;

  /** Writes a non-zero pattern image to the specified file
   *  Zeros have a zeroGray dot whereas non-zeros have a nonZeroGray dot.
   */
  void writeNonZeroPattern ( const char *fileName, const int zoomFactor = 1,
                             const int zeroGray = 255, const int nonZeroGray = 0 ) const ;

  //! Scale row by factor. This is the generic slow implementation.
  virtual void scaleRow ( const int RowNum, const DataType Factor ) {
    std::vector < typename Row < _DataType >::RowEntry > vec;
    this->makeRowEntries ( vec, RowNum );
    for ( typename vector<typename Row<_DataType>::RowEntry >::const_iterator it = vec.begin(); it != vec.end(); ++it ) {
      this->set ( RowNum, it->col, Factor * ( it->value ) );
    }

  }

  //! Multiplies all entries with the scalar alpha. This is the generic slow implementation.
  virtual Matrix<DataType>& operator*= ( const DataType alpha );

  //! Add other matrix to this matrix. This is the generic slow implementation.
  virtual Matrix<DataType>& operator+= ( const Matrix<DataType>& other );

  //! Subtract other matrix from this matrix. This is the generic slow implementation.
  virtual Matrix<DataType>& operator-= ( const Matrix<DataType>& other );

  /** Adds to this matrix a matrix (factor*mat) at the position I,J.
   *
   *  It adds factor * mat to this. It adds either the entire matrix mat or
   *  the upper (this->_numRows - I, this->_numCols - J) part of it, depending on
   *  which is smaller. Thus, the function always remains in the
   *  range of both matrices.
   *
   *  If for I or J a value < 0 is passed, it is set to 0 after computing
   *  the range to copy. That means, you can use this function to copy
   *  only a left upper part of mat to the left upper part of *this:
   *  if mat is a \f$ n times m \f$ matrix, the call
   *  applyMultipleAtPosition(i-n, j-m, mat) copies only the first \f$ i \f$ rows
   *  and \f$ j \f$ columns to the left upper part of *this (or even less,
   *  if *this is smaller than \f$ i \times j \f$).
   *
   *  ATTENTION: This function does not do what you would expect when applied
   *  to a MaskedOp. Internally, a masked matrix stores values that are
   *  overriden when you call MaskedOp::apply, which performs applyMasked
   *  on its underlying operator. These overridden values are simply
   *  copied in this function.
   */
  virtual void addMultipleAtPosition  ( int I, int J, const Matrix<DataType> &mat, DataType factor = 1 );

  //! Sets entries to matrix \f$ (\delta_{ij})_{ij} \f$.
  virtual void setIdentity();

  //! Adds \f$ (\delta_{ij})_{ij} \f$ to the matrix.
  virtual void addIdentity();

  //! \todo overload on FullMatrix
  virtual DataType rowSum ( const int I ) const;

  DataType columnSum ( const int J ) const {
    DataType s = ZOTrait<DataType>::zero;
    for ( int i = 0; i < this->_numRows; ++i ) {
      s += get ( i, J );
    }
    return s;
  }


  //! writes the transpose of this matrix to other matrix
  //! this is a general, slow implementation. can do better on sparse matrices or other types!
  void transposeTo ( Matrix< DataType > &Dest ) const;

  /** Same as makeRowEntries, but the returned vector may only contain nonzero entries.
      \warning This is a generic slow implementation, overload by more efficient method on derived classes if necessary
      \author Schwen
  */
  virtual void makeNonZeroRowEntries ( std::vector< typename Row<DataType>::RowEntry > &vec, const int RowNum ) const;

  //! multiplies two Matrices \f$( this += M_1 \cdot M_2)\f$ and adds the result to the current matrix
  void addMatrixProduct ( const Matrix<DataType> &M1, const Matrix<DataType> &M2 );

  //! change size of matrix, preserving data
  virtual void resize ( const int, const int );

  //! change size of matrix, setting all data to zero
  virtual void reallocate ( const int, const int );

  // increase size of matrix
  void growBy ( const int addl_rows, const int addl_cols ) {
    resize ( this->_numRows + addl_rows, this->_numCols + addl_cols );
  }

  //! Set all entries of Matrix to same value. This is a generic slow implementation! It fails for matrices that have a fixed sparsity structure.
  virtual void setAll ( const DataType value ) {
    for ( int i = 0; i < this->_numRows; ++i )
      for ( int j = 0; j < this->_numCols; ++j )
        set ( i, j, value );
  }

  //! return number of nonzero entries in the matrix (which is NOT the number of stored entries)
  virtual int numNonZeroes ( int I ) const {
    std::vector<typename Row<_DataType>::RowEntry > vec;
    int nNonZeroes = 0;
    this->makeRowEntries ( vec, I );
    for ( typename vector<typename Row<_DataType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
      if ( it->value != 0 )
        ++nNonZeroes;
    }
    return nNonZeroes;
  }

  virtual int maxNumNonZeroesPerRow () const {
    int maxNum = 0;
    for ( int i = 0; i < this->getNumRows(); i++ )
      maxNum = aol::Max ( maxNum, numNonZeroes ( i ) );
    return maxNum;
  }

  virtual int numStoredEntries ( ) const {
    throw aol::UnimplementedCodeException ( "aol::Matrix::numStoredEntries() not implemented, implement on derived subclasses!", __FILE__, __LINE__ );
    return ( -1 );
  }

  //! return number of stored entries in the matrix
  virtual int numStoredEntries ( const int I ) const {
    std::vector<typename Row<_DataType>::RowEntry > vec;
    this->makeRowEntries ( vec, I );

    return static_cast<int> ( vec.size() );
  }

  virtual bool checkForNANsAndINFs() const {
    throw aol::UnimplementedCodeException ( "aol::Matrix::checkForNANsAndINFs() not implemented, implement on derived subclasses!", __FILE__, __LINE__ );
    return false;
  }

  virtual void addTensorProductMultiple ( const Vector < DataType > &, const Vector < DataType > &, const DataType ) {
    throw aol::UnimplementedCodeException ( "aol::Matrix::addTensorProductMultiple() not implemented, implement on derived subclasses!", __FILE__, __LINE__ );
  }

  inline virtual void addTensorProduct ( const Vector < DataType > &, const Vector < DataType > & ) {
    throw aol::UnimplementedCodeException ( "aol::Matrix::addTensorProduct() not implemented, implement on derived subclasses!", __FILE__, __LINE__ );
  }

protected:
  virtual void matrixInit ( const int Rows, const int Columns ) {
    this->_numRows = Rows;
    this->_numCols = Columns;
  }

#ifdef BOUNDS_CHECK
  inline bool boundsCheck ( const int i, const int j, const char* msg, const char* fi, const int li ) const {
    const bool isIn = ( i >= 0 && i < static_cast<int> ( this->_numRows ) && j >= 0 && j < static_cast<int> ( this->_numCols ) );
    if ( !isIn ) {
      char errmsg[1024];
      sprintf ( errmsg, "%s %d %d, bounds %d %d", msg, i, j, this->_numRows, this->_numCols );
      throw OutOfBoundsException ( errmsg, fi, li );
    }
    return ( isIn );
  }
#endif

private:
  aol::Matrix<DataType>& operator= ( const aol::Matrix<DataType>& ); // do not implement (derived classes check whether size matches and copy contents)

protected:

  static const Format& format;
  static bool prettyFormat;
};

//! Should be replaced by template from surfmesh/generalstuff.h
template <class _DataType> ostream& operator<< ( ostream& os, const Matrix<_DataType> &m ) {
  return m.print ( os );
}
template <class _DataType> istream& operator>> ( istream& is, Matrix<_DataType> &m ) {
  return m.read ( is );
}

// ============================================================================================

/**  \brief A Permutation-Matrix class.
  *
  *   \author Droske, von Deylen
  *  \ingroup Matrix
  */
template <class ArgDataType>
class PermutationMatrix : public Matrix<ArgDataType> {

private:
  // hide default constructor.
  PermutationMatrix() {}
  //! do not allow set(...) inherited from Matrix.
  void set ( int, int, ArgDataType ) {}
  //! do not allow add(...) inherited from Matrix.
  void add ( int, int, ArgDataType ) {}

public:
  explicit PermutationMatrix ( int Dimension ) {
    reallocate ( Dimension );
  }

  explicit PermutationMatrix ( const Vector<int> & rhs ) {
    reallocate ( rhs.size () );
    setPermutationVector ( rhs );
  }

  //! Resize Matrix, deleting old contents
  void reallocate ( const int Dimension ) {
    this->_numRows = this->_numCols = Dimension;
    data.reallocate ( Dimension );
    for ( int i = 0; i < data.size(); ++i )
      data[i] = i;
  }

  //! Resize Matrix, deleting old contents
  void reallocate ( const int Dim1, const int Dim2 ) {
    if ( Dim1 == Dim2 )
      reallocate ( Dim1 );
    else
      throw Exception ( "PermutationMatrix::reallocate called for non-square format" );
  }

  ArgDataType get ( int I, int J ) const {
    if ( ( I >= 0 ) && ( J >= 0 ) && ( I < this->_numRows ) && ( J < this->_numCols ) )
      return this->data[ I ] == J ? ZOTrait<ArgDataType>::one : ZOTrait<ArgDataType>::zero;
    else
      throw ( "PermutationMatrix::get(...): indices out of range." );
  }

  //! adjusts the matrix, so that x[i] and x[j] are swapped on multiplication.
  //! i.e. p := p * tau, where p is the permutation defined by the matrix and tau the transposition of indices i and j.
  void makeSwap ( int I, int J ) {
    swap ( data[I], data[J] );
  }

  //! Can permute any vector class
  void apply ( const Vector<ArgDataType> &Arg, Vector<ArgDataType> &Dest ) const {
    if ( Arg.size() != this->_numRows || Dest.size() != this->_numRows ) {
      cerr << "PermutationMatrix::apply(...): Source or Destination dimensions do not match with matrix dimension.";
      return;
    }
    for ( int i = 0; i < this->_numRows; ++i )
      Dest.set ( i, Arg.get ( data[i] ) );
  }

  //! Can permute any vector class
  void applyTransposed ( const Vector<ArgDataType> &Arg, Vector<ArgDataType> &Dest ) const {
    if ( Arg.size() != this->_numRows || Dest.size() != this->_numRows ) {
      cerr << "PermutationMatrix::apply(...): Source or Destination dimensions do not match with matrix dimension.";
      return;
    }
    for ( int i = 0; i < this->_numRows; ++i )
      Dest.set ( data[i], Arg.get ( i ) );
  }

  void setPermutationVector ( const Vector<int> & rhs ) {
    data = rhs;
  }

  const Vector<int> & getPermutationVectorRef () const {
    return data;
  }

  void setZero ( ) {
    throw UnimplementedCodeException ( "aol::PermutationMatrix::setZero not implemented (how? ;-)", __FILE__, __LINE__ );
  }

  void makeRowEntries ( vector<typename Row<ArgDataType>::RowEntry > &vec, const int RowNum ) const {
    vec.resize ( 1 );
    vec[0].col = data[RowNum];;
    vec[0].value = 1;
  }

  aol::PermutationMatrix<ArgDataType>& operator= ( const aol::PermutationMatrix<ArgDataType> &other ) {
#ifdef BOUNDS_CHECK
    if ( ( other.getNumRows() != this->_numRows ) || ( other.getNumCols() != this->_numCols ) ) {
      cerr << "ERROR in aol::PermutationMatrix<_DataType>::operator=() : check dimensions.." << endl;
      return *this;
    }
#endif
    this->data = other.data;
    return ( *this );
  }

protected:
  Vector<int> data;
};

// ============================================================================================

/** \brief In a FullMatrix, all elements are stored.
 *  \author Droske
 *  \ingroup Matrix
 */
template < class _DataType, bool Parallelize = true >
class FullMatrix : public Matrix <_DataType> {
protected:
  Vector<_DataType> data;

public:

  FullMatrix();

  //! MxN matrix
  FullMatrix ( int m, int n );

  //! Create a 1xN matrix from a vector
  explicit FullMatrix ( const Vector<_DataType>& V, bool Transpose = false );

  //! create m by n matrix with entries copied from Vec. The CopyFlag value is passed to the constructor of the internal data vector.
  FullMatrix ( const aol::Vector<_DataType> &Vec, const int m, const int n, CopyFlag copyFlag = aol::DEEP_COPY );

  //! create FullMatrix for GridStructure (you probably do not want this for not-very-small grids)
  explicit FullMatrix ( const qc::GridStructure &grid ) : Matrix<_DataType> ( grid.getNumberOfNodes(), grid.getNumberOfNodes() ) {
    data.reallocate ( this->_numRows * this->_numCols );
  }

  //! Convert linear operator to full matrix, given both dimensions
  FullMatrix ( const Op<Vector<_DataType> >& op, int n, int m );

  explicit FullMatrix ( const Matrix<_DataType>& m );

  FullMatrix ( const FullMatrix<_DataType, Parallelize>& m, CopyFlag copyFlag = aol::DEEP_COPY ); // FullMatrix is used as return value elsewhere, thus cannot be explicit.

  virtual ~FullMatrix() {
    // Destructor of Vector will be called automatically
  }

  //! resize Matrix, deleting old contents
  void reallocate ( const int m, const int n ) {
    this->_numRows = m;
    this->_numCols = n;
    data.reallocate ( m * n );
  }

  //! Writes the matrix header with memory allocation data
  virtual void printHead ( ostream& out ) const;

  //! Reads the matrix Header, reallocating memory
  virtual void readHead ( istream& in );

  //! Get entry (i,j).
  //! Slow due to being virtual
  inline _DataType get ( int i, int j ) const {
    return ref ( i, j );
  }

  //! Returns reference to entry (i,j).
  //! Use this method for access within FullMatrix
  _DataType& ref ( int i, int j ) {
#ifdef BOUNDS_CHECK
    this->boundsCheck ( i, j, "FullMatrix::ref: Index out of bounds", __FILE__, __LINE__ );
#endif
    return data [i*this->_numCols + j];
  }

  //! Returns const reference to entry (i,j).
  //! Use this method for access within FullMatrix
  const _DataType& ref ( int i, int j ) const {
#ifdef BOUNDS_CHECK
    this->boundsCheck ( i, j, "FullMatrix::get: Index out of bounds", __FILE__, __LINE__ );
#endif
    return data [i*this->_numCols + j];
  }

  //! Multiindex access
  const _DataType& operator [] ( const Vec<2, int>& i ) const {
    return ref ( i [0], i [1] );
  }

  //! Multiindex access
  _DataType& operator [] ( const Vec<2, int>& i ) {
    return ref ( i [0], i [1] );
  }

  //! Set entry (i,j).
  //! Slow due to being virtual
  inline void set ( int i, int j, _DataType value ) {
    ref ( i, j ) = value ;
  }

  //! Add to entry (i,j).
  //! Slow due to being virtual
  inline void add ( int i, int j, _DataType value ) {
    ref ( i, j ) += value ;
  }

  void setZero ( );

  FullMatrix<_DataType, Parallelize>& operator*= ( const _DataType alpha );

  FullMatrix<_DataType, Parallelize>& operator*= ( const FullMatrix<_DataType, Parallelize> &mat );

  FullMatrix<_DataType, Parallelize>& operator= ( const Matrix<_DataType> &Mat );

  FullMatrix<_DataType, Parallelize>& operator= ( const FullMatrix<_DataType, Parallelize> &Mat );

  using Matrix<_DataType>::operator+=;
  FullMatrix<_DataType, Parallelize>& operator+= ( const FullMatrix<_DataType, Parallelize> &Mat );

  using Matrix<_DataType>::operator-=;
  FullMatrix<_DataType, Parallelize>& operator-= ( const FullMatrix<_DataType, Parallelize> &Mat );

  void mult ( const Vector<_DataType> &src, Vector<_DataType> &dst ) const {
    // Compiler does not remove if inside the omp pragma which might drastically affect performance for small matrices
    if ( Parallelize ) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 0; i < static_cast<int> ( this->_numRows ); ++i ) {
        dst[ i ] = 0;
        for ( int j = 0; j < static_cast<int> ( this->_numCols ); ++j ) {
          dst[ i ] += ref ( i, j ) * src[ j ];
        }
      }
    } else {
      for ( int i = 0; i < static_cast<int> ( this->_numRows ); ++i ) {
        dst[ i ] = 0;
        for ( int j = 0; j < static_cast<int> ( this->_numCols ); ++j ) {
          dst[ i ] += ref ( i, j ) * src[ j ];
        }
      }
    }
  }

  void applyAdd ( const Vector<_DataType> &src, Vector<_DataType> &dst ) const {
    // Compiler does not remove if inside the omp pragma which might drastically affect performance for small matrices
    if ( Parallelize ) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int i = 0; i < static_cast<int> ( this->_numRows ); ++i ) {
        for ( int j = 0; j < static_cast<int> ( this->_numCols ); ++j ) {
          dst[ i ] += ref ( i, j ) * src[ j ];
        }
      }
    } else {
      for ( int i = 0; i < static_cast<int> ( this->_numRows ); ++i ) {
        for ( int j = 0; j < static_cast<int> ( this->_numCols ); ++j ) {
          dst[ i ] += ref ( i, j ) * src[ j ];
        }
      }
    }
  }

  void applyAddTranspose ( const Vector<_DataType> &src, Vector<_DataType> &dst ) const;

  void applyAddTransposed ( const Vector<_DataType> &src, Vector<_DataType> &dst ) const {
    applyAddTranspose( src, dst );
  }

  void makeProduct ( const Matrix<_DataType> &A, const Matrix<_DataType> &B ) {
    if ( A.getNumRows() != this->_numRows || A.getNumCols() != B.getNumRows() || B.getNumCols() != this->_numCols ) {
      throw Exception ( "FullMatrix::makeProduct: dimensions not compatible.", __FILE__, __LINE__ );
    }
    for ( int i = 0; i < this->_numRows; ++i ) {
      for ( int j = 0; j < this->_numCols; ++j ) {
        _DataType v = 0;
        for ( int k = 0; k < A.getNumCols(); ++k ) {
          v += A.get ( i, k ) * B.get ( k, j );
        }
        set ( i, j, v );
      }
    }
  }

  void makeRowEntries ( vector<typename Row<_DataType>::RowEntry > &vec, const int RowNum ) const {
    if ( !makeRowEntriesWarningPrinted )  {
      cerr << aol::color::red << "Warning: Using makeRowEntries on FullMatrix is slow! Consider using ref() instead. This warning won't be printed again\n" << aol::color::reset;
      makeRowEntriesWarningPrinted = true;
    }
    int size = this->getNumCols ();
    vec.resize ( size );
    for ( int i = 0; i < size; ++i ) {
      vec [i].col = i;
      vec [i].value = ref ( RowNum, i );
    }
  }

  //! Makes this matrix the inverse of Mat if Mat is invertible.
  void makeInverse ( const Matrix<_DataType> &Mat );

  // makes P such that Mat*P^T=[ A_1 A_2 ] where A_1 is invertible
  // where Mat is this matrix
  void makeLeftInverse ( FullMatrix<_DataType, Parallelize>& A_1Inv,
                         FullMatrix<_DataType, Parallelize>& A_2,
                         PermutationMatrix<_DataType>& P );

  void makeSubstitutionMatrix ( const Vector<_DataType> &b,
                                Vector<_DataType> &q,
                                FullMatrix<_DataType, Parallelize> &/*F*/,
                                PermutationMatrix<_DataType>& P );

  //! Method for LU-decomposition.
  void makeLU ( const Matrix<_DataType>& Mat, PermutationMatrix<_DataType>& P );

  //! Solves the LU-system. P is the Permutation generated by makeLU.
  void LUSolve ( Vector<_DataType>& X, const Vector<_DataType> &RHS, const PermutationMatrix<_DataType>& P ) const;

  // Row[i] -> Row[i+1]  i=0..m-2; Row[m-1]->Row[0]
  void shiftRowsUp();

  void apply ( const Vector<_DataType> &Arg, Vector<_DataType> &Dest ) const {
    mult ( Arg, Dest );
  }

  //! transposes this matrix
  void transpose() {
    if ( this->_numRows != this->_numCols ) {
      throw Exception ( "Transposition of non-quadratic matrices not implemented.", __FILE__, __LINE__ );
    }
    for ( int i = 0; i < this->_numRows; ++i ) {
      for ( int j = i + 1; j < this->_numCols; ++j ) {
        _DataType v = ref ( i, j );
        ref ( i, j ) = ref ( j, i );
        ref ( j, i ) = v;
      }
    }
  }

  //! writes transposed matrix to another one
  void transposeTo ( FullMatrix< _DataType, Parallelize > &Dest ) const  {
    if ( ( this->_numRows != Dest.getNumCols() ) || ( this->_numCols != Dest.getNumRows() ) ) {
      throw Exception ( "aol::FullMatrix<_DataType>::transposeTo: matrices incompatible", __FILE__, __LINE__ );
    } else {
      for ( int i = 0; i < this->_numRows; ++i ) {
        for ( int j = 0; j < this->_numCols; ++j ) {
          Dest.ref ( j, i ) = ref ( i, j );
        }
      }
    }
  }

  //! change makeRowEntriesWarningPrinted
  static void setMakeRowEntriesWarningPrinted ( const bool Arg ) {
    makeRowEntriesWarningPrinted = Arg;
  }

  //! Reads a block from the full matrix
  void getBlock ( int r, int c, FullMatrix<_DataType, Parallelize>& block ) const;

  //! Write a block into the full matrix
  void setBlock ( int r, int c, const Matrix<_DataType>& block );

  //! Reads a partial column from the full matrix
  void getSubColumn ( int r, int c, Vector<_DataType>& subcol ) const;

  //! Write a partial column into the full matrix
  void setSubColumn ( int r, int c, const Vector<_DataType>& subcol );

  //! Write a partial column into the full matrix
  void addSubColumn ( int r, int c, const Vector<_DataType>& subcol );

  //! Reads a partial row from the full matrix
  void getSubRow ( int r, int c, Vector<_DataType>& subrow ) const;

  //! Write a partial row into the full matrix
  void setSubRow ( int r, int c, const Vector<_DataType>& subrow );

  _DataType getFrobeniusNormSqr() const;

  bool checkForNANsAndINFs() const {
    return data.checkForNANsAndINFs();
  }

  //! SelfTest
  //! Currently tests only LU decomposition
  static bool isSelfTestOk ();

  //! Returns a reference to the aol::Vector that stores the entries of the matrix.
  const Vector<_DataType>& getDataVectorReference() const {
    return data;
  }

protected:

  void swapRows ( int R1, int R2 );

  void swapColumns ( int C1, int C2 );

  void multRow ( int R, _DataType alpha, int StartIndex );

  //! What does this method do?
  void addToRow ( int R1, int R2, _DataType alpha, int StartIndex );

  static bool makeRowEntriesWarningPrinted;
};

template < class _DataType, bool Parallelize >
bool FullMatrix< _DataType, Parallelize >::makeRowEntriesWarningPrinted = false;

//! A class for symmetric matrices.
/**
 * \ingroup Matrix
 */
template <typename DataType>
class SymmetricMatrix : public Matrix<DataType> {

  SymmetricMatrix() {}
  void init ( int Dimension );
  void destroy();
public:

  explicit SymmetricMatrix ( int Dimension ) {
    init ( Dimension );
  }

  explicit SymmetricMatrix ( const SymmetricMatrix<DataType> &other ) : Matrix<DataType> ( other ) {
    throw UnimplementedCodeException ( "SymmetricMatrix copy constructor not implemented", __FILE__, __LINE__ );
  }

  virtual ~SymmetricMatrix() {
    destroy();
  }

  //! resize Matrix, deleting old contents
  void reallocate ( const int Dimension ) {
    destroy();
    init ( Dimension );
  }

  void reallocate ( const int Dim1, const int Dim2 ) {
    if ( Dim1 == Dim2 ) {
      reallocate ( Dim1 );
    } else {
      throw Exception ( "SymmetricMatrix<>::reallocate called for non-square format", __FILE__, __LINE__ );
    }
  }


  DataType get ( int I, int J ) const;

  void set ( int I, int J, DataType Value );

  Matrix<DataType>& operator*= ( const DataType alpha );


  void apply ( const Vector<DataType> &src, Vector<DataType> &dst ) const;

  void setZero ( );

protected:
  DataType *data;
};


//! A diagonal matrix class.
/**
 * \ingroup Matrix
 */
template <typename DataType>
class DiagonalMatrix : public Matrix<DataType> {
protected:
  Vector<DataType> diag;

public:
  //! standard constructor creating size-0 matrix?
  DiagonalMatrix() {}

  //! constructor
  explicit DiagonalMatrix ( const int Dimension ) :
      Matrix<DataType> ( Dimension, Dimension ),
      diag ( Dimension ) {
  }

  template< typename Struc >
  explicit DiagonalMatrix ( const Struc &stru ) :
      Matrix<DataType> ( stru.getNumberOfNodes(), stru.getNumberOfNodes() ),
      diag ( stru.getNumberOfNodes() ) {
  }

public:
  //! constructor for compatibility
  DiagonalMatrix ( int Dim1, int Dim2 ) : Matrix<DataType> ( Dim1, Dim1 ), diag ( Dim1 ) {
    if ( Dim1 != Dim2 ) {
      throw Exception ( "Trying to set up Diagonal matrix with different number of rows and columns. That won't work!", __FILE__, __LINE__ );
    }
  }

  DiagonalMatrix ( const DiagonalMatrix &mat, CopyFlag copyFlag = DEEP_COPY )
  : Matrix<DataType> ( mat.getNumRows (), mat.getNumCols () ), diag ( mat.diag, copyFlag ) {
    if ( mat.getNumRows () != mat.getNumCols () )
      throw Exception ( "Trying to set up Diagonal matrix with different number of rows and columns. That won't work!", __FILE__, __LINE__ );

    switch ( copyFlag ) {
    // For DEEP_COPY, STRUCT_COPY and FLAT_COPY everything that needs to be done is already done in the constructor of Matrix and diag.
    case DEEP_COPY:
    case STRUCT_COPY:
    case FLAT_COPY:
      break;

    default:
      throw UnimplementedCodeException ( "This CopyFlag is not implemented yet.", __FILE__, __LINE__ );
      break;
    }
  }

  Matrix<DataType>* clone ( CopyFlag copyFlag = DEEP_COPY ) const {
    DiagonalMatrix *mat = new DiagonalMatrix ( *this, copyFlag );
    return mat;
  }

public:
  virtual ~DiagonalMatrix() { }

  //! resize Matrix, preserving old contents
  void resize ( const int Dimension ) {
    diag.resize ( Dimension );
    this->_numRows = this->_numCols = Dimension;
  }

  void resize ( const int Dim1, const int Dim2 ) {
    if ( Dim1 == Dim2 ) {
      resize ( Dim1 );
    } else {
      throw Exception ( "DiagonalMatrix<>::resize called for non-square format", __FILE__, __LINE__ );
    }
  }

  //! resize Matrix, deleting old contents
  void reallocate ( const int Dimension ) {
    diag.reallocate ( Dimension );
    this->_numRows = this->_numCols = Dimension;
  }

  void reallocate ( const int Dim1, const int Dim2 ) {
    if ( Dim1 == Dim2 ) {
      reallocate ( Dim1 );
    } else {
      throw Exception ( "DiagonalMatrix<>::reallocate called for non-square format", __FILE__, __LINE__ );
    }
  }

  const Vector<DataType> &getDiagVectorReference ( ) const {
    return diag;
  }

  //! Returns a matrix entry.
  virtual DataType get ( int I, int J ) const {
#ifdef BOUNDS_CHECK
    this->boundsCheck ( I, J, "DiagonalMatrix::get: Index out of bounds", __FILE__, __LINE__ );
#endif
    return ( I == J ? diag.get ( I ) : ZOTrait<DataType>::zero );
  }

  const DataType operator[] ( int I ) const {
    return diag[I];
  }

  //! Sets a diagonal matrix entry. Others are fixed to 0.
  void set ( int I, int J, DataType Value ) {
#ifdef BOUNDS_CHECK
    this->boundsCheck ( I, J, "DiagonalMatrix::set Index out of bounds", __FILE__, __LINE__ );
#endif
    if ( I == J ) diag[ I ] = Value;
  }

  //! Sets a diagonal matrix entry. Others are fixed to 0.
  void add ( int I, int J, DataType Value ) {
#ifdef BOUNDS_CHECK
    this->boundsCheck ( I, J, "DiagonalMatrix::add: Index out of bounds", __FILE__, __LINE__ );
#endif
    if ( I == J ) diag[ I ] += Value;
  }

  virtual int numNonZeroes ( int I ) const {
    int nNonZeroes = 0;
    if ( diag[ I ] != 0 )
      ++nNonZeroes;
    return nNonZeroes;
  }

  //! Multiply matrix with scalar
  Matrix<DataType>& operator*= ( const DataType alpha ) {
    diag *= alpha;
    return *this;
  }

  Matrix<DataType>& operator*= ( const DiagonalMatrix<DataType> &Mat ) {
    if ( diag.size() != Mat.diag.size() ) {
      throw Exception ( "Matrix dimensions don't match.", __FILE__, __LINE__ );
    }
    for ( int i = 0; i < diag.size(); ++i ) {
      diag[i] *= Mat.diag[i];
    }
    return *this;
  }

  using Matrix<DataType>::operator+=;
  Matrix<DataType>& operator+= ( const DiagonalMatrix<DataType> &Mat ) {
    diag += Mat.diag;
    return *this;
  }

  Matrix<DataType>& operator= ( const DiagonalMatrix<DataType> &Mat ) {
    if ( Mat._numRows != this->_numRows ) {
      throw Exception ( "DiagonalMatrix::operator= :  dimensions don't match.", __FILE__, __LINE__ );
    }
    diag = Mat.diag;
    return *this;
  }

  //! \todo Rename to addMultipleTo
  Matrix<DataType>& addMultiple ( const DiagonalMatrix<DataType> &Mat, DataType Scalar ) {
    diag.addMultiple ( Mat.diag, Scalar );
    return *this;
  }

  void applyAdd ( const Vector<DataType> &src, Vector<DataType> &dst ) const {
    for ( int i = 0; i < this->_numRows; ++i ) {
      dst[ i ] += src.get ( i ) * diag[ i ];
    }
  }

  //! Matrix-vector multiplication.
  void apply ( const Vector<DataType> &src, Vector<DataType> &dst ) const {
    for ( int i = 0; i < this->_numRows; ++i ) {
      dst[ i ] = src.get ( i ) * diag[ i ];
    }
  }

  void applySingle ( Vector<DataType> &Arg ) const {
    for ( int i = 0; i < this->_numRows; ++i ) {
      Arg[ i ] = Arg[ i ] * diag[ i ];
    }
  }

  //! Inverts the matrix.
  void invert ( bool saveZeros = false ) {
    for ( int i = 0; i < this->_numRows; ++i ) {
      if ( saveZeros && diag[i] == 0 )
        diag[ i ] = 0.;
      else
        diag[ i ] = ZOTrait<DataType>::one / diag[ i ];
    }
  }

  //! Sets all entries to zero
  void setZero ( ) {
    diag.setZero();
  }

  //! Sets all entries to the specified Data
  void setAll ( const DataType Data ) {
    diag.setAll ( Data );
  }

  //! Makes this matrix the inverse of the diagonal part of other_mat. Useful e. g. for preconditioning.
  template <typename MatrixType>
  void setToInverseDiagonalOf ( const MatrixType &other_mat ) {
    // missing: size checking.
    for ( int i = 0; i < this->_numRows; ++i ) {
      const DataType omatii ( other_mat.getDiag ( i ) );
      if ( omatii != ZOTrait<DataType>::zero ) {
        diag[i] = ZOTrait<DataType>::one / omatii;
      } else {
        diag[i] = ZOTrait<DataType>::one;
#ifdef VERBOSE
        cerr << "DiagonalMatrix::setToInverseDiagonalOf: got zero entry, setting 1.0" << endl;
#endif
      }
    }
  }

  template <typename MatrixType>
  void setToAbsInverseDiagonalOf ( const MatrixType &other_mat ) {
    // missing: size checking.
    for ( int i = 0; i < this->_numRows; ++i ) {
      const DataType omatii ( other_mat.getDiag ( i ) );
      if ( omatii != ZOTrait<DataType>::zero ) {
        diag[i] = aol::Abs( ZOTrait<DataType>::one / omatii );
      } else {
        diag[i] = ZOTrait<DataType>::one;
#ifdef VERBOSE
        cerr << "DiagonalMatrix::setToInverseDiagonalOf: got zero entry, setting 1.0" << endl;
#endif
      }
    }
  }

  //! Lumps the rows of otherMat to this diagonal matrix
  template< typename MatrixType >
  void setToLumpedCopyOf ( const MatrixType &otherMat ) {
    if ( ( this->getNumRows() != otherMat.getNumRows() ) || ( this->getNumCols() != otherMat.getNumCols() ) )
      throw aol::Exception ( "aol::DiagonalMatrix::setToLumpedCopyOf: matrix sizes do not match", __FILE__, __LINE__ );

    vector<typename Row<DataType>::RowEntry > vec;
    for ( int i = 0; i < otherMat.getNumRows(); ++i ) {
      otherMat.makeRowEntries ( vec, i );
      for ( typename vector<typename Row<DataType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
        this->add ( i, i, it->value );
      }
    }
  }

  virtual void makeRowEntries ( vector<typename Row<DataType>::RowEntry > &vec, const int RowNum ) const {
    vec.resize ( 1 );
    vec[0].col = RowNum;
    vec[0].value = get ( RowNum, RowNum );
  }

};

/**
 * \brief A class for block matrices that does its own memory management.
 *        Every block entry may contain a matrix or just NULL.
 *
 * \author Berkels
 * \ingroup Matrix
 */
template <typename _MatrixType>
class SparseBlockMatrix : public BlockOpBase<typename _MatrixType::DataType, _MatrixType> {
protected:
  typedef typename _MatrixType::DataType RealType;

private:
  typedef int (_MatrixType::*MatrixMemberFunctionPointer)(int x) const;         //typedef for a const member function pointer. Used in maxNumStoredEntriesPerRow and maxNumNonZeroesPerRow.

  //! Vector for delete flags of blocks
  vector<bool> _deleteFlag;

protected:
  void setDeleteFlag ( int i, int j, bool flag ) {
    _deleteFlag[qc::ILexCombine2 ( j, i, this->_n )] = flag;
  }

  bool getDeleteFlag ( int i, int j ) {
    return _deleteFlag[qc::ILexCombine2 ( j, i, this->_n )];
  }

public:
  typedef _MatrixType MatrixType;
  typedef typename MatrixType::DataType DataType;

  SparseBlockMatrix () : BlockOpBase<DataType, MatrixType> ( 0, 0 ), _deleteFlag ( 0 ) {
  }

  SparseBlockMatrix ( const int M, const int N ) : BlockOpBase<DataType, MatrixType> ( M, N ), _deleteFlag ( M * N ) {
    typename vector< bool >::iterator it;
    for ( it = _deleteFlag.begin (); it != _deleteFlag.end (); ++it )
      *it = false;
  }

  //! \note Does not initialize the blocks, i.e. only depends on NumRows.size() and NumCols.size().
  SparseBlockMatrix ( const aol::Vector<int> &NumRows, const aol::Vector<int> &NumCols ) : BlockOpBase<DataType, MatrixType> ( NumRows.size(), NumCols.size() ), _deleteFlag ( this->_m * this->_n ) {
    typename vector< bool >::iterator it;
    for ( it = _deleteFlag.begin (); it != _deleteFlag.end (); ++it )
      *it = false;
  }

  //! \warning Using STRUCT_COPY copies the "block structure", but does not copy the structure of the blocks. Specifically the number of blocks per row and column is copied,
  //!          but not the structure of the blocks themselves.
  explicit SparseBlockMatrix ( const SparseBlockMatrix<MatrixType> &Matrix, CopyFlag copyFlag = DEEP_COPY, bool structCopySubMatrices = false )
  : BlockOpBase<DataType, MatrixType> ( Matrix.getNumRows(), Matrix.getNumCols() ), _deleteFlag ( Matrix.getNumRows() * Matrix.getNumCols() ) {
    switch ( copyFlag ) {
      case STRUCT_COPY: {
        // Treat STRUCT_COPY seperately because of the flag structCopySubMatrices.
        if ( structCopySubMatrices ) {
          for ( int i = 0; i < this->getNumRows (); ++i ) {
            for ( int j = 0; j < this->getNumCols (); ++j ) {
              if ( Matrix.getPointer ( i, j ) ) {
                allocateMatrixAndCopyOtherMatrix ( i, j, Matrix.getReference ( i, j ), aol::STRUCT_COPY );
                this->setDeleteFlag ( i, j, true );
              }
              else
                this->setDeleteFlag ( i, j, false );
            }
          }
        }
      }
      break;
      case FLAT_COPY:
      case DEEP_COPY: {
        for ( int i = 0; i < this->getNumRows(); ++i ) {
          for ( int j = 0; j < this->getNumCols(); ++j ) {
            if ( Matrix.getPointer ( i, j ) ) {
              allocateMatrixAndCopyOtherMatrix ( i, j, * ( Matrix.getPointer ( i, j ) ), copyFlag );
              this->setDeleteFlag ( i, j, true );
            }
            else
              this->setDeleteFlag ( i, j, false );
          }
        }
      }
      break;
      default: {
        throw UnimplementedCodeException ( "This CopyFlag is not implemented yet.", __FILE__, __LINE__ );
      }
      break;
    }
  }

  //! \warning clone using STRUCT_COPY behaves different from this classes' copy constructor with STRUCT_COPY!
  //! \author toelkes
  SparseBlockMatrix<MatrixType>* clone ( CopyFlag copyFlag = DEEP_COPY ) const {
    SparseBlockMatrix<MatrixType> *mat = new SparseBlockMatrix<MatrixType> ( this->getNumRows (), this->getNumCols () );
    for ( int i = 0; i < this->getNumRows (); ++i )
      for ( int j = 0; j < this->getNumCols (); ++j ) {
        if ( this->getPointer ( i, j ) ) {
          mat->setPointer ( i, j, static_cast < MatrixType* > ( this->getPointer ( i, j )->clone ( copyFlag ) ) );
          mat->setDeleteFlag ( i, j, true );
        }
        else
          mat->setDeleteFlag ( i, j, false );
      }

    return mat;
  }

  virtual ~SparseBlockMatrix() {
    deleteBlockEntries();
  }

  void clear( ) {
    throw UnimplementedCodeException ( "clear can't be implemented on SparseBlockMatrix as long as the number of blocks is fixed.", __FILE__, __LINE__ );
  }

  //! Delete all contents, freeing memory, but keeping the number of blocks unchanged.
  void deleteBlockEntries( ) {
    typename vector< MatrixType* >::iterator it;
    typename vector< bool >::iterator deleteFlagIterator;
    for ( it = this->_blockEntryPointers.begin(), deleteFlagIterator = this->_deleteFlag.begin (); it != this->_blockEntryPointers.end(); ++it, ++deleteFlagIterator ) {
      if ( *it && *deleteFlagIterator ) {
        delete ( *it );
        ( *it ) = NULL;
      }
    }
  }

  //! Call setZero() on all non-NULL blocks.
  void setZero( ) {
    typename vector< MatrixType* >::iterator it;
    for ( it = this->_blockEntryPointers.begin(); it != this->_blockEntryPointers.end(); ++it ) {
      if ( *it ) {
        ( *it )->setZero();
      }
    }
  }

  //! Set matrix (i.e. diagonal blocks) to identity.
  void setIdentity () {
    if ( this->getNumRows () != this->getNumCols () ) {
      throw Exception( "Can not set non quadratic matrix to identity!", __FILE__, __LINE__);
    }

    for ( int i = 0; i < this->getNumRows (); ++i )
      for ( int j = 0; j < this->getNumCols (); ++j ) {
        if ( i == j ) {
          if ( !this->getPointer ( i, j ) ) {
            throw Exception( "Diagonal blocks have to be set before the matrix is set to identity!", __FILE__, __LINE__);
          }
          this->getReference ( i, j ).setIdentity ();
        }
        else {
          if ( this->getPointer ( i, j ) ) {
            this->getReference ( i, j ).setZero ();
          }
        }
      }
  }

  void setMatrixPointer ( const int I, const int J, _MatrixType *PMat, const bool DeleteFlag ) {
    if ( this->getPointer ( I, J ) && this->getDeleteFlag ( I, J )  ) {
      delete ( this->getPointer ( I, J ) );
    }
    this->setPointer ( I, J, PMat );
    setDeleteFlag ( I, J, DeleteFlag );
  }

  template <typename InitType>
  MatrixType& allocateMatrix ( const int I, const int J, const InitType &Grid ) {
    if ( this->getPointer ( I, J ) && this->getDeleteFlag ( I, J )  ) {
      delete ( this->getPointer ( I, J ) );
    }
    MatrixType* tempPointer = new MatrixType ( Grid );
    this->setPointer ( I, J, tempPointer );
    this->setDeleteFlag ( I, J, true );
    return ( *tempPointer );
  }

  template <typename InitType, typename NewMatrixType>
  NewMatrixType& allocateMatrix ( const int I, const int J, const InitType &Grid, const NewMatrixType & ) {
    if ( this->getPointer ( I, J ) && this->getDeleteFlag ( I, J ) ) {
      delete ( this->getPointer ( I, J ) );
    }
    NewMatrixType* tempPointer = new NewMatrixType ( Grid );
    this->setPointer ( I, J, static_cast<MatrixType*>(tempPointer) );
    this->setDeleteFlag ( I, J, true );
    return ( *tempPointer );
  }

  template <typename NewMatrixType>
  NewMatrixType& allocateMatrix ( const int I, const int J, const int numRows, const int numCols, const NewMatrixType & ) {
    return ( this->template allocateMatrixOfType<NewMatrixType> ( I, J, numRows, numCols ) );
  }

  template <typename NewMatrixType>
  NewMatrixType& allocateMatrixOfType ( const int I, const int J, const int numRows, const int numCols ) {
    if ( this->getPointer ( I, J ) && this->getDeleteFlag ( I, J ) ) {
      delete ( this->getPointer ( I, J ) );
    }
    NewMatrixType* tempPointer = new NewMatrixType ( numRows, numCols );
    this->setPointer ( I, J, static_cast<MatrixType*>(tempPointer) );
    this->setDeleteFlag ( I, J, true );
    return ( *tempPointer );
  }

  MatrixType& allocateMatrix ( const int I, const int J, const int NumRows, const int NumCols ) {
    if ( this->getPointer ( I, J ) && this->getDeleteFlag ( I, J ) ) {
      delete ( this->getPointer ( I, J ) );
    }
    MatrixType* tempPointer = new MatrixType ( NumRows, NumCols );
    this->setPointer ( I, J, tempPointer );
    this->setDeleteFlag ( I, J, true );
    return ( *tempPointer );
  }

  template <typename NewMatrixType>
  NewMatrixType& allocateMatrixAndCopyOtherMatrix ( const int I, const int J, const NewMatrixType &OtherMatrix, CopyFlag copyFlag = DEEP_COPY ) {
    if ( this->getPointer ( I, J ) && this->getDeleteFlag ( I, J ) ) {
      delete ( this->getPointer ( I, J ) );
    }
    NewMatrixType* tempPointer = new NewMatrixType ( OtherMatrix, copyFlag );
    this->setPointer ( I, J, static_cast<MatrixType*>(tempPointer) );
    this->setDeleteFlag ( I, J, true );
    return ( *tempPointer );
  }

  SparseBlockMatrix<MatrixType>& operator= ( const SparseBlockMatrix<MatrixType> &Mat ) {
#ifdef BOUNDS_CHECK
    if ( (this->getNumCols() != Mat.getNumCols()) || (this->getNumRows() != Mat.getNumRows()) )
      throw Exception( "Size mismatch.", __FILE__, __LINE__);
#endif
    for(int i=0; i<this->getNumRows(); ++i) {
      for(int j=0; j<this->getNumCols(); ++j) {
        if( Mat.getPointer( i, j ) ) {
          allocateMatrixAndCopyOtherMatrix( i, j, *(Mat.getPointer( i, j )) );
        }
        else {
          if ( this->getPointer ( i, j ) && this->getDeleteFlag ( i, j ) ) {
            delete ( this->getPointer ( i, j ) );
            this->setPointer ( i, j, NULL );
	    this->setDeleteFlag ( i, j, false );
          }
        }
      }
    }
    return *this;
  }

  SparseBlockMatrix<MatrixType>& operator*= ( const DataType alpha ) {
    for ( int i = 0; i < this->getNumRows (); ++i ) {
      for ( int j = 0; j < this->getNumCols (); ++j ) {
        if ( this->getPointer ( i, j ) )
          this->getReference ( i, j ) *= alpha;
      }
    }

    return *this;
  }

  SparseBlockMatrix<MatrixType>& operator+= ( const SparseBlockMatrix<MatrixType> &Mat ) {
    this->addMultiple ( Mat, aol::ZOTrait<RealType>::one );
    return *this;
  }

  //! \brief Adds vec1 \f$ \otimes \f$ vec2 to this.
  //! \author Toelkes
  inline void addTensorProduct ( const MultiVector< RealType > &vec1, const MultiVector< RealType > &vec2 ) {
    addTensorProductMultiple ( vec1, vec2, 1.0 );
  }

  //! \brief Adds factor * vec1 \f$ \otimes \f$ vec2 to this.
  //! \author Toelkes
  void addTensorProductMultiple ( const MultiVector< RealType > &vec1, const MultiVector< RealType > &vec2, const RealType factor ) {
    for ( int i = 0; i < this->getNumRows (); ++i ) {
      for ( int j = 0; j < this->getNumCols (); ++j ) {
        if ( this->getPointer ( i, j ) )
          this->getReference ( i, j ).addTensorProductMultiple ( vec1[i], vec2[j], factor );
        else
          throw Exception ( "Can't store tensor product.", __FILE__, __LINE__ );
      }
    }
  }

  void addMultiple ( const SparseBlockMatrix &matrix, const RealType factor ) {
    for ( int i = 0; i < this->getNumRows (); ++i ) {
      for ( int j = 0; j < this->getNumCols (); ++j ) {
        if ( matrix.getPointer ( i, j ) ) {
          if ( this->getPointer ( i, j ) )
            this->getReference ( i, j ).addMultiple ( matrix.getReference ( i, j ), factor );
          else {
            allocateMatrixAndCopyOtherMatrix ( i, j, *( matrix.getPointer ( i, j ) ) );
            this->getReference ( i, j ) *= factor;
          }
        }
      }
    }
  }

  void applyAddTransposed ( const MultiVector<RealType> &MArg, MultiVector<RealType> &MDest ) const {
#ifdef BOUNDS_CHECK
    if ( MArg.numComponents() != this->getNumRows () || MDest.numComponents() != this->getNumCols() ) {
      throw aol::Exception("MultiVectors do not have correct numbers of components", __FILE__, __LINE__ );
    }
#endif

    for ( int i = 0; i < this->getNumRows (); ++i ) {
      for ( int j = 0; j < this->getNumCols (); ++j ) {
        if ( this->getPointer ( i, j ) ) {
          this->getReference ( i, j ).applyAddTransposed ( MArg[i], MDest[j] );
        }
      }
    }
  }

  void transpose () {
    SparseBlockMatrix<MatrixType> tempMat ( *this, DEEP_COPY );
    tempMat.transposeTo ( *this );
  }

  template <typename OutputMatrixType>
  void transposeTo ( SparseBlockMatrix<OutputMatrixType> &OtherMatrix ) const {
#ifdef BOUNDS_CHECK
    if ( ( this->getNumCols() != OtherMatrix.getNumRows() ) || ( this->getNumRows() != OtherMatrix.getNumCols() ) )
      throw Exception ( "Size mismatch.", __FILE__, __LINE__ );
#endif
    OtherMatrix.deleteBlockEntries();
    for ( int i = 0; i < this->getNumRows(); ++i ) {
      for ( int j = 0; j < this->getNumCols(); ++j ) {
        if ( this->getPointer ( i, j ) ) {
          OutputMatrixType &newBlock = OtherMatrix.allocateMatrix ( j, i, this->getReference ( i, j ).getNumCols(), this->getReference ( i, j ).getNumRows() );
          this->getPointer ( i, j )->transposeTo ( newBlock );
        }
      }
    }
  }

  void scaleRow ( const int Row, const RealType ScaleFactor ) {
    for ( int j = 0; j < this->getNumCols(); ++j ) {
      if ( this->getPointer ( Row, j ) ) {
        this->getReference ( Row, j ) *= ScaleFactor;
      }
    }
  }

  void scaleColumn ( const int Column, const RealType ScaleFactor ) {
    for ( int i = 0; i < this->getNumRows(); ++i ) {
      if ( this->getPointer ( i, Column ) ) {
        this->getReference ( i, Column ) *= ScaleFactor;
      }
    }
  }

  //! multiplies two Matrices \f$( this += M_1 \cdot M_2)\f$ and adds the result to the current matrix
  template <typename InputMatrixType1, typename InputMatrixType2>
  void addMatrixProduct ( const SparseBlockMatrix<InputMatrixType1> &M1, const SparseBlockMatrix<InputMatrixType2> &M2 ){
#ifdef BOUNDS_CHECK
    if ( ( this->getNumRows() != M1.getNumRows() ) || ( this->getNumCols() != M2.getNumCols() )
         || ( M1.getNumCols() != M2.getNumRows() ) )
      throw Exception ( "Size mismatch.", __FILE__, __LINE__ );
#endif
    for ( int i = 0; i < this->getNumRows(); ++i ) {
      for ( int j = 0; j < this->getNumCols(); ++j ) {
        for ( int k = 0; k < M2.getNumRows(); ++k )
          // We can only add something if the corresponding blocks in M1 and M2 are both not NULL
          if ( M1.getPointer ( i, k ) && M2.getPointer ( k, j ) ) {
            if ( this->getPointer ( i, j ) == NULL )
              allocateMatrix ( i, j, M1.getPointer ( i, k )->getNumRows(), M2.getPointer ( k, j )->getNumCols() );
            this->getPointer ( i, j )->addMatrixProduct ( *M1.getPointer ( i, k ), *M2.getPointer ( k, j ) );
          }
      }
    }
  }

private:
  const MatrixType* getFirstStoredPointer ( ) const {
    const MatrixType* tempPointer = NULL;
    for ( int blockRow = 0; blockRow < this->getNumRows(); blockRow++ ) {
      for ( int blockCol = 0; blockCol < this->getNumCols(); blockCol++ ) {
        if ( this->getPointer ( blockRow, blockCol ) ) {
          tempPointer = this->getPointer ( blockRow, blockCol );
          break;
        }
      }
      if ( tempPointer )
        break;
    }
    return tempPointer;
  }

public:
  void getBlockRowNums ( Vector<int> &RowNums ) const {
    const MatrixType* tempPointer = getFirstStoredPointer();

    // If no block is set, we cannot determine the size.
    if ( tempPointer == NULL )
      throw aol::InconsistentDataException ( "Cannot determine row size.", __FILE__, __LINE__ );

    const int defaultNumRows = tempPointer->getNumRows();

    RowNums.reallocate ( this->getNumRows() );

    for ( int blockRow = 0; blockRow < this->getNumRows(); blockRow++ ) {
      aol::Vector<int> sizes;
      sizes.reserve ( this->getNumRows() );
      for ( int blockCol = 0; blockCol < this->getNumCols(); blockCol++ ) {
        if ( this->getPointer ( blockRow, blockCol ) )
          sizes.pushBack( this->getReference ( blockRow, blockCol ).getNumRows() );
      }
      if ( sizes.size() == 0 )
        RowNums[blockRow] = defaultNumRows;
      else {
        if ( sizes.numOccurence ( sizes[0] ) != sizes.size() )
          throw aol::DimensionMismatchException ( "Block size mismatch.", __FILE__, __LINE__ );

        RowNums[blockRow] = sizes[0];
      }
    }
  }

  void getBlockColNums ( Vector<int> &ColNums ) const {
    const MatrixType* tempPointer = getFirstStoredPointer();

    // If no block is set, we cannot determine the size.
    if ( tempPointer == NULL )
      throw aol::InconsistentDataException ( "Cannot determine row size.", __FILE__, __LINE__ );

    const int defaultNumCols = tempPointer->getNumCols();

    ColNums.reallocate ( this->getNumCols() );
    for ( int blockCol = 0; blockCol < this->getNumCols(); blockCol++ ) {
      aol::Vector<int> sizes;
      sizes.reserve ( this->getNumCols() );
      for ( int blockRow = 0; blockRow < this->getNumRows(); blockRow++ ) {
        if ( this->getPointer ( blockRow, blockCol ) )
          sizes.pushBack( this->getReference ( blockRow, blockCol ).getNumCols() );
      }
      if ( sizes.size() == 0 )
        ColNums[blockCol] = defaultNumCols;
      else {
        if ( sizes.numOccurence ( sizes[0] ) != sizes.size() )
          throw aol::DimensionMismatchException ( "Block size mismatch.", __FILE__, __LINE__ );

        ColNums[blockCol] = sizes[0];
      }
    }
  }

  int getTotalNumRows ( ) const {
    Vector<int> rowSizes;
    this->getBlockRowNums ( rowSizes );
    return rowSizes.sum();
  }

  int getTotalNumCols ( ) const {
    Vector<int> colSizes;
    this->getBlockColNums ( colSizes );
    return colSizes.sum();
  }

  /**
   * add this blockMatrix to BigMatrix as if it were one big matrix without block structure.
   *
   * Assumes that all matrices in one row have same number of columns and vice versa.
   */
  template <typename BigMatrixType>
  void addUnblockedMatrixTo ( BigMatrixType &BigMatrix, const bool autoResize = true ) const {
    // If no block is set, there is nothing to add.
    if ( getFirstStoredPointer() == NULL )
      return;

    Vector<int> rowSizes, colSizes;
    this->getBlockRowNums ( rowSizes );
    this->getBlockColNums ( colSizes );

    if ( autoResize ) {
      BigMatrix.resize ( rowSizes.sum(), colSizes.sum() );
    }

    vector<typename Row<DataType>::RowEntry > vec;
    int bigIOffset = 0;
    for ( int blockRow = 0; blockRow < this->getNumRows(); blockRow++ ) {
      int bigJOffset = 0;
      for ( int blockCol = 0; blockCol < this->getNumCols(); blockCol++ ) {
        const MatrixType* currentBlock = this->getPointer ( blockRow, blockCol );
        if ( currentBlock ) {
          for ( int i = 0; i < currentBlock->getNumRows(); ++i ) {
            currentBlock->makeRowEntries ( vec, i );
            for ( typename vector<typename Row<RealType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
              BigMatrix.add ( bigIOffset + i,  bigJOffset + it->col,  it->value );
            }
          }
        }
        bigJOffset += colSizes[blockCol];
      }
      bigIOffset += rowSizes[blockRow];
    }
  }

  //! makeRowEntries for row in sparseBlockRow as if the matrix were one big matrix without block structure
  void makeUnblockedRowEntries ( vector<typename Row<RealType>::RowEntry > &vec, const int blockRow, const int Row ) const {
    vec.clear();

    aol::RandomAccessContainer< std::vector< typename aol::Row<RealType>::RowEntry > > REVecVec ( this->getNumCols() );
    for ( int blockj = 0; blockj < this->getNumCols(); ++blockj ) {
      if ( this->getPointer ( blockRow, blockj ) ) {
        this->getReference ( blockRow, blockj ).makeRowEntries ( REVecVec[blockj], Row );
      }
    }

    int offset = 0;
    for ( int j = 0; j < this->getNumCols(); ++j ) {
      vec.reserve ( vec.size() + REVecVec[j].size() );
      for ( unsigned int jj = 0; jj < REVecVec[j].size(); ++jj ) {
        typename aol::Row<RealType>::RowEntry ent ( REVecVec[j][jj] );
        ent.col += offset;
        vec.push_back( ent );
      }
      int i = 0;
      while ( !this->getPointer(i, j) && i < this->getNumRows() ) i++;
      offset += this->getReference(i,j).getNumCols();
    }
  }

    //! makeRowEntries for row in sparseBlockRow as if the matrix were one big matrix without block structure
  void makeUnblockedSortedRowEntries ( vector<typename Row<RealType>::RowEntry > &vec, const int blockRow, const int Row ) const {
    vec.clear();

    std::vector< std::vector< typename aol::Row<RealType>::RowEntry > > REVecVec ( this->getNumCols() );
    for ( int blockj = 0; blockj < this->getNumCols(); ++blockj ) {
      if ( this->getPointer ( blockRow, blockj ) ) {
        this->getReference ( blockRow, blockj ).makeSortedRowEntries ( REVecVec[blockj], Row );
      }
    }

    int offset = 0;
    for ( int j = 0; j < this->getNumCols(); ++j ) {
      vec.reserve ( vec.size() + REVecVec[j].size() );
      for ( unsigned int jj = 0; jj < REVecVec[j].size(); ++jj ) {
        typename aol::Row<RealType>::RowEntry ent ( REVecVec[j][jj] );
        ent.col += offset;
        vec.push_back( ent );
      }
      int i = 0;
      while ( !this->getPointer(i, j) && i < this->getNumRows() ) i++;
      offset += this->getReference(i,j).getNumCols();
    }
  }

  void computeConditionNumberViaOctave() const {
    SparseMatrix<RealType> bigMat ( 1, 1 );
    addUnblockedMatrixTo ( bigMat );
    aol::computeConditionNumberViaOctave<SparseMatrix<RealType>, RealType> ( bigMat );
  }

  void computeConditionNumbersOfBlocksViaOctave() const {
    for ( int blockRow = 0; blockRow < this->getNumRows(); blockRow++ ) {
      for ( int blockCol = 0; blockCol < this->getNumCols(); blockCol++ ) {
        cerr << "Condition number of block " << blockRow << " " << blockCol << ":" << endl;
        if ( this->getPointer ( blockRow, blockCol ) )
          aol::computeConditionNumberViaOctave<MatrixType, RealType> ( this->getReference ( blockRow, blockCol ) );
        else
          cerr << "Block is empty\n";
      }
    }
  }

  bool checkForNANsAndINFs() const {
    for ( int blockRow = 0; blockRow < this->getNumRows(); blockRow++ ) {
      for ( int blockCol = 0; blockCol < this->getNumCols(); blockCol++ ) {
        if ( this->getPointer ( blockRow, blockCol ) && aol::checkForNANsAndINFs<MatrixType, RealType> ( this->getReference ( blockRow, blockCol ) ) )
          return true;
      }
    }
    return false;
  }

  //! \brief Returns the maximum number of x per row, where x is defined by a MatrixMemberFunctionPointer
  int maxNumXPerRow ( MatrixMemberFunctionPointer x ) const {
    int maxNum = 0;

    for ( int blockRow = 0; blockRow < this->getNumRows(); blockRow++ ) {
      int numLocalRows = 0;
      for ( int blockCol = 0; blockCol < this->getNumCols(); blockCol++ ) {
        if ( this->getPointer ( blockRow, blockCol ) ) {
          numLocalRows = this->getPointer ( blockRow, blockCol )->getNumRows();
          break;
        }
      }
      for ( int locRow = 0; locRow < numLocalRows; locRow++ ) {
        int maxNumThisRow = 0;
        for ( int blockCol = 0; blockCol < this->getNumCols(); blockCol++ ) {
          if ( this->getPointer ( blockRow, blockCol ) ) {
      //Call the member function given by x (if this is done more than once, using a #define macro may be a good idea). Add the return value to maxNumThisRow.
            maxNumThisRow += ( this->getPointer ( blockRow, blockCol )->*x )( locRow );
          }
        }
        maxNum = aol::Max ( maxNum, maxNumThisRow );
      }
    }

    return maxNum;
  }

  //! \brief Returns maximum number of non zeroes per row.
  //! \author Toelkes
  int maxNumNonZeroesPerRow ( ) const {
    return maxNumXPerRow( &MatrixType::numNonZeroes );
  }

  //! \brief Returns maximum number of stored entries per row.
  //! \author Toelkes
  int maxNumStoredEntriesPerRow ( ) const {
    return maxNumXPerRow( &MatrixType::numStoredEntries );
  }

  //! \brief Returns maximum number of non zero entries per row in a block computed over all blocks.
  //! \author Tatano
  int maxNumNonZeroesPerBlockRow ( ) const {
    aol::Vector<int> numOfNonZeroesPerRow;
    for (int blockCol=0; blockCol<this->getNumCols(); blockCol++) {
      for (int blockRow=0; blockRow<this->getNumRows(); blockRow++) {
        if ( this->getPointer ( blockRow, blockCol ) ) {
          numOfNonZeroesPerRow.pushBack( this->getPointer ( blockRow, blockCol )->maxNumNonZeroesPerRow() );
        }
      }
    }
    return numOfNonZeroesPerRow.getMaxValue();
  }

  void reallocate ( const int, const int ) {
    throw aol::Exception ( "aol::SparseBlockMatrix::reallocate ( int, int ) not implemented", __FILE__, __LINE__ );
  }

  /*!
   * \brief Returns a numRows \f$ \times \f$ numCols sub block matrix, where the sub block matrix consists of blocks given by blockPositions.
   * \author Toelkes
   * \param[in] numRows Number blocks per row in the sub block matrix.
   * \param[in] numCols Number blocks per column in the sub block matrix.
   * \param[in] blockPositions blockPositions[i][j] contains the position (in the original block matrix) of the block that will be the (i, j) block of the sub block matrix ((-1, -1) if the block should be left empty).
   * \param[out] block The sub block matrix.
   * \warning If the original matrix is destroyed the sub block matrix will use unallocated space!
   * \warning Do not use this function in combination with transposeTo on a SparseBlockMatrix or any other function that reallocates (sub) matrices!
  */
  void getSubBlockMatrix ( const unsigned int numRows, const unsigned int numCols,
      const std::vector < std::vector < std::vector < int > > > &blockPositions, SparseBlockMatrix < MatrixType > &block ) {
    block._m = numRows;
    block._n = numCols;

    for ( unsigned int i = 0; i < blockPositions.size (); ++i ) {
      for ( unsigned int j = 0; j < blockPositions[i].size (); ++j ) {
        if ( blockPositions[i][j][0] != -1 && blockPositions[i][j][1] != -1 ) {
          block.setPointer ( i, j, this->getPointer ( blockPositions[i][j][0], blockPositions[i][j][1] ) );
        }
        else {
          block.setPointer ( i, j, NULL );
        }

	// Never delete any sub blocks of "block", because they are managed by this matrix!
	block.setDeleteFlag ( i, j, false );
      }
    }
  }

  /*!
  * \brief sets some pointers to the pointers of block, managed by blockPositions.
  * \author Heeren
  */
  void setSubBlockMatrix ( const std::vector < std::vector < std::vector < int > > > &blockPositions, SparseBlockMatrix < MatrixType > &block ) {
    for ( unsigned int i = 0; i < blockPositions.size (); ++i )
      for ( unsigned int j = 0; j < blockPositions[i].size (); ++j ){
	int I = blockPositions[i][j][0];
	int J = blockPositions[i][j][1];
        if ( I < this->_m && J < this->_n ){
	  if ( this->getPointer ( I, J ) && this->getDeleteFlag ( I, J )  )
            delete ( this->getPointer ( I, J ) );
          this->setPointer ( I, J, block.getPointer ( i, j ) );
	  // Never delete any sub blocks of this, because they are managed by "block" now!
	  this->setDeleteFlag ( I, J, false );
	}
      }
  }

  //! scaleUnblockedRowEntries for row in sparseBlockRow as if the matrix were one big matrix without block structure
  //! \author Tatano
  void scaleUnblockedRowEntries ( const int blockRow, const int Row, const DataType factor ) {
    for ( int j = 0; j < this->getNumCols(); ++j ) {
      MatrixType* currentBlock = this->getPointer ( blockRow, j );
      if ( currentBlock )
        currentBlock->scaleRow(Row, factor);
    }
  }


  /*!
   * \brief The function scales each row of M with the inverse of the row.norm() and gives back a vector whose elements are the inverse of row.norm(). M has to be the copy of this.
   * \author Tatano
   */
  void normalizeRowsTo(SparseBlockMatrix<MatrixType> &M, aol::Vector<RealType> &vecNorms) const {
    int rowOffset = 0;
    vector<typename Row<RealType>::RowEntry > vec;
    for ( int blockRow = 0; blockRow < this->getNumRows(); blockRow++ ) {
      const MatrixType* currentBlock = this->getPointer ( blockRow, 0 );
      const int numRows = currentBlock -> getNumRows();
        if ( currentBlock ) {
          for ( int i = 0; i < numRows; i++ ) {
            this->makeUnblockedRowEntries(vec, blockRow, i);
            aol::Vector<RealType> vecVals(vec.size());
            for (int j=0; j<vec.size(); j++) {
              vecVals[j] = vec[j].value;
            }
            const RealType norm = vecVals.norm();
            if (norm > 0.) {
              vecNorms[i+rowOffset] = 1./norm;
              M.scaleUnblockedRowEntries ( blockRow, i, 1./norm );
            }
          }
        }
      rowOffset += numRows;
     }
  }

  //! scaleUnblockedRowEntries for row in sparseBlockRow as if the matrix were one big matrix without block structure
  //! \author Tatano
  void scaleUnblockedColEntries ( const int blockCol, const int Col, const DataType factor ) {
    for ( int i = 0; i < this->getNumRows(); ++i ) {
      MatrixType* currentBlock = this->getPointer ( i, blockCol );
      if ( this->getPointer ( i, blockCol )->getColumnIndexReference().size() )
        currentBlock->scaleCol(Col, factor);
    }
  }

  /*!
   * \brief The function scales the matrix M columnwise. M has to be the copy of this.
   * \author Tatano
   */
  void scaleColsWithFactorsTo(SparseBlockMatrix<MatrixType> &M, const aol::Vector<RealType> &vecFactors) const {
    int colOffset = 0;
    for ( int blockCol = 0; blockCol < this->getNumCols(); blockCol++ ) {
      const MatrixType* currentBlock = this->getPointer ( 0, blockCol );
      const int numCols = currentBlock -> getNumCols();
      if ( currentBlock ) {
        for ( int i = 0; i < numCols; i++ ) {
            M.scaleUnblockedColEntries ( blockCol, i, vecFactors[i+colOffset] );
        }
      }
      colOffset += numCols;
    }
  }
};

/** \brief A class for block matrices based on BlockOpBase which does its own memory management
 *  \todo combine this class with SparseBlockMatrix
 *  \author Schwen
 *  \ingroup Matrix
 */
template< class _MatrixType >
class BlockMatrix : public BlockOpBase< typename _MatrixType::DataType, _MatrixType > {
public:
  typedef _MatrixType MatrixType;
  typedef typename MatrixType::DataType   DataType;

public:
  /** Standard constructor creating 0x0 block matrix.
   */
  BlockMatrix ( ) : BlockOpBase< typename MatrixType::DataType, MatrixType > ( 0, 0 ) {
  }

  /** Constructor for BlockMatrix ( M, N, m, n ): M x N block structure with size m x n
   *  There should not be default behavior as there is no canonical order of the arguments.
   */
  BlockMatrix ( const int num_block_rows, const int num_block_cols, const int num_rows, const int num_cols ) : BlockOpBase< typename MatrixType::DataType, MatrixType > ( num_block_rows, num_block_cols ) {
    for ( typename vector< MatrixType* >::iterator it = this->_blockEntryPointers.begin(); it != this->_blockEntryPointers.end(); ++it ) {
      ( *it ) = new MatrixType ( num_rows, num_cols );
    }
  }

  /** Constructor for BlockMatrix ( M, N, m, n ): M x N block structure with size m x n
   *  There should not be default behavior as there is no canonical order of the arguments.
   */
  explicit BlockMatrix ( const BlockMatrix<MatrixType> &Other, const CopyFlag Flag ) : BlockOpBase< typename MatrixType::DataType, MatrixType > ( Other.getNumRows(), Other.getNumCols() ) {
    if ( Flag != STRUCT_COPY )
      throw aol::Exception( "BlockMatrix copy constructor only implemented for STRUCT_COPY" );
    for ( typename vector< MatrixType* >::iterator it = this->_blockEntryPointers.begin(); it != this->_blockEntryPointers.end(); ++it ) {
      ( *it ) = new MatrixType ( Other.getReference( 0, 0 ).getNumRows(), Other.getReference( 0, 0 ).getNumCols() );
    }
  }

  /** Constructor creating BlockMatrix for a grid (e. g. qc::GridDefinition), dim by dim blocks
   */
  template< class GridType>
  explicit BlockMatrix ( const GridType &grid )
      : BlockOpBase< typename MatrixType::DataType, MatrixType > ( grid.getDimOfWorld(), grid.getDimOfWorld() ) {

    for ( typename vector< MatrixType* >::iterator it = this->_blockEntryPointers.begin(); it != this->_blockEntryPointers.end(); ++it ) {
      ( *it ) = new MatrixType ( grid );
    }
  }

  template <qc::Dimension Dim>
  explicit BlockMatrix ( const qc::GridSize<Dim> &gridSize )
      : BlockOpBase< typename MatrixType::DataType, MatrixType > ( gridSize.Dim, gridSize.Dim ) {
    for ( typename vector< MatrixType* >::iterator it = this->_blockEntryPointers.begin(); it != this->_blockEntryPointers.end(); ++it ) {
      ( *it ) = new MatrixType ( gridSize );
    }
  }

  /** Constructor creating BlockMatrix for a grid with specified number of blocks
   */
  template< class GridType>
  BlockMatrix ( const GridType &grid, const int num_block_rows, const int num_block_cols ) : BlockOpBase< typename MatrixType::DataType, MatrixType > ( num_block_rows, num_block_cols ) {
    for ( typename vector< MatrixType* >::iterator it = this->_blockEntryPointers.begin(); it != this->_blockEntryPointers.end(); ++it ) {
      ( *it ) = new MatrixType ( grid );
    }
  }

  /** Constructor for BlockMatrix  with different sizes in the blocks.
   */
  BlockMatrix ( const Vector<int> &BlockRowSizes, const Vector<int> &BlockColSizes ) : BlockOpBase< typename MatrixType::DataType, MatrixType > ( BlockRowSizes.size(), BlockColSizes.size() ) {
    for ( int i = 0; i < BlockRowSizes.size(); ++i )
      for ( int j = 0; j < BlockColSizes.size(); ++j )
        this->setPointer ( i, j, new MatrixType ( BlockRowSizes[i], BlockColSizes[j] ) );
  }

  /** Copy constructor making deep copy */
  explicit BlockMatrix ( const BlockMatrix < MatrixType > & other ) : BlockOpBase< typename MatrixType::DataType, MatrixType > ( other.getNumRows(), other.getNumCols() ) {
    for ( int i = 0; i < this->getNumRows(); ++i )
      for ( int j = 0; j < this->getNumCols(); ++j )
        this->setPointer ( i, j, new MatrixType ( other.getReference ( i, j ) ) );
  }

  /** Assignmnent operator */
  BlockMatrix < MatrixType >& operator= ( const BlockMatrix < MatrixType > & /*other*/ ) {
    throw ( aol::UnimplementedCodeException ( "aol::BlockMatrix::operator= not implemented", __FILE__, __LINE__ ) );
    return ( *this );
  }

  /** Destructor
   */
  virtual ~BlockMatrix ( ) {
    deleteBlockEntries();
  }

  void reallocate ( const int, const int ) {
    throw aol::Exception ( "aol::BlockMatrix::reallocate ( int, int ) not implemented", __FILE__, __LINE__ );
  }

  void reallocate ( const int num_block_rows, const int num_block_cols, const int num_rows, const int num_cols ) {
    deleteBlockEntries();

    this->_m = num_block_rows;
    this->_n = num_block_cols;

    this->_blockEntryPointers.resize ( num_block_rows * num_block_cols );

    // create new matrices
    for ( typename vector< MatrixType* >::iterator it = this->_blockEntryPointers.begin(); it != this->_blockEntryPointers.end(); ++it ) {
      ( *it ) = new MatrixType ( num_rows, num_cols );
    }
  }


  /** Scale by scalar factor
  */
  BlockMatrix<MatrixType >& operator*= ( const DataType alpha ) {
    for ( typename vector< MatrixType* >::iterator it = this->_blockEntryPointers.begin(); it != this->_blockEntryPointers.end(); ++it ) {
      * ( *it ) *= alpha;
    }
    return *this;
  }

  /** Add other block matrix (which must have the correct structure) to this matrix.
   *  \todo structure comparison
   */
  BlockMatrix<MatrixType>& operator+= ( const BlockMatrix<MatrixType> &other ) {
    typename vector< MatrixType* >::iterator it;
    typename vector< MatrixType* >::const_iterator oit;
    for ( it = this->_blockEntryPointers.begin(),  oit = other._blockEntryPointers.begin() ; it != this->_blockEntryPointers.end(); it++, oit++ ) {
      * ( *it ) += * ( *oit );
    }
    return *this;
  }

public:
  void addMatrix ( const int i, const int j, const MatrixType &Matrix ) {
    this->getReference ( i, j ) += Matrix;
  }

  void addMultiple ( const  BlockMatrix<MatrixType> &Matrix, typename MatrixType::DataType Factor ) {
    for ( int i = 0; i < this->getNumRows(); i++ )
      for ( int j = 0; j < this->getNumCols(); j++ )
        this->getReference ( i, j ).addMultiple ( Matrix.getReference ( i, j ), Factor );
  }


  // add this blockMatrix to bigMatrix as if it were one big matrix without block structure. This assumes that all matrices in one row have same number of columns and vice versa.
  void addUnblockedMatrixTo ( Matrix<DataType> &bigMatrix, const bool autoResize = true ) const {

    Vector<int> rowSizes ( this->getNumRows() ), colSizes ( this->getNumCols() );
    for ( int blockRow = 0; blockRow < this->getNumRows(); blockRow++ ) {
      rowSizes[blockRow] = this->getReference ( blockRow, 0 ).getNumRows();
    }
    for ( int blockCol = 0; blockCol < this->getNumCols(); blockCol++ ) {
      colSizes[blockCol] = this->getReference ( 0, blockCol ).getNumCols();
    }

    if ( autoResize ) {
      bigMatrix.resize ( rowSizes.sum(), colSizes.sum() );
    }

    vector<typename Row<DataType>::RowEntry > vec;
    int bigIOffset = 0;
    for ( int blockRow = 0; blockRow < this->getNumRows(); blockRow++ ) {
      int bigJOffset = 0;
      for ( int blockCol = 0; blockCol < this->getNumCols(); blockCol++ ) {
        const MatrixType& currentBlock = this->getReference ( blockRow, blockCol );
        for ( int i = 0; i < currentBlock.getNumRows(); ++i ) {
          currentBlock.makeRowEntries ( vec, i );
          for ( typename vector<typename Row<DataType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
            bigMatrix.add ( bigIOffset + i,  bigJOffset + it->col,  it->value );
          }
        }
        bigJOffset += this->getReference ( 0, blockCol ).getNumCols();
      }
      bigIOffset += this->getReference ( blockRow, 0 ).getNumRows();
    }


  }

  //! Transpose each block of this BlockMatrix to the block of other BlockMatrix at transposed position
  void transposeTo ( BlockMatrix< MatrixType > &other ) const {
    for ( int i = 0; i < this->getNumRows(); ++i )
      for ( int j = 0; j < this->getNumCols(); ++j )
        this->getReference ( i, j ).transposeTo ( other.getReference ( j, i ) );
  }

  void transpose( ) {
    BlockMatrix<MatrixType> tempMat ( *this );
    tempMat.transposeTo ( *this );
  }

  void setBlockRow ( const int blockI, const int i, const MultiVector<DataType> &blockRow ) {
    for ( int j = 0; j < this->getNumCols(); ++j )
      this->getReference ( blockI, j ).setRow ( i, blockRow[j] );
  }

  // obvious similar method setBlockCol not implemented yet.

  //! Add product of two block matrices to this matrix. Left and Right type may be different.
  template < typename BMTypeL, typename BMTypeR >
  void addMatrixProduct ( const BMTypeL &Ml, const BMTypeR &Mr ) {
    if ( ( this->_m != Ml.getNumRows() ) || ( this->_n != Mr.getNumCols() ) || ( Ml.getNumCols() != Mr.getNumRows() ) ) {
      throw Exception ( "BlockMatrix::addMatrixProduct: incompatible matrices", __FILE__, __LINE__ );
    }

    for ( int i = 0; i < this->_m; ++i ) {
      for ( int j = 0; j < this->_n; ++j ) {
        for ( int k = 0; k < Ml.getNumCols(); ++k ) {
          this->getReference ( i, j ).addMatrixProduct ( Ml.getReference ( i, k ), Mr.getReference ( k, j ) );
        }
      }
    }
  }

  //! makeRowEntries for row in blockRow as if the matrix were one big matrix without block structure
  void makeUnblockedRowEntries ( vector<typename Row<DataType>::RowEntry > &vec, const int blockRow, const int Row ) const {
    vec.clear();

    std::vector< std::vector< typename aol::Row<DataType>::RowEntry > > REVecVec ( this->getNumCols() );
    for ( int blockj = 0; blockj < this->getNumCols(); ++blockj ) {
      this->getReference ( blockRow, blockj ).makeRowEntries ( REVecVec[blockj], Row );
    }

    int offset = 0;
    for ( int j = 0; j < this->getNumCols(); ++j ) {
      vec.reserve ( vec.size() + REVecVec[j].size() );
      for ( unsigned int jj = 0; jj < REVecVec[j].size(); ++jj ) {
        typename aol::Row<DataType>::RowEntry ent ( REVecVec[j][jj] );
        ent.col += offset;
        vec.push_back ( ent );
      }
      offset += this->getReference ( 0, j ).getNumCols();
    }
  }

  //! returns number of nonzero entries in specific blockRow in speciic row
  int numNonZeroes ( const int blockRow, const int Row ) const {
    int numNonZeroes = 0;
    for ( int blockj = 0; blockj < this->getNumCols(); ++blockj )
      numNonZeroes += this->getReference ( blockRow, blockj ).numNonZeroes( Row );
    return numNonZeroes;
  }

  int getTotalNumRows ( ) const {
    int counter = 0;
    for ( int i = 0; i < this->getNumRows(); ++i ) {
      counter += this->getReference ( i, 0 ).getNumRows();
    }
    return ( counter );
  }

  /*
   * \author Tatano
   */
  int getTotalNumCols ( ) const {
    int counter = 0;
    for ( int i = 0; i < this->getNumCols(); ++i ) {
      counter += this->getReference ( 0, i ).getNumCols();
    }
    return ( counter );
  }

private:
  //! destroy old contents
  void deleteBlockEntries ( ) {
    for ( typename vector< MatrixType* >::iterator it = this->_blockEntryPointers.begin(); it != this->_blockEntryPointers.end(); ++it ) {
      if ( ( *it ) != NULL ) {
        delete ( *it );
      }
      ( *it ) = NULL;
    }
  }
}; // end class BlockMatrix

}

#endif
