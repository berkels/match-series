#ifndef __CSRMATRIX_H
#define __CSRMATRIX_H

#include <rows.h>
#include <sparseMatrices.h>
#include <rectangularGrid.h>

// Deprecated typedefs still necessary for some finished projects.
typedef double _DOUBLE_PRECISION_t;
typedef int _INTEGER_t;

namespace aol {

/** \brief Matrix in compressed sparse row storage.
 *
 *  Unlike aol::CSRMatrix, this matrix can be contructed from any matrix type that supports makeSortedRowEntries.
 *
 *  \todo Merge into aol::CSRMatrix?
 *  \ingroup Matrix
 */
template <typename BaseClass = aol::GenSparseOp<double>, typename DataType = double, typename IndexType = int >
class CSR_Matrix : public aol::CSRMatrix<DataType, IndexType> {
public:
protected:
  CSR_Matrix( ) : aol::CSRMatrix<DataType, IndexType> () {}

  using aol::CSRMatrix<DataType, IndexType>::_numRows;
  using aol::CSRMatrix<DataType, IndexType>::_numCols;

public:
  template <typename MatrixType>
  explicit CSR_Matrix ( const MatrixType &Mat, CopyFlag copyFlag = DEEP_COPY ) : aol::CSRMatrix<DataType, IndexType> ( Mat.getNumRows(), Mat.getNumCols() ) {
    if ( copyFlag != DEEP_COPY )
      throw aol::Exception ( "Unsupported CopyFlag", __FILE__, __LINE__);
    init ( Mat );
  }

  explicit CSR_Matrix ( const CSR_Matrix &Mat ) : aol::CSRMatrix<DataType, IndexType> ( Mat.getNumRows(), Mat.getNumCols() ) {
    init ( Mat );
  }

  ~CSR_Matrix( ) {}

  int numNonZeroes( ) const {
    throw aol::Exception ( "CSR_Matrix::numNonZeroes() does not return what the method name suggests", __FILE__, __LINE__ );
    // attention: numNonZeroes is used within this matrix with current return value!
    return this->_indPointer[ this->_numRows ] - this->_indPointer[ 0 ];
  }

  int numNonZeroes ( IndexType RowNum ) const {
    throw aol::Exception ( "CSR_Matrix::numNonZeroes ( int RowNum ) does not return what the method name suggests", __FILE__, __LINE__ );
    return this->_indPointer[ RowNum+1 ] - this->_indPointer[ RowNum ];
  }

  int numStoredEntries( ) const {
    return this->_indPointer[ this->_numRows ] - this->_indPointer[ 0 ];
  }

  int numStoredEntries ( const IndexType RowNum ) const {
    return this->_indPointer[ RowNum+1 ] - this->_indPointer[ RowNum ];
  }

  //! Return vector of row entries. Entries need not be sorted with respect to column index and zeros may be contained.
  void makeRowEntries ( vector<typename aol::Row<DataType>::RowEntry > &vec, const IndexType RowNum ) const {
    const int l = this->_indPointer[RowNum+1] - this->_indPointer[RowNum];
    vec.resize ( l );
    const IndexType *c = & ( this->_index[this->_indPointer[RowNum]] );
    const DataType *v = & ( this->_value[this->_indPointer[RowNum]] );
    for ( int i = 0; i < l; i++ ) {
      vec[i].value = *v++;
      vec[i].col = *c++;
    }
  }

  //! Same as makeRowEntries, only that entries have to be sorted wrt column index.
  // CSRs are always sorted.
  void makeSortedRowEntries ( vector<typename aol::Row<DataType>::RowEntry > &vec, const IndexType RowNum ) const {
    makeRowEntries ( vec, RowNum );
  }

  // this only sets all elements to zero but keeps the matrix structure
  void setZero( ) {
    this->_value.setZero ();
  }

  CSR_Matrix& operator*= ( const DataType alpha ) {
    for ( int i = 0; i < numStoredEntries(); ++i ) {
      this->_value[i] *= alpha;
    }
    return ( *this );
  }

  // assumes both have the same structure, be careful.
  void addMultiple ( const CSR_Matrix &other, const DataType factor = 1.0 ) {
    const int num = numStoredEntries();
    for ( int i = 0; i < num; i++ ) {
      this->_value[i] += factor * other._values[i];
    }
  }

  void set ( IndexType I, IndexType J,  DataType value ) {
    CSRMatrix<DataType, IndexType>::set ( I, J, value );
  }

  void add ( IndexType I, IndexType J,  DataType value ) {
    CSRMatrix<DataType, IndexType>::add ( I, J, value );
  }

  DataType get ( IndexType I, IndexType J ) const {
    return CSRMatrix<DataType, IndexType>::get ( I, J );
  }

  DataType getDiag ( IndexType I ) const  {
    return CSRMatrix<DataType, IndexType>::get ( I, I );
  }

  void applyAdd ( const aol::Vector<DataType> &Arg, aol::Vector<DataType> &Dest ) const {
    CSRMatrix<DataType, IndexType>::applyAdd ( Arg, Dest );
  }

  void apply ( const aol::Vector<DataType> &Arg, aol::Vector<DataType> &Dest ) const {
    CSRMatrix<DataType, IndexType>::apply ( Arg, Dest );
  }

protected:

  template <typename MatrixType>
  void init ( const MatrixType &Mat ) {
    const int num = Mat.numStoredEntries();
    this->_numRows = Mat.getNumRows();
    this->_numCols = Mat.getNumCols();
    this->_value.resize ( num );
    this->_indPointer.resize ( this->_numRows+1 );
    this->_index.resize ( num );

    IndexType r = static_cast<IndexType> ( 0 );
    do {
      if ( r >= Mat.getNumRows () )
        break;

      this->_indPointer[r] = static_cast<IndexType> ( 0 );
    } while ( Mat.numStoredEntries ( r++ ) == 0 );
    --r;

    IndexType counter = static_cast<IndexType> ( 0 );
    for ( ; r < Mat.getNumRows (); ++r ) {
      this->_indPointer[r + 1] = this->_indPointer[r] + Mat.numStoredEntries ( r );

      vector<typename aol::Row<typename MatrixType::DataType>::RowEntry > vec;
      Mat.makeSortedRowEntries ( vec, r );

      for ( typename vector<typename aol::Row<typename MatrixType::DataType>::RowEntry >::iterator it = vec.begin();
          it != vec.end(); ++it ) {
        this->_index[counter] = it->col;
        this->_value[counter++] = it->value;
      }
    }
  }
};


} // end namespace aol

#endif
