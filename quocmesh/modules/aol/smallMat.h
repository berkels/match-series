#ifndef __SMALLMAT_H
#define __SMALLMAT_H

#include <aol.h>
#include <vec.h>
#include <smallVec.h>

namespace aol {

//! \brief A Matrix of any x- and y- dimension (template-parameters)
//!        derived from the vec-class (no use of STL-vectors) (ON)
//! @ingroup Matrix
template <int numRows, int numCols, typename _DataType>
class Mat {
protected:
  // Allow dummy matrices with numRows == 0, but avoid arrays of size 0.
  Vec<numCols, _DataType> _row[numRows > 0 ? numRows : 1];

public:
  typedef _DataType DataType;

  //! Constructor
  Mat() {
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        this->_row[i][j] = ZOTrait<_DataType>::zero;
  }

  //! Copy-constructor
  Mat ( const Mat<numRows, numCols, _DataType> &rhs ) {
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        this->_row[i][j] = rhs._row[i][j];
  }

  //! Copy-constructor for structure or data copying
  Mat ( const Mat<numRows, numCols, _DataType> &rhs, CopyFlag copyFlag ) {
    switch ( copyFlag ) {
      case DEEP_COPY:
        for ( int i = 0; i < numRows; ++i )
          for ( int j = 0; j < numCols; ++j )
            this->_row[i][j] = rhs._row[i][j];
        break;
      case STRUCT_COPY:
        for ( int i = 0; i < numRows; ++i )
          for ( int j = 0; j < numCols; ++j )
            this->_row[i][j] = ZOTrait<DataType>::zero;
        break;
      default:
        string errorMessage = strprintf( "Copying a Mat is not possible with copyFlag=%d", copyFlag );
        throw Exception( errorMessage, __FILE__, __LINE__);
    }
  }

  // ---------------------------------------------

  //! operator=
  Mat<numRows, numCols, _DataType>& operator= ( const Mat<numRows, numCols, _DataType> &rhs ) {
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        this->_row[i][j] = rhs._row[i][j];
    return *this;
  }

  Vec<numCols, _DataType>&         operator[] ( const int i )       {
#ifdef BOUNDS_CHECK
    rowBoundsCheck ( i, __FILE__, __LINE__ );
#endif
    return this->_row[i];
  }

  const Vec<numCols, _DataType>&   operator[] ( const int i ) const {
#ifdef BOUNDS_CHECK
    rowBoundsCheck ( i, __FILE__, __LINE__ );
#endif
    return this->_row[i];
  }

  _DataType get ( const int I, const int J ) const           {
#ifdef BOUNDS_CHECK
    rowBoundsCheck ( I, __FILE__, __LINE__ );
#endif
    return this->_row[I][J];
  }
  void set ( const int I, const int J, _DataType Value ) {
#ifdef BOUNDS_CHECK
    rowBoundsCheck ( I, __FILE__, __LINE__ );
#endif
    this->_row[I][J] = Value;
  }

  _DataType getDiag ( const int I ) {
#ifdef BOUNDS_CHECK
    rowBoundsCheck ( I, __FILE__, __LINE__ );
#endif
    return get(I, I);
  }

  int getNumCols() const { return numCols; }
  int getNumRows() const { return numRows; }

  //! equality comparison of two Mats
  bool operator== ( const Mat<numRows, numCols, _DataType> &other ) const {
    bool ret = true;
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        ret &= (*this)[i][j] == other[i][j];
    return ( ret );
  }

  //! inequality comparison of two Mats
  bool operator!= ( const Mat<numRows, numCols, _DataType> &other ) const {
    return ( ! ( *this == other ) );
  }

  void add ( const int I, const int J, const _DataType Value ) {
#ifdef BOUNDS_CHECK
    rowBoundsCheck ( I, __FILE__, __LINE__ );
#endif

    this->_row[I][J] += Value;
  }

  void setZero( ) {
    for ( int i = 0; i < numRows; ++i )
      this->_row[i].setZero();
  }

  void setIdentity() {
    QUOC_ASSERT(numRows == numCols);
    setZero();
    for (int i = 0; i < numRows; ++i)
      (*this)[i][i] = 1.;
  }

  void addToDiagonal( const _DataType Value ) {
    QUOC_ASSERT(numRows <= numCols);
    for (int i = 0; i < numRows; ++i)
      (*this)[i][i] += Value;
  }

  template <class T>
  void setRow ( const int Index, const Vec<numCols, T> &Vec ) {
#ifdef BOUNDS_CHECK
    rowBoundsCheck ( Index, __FILE__, __LINE__ );
#endif
    for ( int j = 0; j < numCols; ++j )
      this->_row[Index][j] = static_cast<T> ( Vec.get ( j ) );
  }

  template <class T>
  void addToRow ( const int Index, const Vec<numCols, T> &Vec ) {
#ifdef BOUNDS_CHECK
    rowBoundsCheck ( Index, __FILE__, __LINE__ );
#endif
    for ( int j = 0; j < numCols; ++j )
      this->_row[Index][j] += static_cast<T> ( Vec.get ( j ) );
  }

  template <class T>
  void getRow ( const int Index, Vec<numCols, T> &Dest ) const {
#ifdef BOUNDS_CHECK
    rowBoundsCheck ( Index, __FILE__, __LINE__ );
#endif
    for ( int j = 0; j < numCols; ++j )
      Dest[j] = static_cast<T> ( this->_row[Index][j] );
  }

  template <class T>
  void getRow ( const int Index, Vector<T> &Dest ) const {
#ifdef BOUNDS_CHECK
    rowBoundsCheck ( Index, __FILE__, __LINE__ );
#endif
    for ( int j = 0; j < numCols; ++j )
      Dest[j] = static_cast<T> ( this->_row[Index][j] );
  }

  void scaleRow ( const int Index, const DataType Factor ) {
#ifdef BOUNDS_CHECK
    rowBoundsCheck ( Index, __FILE__, __LINE__ );
#endif
    for ( int j = 0; j < numCols; ++j )
      this->_row[Index][j] *= Factor;
  }

  template <class T>
  void setCol ( const int Index, const Vec<numRows, T> &Vec ) {
    for ( int j = 0; j < numRows; ++j )
      this->_row[j][Index] = static_cast<T> ( Vec.get ( j ) );
  }

  template <class T>
  void addToCol ( const int Index, const Vec<numRows, T> &Vec ) {
    for ( int j = 0; j < numRows; ++j )
      this->_row[j][Index] += static_cast<T> ( Vec.get ( j ) );
  }

  template <class T>
  void getCol ( const int Index, Vec<numRows, T> &Dest ) const {
    for ( int j = 0; j < numRows; ++j )
      Dest[j] = static_cast<T> ( this->_row[j][Index] );
  }

  template <class T>
  void getCol ( const int Index, Vector<T> &Dest ) const {
    for ( int j = 0; j < numRows; ++j )
      Dest[j] = static_cast<T> ( this->_row[j][Index] );
  }

  template <class T>
  void getColumn ( const int Index, Vec<numRows, T> &Dest ) const {
    getCol(Index, Dest);
  }

  template <class T>
  void getColumn ( const int Index, Vector<T> &Dest ) const {
    getCol(Index, Dest);
  }

  template <class T>
  void getSubColumn ( int startRow, int col, Vector<T> & dest) const {
    QUOC_ASSERT(startRow + dest.size() <= numRows);
    for (int i = 0; i < dest.size(); ++i)
      dest[i] = (*this)[startRow + i][col];
  }

  template <class T, int size>
  void getSubColumn ( int startRow, int col, Vec<size, T> & dest) const {
    QUOC_ASSERT(startRow + size <= numRows);
    for (int i = 0; i < size; ++i)
      dest[i] = (*this)[startRow + i][col];
  }

  template <class T>
  void getSubRow ( int row, int startCol, Vector<T> & dest) const {
    QUOC_ASSERT(startCol + dest.size() <= numCols);
    for (int i = 0; i < dest.size(); ++i)
      dest[i] = (*this)[row][startCol + i];
  }

  template <class T, int size>
  void getSubRow ( int row, int startCol, Vec<size, T> & dest) const {
    QUOC_ASSERT(startCol + size <= numCols);
    for (int i = 0; i < size; ++i)
      dest[i] = (*this)[row][startCol + i];
  }

  template <class T>
  void setSubColumn ( int startRow, int col, const Vector<T> & subCol) {
    QUOC_ASSERT(startRow + subCol.size() <= numRows);
    for (int i = 0; i < subCol.size(); ++i)
      (*this)[startRow + i][col] = subCol[i];
  }

  template <class T, int size>
  void setSubColumn ( int startRow, int col, const Vec<size, T> & subCol) {
    QUOC_ASSERT(startRow + size <= numRows);
    for (int i = 0; i < size; ++i)
      (*this)[startRow + i][col] = subCol[i];
  }

  template <class T>
  void setSubRow ( int row, int startCol, const Vector<T> & subRow) {
    QUOC_ASSERT(startCol + subRow.size() <= numCols);
    for (int i = 0; i < subRow.size(); ++i)
      (*this)[row][startCol + i] = subRow[i];
  }

  template <class T, int size>
  void setSubRow ( int row, int startCol, const Vec<size, T> & subRow) {
    QUOC_ASSERT(startCol + size <= numCols);
    for (int i = 0; i < size; ++i)
      (*this)[row][startCol + i] = subRow[i];
  }

  template <class T, int rowsSubMatrix, int colsSubMatrix>
  void getSubMatrix ( int startRow, int startCol, Mat<rowsSubMatrix,colsSubMatrix,T> & subMatrix) const {
    QUOC_ASSERT(startRow + rowsSubMatrix <= numRows);
    QUOC_ASSERT(startCol + colsSubMatrix <= numCols);
    for (int i = 0; i < rowsSubMatrix; ++i)
      for (int j = 0; j < colsSubMatrix; ++j)
        subMatrix[i][j] = (*this)[startRow + i][startCol + j];
  }

  template <class T, int rowsSubMatrix, int colsSubMatrix>
  void setSubMatrix ( int startRow, int startCol, const Mat<rowsSubMatrix,colsSubMatrix,T> & subMatrix) {
    QUOC_ASSERT(startRow + rowsSubMatrix <= numRows);
    QUOC_ASSERT(startCol + colsSubMatrix <= numCols);
    for (int i = 0; i < rowsSubMatrix; ++i)
      for (int j = 0; j < colsSubMatrix; ++j)
        (*this)[startRow + i][startCol + j] = subMatrix[i][j];
  }

  template <class T>
  void multAdd ( const Vec<numCols, T> &Arg, Vec<numRows, T> &Dest ) const {
    for ( int i = 0; i < numRows; ++i ) {
      for ( int j = 0; j < numCols; ++j )
        Dest[i] += this->_row[i][j] * Arg.get ( j );
    }
  }

  template <class T>
  void multAdd ( const Vector<T> &Arg, Vec<numRows, T> &Dest ) const {
    for ( int i = 0; i < numRows; ++i ) {
      for ( int j = 0; j < numCols; ++j )
        Dest[i] += this->_row[i][j] * Arg.get ( j );
    }
  }

  template <class T>
  inline void applyAdd ( const Vec<numCols, T> &Arg, Vec<numRows, T> &Dest ) const {
    multAdd( Arg, Dest );
  }

  //! \f$ A \mapsto \alpha (A + A^T) \f$
  void symmetrize ( const _DataType Alpha = 0.5 ) {
    if( numRows != numCols )
      throw aol::Exception( "symmetrize needs numRows == numCols!", __FILE__, __LINE__ );
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j <= i; ++j )
        (*this) [i][j] = (*this) [j][i] = Alpha * ( (*this) [i][j] + (*this) [j][i] );
  }

  bool isExactlyDiagonal () const {
    bool ret = true;
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        if ( i != j && this->_row[i][j] != ZOTrait<_DataType>::zero )
          ret = false;
    return ret;
  }

  bool checkForNANsAndINFs() const {
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        if ( !aol::isFinite ( this->_row[i][j] ) ) return true;
    return false;
  }

  //! \f$ x \mapsto Ax \f$
  //
  template <class T>
  void mult ( const Vec<numCols, T> &Arg, Vec<numRows, T> &Dest ) const {
    Dest.setZero();
    multAdd ( Arg, Dest );
  }

  //! \f$ x \mapsto Ax \f$
  template <class T>
  void apply ( const Vec<numCols, T> &Arg, Vec<numRows, T> &Dest ) const {
    Dest.setZero();
    multAdd ( Arg, Dest );
  }

  //! \f$ x \mapsto A^Tx \f$
  //
  template <class T>
  void multTransposed ( const Vec<numRows, T> &Arg, Vec<numCols, T> &Dest ) const {
    Dest.setZero();
    for ( int i = 0; i < numCols; ++i ) {
      for ( int j = 0; j < numRows; ++j )
        Dest[i] += this->_row[j][i] * Arg.get ( j );
    }
  }

  void elimRowColTo ( const int Row, const int Col, Mat < numRows - 1, numCols - 1, _DataType > &subMat ) {
    for ( int row = 0; row < Row; ++row ) {
      for ( int col = 0; col < Col; ++col ) {
        subMat[row][col] = ( *this ) [row][col];
      }
      for ( int col = Col; col < numCols - 1; ++col ) {
        subMat[row][col] = ( *this ) [row][col+1];
      }
    }
    for ( int row = Row; row < numRows - 1; ++row ) {
      for ( int col = 0; col < Col; ++col ) {
        subMat[row][col] = ( *this ) [row+1][col];
      }
      for ( int col = Col; col < numCols - 1; ++col ) {
        subMat[row][col] = ( *this ) [row+1][col+1];
      }
    }
  }

  //! \f$ A \mapsto A + B \f$
  Mat<numRows, numCols, _DataType> &operator+= ( const Mat<numRows, numCols, _DataType> &Other ) {
    for ( int i = 0; i < numRows; ++i ) {
      for ( int j = 0; j < numCols; ++j ) {
        this->_row[i][j] += Other._row[i][j];
      }
    }
    return *this;
  }

  //! \f$ A \mapsto A - B \f$
  Mat<numRows, numCols, _DataType> &operator-= ( const Mat<numRows, numCols, _DataType> &Other ) {
    for ( int i = 0; i < numRows; ++i ) {
      for ( int j = 0; j < numCols; ++j ) {
        this->_row[i][j] -= Other._row[i][j];
      }
    }
    return *this;
  }

  //! ":"-Product: returns *this : mat
  _DataType ddprod ( const Mat<numRows, numCols, _DataType>& mat ) const {
    _DataType result = 0;
    for ( int i = 0; i < numRows; ++i ) {
      for ( int j = 0; j < numCols; ++j ) {
        result += this->_row[i][j] * mat._row[i][j];
      }
    }
    return result;
  }

  //! Frobenius dot product
  _DataType dotProduct ( const Mat<numRows, numCols, _DataType>& mat ) const {
    return this->ddprod( mat );
  }

  //! ":^T"-Product : returns *this : mat^T
  _DataType ddprodT ( const Mat<numRows, numCols, _DataType>& mat ) const {
    _DataType result = 0;
    for ( int i = 0; i < numRows; ++i ) {
      for ( int j = 0; j < numCols; ++j ) {
        result += this->_row[i][j] * mat._row[j][i];
      }
    }
    return result;
  }

  Vec<numRows, _DataType> operator * ( const Vec<numCols, _DataType>& vec ) const {
    Vec<numRows, _DataType> Res;
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        Res[i] += this->_row[i][j] * vec[j];

    return Res;
  }

  Vector<_DataType> operator * ( const Vector<_DataType>& vec ) const {
    Vector<_DataType> Res(numRows);
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        Res[i] += this->_row[i][j] * vec[j];

    return Res;
  }

  //! Writes the transpose of this to itself. See the other transpose* methods
  //! for the naming convention of data im/export
  void transpose( ) {
    inplaceTranspose ( *this );
  }


  // writes the transpose of mat to this Mat
  void transposeFrom ( const Mat<numCols, numRows, _DataType> &mat ) {
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        (*this)[i][j] = mat[j][i];
  }

  // writes the transpose of this Mat to mat
  void transposeTo ( Mat<numCols, numRows, _DataType> &mat ) const {
    for ( int i = 0; i < numCols; ++i )
      for ( int j = 0; j < numRows; ++j )
        mat[i][j] = (*this)[j][i];
  }

  // returns a Mat which is the transpose of this Mat
  Mat<numCols, numRows, _DataType> transposed ( ) const {
    Mat<numCols, numRows, _DataType> ret;
    this->transposeTo ( ret );
    return ret;
  }

  //! \f$ A \mapsto \alpha A \f$
  Mat<numRows, numCols, _DataType> &operator*= ( const _DataType Alpha ) {
    for ( int i = 0; i < numRows; ++i )
      this->_row[i] *= Alpha;
    return *this;
  }

  //! \f$ A \mapsto \alpha^{-1} A \f$
  Mat<numRows, numCols, _DataType> &operator/= ( const _DataType Alpha ) {
    for ( int i = 0; i < numRows; ++i )
      this->_row[i] /= Alpha;
    return *this;
  }

  void makeTensorProduct ( const Vec<numRows, _DataType> &A, const Vec<numCols, _DataType> &B ) {
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        this->_row[i][j] = A[i] * B[j];
  }

  void addTensorProduct ( const Vec<numRows, _DataType> &A, const Vec<numCols, _DataType> &B ) {
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        this->_row[i][j] += A[i] * B[j];
  }

  //! \f$ A = I-(a\otimes a) \f$
  void makeProjectionMatrix ( const Vec<numRows, _DataType> &A ) {
#ifdef BOUNDS_CHECK
  if( numRows != numCols )
    throw aol::Exception( "makeProjectionMatrix needs numRows == numCols!", __FILE__, __LINE__ );
#endif
    for ( int i = 0; i < numRows; ++i ){
      this->_row[i][i] = aol::ZOTrait<_DataType>::one - aol::Sqr( A[i] );

      for( int j = i+1; j < numCols; ++j){
        this->_row[i][j] = - A[i] * A[j];
      }

      for( int j = 0; j < i; ++j){
        this->_row[i][j] = this->_row[j][i];
      }
    }
  }


  template <class AnotherType, int dimBoth>
  void makeProduct ( const Mat<numRows, dimBoth, AnotherType> &A,
                     const Mat<dimBoth, numCols, AnotherType> &B ) {
    setZero( );
    for ( int i = 0; i < numRows; ++i ) {
      for ( int j = 0; j < numCols; ++j ) {
        for ( int k = 0; k < dimBoth; ++k ) {
          this->_row[i][j] += static_cast<_DataType> ( A[i][k] * B[k][j] );
        }
      }
    }
  }

  template <class AnotherType, int dimBoth>
  void makeProductAtransposedB ( const Mat<dimBoth, numRows, AnotherType> &A,
                                 const Mat<dimBoth, numCols, AnotherType> &B ) {
    setZero( );
    for ( int i = 0; i < numRows; ++i ) {
      for ( int j = 0; j < numCols; ++j ) {
        for ( int k = 0; k < dimBoth; ++k ) {
          this->_row[i][j] += static_cast<_DataType> ( A[k][i] * B[k][j] );
        }
      }
    }
  }

  template <class AnotherType, int dimBoth>
  void makeProductABtransposed ( const Mat<numRows, dimBoth, AnotherType> &A,
                                 const Mat<numCols, dimBoth, AnotherType> &B ) {
    setZero( );
    for ( int i = 0; i < numRows; ++i ) {
      for ( int j = 0; j < numCols; ++j ) {
        for ( int k = 0; k < dimBoth; ++k ) {
          this->_row[i][j] += static_cast<_DataType> ( A[i][k] * B[j][k] );
        }
      }
    }
  }

  //! \f$ \mbox{*this} \mapsto ABA^{T} \f$
  template <class AnotherType, int dimBoth>
  void makeProductABAtransposed ( const Mat<numRows, dimBoth, AnotherType> &A,
                                  const Mat<dimBoth, dimBoth, AnotherType> &B ) {
    if ( numCols != numRows )
      throw aol::Exception ( "Mat<...>::makeProductABAtransposed works only for quadratic matrices", __FILE__, __LINE__ );

    Mat<numRows, dimBoth, AnotherType> AB;
    AB.makeProduct ( A, B );
    makeProductABtransposed( AB, A );
  }

  //! \f$ \mbox{*this} \mapsto \frac{1}{2}(A^{T}B + B^{T}A) \f$
  template <class AnotherType, int dimBoth>
  void makeSymmetricProduct ( const Mat<dimBoth, numRows, AnotherType> &A,
                                 const Mat<dimBoth, numCols, AnotherType> &B ) {
    if ( numCols != numRows ) {
      throw aol::Exception ( "Mat<...>::makeSymmetricProduct works only for quadratic matrices", __FILE__, __LINE__ );
    }
    setZero( );
    Mat<numRows, numCols, _DataType> aux;
    aux.makeProductAtransposedB( A, B );
    makeProductAtransposedB( B, A );
    (*this)+=aux;
    (*this)*=.5;
  }

  //! dagger product: returns *this dagger A
  //! with M dagger A = (M M^T)^-1 M A^T
  Mat<numRows, numRows, _DataType> daggerProduct(const Mat<numRows, numCols, _DataType> & A) const
  {
    Mat<numRows, numRows, _DataType> MMT;
    MMT.makeProductABtransposed(*this, *this);
    Mat<numRows, numRows, _DataType> MMTinv;
    MMTinv.makeInverse(MMT);

    // abuse MMT for *this * A^T
    MMT.makeProductABtransposed(*this, A);

    Mat<numRows, numRows, _DataType> erg;
    erg.makeProduct(MMTinv, MMT);
    return erg;
  }

  //! sharp product: returns *this sharp A
  //! with M sharp A = (M M^T)^-1 A M^T
  Mat<numRows, numRows, _DataType> sharpProduct(const Mat<numRows, numCols, _DataType> & A) const
  {
    Mat<numRows, numRows, _DataType> MMT;
    MMT.makeProductABtransposed(*this, *this);
    Mat<numRows, numRows, _DataType> MMTinv;
    MMTinv.makeInverse(MMT);

    // abuse MMT for *this * A^T
    MMT.makeProductABtransposed(A, *this);

    Mat<numRows, numRows, _DataType> erg;
    erg.makeProduct(MMTinv, MMT);
    return erg;
  }

  _DataType normSqr() const {
    _DataType res = 0;
    for ( int i = 0; i < numRows; ++i ) {
      res += this->_row[i].normSqr();
    }
    return res;
  }

  //! Frobenius norm.
  _DataType norm() const {
    return sqrt( this->normSqr() );
  }

  void clamp( const _DataType Min, const _DataType Max ) {
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        this->_row[i][j] = aol::Clamp( this->_row[i][j], Min, Max );
  }

  _DataType infinityNorm() const {
    _DataType res = 0;
    for ( int i = 0; i < numRows; ++i ) {
      _DataType temp = 0;
      for ( int j = 0; j < numCols; ++j ) {
        temp += Abs ( this->_row[i][j] );
      }
      if ( temp > res ) res = temp;
    }
    return ( res );
  }


  ostream& print ( ostream& out ) const {
    for ( int i = 0; i < numRows - 1; ++i )
      out << this->_row[i] << endl;
    out << this->_row[numRows-1];
    return ( out );
  }

  istream& read ( istream& in ) {
    for ( int i = 0; i < numRows; ++i )
      in >> this->_row [i];
    return in;
  }

  //! swap two rows
  void swapRows ( const int i, const int j ) {
    _DataType tmp;
    for ( int l = 0; l < numCols; ++l ) {
      tmp = this->_row[i][l];
      this->_row[i][l] = this->_row[j][l];
      this->_row[j][l] = tmp;
    }
  }

  //! operator *=
  Mat<numRows, numCols, _DataType> &operator*= ( const Mat<numCols, numCols, _DataType>& mat ) {
    if ( static_cast<const void*>(this) == static_cast<const void*>(&mat) ) {
      cerr << "Matrix<_DataType>::operator*= :  don't multiply with the same matrix!" << endl;
    } else {
      aol::Vec<numCols, _DataType> vec;
      int i, j, k ;
      _DataType val;
      for ( i = 0; i < numRows ; ++i ) {
        vec = _row[i];

        for ( j = 0; j < numCols ; ++j ) {
          val = 0;
          for ( k = 0; k < numCols; ++k ) {
            // Note: vec[k] is a copy of this->_row[i][k]
            val += mat[k][j] * vec[k];
          }
          _row[i][j] = val;
        }
      }
    }
    return *this;
  }

  //! \f$ A \mapsto BA \f$
  void leftMult ( const Mat<numRows, numRows, _DataType>& mat ) {
    if ( static_cast<const void*>(this) == static_cast<const void*>(&mat) )
      cerr << "Matrix<_DataType>::leftMult :  don't multiply with the same matrix!" << endl;
    else {
      aol::Vec<numRows, _DataType> vec;
      for ( int j = 0; j < numCols; ++j ) {
        // save j-th column in vec
        for ( int k = 0; k < numRows; ++k )
          vec[k] = _row[k][j];
        for ( int i = 0; i < numRows ; ++i ) {
          _row[i][j] = 0;
          for ( int k = 0; k < numRows; ++k )
            _row[i][j] += mat[i][k] * vec[k];
        }
      }
    }
  }

  //! \f$ A \mapsto AB \f$
  void rightMult ( const Mat<numRows, numRows, _DataType>& mat ) {
    if ( static_cast<const void*>(this) == static_cast<const void*>(&mat) )
      cerr << "Matrix<_DataType>::rightMult :  don't multiply with the same matrix!" << endl;
    else {
      aol::Vec<numRows, _DataType> vec;
      for ( int j = 0; j < numRows; ++j ) {
        // save j-th row in vec
        for ( int k = 0; k < numRows; ++k )
          vec[k] = _row[j][k];
        for ( int i = 0; i < numCols; ++i ) {
          _row[j][i] = 0;
          for ( int k = 0; k < numCols; ++k )
            _row[j][i] += vec[k] * mat[k][i];
        }
      }
    }
  }

  //! \f$ A \mapsto A + \alpha B \f$
  Mat<numRows, numCols, _DataType> & addMultiple ( const Mat<numRows, numCols, _DataType>& mat, const _DataType& alpha ) {
    if ( static_cast<const void*>(this) == static_cast<const void*>(&mat) )
      cerr << "Matrix<_DataType>::addMultiple :  don't add the same matrix!" << endl;
    else {
      for ( int j = 0; j < numCols; ++j )
        for ( int i = 0; i < numRows ; ++i )
            _row[i][j] += alpha * mat[i][j];
    }
    return *this;
  }


  /** computes the inverse of the given argument matrix and stores it into this matrix.
   *  This method uses the Gauss-Jacobi Algorithm and throws an exception if the matrix
   *  is not invertible
   */
  void makeInverse ( const Mat<numRows, numCols, _DataType> &a ) {
    if ( numCols != numRows ) {
      throw aol::Exception ( "Mat<...>::makeInverse works only for quadratic matrices", __FILE__, __LINE__ );
    }

    const int N = numCols;
    int  indxc[numCols], indxr[numCols], ipiv[numCols] = { 0 };
    int  i, icol = 0, irow = 0, j, k, l, ll;
    _DataType big, dum, pivinv;

    ( *this ) = a;

    for ( i = 0; i < N; ++i ) {
      big = 0.0;
      for ( j = 0; j < N; ++j )
        if ( ipiv[j] != 1 ) {
          for ( k = 0; k < N; ++k ) {
            if ( ipiv[k] == 0 ) {
              if ( fabs ( this->_row[j][k] ) >= big ) {
                big = fabs ( this->_row[j][k] );
                irow = j;
                icol = k;
              }
            } else
              if ( ipiv[k] > 1 ) throw aol::Exception ( "Mat<...>::makeInverse: Singular Matrix-1", __FILE__, __LINE__ );
          }
        }
      ++ ( ipiv[icol] );
      if ( irow != icol ) {
        for ( l = 0; l < N; ++l ) std::swap ( this->_row[irow][l], this->_row[icol][l] );
      }
      indxr[i] = irow;
      indxc[i] = icol;
      if ( this->_row[icol][icol] == 0.0 ) throw aol::Exception ( "Mat<...>::makeInverse: Singular Matrix-2", __FILE__, __LINE__ );
      pivinv = aol::NumberTrait<_DataType>::one / this->_row[icol][icol];
      this->_row[icol][icol] = 1.0;

      for ( l = 0; l < N; ++l ) this->_row[icol][l] *= pivinv;
      for ( ll = 0; ll < N; ++ll )
        if ( ll != icol ) {
          dum = this->_row[ll][icol];
          this->_row[ll][icol] = 0.0;
          for ( l = 0; l < N; ++l ) this->_row[ll][l] -= this->_row[icol][l] * dum;
        }
    }
    for ( l = N - 1; l >= 0; l-- ) {
      if ( indxr[l] != indxc[l] )
        for ( k = 0; k < N; ++k )
          std::swap ( this->_row[k][indxr[l]], this->_row[k][indxc[l]] );
    }
  }

  //! computes the trace of this matrix.
  _DataType tr() const {
    _DataType tr = static_cast<DataType> ( 0 );
    for ( int i = 0; i < numRows; ++i )
      tr += this->_row[i][i];

    return tr;
  }

protected:
#ifdef BOUNDS_CHECK
  inline void rowBoundsCheck ( const int i, const char* fi, const int li ) const {
    if ( i < 0 || i >= numRows ) {
      char errmsg[1024];
      sprintf( errmsg, "aol::Mat: index %d out of bounds (%d)", i, numRows );
      throw aol::OutOfBoundsException ( errmsg, fi, li );
    }
  }
#endif
};

//! transpose symmetrix Mat<N,N, _DataType> to itself (in place)
template <int N, typename _DataType>
void inplaceTranspose ( aol::Mat<N, N, _DataType> &mat ) {
  for ( int i = 0; i < N; ++i ) {
    for ( int j = i + 1; j < N; ++j ) {
      const _DataType t = mat[i][j];
      mat[i][j] = mat[j][i];
      mat[j][i] = t;
    }
  }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// THE CLASS MATRIX22 DERIVED FROM MAT
// / /////////////////////////////////////////////////////////

//! \brief A simple 2x2 matrix.
//!        Including all methods which aren't in class Mat
//! @ingroup Matrix
template <typename _DataType>
class Matrix22 : public Mat<2, 2, _DataType> {


public:
  typedef _DataType DataType;

  //! Constructor
  Matrix22()
      : Mat<2, 2, _DataType>() { }

  //! Constructor (for return value optimization)
  //! fills this with:
  //! [ a b ]
  //! [ c d ]
  Matrix22 ( const _DataType a, const _DataType b,
             const _DataType c, const _DataType d ) {
    this->_row[0][0] = a;     this->_row[0][1] = b;
    this->_row[1][0] = c;     this->_row[1][1] = d;
  }
  
  //! fills this with:
  //! [ a a ]
  //! [ a a ]
  Matrix22 ( const _DataType a ) {
    this->_row[0][0] = a;     this->_row[0][1] = a;
    this->_row[1][0] = a;     this->_row[1][1] = a;
  }

  //! Copy-constructor
  Matrix22 ( const Matrix22<_DataType> &rhs )
      : Mat<2, 2, _DataType> ( rhs ) {}

  Matrix22 ( const Mat<2, 2, _DataType> & rhs )
    : Mat<2, 2, _DataType> ( rhs ) {}

  //! fills matrix with parameters passed
  void fill ( const _DataType a, const _DataType b,
              const _DataType c, const _DataType d ) {

    this->_row[0][0] = a;
    this->_row[0][1] = b;
    this->_row[1][0] = c;
    this->_row[1][1] = d;
  }

  void fill ( const aol::Vector<_DataType> &Vec ) {
    fill ( Vec[0], Vec[1], Vec[2], Vec[3] );
  }

  //! writes the inverse of this matrix to the argument
  //! \attention opposite of aol::Mat<>::makeInverse
  void make_inverse ( Matrix22<_DataType> &inv ) const {
    _DataType determinante = det();
    if ( determinante == _DataType ( 0 ) ) {
      throw ( Exception ( "matrix not invertible.", __FILE__, __LINE__ ) );
    }

    inv[0][0] =  ( *this ) [1][1];
    inv[0][1] = - ( *this ) [0][1];
    inv[1][0] = - ( *this ) [1][0];
    inv[1][1] =  ( *this ) [0][0];
    inv /= determinante;
  }

  //! return the inverse if it exists (otherwise throw...)
  Matrix22<_DataType> inverse() const {
    Matrix22<_DataType> inv;
    make_inverse ( inv );
    return inv;
  }
  
  //! inverts *this (otherwise throw...)
  void invert() {
    Matrix22<_DataType> inv(*this);
    inv.make_inverse ( *this );
  }

  void setIdentity( ) {
    for ( int i = 0; i < 2; ++i ) {
      for ( int j = 0; j < 2; ++j ) {
        this->_row[i][j] = static_cast<_DataType> ( i == j );
      }
    }
  }

  //! The cofactor matrix is defined as \f$ \mbox{Cof} A := \mbox{det} A \cdot A^{-T} \f$
  //
  template <class T>
  void makeCofactorMatrix ( const Matrix22<T> &Mat ) {
    this->_row[0][0] =   Mat._row[1][1];
    this->_row[1][1] =   Mat._row[0][0];
    this->_row[1][0] = - Mat._row[0][1];
    this->_row[0][1] = - Mat._row[1][0];
  }

  //! qc::Computes the matrix \f$ \partial_{a_{ij}} \mbox{Cof} A \f$.
  //
  template <typename T>
  void makeDerivedCofactorMatrix ( const Matrix22<T> &, const int I, const int J ) {
    this->setZero();
    if ( I == J ) this->_row[ ( I+1 ) %2 ][ ( J+1 ) %2 ] = 1.;
    else this->_row[ ( I+1 ) %2 ][ ( J+1 ) %2 ] = -1.;
  }

  template <class T>
  static void multWithDerivedCofactorMatrix ( const int I, const int J,
                                              const Vec2<T> &Arg, Vec2<T> &Dest ) {
    if ( I == J ) {
      Dest[ I ] = 0.;
      Dest[ ( I+1 ) %2 ] = Arg.get ( ( I + 1 ) % 2 );
    } else {
      Dest[ J ] =  - Arg.get ( I );
      Dest[ I ] = 0.;
    }
  }

  //! computes the determinant of this matrix.
  _DataType det( ) const {
    return this->_row[0][0]*this->_row[1][1] - this->_row[1][0]*this->_row[0][1];
  }

  //! computes the trace of this matrix.
  _DataType tr( ) const {
    return this->_row[0][0] + this->_row[1][1];
  }

  using Mat<2, 2, _DataType>::operator*=;

  /**
   * \f$ A\mapsto A*B \f$
   */
  Matrix22<_DataType> &operator*= ( const Matrix22<_DataType> &Other ) {
    Matrix22<_DataType> tmp;
    tmp = *this;
    this->setZero( );
    for ( int i = 0; i < 2; ++i ) {
      for ( int j = 0; j < 2; ++j ) {
        for ( int k = 0; k < 2; ++k ) {
          this->_row[i][j] += tmp._row[i][k] * Other._row[k][j];
        }
      }
    }
    return *this;
  }

  void eigenValues ( Vec2<_DataType> &Eig ) const;

  void eigenVector ( Vec2<_DataType> &EigV, _DataType lambda ) const {
    EigV[0] = -this->_row[0][1];
    EigV[1] = this->_row[0][0] - lambda;

    if ( EigV [0] == 0 && EigV [1] == 0 ) {

      EigV[0] = this->_row[1][1] - lambda;
      EigV[1] = -this->_row[1][0];
    }

    EigV.normalize ();
  }

  //! Stores eigenvectors in columns
  void eigenVectors ( Matrix22 &Eig, Vec2<_DataType> lambda ) const {
    if ( lambda [0] == lambda [1] ) {
      Eig.setIdentity ();
      return;
    }

    Vec2<_DataType> eig [2];
    for ( int i = 0; i < 2; ++i )
      eigenVector ( eig [i], lambda [i] );

    Eig.fill ( eig [0][0], eig [1][0], eig [0][1], eig [1][1] );
  }

  //! Stores eigenvectors in columns
  void eigenVectors ( Matrix22 &Eig ) const {
    Vec2<_DataType> lambda;
    eigenValues ( lambda );

    if ( lambda [0] == lambda [1] ) {
      Eig.setIdentity ();
      return;
    }

    Vec2<_DataType> eig [2];
    for ( int i = 0; i < 2; ++i )
      eigenVector ( eig [i], lambda [i] );

    Eig.fill ( eig [0][0], eig [1][0], eig [0][1], eig [1][1] );
  }

  Vec2<_DataType> operator* ( const Vec2<_DataType>& vec ) const {
    Vec2<_DataType> Res;
    for ( int i = 0; i < 2; ++i )
      for ( int j = 0; j < 2; ++j )
        Res[i] += this->_row[i][j] * vec[j];

    return( Res );
  }

  //! makes this matrix a rotation matrix with angle Alpha
  void makeRotation ( _DataType Alpha ) {
    const _DataType co = cos ( Alpha ), si = sin ( Alpha );
    this->_row[0][0] =  co;
    this->_row[0][1] = -si;
    this->_row[1][0] =  si;
    this->_row[1][1] =  co;
  }
  
  //! rotates this matrix to a different frame, i.e. \f$ A \mapsto R A R^T \f$
  void rotate ( _DataType Alpha ) {
    Matrix22 rot;
    rot.makeRotation( Alpha );
    Matrix22 tmp( *this );
    this->setZero();
    
    for ( int i = 0; i < 2; ++i )
      for ( int l = 0; l < 2; ++l )
        for ( int k = 0; k < 2; ++k )
          for ( int j = 0; j < 2; ++j )
            this->_row[i][l] += tmp._row[j][k] * rot._row[i][j] * rot._row[l][k];
  }
};


// forward declaration for conversion constructor
template <typename _DataType> class Matrix33Symm;

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// THE CLASS MATRIX33 DERIVEDED FROM MAT
// / /////////////////////////////////////////////////////////

//! \brief A simple 3x3 matrix.
//!        Including all methods which aren't in class Mat
//! @ingroup Matrix
template <typename _DataType>
class Matrix33 : public Mat<3, 3, _DataType> {

public:
  typedef _DataType DataType;

  //! Constructor
  Matrix33()
      : Mat<3, 3, _DataType>() { }

  //! Constructor (for return value optimization)
  //! fills this with:
  //! [ a b c ]
  //! [ d e f ]
  //! [ g h i ]
  Matrix33 ( const _DataType a, const _DataType b, const _DataType c,
             const _DataType d, const _DataType e, const _DataType f,
             const _DataType g, const _DataType h, const _DataType i ) {
    this->_row[0][0] = a;     this->_row[0][1] = b;    this->_row[0][2] = c;
    this->_row[1][0] = d;     this->_row[1][1] = e;    this->_row[1][2] = f;
    this->_row[2][0] = g;     this->_row[2][1] = h;    this->_row[2][2] = i;
  }

  //! Copy-constructor
  Matrix33 ( const Matrix33<_DataType> &rhs )
      : Mat<3, 3, _DataType> ( rhs ) {}

  //! Copy-constructor
  Matrix33 ( const Mat<3, 3, _DataType> &rhs )
      : Mat<3, 3, _DataType> ( rhs ) {}

  //! Copy from symmetric matrix
  explicit Matrix33 ( const Matrix33Symm<_DataType>& mat );

  //! fills matrix with parameters passed
  void fill (  const _DataType a, const _DataType b, const _DataType c,
               const _DataType d, const _DataType e, const _DataType f,
               const _DataType g, const _DataType h, const _DataType i ) {

    this->_row[0][0] = a;     this->_row[0][1] = b;    this->_row[0][2] = c;
    this->_row[1][0] = d;     this->_row[1][1] = e;    this->_row[1][2] = f;
    this->_row[2][0] = g;     this->_row[2][1] = h;    this->_row[2][2] = i;
  }

  //! calc the inverse matrix with adjoints
  Matrix33<_DataType> inverse() const {
    _DataType determinante ( det() );
    if ( determinante == _DataType ( 0 ) ) {
      throw ( Exception ( "matrix not invertible.", __FILE__, __LINE__ ) );
    }

    Matrix33 res;
    _DataType A;

    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        // calc the determinant, the sign is automatically ok
        A = this->_row[ ( i+1 ) % 3][ ( j+1 ) % 3] * this->_row[ ( i+2 ) % 3][ ( j+2 ) % 3]
            - this->_row[ ( i+1 ) % 3][ ( j+2 ) % 3] * this->_row[ ( i+2 ) % 3][ ( j+1 ) % 3];
        // consider the transposed matrix:
        res.set ( j, i, A / determinante );
      }
    }

    return res;
  }

  void setIdentity( ) {
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        ( this->_row[i] ) [j] = static_cast<_DataType> ( i == j );
      }
    }
  }

  using Mat<3, 3, _DataType>::setRow;

  void setRow ( const int Index, const _DataType V1, _DataType V2, _DataType V3 ) {
    ( this->_row[Index] ) [0] = V1;
    ( this->_row[Index] ) [1] = V2;
    ( this->_row[Index] ) [2] = V3;
  }

  //! The Cofactormatrix is defined as \f$ \mbox{Cof} A := \mbox{det} A \cdot A^{-T} \f$
  //
  template <class T>
  void makeCofactorMatrix ( const Matrix33<T> &Mat ) {
    this->_row[0][0] =   Mat.get ( 1, 1 ) * Mat.get ( 2, 2 ) - Mat.get ( 1, 2 ) * Mat.get ( 2, 1 );
    this->_row[1][0] = - Mat.get ( 0, 1 ) * Mat.get ( 2, 2 ) + Mat.get ( 2, 1 ) * Mat.get ( 0, 2 );
    this->_row[2][0] =   Mat.get ( 0, 1 ) * Mat.get ( 1, 2 ) - Mat.get ( 1, 1 ) * Mat.get ( 0, 2 );

    this->_row[0][1] = - Mat.get ( 1, 0 ) * Mat.get ( 2, 2 ) + Mat.get ( 2, 0 ) * Mat.get ( 1, 2 );
    this->_row[1][1] =   Mat.get ( 0, 0 ) * Mat.get ( 2, 2 ) - Mat.get ( 0, 2 ) * Mat.get ( 2, 0 );
    this->_row[2][1] = - Mat.get ( 0, 0 ) * Mat.get ( 1, 2 ) + Mat.get ( 1, 0 ) * Mat.get ( 0, 2 );

    this->_row[0][2] =   Mat.get ( 1, 0 ) * Mat.get ( 2, 1 ) - Mat.get ( 1, 1 ) * Mat.get ( 2, 0 );
    this->_row[1][2] = - Mat.get ( 0, 0 ) * Mat.get ( 2, 1 ) + Mat.get ( 0, 1 ) * Mat.get ( 2, 0 );
    this->_row[2][2] =   Mat.get ( 0, 0 ) * Mat.get ( 1, 1 ) - Mat.get ( 1, 0 ) * Mat.get ( 0, 1 );
  }



  //! qc::Computes the matrix \f$ \partial_{a_{ij}} \mbox{Cof} A \f$.
  //
  template <class T>
  void makeDerivedCofactorMatrix ( const Matrix33<T> &Mat, int I, int J ) {
    this->setZero();

    this->_row[ (I + 1) % 3 ][ (J + 1) % 3 ] =   Mat.get( (I + 2) % 3, (J + 2) % 3 );
    this->_row[ (I + 2) % 3 ][ (J + 2) % 3 ] =   Mat.get( (I + 1) % 3, (J + 1) % 3 );
    this->_row[ (I + 1) % 3 ][ (J + 2) % 3 ] = - Mat.get( (I + 2) % 3, (J + 1) % 3 );
    this->_row[ (I + 2) % 3 ][ (J + 1) % 3 ] = - Mat.get( (I + 1) % 3, (J + 2) % 3 );

  }

  //! computes the determinant of this matrix.
  //
  _DataType det() const {
    return this->_row[0][0]*this->_row[1][1]*this->_row[2][2] + this->_row[0][1]*this->_row[1][2]*this->_row[2][0] +
           this->_row[0][2]*this->_row[1][0]*this->_row[2][1] - this->_row[0][2]*this->_row[1][1]*this->_row[2][0] -
           this->_row[0][1]*this->_row[1][0]*this->_row[2][2] - this->_row[0][0]*this->_row[1][2]*this->_row[2][1];
  }

  //! computes the trace of this matrix.
  //
  _DataType tr() const {
    return this->_row[0][0] + this->_row[1][1] + this->_row[2][2];
  }

  using Mat<3, 3, _DataType>::operator*=;

  //! \f$ A\mapsto A*B \f$
  //
  Matrix33<_DataType> &operator*= ( const Matrix33<_DataType> &Other ) {
    Matrix33<_DataType> tmp;
    tmp = *this;
    this->setZero( );
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        for ( int k = 0; k < 3; ++k ) {
          this->_row[i][j] += tmp._row[i][k] * Other._row[k][j];
        }
      }
    }
    return *this;
  }

  //! \f$ A \mapsto A^T \f$
  //
  void transpose() {
    inplaceTranspose ( *this );
  }

  //! return (this matrix) transposed (two copies necessary)
  aol::Matrix33<_DataType> transposed () const {
    aol::Matrix33<_DataType> tmp = *this;
    tmp.transpose();
    return ( tmp );
  }

  //! Set to rotation matrix in the yz plane (about x axis)
  void setRotationAboutX ( _DataType alpha ) {
    const _DataType co = cos ( alpha ), si = sin ( alpha );
    this->_row[0][0] = 1.0;    this->_row[0][1] = 0.0;    this->_row[0][2] = 0.0;
    this->_row[1][0] = 0.0;    this->_row[1][1] = co ;    this->_row[1][2] = -si;
    this->_row[2][0] = 0.0;    this->_row[2][1] = si ;    this->_row[2][2] = co ;
  }

  //! Set to rotation matrix in the xz plane (about y axis)
  void setRotationAboutY ( _DataType alpha ) {
    const _DataType co = cos ( alpha ), si = sin ( alpha );
    this->_row[0][0] = co ;    this->_row[0][1] = 0.0;    this->_row[0][2] = -si;
    this->_row[1][0] = 0.0;    this->_row[1][1] = 1.0;    this->_row[1][2] = 0.0;
    this->_row[2][0] = si ;    this->_row[2][1] = 0.0;    this->_row[2][2] = co ;
  }

  //! Set to rotation matrix in the xy plane (about z axis)
  void setRotationAboutZ ( _DataType alpha ) {
    const _DataType co = cos ( alpha ), si = sin ( alpha );
    this->_row[0][0] = co ;    this->_row[0][1] =-si ;    this->_row[0][2] = 0.0;
    this->_row[1][0] = si ;    this->_row[1][1] = co ;    this->_row[1][2] = 0.0;
    this->_row[2][0] = 0.0;    this->_row[2][1] = 0.0;    this->_row[2][2] = 1.0;
  }

  //! Set to rotation matrix about the ith axis (angle=``alpha'',axis=``Axis'')
  void setRotationAboutAxis ( _DataType Alpha, short Axis ) {
    switch ( Axis ) {
    case 0: this->setRotationAboutX( Alpha ); break;
    case 1: this->setRotationAboutY( Alpha ); break;
    case 2: this->setRotationAboutZ( Alpha ); break;
    default: this->setRotationAboutZ( Alpha );
    }
  }

  //! Set to rotation matrix about the axis specified with the qc::Comp argument.
  void setRotationAboutAxis ( _DataType Alpha, qc::Comp Component ) {
    switch ( Component ) {
    case qc::QC_X: this->setRotationAboutX( Alpha ); break;
    case qc::QC_Y: this->setRotationAboutY( Alpha ); break;
    case qc::QC_Z: this->setRotationAboutZ( Alpha ); break;
    default: throw Exception( "Invalid qc::Comp value", __FILE__, __LINE__);
    }
  }

  Vec3<_DataType> operator* ( const Vec3<_DataType>& vec ) const {
    Vec3<_DataType> Res;
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        Res[i] += this->_row[i][j] * vec[j];

    return( Res );
  }

};


/** a simple 4x4 matrix
 */
template <typename _DataType>
class Matrix44 : public Mat<4, 4, _DataType> {

public:
  typedef _DataType DataType;

  /** a default constructor
   */
  Matrix44<_DataType> ( void ) : Mat<4, 4, _DataType>() {}

  //! Constructor
  //! fills this with:
  //! [ a b c d ]
  //! [ e f g h ]
  //! [ i j k l ]
  //! [ m n o p ]
  Matrix44 ( const _DataType a, const _DataType b, const _DataType c, const _DataType d,
             const _DataType e, const _DataType f, const _DataType g, const _DataType h,
             const _DataType i, const _DataType j, const _DataType k, const _DataType l,
             const _DataType m, const _DataType n, const _DataType o, const _DataType p ) {
    this->_row[0][0] = a;    this->_row[0][1] = b;    this->_row[0][2] = c;    this->_row[0][3] = d;
    this->_row[1][0] = e;    this->_row[1][1] = f;    this->_row[1][2] = g;    this->_row[1][3] = h;
    this->_row[2][0] = i;    this->_row[2][1] = j;    this->_row[2][2] = k;    this->_row[2][3] = l;
    this->_row[3][0] = m;    this->_row[3][1] = n;    this->_row[3][2] = o;    this->_row[3][3] = p;
  }

  Matrix44(const Mat<4, 4, _DataType> & rhs) : Mat<4, 4, _DataType>(rhs) {
  }

  /** compute the determinant of this matrix
   */
  _DataType det ( void ) const {
    // We Laplace-expand the determinant along the first row
    _DataType determinant = aol::NumberTrait<_DataType>::zero;
    aol::Matrix33<_DataType> dummyMatrix;

    dummyMatrix.fill ( this->_row[1][1], this->_row[1][2], this->_row[1][3],
                       this->_row[2][1], this->_row[2][2], this->_row[2][3],
                       this->_row[3][1], this->_row[3][2], this->_row[3][3] );
    determinant += this->_row[0][0] * dummyMatrix.det();

    dummyMatrix.fill ( this->_row[1][0], this->_row[1][2], this->_row[1][3],
                       this->_row[2][0], this->_row[2][2], this->_row[2][3],
                       this->_row[3][0], this->_row[3][2], this->_row[3][3] );
    determinant -= this->_row[0][1] * dummyMatrix.det();

    dummyMatrix.fill ( this->_row[1][0], this->_row[1][1], this->_row[1][3],
                       this->_row[2][0], this->_row[2][1], this->_row[2][3],
                       this->_row[3][0], this->_row[3][1], this->_row[3][3] );
    determinant += this->_row[0][2] * dummyMatrix.det();

    dummyMatrix.fill ( this->_row[1][0], this->_row[1][1], this->_row[1][2],
                       this->_row[2][0], this->_row[2][1], this->_row[2][2],
                       this->_row[3][0], this->_row[3][1], this->_row[3][2] );
    determinant -= this->_row[0][3] * dummyMatrix.det();

    return determinant;
  }
};


template <int numRows, int numCols, class T>
inline ostream &operator<< ( ostream &os, const Mat<numRows, numCols, T> &m ) {
  return m.print ( os );
}

template <int numRows, int numCols, class T>
inline istream &operator>> ( istream &is, Mat<numRows, numCols, T> &m ) {
  return m.read ( is );
}

template <class _DataType>
class Matrix33Symm {
public:
  typedef _DataType DataType;

  //! Constructor
  Matrix33Symm()
      : _row ( 3 ) { }

  //! Copy-constructor
  Matrix33Symm ( const Matrix33Symm<_DataType> &rhs )
      : _row ( 3 ) {
    _row[0] = rhs._row[0];
    _row[1] = rhs._row[1];
    _row[2] = rhs._row[2];
  }

  //! Constructor (for return value optimization)
  //! fills this with:
  //! [ a b c ]
  //! [ d e f ]
  //! [ g h i ]
  Matrix33Symm ( const _DataType a, const _DataType b, const _DataType c,
                 const _DataType d, const _DataType e, const _DataType f,
                 const _DataType g, const _DataType h, const _DataType i )
      : _row ( 3 ) {

    if ( b != d || c != g || f != h ) throw aol::InconsistentDataException ( "Matrix33Symm::Matrix33Symm : must be symmetric", __FILE__, __LINE__ );

    _row[0][0] = a;
    _row[0][1] = b;
    _row[0][2] = c;
    _row[1][0] = d;
    _row[1][1] = e;
    _row[1][2] = f;
    _row[2][0] = g;
    _row[2][1] = h;
    _row[2][2] = i;
  }

  //! Constructor (for return value optimization)
  //! fills this with:
  //! [ a b c ]
  //! [ b d e ]
  //! [ c e f ]
  Matrix33Symm ( const _DataType a, const _DataType b, const _DataType c,
                 const _DataType d, const _DataType e, const _DataType f )
      : _row ( 3 ) {

    _row[0][0] = a;
    _row[0][1] = b;
    _row[0][2] = c;
    _row[1][0] = b;
    _row[1][1] = d;
    _row[1][2] = e;
    _row[2][0] = c;
    _row[2][1] = e;
    _row[2][2] = f;
  }

  //! operator=
  Matrix33Symm<_DataType>& operator= ( const Matrix33Symm<_DataType> &rhs ) {
    _row[0] = rhs._row[0];
    _row[1] = rhs._row[1];
    _row[2] = rhs._row[2];
    return *this;
  }


  //! Calulates the Eigenvalues and Eigenvectors of the matrix using Jacobi-iteration.
  //! eigenVals[i] corresponds to the i-th column of eigenVecs
  void eigenVectors ( Vec3<_DataType> &eigenVals, Matrix33<_DataType> &eigenVecs );


  ostream& print ( ostream& out ) const {
    out << _row[0] << endl;
    out << _row[1] << endl;
    out << _row[2];
    return out;
  }

  const Vec3<_DataType>& operator[] ( const int i ) const {
    return _row[i];
  }

  _DataType get ( const int i, const int j ) const {
      return ( _row[i] ) [j];
    }
  void     set ( const int i, const int j, const _DataType Value ) {
    ( _row[i] ) [j] = Value;
    ( _row[j] ) [i] = Value;
  }

  void add ( const int I, const int J, const _DataType Value ) {
    _row[I][J] += Value;
    _row[J][I] += Value;
  }

  void setZero() {
    _row[0].setZero();
    _row[1].setZero();
    _row[2].setZero();
  }

  template<class T>
  void setIdentity() {
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        ( _row[i] ) [j] = static_cast<T> ( i == j );
      }
    }
  }

  void setRow ( const int index, const _DataType v1, _DataType v2, _DataType v3 ) {
    ( _row[index] ) [0] = v1;
    ( _row[index] ) [1] = v2;
    ( _row[index] ) [2] = v3;
    ( _row[0] ) [index] = v1;
    ( _row[1] ) [index] = v2;
    ( _row[2] ) [index] = v3;
  }

  template <class T>
  void setRow ( const int /*index*/, const Vec3<T> &vec ) {
    ( _row[this->Index] ) [0] = static_cast<T> ( vec.get ( 0 ) );
    ( _row[0] ) [this->Index] = static_cast<T> ( vec.get ( 0 ) );

    ( _row[this->Index] ) [1] = static_cast<T> ( vec.get ( 1 ) );
    ( _row[1] ) [this->Index] = static_cast<T> ( vec.get ( 1 ) );

    ( _row[this->Index] ) [2] = static_cast<T> ( vec.get ( 2 ) );
    ( _row[2] ) [this->Index] = static_cast<T> ( vec.get ( 2 ) );
  }

  template <class T>
  void getRow ( const int index, Vec3<T> &dest ) const {
    dest[0] = static_cast<T> ( ( _row ) [index][0] );
    dest[1] = static_cast<T> ( ( _row ) [this->Index][1] );
    dest[2] = static_cast<T> ( ( _row ) [this->Index][2] );
  }

  template <class T>
  void getColumn ( const int /*index*/, Vec3<T> &dest ) const {
    dest[0] = static_cast<T> ( ( _row ) [0][this->Index] );
    dest[1] = static_cast<T> ( ( _row ) [1][this->Index] );
    dest[2] = static_cast<T> ( ( _row ) [2][this->Index] );
  }


  //! \f$ x \mapsto Ax \f$
  //
  template <class T>
  void mult ( const Vec3<T> &Arg, Vec3<T> &Dest ) const {
    Dest[0] = _row[0][0] * Arg.get ( 0 ) + _row[0][1] * Arg.get ( 1 ) + _row[0][2] * Arg.get ( 2 );
    Dest[1] = _row[1][0] * Arg.get ( 0 ) + _row[1][1] * Arg.get ( 1 ) + _row[1][2] * Arg.get ( 2 );
    Dest[2] = _row[2][0] * Arg.get ( 0 ) + _row[2][1] * Arg.get ( 1 ) + _row[2][2] * Arg.get ( 2 );
  }

  //! computes the determinant of this matrix.
  //
  _DataType det() const {
    return
      _row[0][0]*_row[1][1]*_row[2][2] +
      _row[0][1]*_row[1][2]*_row[2][0] +
      _row[0][2]*_row[1][0]*_row[2][1]
      - _row[0][2]*_row[1][1]*_row[2][0]
      - _row[0][1]*_row[1][0]*_row[2][2]
      - _row[0][0]*_row[1][2]*_row[2][1];
  }

  //! computes the trace of this matrix.
  //
  _DataType tr() const {
    return _row[0][0] + _row[1][1] + _row[2][2];
  }

  //! \f$ A \mapsto A + B \f$
  //
  Matrix33Symm<_DataType> &operator+= ( const Matrix33Symm<_DataType> &Other ) {
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        _row[i][j] += Other._row[i][j];
      }
    }
    return *this;
  }

  //! \f$ A \mapsto A - B \f$
  //
  Matrix33Symm<_DataType> &operator-= ( const Matrix33Symm<_DataType> &Other ) {
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        _row[i][j] -= Other._row[i][j];
      }
    }
    return *this;
  }

  //! \f$ A\mapsto A*B \f$
  //
  Matrix33Symm<_DataType> &operator*= ( const Matrix33Symm<_DataType> &Other ) {
    Matrix33Symm<_DataType> tmp;
    tmp = *this;
    setZero();
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        for ( int k = 0; k < 3; ++k ) {
          _row[i][j] += tmp._row[i][k] * Other._row[k][j];
        }
      }
    }
    return *this;
  }

  template< class T >
  void makeProduct ( const Matrix33Symm<T> &Mat1, const Matrix33Symm<T> &Mat2 ) {
    setZero();
    for ( int i = 0;i < 3; ++i ) for ( int j = 0;j < 3; ++j ) for ( int k = 0;k < 3; ++k ) {
          this->_v[i][j] += Mat1.get ( i, k ) * Mat2.get ( k, j );
        }
  }

  void makeProduct ( const Matrix33Symm<_DataType> &A, const Matrix33Symm<_DataType> &B ) {
    setZero();
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        for ( int k = 0; k < 3; ++k ) {
          _row[i][j] += A._row[i][k] * B._row[k][j];
        }
      }
    }
  }


  //! \f$ A \mapsto \alpha A \f$
  //
  Matrix33Symm<_DataType> &operator*= ( _DataType Alpha ) {
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        _row[i][j] *= Alpha;
    return *this;
  }


  //! \f$ A \mapsto \alpha^{-1} A \f$
  //
  Matrix33Symm<_DataType> &operator/= ( _DataType Alpha ) {
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        _row[i][j] /= Alpha;
    return *this;
  }

  //! \f$ A \mapsto || A ||_{\infty} \f$
  //
  _DataType infinityNorm() {
    _DataType res = 0;
    for ( int i = 0; i < 3; ++i ) {
      _DataType temp = 0;
      for ( int j = 0; j < 3; ++j ) {
        temp += Abs ( _row[i][j] );
      }
      if ( temp > res ) res = temp;
    }
    return ( res );
  }


  //! \f$ A \mapsto A^T \f$
  //
  void transpose() {
    // Do nothing
  }

  //! \f$ Solves Ax = b \f$
  void applyInverseTo( const Vec3<_DataType> &b, Vec3<_DataType> &x  ) const {
    _DataType determinante = det();
    if ( determinante == _DataType ( 0 ) )
      throw ( Exception ( "matrix not invertible.", __FILE__, __LINE__ ) );

    Matrix33Symm<_DataType> AInverse( _row[1][1]*_row[2][2]  - _row[2][1]*_row[2][1], _row[0][2]*_row[1][2] - _row[0][1]*_row[2][2],
                                      _row[0][1]*_row[1][2]  - _row[1][1]*_row[0][2], _row[0][0]*_row[2][2] - _row[0][2]*_row[0][2],
                                      _row[0][1]*_row[0][2]  - _row[0][0]*_row[1][2], _row[0][0]*_row[1][1] - _row[0][1]*_row[0][1]);
    AInverse.mult( b, x );
    x /= determinante;
  }

  //! \f$ Solves Ax = b \f$
  void applyInverseTo( const Vec3<_DataType> &b, aol::Vector<_DataType> &x  ) const {
    if ( x.size() != 3 )
      throw ( Exception ( "Wrong vector dimension! Should be 3!", __FILE__, __LINE__ ) );
    Vec3<_DataType> v;
    applyInverseTo( b, v );
    for( int i = 0; i < 3; i++ )
      x[i] = v[i];
  }

private:
  vector< Vec3<_DataType> >  _row;
};

template <typename RealType, int Dim> struct MatDimTrait {
  typedef Mat< Dim, Dim, RealType> MatType;
};

template <typename RealType> struct MatDimTrait<RealType, 2> {
  typedef Matrix22<RealType> MatType;
};

template <typename RealType> struct MatDimTrait<RealType, 3> {
  typedef Matrix33<RealType> MatType;
};

template <typename RealType> struct MatDimTrait<RealType, 4> {
  typedef Matrix44<RealType> MatType;
};

template <class T>
inline ostream &operator<< ( ostream &os, const Matrix33Symm<T> &m ) {
  return m.print ( os );
}

} // end namespace

#endif
