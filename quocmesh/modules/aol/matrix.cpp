#include <matrix.h>
#include <sparseMatrices.h>
#include <op.h>

#include <scalarArray.h>

template <class T>
qc::ScalarArray<unsigned char, qc::QC_2D>* aol::Matrix<T>::createNonZeroPattern ( const int zoomFactor, const int zeroGray, const int nonZeroGray ) const {
  typedef qc::ScalarArray<unsigned char, qc::QC_2D> IMG_TYPE;

  const int columns = this->getNumCols();
  const int rows = this->getNumRows();
  const int dimX = columns * zoomFactor;
  const int dimY = rows * zoomFactor;

  IMG_TYPE *img = new IMG_TYPE ( dimX, dimY );

  img->setAll ( zeroGray );
  for ( int l = 0; l < columns; ++l ) {
    for ( int i = 0; i < rows; ++i ) {
      if ( get ( i, l ) != 0 ) {
        const int x = l * zoomFactor;
        const int y = i * zoomFactor;
        img->set ( x, y, nonZeroGray );
      }
    }
  }
  return img;
}

template <class _DataType>
void aol::Matrix<_DataType>::writeNonZeroPattern ( const char *fileName, const int zoomFactor, const int zeroGray, const int nonZeroGray ) const {
  ofstream file ( fileName );
  qc::ScalarArray<unsigned char, qc::QC_2D> *img = createNonZeroPattern ( zoomFactor, zeroGray, nonZeroGray );
  img->save ( file, qc::PGM_UNSIGNED_CHAR_BINARY );
  delete img;
}


template <class _DataType>
void aol::Matrix<_DataType>::setRowToZero ( const int row ) {
  for ( int i = 0; i < this->_numCols; ++i ) {
    set ( row, i, 0 );
  }
}

template <class _DataType>
void aol::Matrix<_DataType>::setColToZero ( const int col ) {
  for ( int i = 0; i < this->_numRows; ++i ) {
    set ( i, col, 0 );
  }
}

template <class _DataType>
void aol::Matrix<_DataType>::setRowColToZero ( const int rowCol ) {
  setRowToZero ( rowCol );
  setColToZero ( rowCol );
}

template <class _DataType>
void aol::Matrix<_DataType>::applyAdd ( const aol::Vector<_DataType> &arg, aol::Vector<_DataType> &dest ) const {
  if ( arg.size()  != this->getNumCols() ||
       dest.size() != this->getNumRows() ) {
    cerr << "arg.size = " << arg.size() << " dest.size = " << dest.size() << endl;
    throw ( aol::Exception ( "Wrong size of Matrix, arg and dest.", __FILE__, __LINE__ ) );
  }
  cerr << "WARNING!!!!!! Generic, slow aol::Matrix<T>::applyAdd called!\n";
  if ( arg.size() != this->_numCols ) {
    cerr << "ERROR in aol::Matrix<T>::applyAdd: number of columns doesn't match aol::Vectorlength of src!" << endl;
  } else if ( dest.size() != this->_numRows ) {
    cerr << "ERROR in aol::Matrix<T>::applyAdd: number of rows doesn't match aol::Vectorlength of dst!" << endl;
  } else {
    int deleteFlag = 0;

    const aol::Vector<_DataType>* Src;

    if ( &arg == &dest ) {
      Src = new aol::Vector<_DataType> ( arg );
      deleteFlag = 1;
    } else {
      Src = &arg;
    }
    int i, j;
    for ( i = 0; i < this->_numRows; ++i ) {
      for ( j = 0 ; j < this->_numCols ; ++j ) {
        dest[i] += get ( i, j ) * Src->get ( j );
      }
    }
    if ( deleteFlag ) delete Src;
  }
}

template <class _DataType>
_DataType aol::Matrix<_DataType>::getMaxValue() const {
  DataType m = get ( 0, 0 ), v;
  for ( int i = 0; i < this->getNumRows(); ++i ) {
    vector< typename Row<DataType>::RowEntry > vec;
    this->makeRowEntries ( vec, i );
    for ( typename vector<typename Row<DataType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
      v = it->value;
      if ( m < v ) {
        m = v;
      }
    }
  }

  return m;
}


template <class _DataType>
_DataType aol::Matrix<_DataType>::getFrobeniusNormSqr() const {
  DataType frobeniusNormSqr = ZOTrait<DataType>::zero;

  vector< typename Row<DataType>::RowEntry > vec;
  for ( int i = 0; i < this->getNumRows(); ++i ) {
    this->makeRowEntries ( vec, i );
    for ( typename vector<typename Row<DataType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
      frobeniusNormSqr += Sqr ( it->value );
    }
  }

  return frobeniusNormSqr;
}


template <class _DataType>
void aol::Matrix<_DataType>::getColumn ( int j, aol::Vector<_DataType> &v ) const {
  if ( j >= this->_numCols  || j < 0 ) {
    cerr << "ERROR in aol::Matrix<T>::getColumn: index j not in range!" << endl;
  } else if ( v.size() != this->_numRows ) {
    cerr << "ERROR in aol::Matrix<T>::getColumn: aol::VectorLength not equal to number of rows of aol::Matrix!" << endl;
  } else {
    for ( int i = 0 ; i < this->_numRows ; ++i ) {
      v[i] = get ( i, j );
    }
  }
}

template <class _DataType>
void aol::Matrix<_DataType>::getColumn ( int c, aol::MultiVector<_DataType>& Column ) const {
  if ( c >= this->getNumCols () || c < 0 ) {
    throw aol::Exception ( "aol::Matrix<T>::getColumn: index c not in range!", __FILE__, __LINE__ );
  } else if ( Column.getTotalSize() != this->getNumRows() ) {
    throw aol::Exception ( "aol::Matrix<T>::getColumn: aol::MultiVector total length not equal to number of rows of aol::Matrix!", __FILE__, __LINE__ );
  }
  int number = 0;
  for ( int d = 0; d < Column.numComponents(); ++d ) {
    for ( int i = 0; i < Column[d].size(); ++i ) {
      Column[d][i] = get ( number, c );
      number++;
    }
  }
}

template <class _DataType>
void aol::Matrix<_DataType>::getRow ( int i, aol::Vector<_DataType> &v ) const {
  if ( i >= this->_numRows  || i < 0 ) {
    cerr << "ERROR in aol::Matrix<T>::getRow: index i not in range!" << endl;
  } else if ( v.size() != this->_numCols ) {
    cerr << "ERROR in aol::Matrix<T>::getRow: aol::VectorLength not equal to number of columns of aol::Matrix!" << endl;
  } else {
    for ( int j = 0 ; j < this->_numCols ; ++j ) {
      v[j] = get ( i, j );
    }
  }
}

template <class _DataType>
void aol::Matrix<_DataType>::setColumn ( int j, const aol::Vector<_DataType> &v ) {
  if ( j >= this->_numCols  || j < 0 ) {
    cerr << "ERROR in aol::Matrix<T>::setColumn: index j not in range!" << endl;
  } else if ( v.size() != this->_numRows ) {
    cerr << "ERROR in aol::Matrix<T>::setColumn: aol::VectorLength not equal to number of rows of aol::Matrix!" << endl;
  } else {
    for ( int i = 0 ; i < this->_numRows ; ++i ) {
      set ( i, j, v.get ( i ) );
    }
  }
}

template <class _DataType>
void aol::Matrix<_DataType>::setColumn ( int c, const aol::MultiVector<_DataType>& Column ) {
  if ( c >= this->getNumCols () || c < 0 ) {
    throw aol::Exception ( "aol::Matrix<T>::setColumn: index c not in range!", __FILE__, __LINE__ );
  } else if ( Column.getTotalSize() != this->getNumRows() ) {
    throw aol::Exception ( "aol::Matrix<T>::setColumn: aol::MultiVector total length not equal to number of rows of aol::Matrix!", __FILE__, __LINE__ );
  }
  int number = 0;
  for ( int d = 0; d < Column.numComponents(); ++d ) {
    for ( int i = 0; i < Column[d].size(); ++i ) {
      set ( number, c, Column[d][i] );
      number++;
    }
  }
}

template <class _DataType>
void aol::Matrix<_DataType>::addColumn ( int j, const Vector<_DataType> &v, const _DataType Factor ) {
  if ( j >= this->_numCols  || j < 0 ) {
    cerr << "ERROR in aol::Matrix<T>::addColumn: index j not in range!" << endl;
  } else if ( v.size() != this->_numRows ) {
    cerr << "ERROR in aol::Matrix<_DataType>::addColumn: aol::VectorLength not equal to number of rows of aol::Matrix!" << endl;
  } else {
    for ( int i = 0 ; i < this->_numRows ; ++i ) {
      add ( i, j, Factor*v.get ( i ) );
    }
  }
}

template <class _DataType>
void aol::Matrix<_DataType>::addColumn ( int c, const aol::MultiVector<_DataType>& Column, const _DataType Factor ) {
  if ( c >= this->getNumCols () || c < 0 ) {
    throw aol::Exception ( "aol::Matrix<_DataType>::addColumn: index c not in range!", __FILE__, __LINE__ );
  } else if ( Column.getTotalSize() != this->getNumRows() ) {
    throw aol::Exception ( "aol::Matrix<_DataType>::addColumn: aol::MultiVector total length not equal to number of rows of aol::Matrix!", __FILE__, __LINE__ );
  }
  int number = 0;
  for ( int d = 0; d < Column.numComponents(); ++d ) {
    for ( int i = 0; i < Column[d].size(); ++i ) {
      add ( number, c, Factor*Column[d][i] );
      number++;
    }
  }
}

template <class _DataType>
void aol::Matrix<_DataType>::setRow ( int i, const aol::Vector<_DataType> &v ) {
  if ( i >= this->_numRows  || i < 0 ) {
    cerr << "ERROR in aol::Matrix<_DataType>::getColumn: index i not in range!" << endl;
  } else if ( v.size() != this->_numCols ) {
    cerr << "ERROR in aol::Matrix<_DataType>::getColumn: aol::VectorLength not equal"
         << " to number of columns of aol::Matrix!" << endl;
  } else {
    for ( int j = 0 ; j < this->_numCols ; ++j ) {
      set ( i, j, v.get ( j ) );
    }
  }
}


//! untested
template <class _DataType>
void aol::Matrix<_DataType>::transposeTo ( aol::Matrix<_DataType> &Dest ) const {
  // attention: derived classes must set Dest to zero!
  if ( ( this->_numRows != Dest.getNumCols() ) || ( this->_numCols != Dest.getNumRows() ) ) {
    throw Exception ( "aol::Matrix<_DataType>::transposeTo: matrices incompatible", __FILE__, __LINE__ );
  } else {
    for ( int i = 0; i < this->getNumRows(); ++i ) {
      vector< typename Row<DataType>::RowEntry > vec;
      this->makeRowEntries ( vec, i );
      for ( typename vector<typename Row<DataType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
        Dest.set ( it->col, i, it->value );
      }
    }
  }
}


template <class _DataType, bool Parallelize>
aol::FullMatrix<_DataType, Parallelize>& aol::FullMatrix<_DataType, Parallelize>::operator*= ( const aol::FullMatrix<_DataType, Parallelize> &mat ) {

  if ( mat.getNumRows() != mat.getNumCols() )
    throw aol::Exception ( "Matrix<_DataType>::operator*= :  rows != columns for right side!", __FILE__, __LINE__ );
  if ( this->_numRows != this->_numCols )
    throw aol::Exception ( "Matrix<_DataType>::operator*= :  rows != columns for left side!", __FILE__, __LINE__ );
  if ( this->_numCols != mat.getNumRows() )
    throw aol::Exception ( "Matrix<_DataType>::operator*= :  dimensions don't match!", __FILE__, __LINE__ );
  if ( this == &mat )
    throw aol::Exception ( "Matrix<_DataType>::operator*= :  don't multiply with the same matrix!", __FILE__, __LINE__ );


  aol::Vector<_DataType> vec ( this->getNumCols() );
  for ( int i = 0; i < this->getNumRows(); ++i ) {
    this->getRow ( i, vec );
    for ( int j = 0; j < this->getNumCols(); ++j ) {
      _DataType val = 0;
      for ( int k = 0; k < this->getNumCols(); ++k ) {
        val += mat.get ( k, j ) * vec[k];
      }
      ref ( i, j ) = val;
    }
  }

  return *this;
}

template <class _DataType>
aol::Matrix<_DataType>& aol::Matrix<_DataType>::operator+= ( const aol::Matrix<_DataType> &other ) {
  if ( this == &other ) {
    throw aol::Exception ( "aol::Matrix<_DataType>::operator+=: Don't subtract matrix from itself", __FILE__, __LINE__ );
    // maybe this works and the exception is not necessary
  }

  if ( this->getNumRows() != other.getNumRows() || this->getNumCols() != other.getNumCols() ) {
    throw aol::Exception ( "aol::Matrix<_DataType>::operator+= : Dimensions don't match", __FILE__, __LINE__ );
  }

  for ( int i = 0; i < this->getNumRows(); ++i ) {
    std::vector<typename Row<_DataType>::RowEntry > vec;
    other.makeRowEntries ( vec, i );
    for ( typename vector<typename Row<_DataType>::RowEntry >::const_iterator it = vec.begin(); it != vec.end(); ++it ) {
      this->add ( i, it->col, it->value );
    }
  }

  return ( *this );
}

template <class _DataType>
aol::Matrix<_DataType>& aol::Matrix<_DataType>::operator-= ( const aol::Matrix<_DataType> &other ) {
  if ( this == &other ) {
    throw aol::Exception ( "aol::Matrix<_DataType>::operator-=: Don't subtract matrix from itself", __FILE__, __LINE__ );
    // maybe this works and the exception is not necessary
  }

  if ( this->getNumRows() != other.getNumRows() || this->getNumCols() != other.getNumCols() ) {
    throw aol::Exception ( "aol::Matrix<_DataType>::operator-= : Dimensions don't match", __FILE__, __LINE__ );
  }

  for ( int i = 0; i < this->getNumRows(); ++i ) {
    std::vector<typename Row<_DataType>::RowEntry > vec;
    other.makeRowEntries ( vec, i );
    for ( typename vector<typename Row<_DataType>::RowEntry >::const_iterator it = vec.begin(); it != vec.end(); ++it ) {
      this->add ( i, it->col, - ( it->value ) );
    }
  }

  return ( *this );
}

template <class _DataType>
aol::Matrix<_DataType>& aol::Matrix<_DataType>::operator*= ( const _DataType alpha ) {
  for ( int i = 0; i < this->getNumRows(); ++i ) {
    std::vector<typename Row<_DataType>::RowEntry > vec;
    this->makeRowEntries ( vec, i );
    for ( typename vector<typename Row<_DataType>::RowEntry >::const_iterator it = vec.begin(); it != vec.end(); ++it ) {
      this->set ( i, it->col, alpha * ( it->value ) );
    }
  }

  return ( *this );
}


template <class _DataType>
void aol::Matrix<_DataType>::addMultipleAtPosition ( int I, int J , const aol::Matrix<_DataType> &mat, _DataType factor ) {
  int i, j, iend, jend;

  iend = I + mat.getNumRows();
  if ( iend >= this->_numRows ) iend = this->_numRows;
  jend = J + mat.getNumCols();
  if ( jend >= this->_numCols ) jend = this->_numCols;

  if ( I < 0 ) I = 0;
  if ( J < 0 ) J = 0;

  for ( i = I; i < iend ; ++i )
    for ( j = J; j < jend ; ++j ) {
      _DataType a_ij ( mat.get ( i - I, j - J ) );
      if ( a_ij != ZOTrait<_DataType>::zero )
        set ( i, j, get ( i, j ) + factor * a_ij );
    }
}


template <class _DataType>
void aol::Matrix<_DataType>::setIdentity() {
  if ( this->_numRows == this->_numCols ) {
    this->setZero();
    for ( int i = 0 ; i < this->_numRows; ++i )
      set ( i, i, 1 );
  } else {
    throw aol::Exception ( "Cannot set nonquadratic matrix to identity.", __FILE__, __LINE__ );
  }
}

template <class _DataType>
void aol::Matrix<_DataType>::addIdentity() {
  if ( this->_numRows == this->_numCols ) {
    for ( int i = 0 ; i < this->_numRows; ++i )
      add ( i, i, 1 );
  } else {
    throw aol::Exception ( "Cannot add identity to a nonquadratic matrix.", __FILE__, __LINE__ );
  }
}

template <class _DataType>
_DataType aol::Matrix<_DataType>::rowSum ( const int I ) const {
  DataType s = ZOTrait<DataType>::zero;
  vector< typename Row<DataType>::RowEntry > vec;
  this->makeRowEntries ( vec, I );
  for ( typename vector<typename Row<DataType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
    s += it->value;
  }

  return s;
}

template <class _DataType>
void aol::Matrix<_DataType>::makeNonZeroRowEntries ( vector<typename Row<_DataType>::RowEntry > &vec, const int RowNum ) const {
  vector< typename aol::Row<DataType>::RowEntry > allrevec;
  this->makeRowEntries ( allrevec, RowNum );

  vec.reserve ( allrevec.size() ); // depending on whether we have many zeroes or not, this may be inefficient or not ...
  for ( typename std::vector< typename aol::Row<DataType>::RowEntry >::const_iterator it = allrevec.begin(); it != allrevec.end(); ++it ) {
    if ( it->value != aol::NumberTrait<DataType>::zero ) {
      vec.push_back ( *it );
    }
  }
}

template <class _DataType>
void aol::Matrix<_DataType>::addMatrixProduct ( const Matrix<_DataType> &M1, const Matrix<_DataType> &M2 ) {
  if ( this->getNumRows() != M1.getNumRows() || this->getNumCols() != M2.getNumCols() || M1.getNumCols() != M2.getNumRows() ) {
    throw Exception ( "aol::Matrix::addMatrixProduct: incompatible matrix sizes.\n", __FILE__, __LINE__ );
  }

  std::vector<typename Row<_DataType>::RowEntry > vec1, vec2;

  for ( int i = 0; i < this->getNumRows(); ++i ) {
    M1.makeRowEntries ( vec1, i );
    for ( typename vector<typename Row<_DataType>::RowEntry >::iterator it1 = vec1.begin(); it1 != vec1.end(); ++it1 ) {
      M2.makeRowEntries ( vec2, it1->col );
      for ( typename vector<typename Row<_DataType>::RowEntry >::iterator it2 = vec2.begin(); it2 != vec2.end(); ++it2 ) {
        add ( i, it2->col, it1->value * it2->value );
      }
    }
  }
}

template <class _DataType>
void aol::Matrix<_DataType>::resize ( const int, const int ) {
  throw aol::Exception ( "aol::Matrix<>::resize must be overloaded on derived classes where you want to use it.", __FILE__, __LINE__ );
}

template <class _DataType>
void aol::Matrix<_DataType>::reallocate ( const int, const int ) {
  throw aol::Exception ( "aol::Matrix<>::reallocate must be overloaded on derived classes where you want to use it.", __FILE__, __LINE__ );
}

template <class _DataType>
ostream& aol::Matrix<_DataType>::print ( ostream& os ) const {
  int m = this->getNumRows (), n = this->getNumCols ();
  int i, j;

  printHead ( os );

  os << format ( ' ' ) << ( prettyFormat ? "|" : " " );
  for ( j = 0 ; j < n ; ++j ) {
    os << format ( j ) << " ";
  }
  os << endl;

  if ( prettyFormat ) {
    os << format ( '-' ) << "+";
    for ( j = 0 ; j < n ; ++j ) {
      os << format ( '-' ) << "-";
    }
    os << endl;
  }

  for ( i = 0 ; i < m ; ++i ) {
    os << format ( i ) << ( prettyFormat ? "|" : " " );
    for ( j = 0 ; j < n ; ++j )
      os << format ( get ( i, j ) ) << " ";
    os << endl;
  }
  return os;
}

template <class _DataType>
ostream& aol::Matrix<_DataType>::printSparse ( ostream& os, const aol::Format DataFormatter ) const {
  this->printHead ( os );
  for ( int i = 0; i < this->_numRows; ++i ) {
    bool printedHeader = false;
    std::vector<typename Row<_DataType>::RowEntry > vec;
    this->makeRowEntries ( vec, i );
    for ( typename vector<typename Row<_DataType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
      if( it->value != 0. ){
	if( !printedHeader ){
	  os << "Row " << aol::intFormat ( i ) << ": ";
	  printedHeader = true;
	}
        os << "(" << aol::intFormat ( i ) << ", " << aol::intFormat ( it->col ) << "): " << DataFormatter ( it->value ) << "   ";
      }
    }
    if( printedHeader )
      os << endl;
  }
  return os;
}

template <class _DataType>
ostream& aol::Matrix<_DataType>::printSparseOctave ( ostream& os ) const {
  os << "A = sparse(" << this->_numRows << "," <<  this->_numCols << ");" << endl;
  for ( int i = 0; i < this->_numRows; ++i ) {
    std::vector<typename Row<_DataType>::RowEntry > vec;
    this->makeRowEntries ( vec, i );
    for ( typename vector<typename Row<_DataType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
      os << "A (" << i + 1 << ", " << it->col + 1 << ") = " << aol::longScientificFormat ( it->value ) << ";" << endl;
    }
    os << endl;
  }
  return os;
}

template <class _DataType>
void aol::Matrix<_DataType>::printHead ( ostream& os ) const {
  os << "Matrix " << this->getNumRows () << " x " << this->getNumCols () << endl;
}

template <class _DataType>
void aol::Matrix<_DataType>::readHead ( istream& /*is*/ ) {
  throw UnimplementedCodeException ( "Cannot read abstract Matrix.", __FILE__, __LINE__ );
}

template <class _DataType, bool Parallelize>
void aol::FullMatrix<_DataType, Parallelize>::printHead ( ostream& os ) const {
  os << "FullMatrix " << this->getNumRows () << " x " << this->getNumCols () << endl;
}

template <class _DataType, bool Parallelize>
void aol::FullMatrix<_DataType, Parallelize>::readHead ( istream& is ) {
  string s;
  char c = 0;
  int m, n;

  is >> s >> m >> c >> n;
  if ( s != "FullMatrix" || c != 'x' ) throw FileException ( "File aol::Format Mismatch: What is the Matrix?", __FILE__, __LINE__ );

  reallocate ( m, n );
}

template <class _DataType, bool Parallelize>
bool aol::FullMatrix<_DataType, Parallelize>::isSelfTestOk () {
  const int dim = 6;
  int i, j;

  cerr << "--- Testing FullMatrix<double> ... \n";

  FullMatrix<double> m ( dim, dim );
  for ( i = 0; i < dim; ++i )
    for ( j = 0; j < dim; ++j )
      m.set ( i, j, pow ( static_cast <double> ( i ), static_cast <double> ( j ) ) );


  FullMatrix<double> mi ( dim, dim );
  PermutationMatrix<double> p ( dim );

  mi.makeLU ( m, p );

  Vector<double> b ( dim );
  Vector<double> x ( dim );
  for ( i = 0; i < dim; ++i )
    b [i] = sin ( static_cast <double> ( i ) );
  mi.LUSolve ( x, b, p );

  double tmp[] = {0.0, 0.885486, 0.267195, -0.391936, 0.086137, -0.005411};
  // The tmp pointer is not 16 byte aligned, so we can't use it directly in the Vector in case USE_SSE is defined.
  Vector<double> correct ( 6 );
  correct.readFromBuffer ( tmp );

  correct -= x;
  cerr << "Difference is " << correct.norm () << endl;

  // test transposition of non-square matrices
  FullMatrix<double> m2 ( 2, 4 ), m2t ( 4, 2 );
  m2.set ( 1, 2, 3.14 );
  m2t.set ( 2, 1, 4.13 );
  m2.transposeTo ( m2t );

  if ( ( correct.norm () < 1E-4 ) && ( m2.ref ( 1, 2 ) - m2t.ref ( 2, 1 ) < 1e-4 ) ) {
    cerr << "Test of FullMatrix<double> ................................................... OK.\n";
    return true;
  }
  cerr << "Test of FullMatrix<double> ............................................... FAILED.\n";
  return false;
}

namespace {
void eatsep ( istream& is ) {
  char c = 0;
  is >> c;
  if ( c != '|' && c != '-' && c != '+' && c != ':' && c != ';' && c != '/' && c != '(' && c != ')' )
    is.putback ( c );
}
}

template <class _DataType>
istream& aol::Matrix<_DataType>::read ( istream& is ) {
  string s;
  char c = 0;
  int i, j;
  int index;
  _DataType value;

  // Read header
  readHead ( is );

  // Column numbers
  for ( j = 0 ; j < this->getNumCols () ; ++j ) {
    eatsep ( is );
    is >> index;
    if ( index != j ) throw aol::FileException ( "File aol::Format Mismatch: Column Indices Incorrect", __FILE__, __LINE__ );
  }

  // Skip seperator line
  is >> c;
  if ( c == '-' ) is >> s;
  else is.putback ( c );

  // Read data
  for ( i = 0 ; i < this->getNumRows () ; ++i ) {
    eatsep ( is );
    is >> index;
    if ( index != i ) throw aol::FileException ( "File aol::Format Mismatch: Row Index Incorrect", __FILE__, __LINE__ );
    for ( j = 0 ; j < this->getNumCols () ; ++j ) {
      eatsep ( is );
      is >> value;
      set ( i, j, value );
    }
  }
  return is;
}

template <class _DataType, bool Parallelize>
aol::FullMatrix<_DataType, Parallelize>::FullMatrix ()
    : Matrix<_DataType> ( 0, 0 ) {
  data.reallocate ( 0 );
}

template <class _DataType, bool Parallelize>
aol::FullMatrix<_DataType, Parallelize>::FullMatrix ( int m, int n )
    : Matrix<_DataType> ( m, n ) {
  data.reallocate ( this->_numRows * this->_numCols );
}

template <class _DataType, bool Parallelize>
aol::FullMatrix<_DataType, Parallelize>::FullMatrix ( const aol::Vector<_DataType>& V, bool Transpose )
    : Matrix<_DataType> ( Transpose ? 1 : V.size (), Transpose ? V.size () : 1 ) {
  data.reallocate ( this->_numRows * this->_numCols );

  if ( Transpose )
    for ( int i = 0; i < V.size (); ++i )
      ref ( 0, i ) = V.get ( i );
  else
    for ( int i = 0; i < V.size (); ++i )
      ref ( i, 0 ) = V.get ( i );
}

template <class _DataType, bool Parallelize>
aol::FullMatrix<_DataType, Parallelize>::FullMatrix ( const aol::Vector<_DataType>& Vec, const int m, const int n, CopyFlag copyFlag )
  : Matrix<_DataType> ( m, n ), data ( Vec, copyFlag ) {
  if ( Vec.size() != m * n ) {
    cerr << Vec.size() << " " << m << " " << n << endl;
    throw aol::Exception ( "aol::FullMatrix ( Vec, m, n ): incorrect size", __FILE__, __LINE__ );
  }
}



template <class _DataType, bool Parallelize>
aol::FullMatrix<_DataType, Parallelize>::FullMatrix ( const aol::Op<aol::Vector<_DataType> >& op, int m, int n )
    : Matrix<_DataType> ( m, n ) {
  data.reallocate ( m * n );
  aol::Vector<_DataType> ei ( n ), col ( m );
  for ( int i = 0; i < n; ++i ) {
    ei [i] = 1;
    if ( i != 0 ) ei [i-1] = 0;
    op.apply ( ei, col );
    this->setSubColumn ( 0, i, col );
  }
}

template <class _DataType, bool Parallelize>
aol::FullMatrix<_DataType, Parallelize>::FullMatrix ( const aol::Matrix<_DataType>& m )
    : Matrix<_DataType> ( m ) { 
  data.reallocate ( this->_numRows * this->_numCols );

  for ( int i = 0; i < this->_numRows; ++i )
    for ( int j = 0; j < this->_numCols; ++j )
      ref ( i, j ) =  m.get ( i, j );
}

template <class _DataType, bool Parallelize>
aol::FullMatrix<_DataType, Parallelize>::FullMatrix ( const aol::FullMatrix<_DataType, Parallelize>& m, CopyFlag copyFlag )
    : Matrix<_DataType> ( m ) {
  if( copyFlag != aol::DEEP_COPY )
    throw aol::Exception ( "aol::FullMatrix<_DataType, Parallelize>::FullMatrix: unimplemented copyFlag!", __FILE__, __LINE__ );
    
  data.reallocate ( this->_numRows * this->_numCols );

  for ( int i = 0; i < this->_numRows; ++i )
    for ( int j = 0; j < this->_numCols; ++j )
      ref ( i, j ) =  m.ref ( i, j );
}

template <class _DataType, bool Parallelize>
void aol::FullMatrix<_DataType, Parallelize>::setZero ( ) {
  int i, j;
  for ( i = 0 ; i < this->_numRows ; ++i ) {
    for ( j = 0 ; j < this->_numCols ; ++j ) {
      ref ( i, j ) = 0;
    }
  }
}

template <class _DataType, bool Parallelize>
aol::FullMatrix<_DataType, Parallelize>& aol::FullMatrix<_DataType, Parallelize>::operator*= ( const _DataType alpha ) {
  int i, j;
  for ( i = 0 ; i < this->_numRows ; ++i ) {
    for ( j = 0 ; j < this->_numCols ; ++j ) {
      ref ( i, j ) *= alpha;
    }
  }
  return *this;
}

template <class _DataType, bool Parallelize>
aol::FullMatrix<_DataType, Parallelize>& aol::FullMatrix<_DataType, Parallelize>::operator= ( const aol::Matrix<_DataType> &mat ) {
  int i, j;
#ifdef BOUNDS_CHECK
  if ( ( mat.getNumRows() != this->_numRows ) || ( mat.getNumCols() != this->_numCols ) ) {
    cerr << "ERROR in aol::FullMatrix<_DataType, Parallelize>::operator=() : check dimensions.." << endl;
    return *this;
  }
#endif
  for ( i = 0 ; i < this->_numRows ; ++i ) {
    for ( j = 0 ; j < this->_numCols ; ++j ) {
      ref ( i, j ) = mat.get ( i, j );
    }
  }
  return *this;
}

template <class _DataType, bool Parallelize>
aol::FullMatrix<_DataType, Parallelize>& aol::FullMatrix<_DataType, Parallelize>::operator= ( const aol::FullMatrix<_DataType, Parallelize> &mat ) {
#ifdef BOUNDS_CHECK
  if ( ( mat.getNumRows() != this->_numRows ) || ( mat.getNumCols() != this->_numCols ) ) {
    cerr << "ERROR in aol::FullMatrix<_DataType, Parallelize>::operator=() : check dimensions.." << endl;
    return *this;
  }
#endif
  data = mat.data;
  return *this;
}

template <class _DataType, bool Parallelize>
aol::FullMatrix<_DataType, Parallelize>& aol::FullMatrix<_DataType, Parallelize>::operator+= ( const aol::FullMatrix<_DataType, Parallelize> &mat ) {
#ifdef BOUNDS_CHECK
  if ( ( mat.getNumRows() != this->_numRows ) || ( mat.getNumCols() != this->_numCols ) ) {
    cerr << "ERROR in aol::FullMatrix<_DataType, Parallelize>::operator+=() : check dimensions.." << endl;
    return *this;
  }
#endif
  for ( int i = 0 ; i < this->_numRows ; ++i ) {
    for ( int j = 0 ; j < this->_numCols ; ++j ) {
      ref ( i, j ) += mat.ref ( i, j );
    }
  }
  return *this;
}

template <class _DataType, bool Parallelize>
aol::FullMatrix<_DataType, Parallelize>& aol::FullMatrix<_DataType, Parallelize>::operator-= ( const aol::FullMatrix<_DataType, Parallelize> &mat ) {
#ifdef BOUNDS_CHECK
  if ( ( mat.getNumRows() != this->_numRows ) || ( mat.getNumCols() != this->_numCols ) ) {
    cerr << "ERROR in aol::FullMatrix<_DataType, Parallelize>::operator+=() : check dimensions.." << endl;
    return *this;
  }
#endif
  for ( int i = 0 ; i < this->_numRows ; ++i ) {
    for ( int j = 0 ; j < this->_numCols ; ++j ) {
      ref ( i, j ) -= mat.ref ( i, j );
    }
  }
  return *this;
}

template <class _DataType, bool Parallelize>
void aol::FullMatrix<_DataType, Parallelize>::applyAddTranspose ( const Vector<_DataType> &src, Vector<_DataType> &dst ) const {
  for ( int i = 0; i < static_cast<int> ( this->_numRows ); ++i ) {
    for ( int j = 0; j < static_cast<int> ( this->_numCols ); ++j ) {
      dst[ j ] += ref ( i, j ) * src[ i ];
    }
  }
}

template <class _DataType, bool Parallelize>
void aol::FullMatrix<_DataType, Parallelize>::makeInverse ( const aol::Matrix<_DataType> &mat ) {
#ifdef BOUNDS_CHECK
  if ( ( mat.getNumRows() != this->_numRows ) || ( mat.getNumCols() != this->_numCols ) ) {
    cerr << "ERROR in aol::FullMatrix<_DataType, Parallelize>::makeInverse : check dimensions.." << endl;
    return;
  }
  if ( this->_numRows != this->_numCols ) {
    cerr << "ERROR in aol::FullMatrix<_DataType, Parallelize>::makeInverse : only rectangular matrices can be inverted!" << endl;
    return;
  }
#endif
  aol::FullMatrix<_DataType, Parallelize> A ( this->_numRows, this->_numCols );
  A = mat;

  this->setIdentity();
  _DataType pivot;
  int i, j;

  for ( j = 0; j < this->_numCols; ++j ) {
    i = j;
    while ( ( i < this->_numRows ) && ( A.ref ( i, j ) == 0.0 ) ) ++i;
    if ( i == this->_numRows ) {
      cerr << "ERROR in gsMatrix::makeInverse: singular matrix...\n";
      return;
    }
    if ( i != j ) {
      swapRows ( i, j );
      A.swapRows ( i, j );
    }
    pivot = A.ref ( j, j );

    A.multRow ( j, ( static_cast< _DataType > ( 1.0 ) ) / pivot, j );
    multRow ( j, ( static_cast< _DataType > ( 1.0 )  ) / pivot, 0 );

    for ( i = 0; i < this->_numRows; ++i ) {
      if ( i != j ) {
        pivot = ( static_cast< _DataType > ( - 1.0 ) ) * A.ref ( i, j );
        addToRow ( i, j, pivot, 0 );
        A.addToRow ( i, j, pivot, j );
      }
    }
  }
}

template <class _DataType, bool Parallelize>
void aol::FullMatrix<_DataType, Parallelize>::makeLeftInverse ( aol::FullMatrix<_DataType, Parallelize>& A_1Inv,
                                                   aol::FullMatrix<_DataType, Parallelize>& A_2,
                                                   aol::PermutationMatrix<_DataType>& P ) {
#ifdef BOUNDS_CHECK
  if ( this->_numRows >= this->_numCols ) {
    cerr << "ERROR in aol::FullMatrix<_DataType, Parallelize>::makeLeftInverse: doesn't make sense for M>=N!\n";
    return;
  }
  if ( ( A_1Inv.getNumRows() != this->_numRows ) || ( A_1Inv.getNumCols() != this->_numRows )
       || ( A_2.getNumRows() != this->_numRows ) || ( A_2.getNumCols() != this->_numCols - this->_numRows ) ) {
    cerr << "ERROR in aol::FullMatrix<_DataType, Parallelize>::makeLeftInverse : check dimensions.." << endl;
    return;
  }
#endif
  int i, j;
  /* unsigned */
  int SwapCol;
  aol::FullMatrix<_DataType, Parallelize> A ( this->_numRows, this->_numCols );
  aol::FullMatrix<_DataType, Parallelize> Swapped ( this->_numRows, this->_numCols );

  // FIX ERROR in operator=
  A = *this;
  Swapped = *this;

  A_1Inv.setIdentity();
  _DataType pivot;

  SwapCol = - 1;
  for ( j = 0; j < this->_numRows; ++j ) {
    do {
      i = j;
      SwapCol++;
      while ( ( i < this->_numRows ) && ( A.ref ( i, SwapCol ) == 0.0 ) ) ++i;
    } while ( i == this->_numRows && SwapCol <  static_cast< int > ( this->_numCols - 1 ) );

    if ( SwapCol == static_cast<int> ( this->_numCols - 1 ) && i == this->_numRows ) {
      cerr << "ERROR in gsMatrix::makeLeftInverse: rank of matrix < M \n";
      return;
    }

    if ( SwapCol != static_cast<int> ( j ) ) {
      P.makeSwap ( j, SwapCol );
      A.swapColumns ( j, SwapCol );
      Swapped.swapColumns ( j, SwapCol );
    }

    if ( i != j ) {
      A_1Inv.swapRows ( i, j );
      A.swapRows ( i, j );
    }
    pivot = A.ref ( j, j );

    A.multRow ( j, ( static_cast< _DataType > ( 1.0 ) ) / pivot, j );
    A_1Inv.multRow ( j, ( static_cast< _DataType > ( 1.0 )  ) / pivot, 0 );

    for ( i = 0; i < this->_numRows; ++i ) {
      if ( i != j ) {
        pivot = ( static_cast< _DataType > ( - 1.0 )  ) * A.ref ( i, j );
        A_1Inv.addToRow ( i, j, pivot, 0 );
        A.addToRow ( i, j, pivot, j );
      }
    }
  }

  // now set A_2
  for ( i = 0; i < this->_numRows; ++i )
    for ( j = 0; j < this->_numCols - this->_numRows; ++j )
      A_2.ref ( i, j ) = Swapped.ref ( i, this->_numRows + j );
}

template <class _DataType, bool Parallelize>
void aol::FullMatrix<_DataType, Parallelize>::makeSubstitutionMatrix ( const aol::Vector<_DataType> &b,
                                                          aol::Vector<_DataType> &q,
                                                          aol::FullMatrix<_DataType, Parallelize> &/*F*/,
                                                          aol::PermutationMatrix<_DataType>& P ) {
  int i;
  aol::FullMatrix<_DataType, Parallelize> A_1Inv ( this->_numRows, this->_numRows );
  aol::FullMatrix<_DataType, Parallelize> A_2 ( this->_numRows, ( this->_numCols - this->_numRows ) );
  aol::FullMatrix<_DataType, Parallelize> Upper ( this->_numRows, this->_numCols - this->_numRows );

  cerr << " S undeclared\n";
  //  S.setZero( );

  makeLeftInverse ( A_1Inv, A_2, P );

  Upper = A_1Inv;
  Upper *= static_cast< _DataType > ( -1.0 ) ;
  Upper *= A_2;

  //  S.addToMatrix( 0, 0, Upper );
  cerr << "S undeclared\n";

  /* for ( i = 0; i < N - M; ++i ) {
     S.data[ i + M ][ i ] = 1.0f;
   }*/

  aol::Vector<_DataType> v ( this->_numRows );
  A_1Inv.mult ( b, v );

  q.setZero();

  for ( i = 0; i < this->_numRows; ++i ) {
    q[i] = v[i];
  }
}


template <class _DataType, bool Parallelize>
void aol::FullMatrix<_DataType, Parallelize>::makeLU ( const aol::Matrix<_DataType> &mat, aol::PermutationMatrix<_DataType>& P ) {
#ifdef BOUNDS_CHECK
  if ( mat.getNumRows() != this->_numRows || mat.getNumCols() != this->_numCols || this->_numRows != this->_numCols ) {
    cerr << "ERROR in aol::FullMatrix<_DataType, Parallelize>::makeLU: mat has not the same dimensions as this matrix or is not quadratic! \n";
    return;
  }
  if ( P.getNumRows() != this->_numRows ) {
    cerr << "ERROR in aol::FullMatrix<_DataType, Parallelize>::makeLU: PermutationMatrix P has not the correction dimension!\n";
    return;
  }
#endif
  ( *this ) = mat;
  int i, j, k, pivotRow = 0;
  _DataType pivot, sum, max, rel, l;

  for ( i = 0; i < this->_numRows - 1; ++i ) {
    // first step: search pivot element;
    max = 0;
    for ( k = i; k < this->_numRows; ++k ) {
      sum = 0;
      for ( j = i; j < this->_numCols; ++j ) {
        sum += aol::Abs ( ref ( i, j ) );
      }
      rel = aol::Abs ( ref ( k, i ) ) / sum;
      if ( rel > max ) {
        max = rel;
        pivotRow = k;
      }
    }
    if ( max == 0 ) {
      cerr << "makeLU: aol::Matrix not regular. returning...\n";
      return;
    }
    if ( i != pivotRow ) {
      swapRows ( i, pivotRow );
      P.makeSwap ( i, pivotRow );
    }

    pivot = ref ( i, i );
    // second step: create upper right, and store lower left elements
    for ( k = i + 1; k < this->_numRows; ++k ) {
      l = ref ( k, i ) / pivot;
      addToRow ( k, i, ( static_cast<_DataType> ( -1.0 ) ) *l, i );
      ref ( k, i ) = l;
    }
  }
}


template <class _DataType, bool Parallelize>
void aol::FullMatrix<_DataType, Parallelize>::LUSolve ( aol::Vector<_DataType>& X, const aol::Vector<_DataType> &RHS, const aol::PermutationMatrix<_DataType>& P ) const {
#ifdef BOUNDS_CHECK
  if ( X.size() != this->_numRows ) {
    cerr << "ERROR in LUSolve: Length of X compatible not to this matrix... \n";
    return;
  }
  if ( RHS.size() != this->_numRows ) {
    cerr << "ERROR in LUSolve: Length of RightHandSide not compatible to this matrix... \n";
    return;
  }
  if ( P.getNumRows() != this->_numRows ) {
    cerr << "ERROR in LUSolve: Dimension of P not compatible to this matrix... \n";
    return;
  }
#endif
  aol::Vector<_DataType> rhs ( this->_numRows ), Y ( this->_numRows );

  // multiply right hand side with permutation matrix
  P.apply ( RHS, rhs );

  int i, j; // May not be unsigned due to >=0 comparison below
  _DataType val;
  // first step: solve L*y = rhs

  for ( i = 0; i < static_cast<int> ( this->_numRows ); ++i ) {
    val = rhs[ i ];
    for ( j = 0; j < i; ++j ) {
      val -= Y[j] * ref ( i, j );
    }
    Y[i] = val;
  }

  // second step: solve R*x = y
  // wouldn't for ( i = M; i-- != 0; ) { be cool (and work with unsigned, too)?
  for ( i = this->_numRows - 1; i >= 0; --i ) {
    val = Y[i];
    for ( j = i + 1; j < static_cast<int> ( this->_numCols ); ++j ) {
      val -= X[j] * ref ( i, j );
    }
    X[i] = val / ref ( i, i );
  }
}

template <class _DataType, bool Parallelize>
inline void aol::FullMatrix<_DataType, Parallelize>::shiftRowsUp() {
  // Unused, untested
  cerr << "Shifting Rows is slow!" << endl;
  for ( int i = 0; i < this->_numRows; ++i )
    swapRows ( i, ( i + 1 ) % this->_numRows );
}

template <class _DataType, bool Parallelize>
inline void aol::FullMatrix<_DataType, Parallelize>::swapRows ( int R1, int R2 ) {
  for ( int j = 0; j < this->_numCols; ++j ) {
    _DataType tmp = ref ( R1, j );
    ref ( R1, j ) = ref ( R2, j );
    ref ( R2, j ) = tmp;
  }
}

template <class _DataType, bool Parallelize>
inline void aol::FullMatrix<_DataType, Parallelize>::swapColumns ( int C1, int C2 ) {
  _DataType tmp;
  for ( int i = 0; i < this->_numRows; ++i ) {
    tmp = ref ( i, C1 );
    ref ( i, C1 ) = ref ( i, C2 );
    ref ( i, C2 ) = tmp;
  }
}

template <class _DataType, bool Parallelize>
inline void aol::FullMatrix<_DataType, Parallelize>::multRow ( int R1, _DataType alpha, int StartIndex ) {
  for ( int j = StartIndex; j < this->_numCols; ++j ) {
    ref ( R1, j ) *= alpha;
  }
}

template <class _DataType, bool Parallelize>
inline void aol::FullMatrix<_DataType, Parallelize>::addToRow ( int R1, int R2, _DataType alpha, int StartIndex ) {
  for ( int j = StartIndex; j < this->_numCols; ++j ) {
    ref ( R1, j ) += ref ( R2, j ) * alpha;
  }
}

template <class _DataType, bool Parallelize>
void aol::FullMatrix<_DataType, Parallelize>::getBlock ( int r, int c, aol::FullMatrix<_DataType, Parallelize>& block ) const {
  int h = block.getNumRows ();
  int w = block.getNumCols ();
#ifdef BOUNDS_CHECK
  if ( r + h > this->getNumRows () || c + w > this->getNumCols () )
    throw aol::Exception ( "aol::FullMatrix<_DataType, Parallelize>::getBlock: Parameter mismatch: Block does not fit", __FILE__, __LINE__ );
#endif
  for ( int i = 0; i < h; ++i )
    for ( int j = 0; j < w; ++j )
      block.set ( i, j, ref ( r + i, c + j ) );
}

template <class _DataType, bool Parallelize>
void aol::FullMatrix<_DataType, Parallelize>::setBlock ( int r, int c, const aol::Matrix<_DataType>& block ) {
  int h = block.getNumRows ();
  int w = block.getNumCols ();
#ifdef BOUNDS_CHECK
  if ( r + h > this->getNumRows () || c + w > this->getNumCols () )
    throw aol::Exception ( "aol::FullMatrix<_DataType, Parallelize>::setBlock: Parameter mismatch: Block does not fit", __FILE__, __LINE__ );
#endif
  for ( int i = 0; i < h; ++i )
    for ( int j = 0; j < w; ++j )
      ref ( r + i, c + j ) = block.get ( i, j );
}

template <class _DataType, bool Parallelize>
void aol::FullMatrix<_DataType, Parallelize>::getSubColumn ( int r, int c, aol::Vector<_DataType>& subcol ) const {
  int h = subcol.size ();
#ifdef BOUNDS_CHECK
  if ( r + h > this->getNumRows () || c >= this->getNumCols () )
    throw aol::Exception ( "aol::FullMatrix<_DataType, Parallelize>::getSubColumn: Parameter mismatch: Block does not fit", __FILE__, __LINE__ );
#endif
  for ( int i = 0; i < h; ++i )
    subcol.set ( i, ref ( r + i, c ) );
}

template <class _DataType, bool Parallelize>
void aol::FullMatrix<_DataType, Parallelize>::setSubColumn ( int r, int c, const aol::Vector<_DataType>& subcol ) {
  int h = subcol.size ();
#ifdef BOUNDS_CHECK
  if ( r + h > this->getNumRows () || c >= this->getNumCols () )
    throw aol::Exception ( "aol::FullMatrix<_DataType, Parallelize>::setSubColumn: Parameter mismatch: Block does not fit", __FILE__, __LINE__ );
#endif
  for ( int i = 0; i < h; ++i )
    ref ( r + i, c ) = subcol.get ( i );
}

template <class _DataType, bool Parallelize>
void aol::FullMatrix<_DataType, Parallelize>::addSubColumn ( int r, int c, const aol::Vector<_DataType>& subcol ) {
  int h = subcol.size ();
#ifdef BOUNDS_CHECK
  if ( r + h > this->getNumRows () || c >= this->getNumCols () )
    throw aol::Exception ( "aol::FullMatrix<_DataType, Parallelize>::addSubColumn: Parameter mismatch: Block does not fit", __FILE__, __LINE__ );
#endif
  for ( int i = 0; i < h; ++i )
    ref ( r + i, c ) += subcol.get ( i );
}

template <class _DataType, bool Parallelize>
void aol::FullMatrix<_DataType, Parallelize>::getSubRow ( int r, int c, aol::Vector<_DataType>& subrow ) const {
  int w = subrow.size ();
#ifdef BOUNDS_CHECK
  if ( r >= this->getNumRows () || c + w > this->getNumCols () )
    throw aol::Exception ( "aol::FullMatrix<_DataType, Parallelize>::getSubRow: Parameter mismatch: Block does not fit", __FILE__, __LINE__ );
#endif
  for ( int i = 0; i < w; ++i )
    subrow.set ( i, ref ( r, c + i ) );
}

template <class _DataType, bool Parallelize>
void aol::FullMatrix<_DataType, Parallelize>::setSubRow ( int r, int c, const aol::Vector<_DataType>& subrow ) {
  int w = subrow.size ();
#ifdef BOUNDS_CHECK
  if ( r >= this->getNumRows () || c + w > this->getNumCols () )
    throw aol::Exception ( "aol::FullMatrix<_DataType, Parallelize>::setSubRow: Parameter mismatch: Block does not fit", __FILE__, __LINE__ );
#endif
  for ( int i = 0; i < w; ++i )
    ref ( r, c + i ) = subrow.get ( i );
}

template <class _DataType, bool Parallelize>
_DataType aol::FullMatrix<_DataType, Parallelize>::getFrobeniusNormSqr() const {
  _DataType frobeniusNormSqr = ZOTrait<_DataType>::zero;

  for ( int i = 0; i < this->getNumRows(); ++i )
    for ( int j = 0; j < this->getNumCols(); ++j )
      frobeniusNormSqr += aol::Sqr ( ref ( i, j ) );

  return ( frobeniusNormSqr );

}


template <typename DataType>
void aol::SymmetricMatrix<DataType>::init ( int Dimension ) {
  this->_numRows = this->_numCols = Dimension;

  data = new DataType[ ( this->_numCols* ( this->_numCols+1 ) ) >> 1 ];
  if ( data == NULL ) {
    cerr << "ERROR in aol::SymmetricMatrix<DataType>::aol::SymmetricMatrix: couldn't allocate memory!" << endl;
  }

  setZero();
}

template <typename DataType>
void aol::SymmetricMatrix<DataType>::destroy ( ) {
  if ( data != NULL )
    delete[] data;
}

template <typename DataType>
DataType aol::SymmetricMatrix<DataType>::get ( int I, int J ) const {
#ifdef BOUNDS_CHECK
  if ( I >= this->_numCols ) {
    cerr << "aol::SymmetricMatrix<DataType>::get() : I = " << I << " >= N! Using N-1..." << endl;
    I = this->_numCols - 1;
  } else if ( I < 0 ) {
    cerr << "aol::SymmetricMatrix<DataType>::get() : I = " << I << " < 0! Using 0..." << endl;
    I = 0;
  } else if ( J >= this->_numCols ) {
    cerr << "aol::SymmetricMatrix<DataType>::get() : J = " << J << " >= N! Using N-1..." << endl;
    J = this->_numCols - 1;
  } else if ( J < 0 ) {
    cerr << "aol::SymmetricMatrix<DataType>::get() : J = " << J << " < 0! Using 0..." << endl;
    J = 0;
  }
#endif
  if ( J >= I )
    return data[ ( ( this->_numCols* ( this->_numCols+1 )- ( this->_numCols-I ) * ( this->_numCols-I+1 ) ) >>1 ) + J-I];
  else
    return data[ ( ( this->_numCols* ( this->_numCols+1 )- ( this->_numCols-J ) * ( this->_numCols-J+1 ) ) >>1 ) + I-J];
}

template <typename DataType>
void aol::SymmetricMatrix<DataType>::set ( int I, int J, DataType Value ) {
#ifdef BOUNDS_CHECK
  if ( I >= this->_numCols ) {
    cerr << "aol::SymmetricMatrix<DataType>::set() : I = " << I << " >= N! Using N-1..." << endl;
    I = this->_numCols - 1;
  } else if ( I < 0 ) {
    cerr << "aol::SymmetricMatrix<DataType>::set() : I = " << I << " < 0! Using 0..." << endl;
    I = 0;
  } else if ( J >= this->_numCols ) {
    cerr << "aol::SymmetricMatrix<DataType>::set() : J = " << J << " >= N! Using N-1..." << endl;
    J = this->_numCols - 1;
  } else if ( J < 0 ) {
    cerr << "aol::SymmetricMatrix<DataType>::set() : J = " << J << " < 0! Using 0..." << endl;
    J = 0;
  }
#endif
  if ( J >= I )
    data[ ( ( this->_numCols* ( this->_numCols+1 )- ( this->_numCols-I ) * ( this->_numCols-I+1 ) ) >>1 ) + J-I] = Value;
  else
    data[ ( ( this->_numCols* ( this->_numCols+1 )- ( this->_numCols-J ) * ( this->_numCols-J+1 ) ) >>1 ) + I-J] = Value;
}

template <typename DataType>
void aol::SymmetricMatrix<DataType>::setZero ( ) {
  for ( int i = 0; i < this->_numCols* ( this->_numCols + 1 ) / 2; ++i )
    data[i] = 0;
}

template <typename DataType>
aol::Matrix<DataType>& aol::SymmetricMatrix<DataType>::operator*= ( const DataType alpha ) {
  for ( int i = 0; i < this->_numCols* ( this->_numCols + 1 ) / 2; ++i )
    data[i] *= alpha;
  return *this;
}

template <typename DataType>
void aol::SymmetricMatrix<DataType>::apply ( const aol::Vector<DataType> &Src, aol::Vector<DataType> &Dst ) const {
#ifdef BOUNDS_CHECK
  if ( ( Src.size() != this->_numCols ) || ( Dst.size() != this->_numCols ) ) {
    std::cout << "ERROR in aol::SymmetricMatrix<DataType>::mult: incompatible dimensions!" << endl;
    return;
  }
#endif
  Dst.setZero();
  int i, j;
  int offset;

  // upper triangle with diagonals
  for ( i = 0; i < this->_numCols ; ++i ) {
    offset = ( this->_numCols * ( this->_numCols + 1 ) - ( this->_numCols - i ) * ( this->_numCols - i + 1 ) ) / 2;
    for ( j = i; j < this->_numCols; ++j ) {
      Dst[i] += data[ offset + j-i] * Src.get ( j );
    }
  }
  // lower triangle
  for ( i = 1; i < this->_numCols ; ++i )
    for ( j = 0; j < i; ++j ) {
      Dst[i] += data[ ( this->_numCols* ( this->_numCols+1 )- ( this->_numCols-j ) * ( this->_numCols-j+1 ) ) /2 + i-j ] * Src.get ( j );
    }
}

// ============================================================================================
namespace aol {
template <> const aol::Format& aol::Matrix<int>::format = aol::mixedFormat;
template <> const aol::Format& aol::Matrix<float>::format = aol::mixedFormat;
template <> const aol::Format& aol::Matrix<double>::format = aol::mixedFormat;
template <> const aol::Format& aol::Matrix<long double>::format = aol::mixedFormat;
template <> bool aol::Matrix<int>::prettyFormat = true;
template <> bool aol::Matrix<float>::prettyFormat = true;
template <> bool aol::Matrix<double>::prettyFormat = true;
template <> bool aol::Matrix<long double>::prettyFormat = true;
}

template class aol::FullMatrix<int, true>;
template class aol::FullMatrix<double, true>;
template class aol::FullMatrix<float, true>;
template class aol::FullMatrix<long double, true>;
template class aol::FullMatrix<int, false>;
template class aol::FullMatrix<double, false>;
template class aol::FullMatrix<float, false>;
template class aol::FullMatrix<long double, false>;
template class aol::DiagonalMatrix<float>;
template class aol::DiagonalMatrix<double>;
template class aol::DiagonalMatrix<long double>;
template class aol::SymmetricMatrix<float>;
template class aol::SymmetricMatrix<double>;
template class aol::SymmetricMatrix<long double>;
template class aol::SparseMatrix<float>;
template class aol::SparseMatrix<double>;
template class aol::SparseMatrix<long double>;
template class aol::Matrix<int>;
template class aol::Matrix<float>;
template class aol::Matrix<double>;
template class aol::Matrix<long double>;
template class aol::PermutationMatrix<float>;
template class aol::PermutationMatrix<double>;
template class aol::PermutationMatrix<long double>;
