#include "diagBandMatrix.h"

namespace aol {

template < typename DataType, int BLower, int BUpper >
DiagBandMatrix<DataType, BLower, BUpper>::DiagBandMatrix ( ) : aol::Matrix<DataType> ( 0, 0 ), _pData ( NULL ), _nBands ( 1 + BLower + BUpper ) {
  this->reallocate ( 0 );
}


template < typename DataType, int BLower, int BUpper >
DiagBandMatrix<DataType, BLower, BUpper>::~DiagBandMatrix ( ) {
  if ( _pData != NULL ) {
    delete[] ( _pData );
  }
}


template < typename DataType, int BLower, int BUpper >
DiagBandMatrix<DataType, BLower, BUpper>::DiagBandMatrix ( const DiagBandMatrix<DataType,BLower,BUpper> &Other ) : aol::Matrix<DataType> ( Other ), _pData ( NULL ), _nBands ( 1 + BLower + BUpper ) {
  if ( _pData != NULL ) {
    delete[] ( _pData );
  }
  _pData = new DataType[ this->getNumRows() * _nBands ];
  memcpy ( _pData, Other._pData, this->getNumRows() * _nBands * sizeof ( DataType ) );
}


template < typename DataType, int BLower, int BUpper >
DiagBandMatrix<DataType,BLower,BUpper>&  DiagBandMatrix<DataType, BLower, BUpper>::operator= ( const DiagBandMatrix<DataType,BLower,BUpper> &Other ) {
  if ( this->getNumRows() != Other.getNumRows() ) {
    throw aol::Exception ( "aol::DiagBandMatrix::operator= trying to assign an incompatible matrix.", __FILE__, __LINE__ );
  }
  // _nBands is then correct
  if ( this != &Other ) { // beware of self-assignment
    memcpy ( _pData, Other._pData, this->getNumRows() * _nBands * sizeof ( DataType ) );
  }
  return ( *this );
}


template < typename DataType, int BLower, int BUpper >
DiagBandMatrix<DataType, BLower, BUpper>::DiagBandMatrix ( const int NRowsCols ) : aol::Matrix<DataType> ( NRowsCols, NRowsCols ), _pData ( NULL ), _nBands ( 1 + BLower + BUpper ) {
  this->reallocate ( NRowsCols );
}

template < typename DataType, int BLower, int BUpper >
DiagBandMatrix<DataType, BLower, BUpper>::DiagBandMatrix ( const qc::GridStructure &Size ) : aol::Matrix<DataType> ( Size.getNumberOfNodes(), Size.getNumberOfNodes() ), _pData ( NULL ), _nBands ( 1 + BLower + BUpper ) {
  this->reallocate ( Size.getNumberOfNodes() );
}

template < typename DataType, int BLower, int BUpper >
DiagBandMatrix<DataType, BLower, BUpper>::DiagBandMatrix ( const int NRows, const int NCols ) : aol::Matrix<DataType> ( NRows, NCols ), _pData ( NULL ), _nBands ( 1 + BLower + BUpper ) {
  if ( NRows != NCols ) {
    throw aol::Exception ( "aol::DiagBandMatrix must be square" );
  }
  this->reallocate ( NRows );
}


template < typename DataType, int BLower, int BUpper >
void DiagBandMatrix<DataType, BLower, BUpper>::reallocate ( const int nRowsCols ) {
  this->matrixInit ( nRowsCols, nRowsCols );
  if ( _pData != NULL ) {
    delete[] ( _pData );
  }
  _pData = new DataType[ this->getNumRows() * _nBands ];
  setZero();
}


template < typename DataType, int BLower, int BUpper >
void DiagBandMatrix<DataType, BLower, BUpper>::reallocate ( const int nRows, const int nCols ) {
  if ( nRows == nCols ) {
    this->reallocate ( nRows );
  } else {
    throw aol::Exception ( "aol::DiagBandMatrix::reallocate: non-sqare size provided.", __FILE__, __LINE__ );
  }
}


template < typename DataType, int BLower, int BUpper >
DataType DiagBandMatrix<DataType, BLower, BUpper>::get ( int I, int J ) const {
#ifdef BOUNDS_CHECK
  this->boundsCheck ( I, J, "aol::DiagBandMatrix::get: Index out of bounds", __FILE__, __LINE__ );
#endif
  return ( validEntry ( I, J ) ? _pData[ dataIndex ( I, J ) ] : aol::NumberTrait<DataType>::zero );
}


template < typename DataType, int BLower, int BUpper >
void DiagBandMatrix<DataType, BLower, BUpper>::set ( int I, int J, DataType Val ) {
#ifdef BOUNDS_CHECK
  this->boundsCheck ( I, J, "aol::DiagBandMatrix::set: Index out of bounds", __FILE__, __LINE__ );
  if ( ! ( validEntry ( I, J ) ) && ( Val != aol::NumberTrait<DataType>::zero ) ) {
    cerr << I << " " << J << " " << BLower << " " << BUpper << endl;
    throw aol::Exception ( "aol::DiagBandMatrix::set: Illegal position", __FILE__, __LINE__ );
  }
#endif
  _pData[ dataIndex ( I, J ) ] = Val;
}


template < typename DataType, int BLower, int BUpper >
void DiagBandMatrix<DataType, BLower, BUpper>::add ( int I, int J, DataType Val ) {
#ifdef BOUNDS_CHECK
  this->boundsCheck ( I, J, "aol::DiagBandMatrix::add: Index out of bounds", __FILE__, __LINE__ );
  if ( ! ( validEntry ( I, J ) ) && ( Val != aol::NumberTrait<DataType>::zero ) ) {
    throw aol::Exception ( "aol::DiagBandMatrix::add: Illegal position", __FILE__, __LINE__ );
  }
#endif
  _pData[ dataIndex ( I, J ) ] += Val;
}


template < typename DataType, int BLower, int BUpper >
void DiagBandMatrix<DataType, BLower, BUpper>::setZero ( ) {
  memset ( _pData, 0, _nBands * this->getNumRows() * sizeof ( DataType ) );
}


template < typename DataType, int BLower, int BUpper >
void DiagBandMatrix<DataType, BLower, BUpper>::apply ( const Vector<DataType> &Arg, Vector<DataType> &Dest ) const {
  const int N = this->getNumRows();

#ifdef BOUNDS_CHECK
  if ( this->getNumRows() != Dest.size() || this->getNumCols() != Arg.size() ) {
    throw ( Exception ( "aol::DiagBandMatrix::apply: incorrect sizes.", __FILE__, __LINE__ ) );
  }
  if ( N < _nBands ) {
    throw aol::Exception ( "aol::DiagBandMatrix::apply: matrix too small", __FILE__, __LINE__ );
  }
#endif

  for ( int i = 0; i < BLower; ++i ) {
    DataType result = aol::ZOTrait<DataType>::zero;
    for ( int j = 0; j <= i + BUpper; ++j ) {
      result += _pData [ dataIndex ( i, j ) ] * Arg[ j ];
    }
    Dest[i] = result;
  }

  for ( int i = BLower; i < N - BUpper ; ++i ) {
    DataType result = aol::ZOTrait<DataType>::zero;
    for ( int offs = -BLower; offs <= BUpper; ++offs ) { // this loop will hopefully be unrolled
      result += _pData [ dataIndexRO ( i, offs ) ] * Arg[ i + offs ];
    }
    Dest[i] = result;
  }

  for ( int i = N - BUpper; i < N ; ++i ) {
    DataType result = aol::ZOTrait<DataType>::zero;
    for ( int j = i - BLower; j < N; ++j ) {
      result += _pData [ dataIndex ( i, j ) ] * Arg[ j ];
    }
    Dest[i] = result;
  }
}


template < typename DataType, int BLower, int BUpper >
void DiagBandMatrix<DataType, BLower, BUpper>::applyAdd ( const Vector<DataType> &Arg, Vector<DataType> &Dest ) const {
  const int N = this->getNumRows();

#ifdef BOUNDS_CHECK
  if ( this->getNumRows() != Dest.size() || this->getNumCols() != Arg.size() ) {
    throw ( Exception ( "aol::DiagBandMatrix::apply: incorrect sizes.", __FILE__, __LINE__ ) );
  }
  if ( N < _nBands ) {
    throw aol::Exception ( "aol::DiagBandMatrix::apply: matrix too small", __FILE__, __LINE__ );
  }
#endif

  for ( int i = 0; i < BLower; ++i ) {
    DataType result = aol::ZOTrait<DataType>::zero;
    for ( int j = 0; j <= i + BUpper; ++j ) {
      result += _pData [ dataIndex ( i, j ) ] * Arg[ j ];
    }
    Dest[i] += result;
  }

  for ( int i = BLower; i < N - BUpper ; ++i ) {
    DataType result = aol::ZOTrait<DataType>::zero;
    for ( int offs = -BLower; offs <= BUpper; ++offs ) { // this loop will hopefully be unrolled
      result += _pData [ dataIndexRO ( i, offs ) ] * Arg[ i + offs ];
    }
    Dest[i] += result;
  }

  for ( int i = N - BUpper; i < N ; ++i ) {
    DataType result = aol::ZOTrait<DataType>::zero;
    for ( int j = i - BLower; j < N; ++j ) {
      result += _pData [ dataIndex ( i, j ) ] * Arg[ j ];
    }
    Dest[i] += result;
  }
}


template < typename DataType, int BLower, int BUpper >
void DiagBandMatrix<DataType, BLower, BUpper>::makeRowEntries ( vector<typename Row<DataType>::RowEntry > &Vec, const int RowNum ) const {
  if ( RowNum < BLower ) {
    Vec.resize ( RowNum + BUpper + 1 );
    for ( int j = 0; j <= RowNum + BUpper; ++j ) {
      Vec[j].col = j;
      Vec[j].value = _pData [ dataIndex ( RowNum, j ) ];
    }
  } else if ( RowNum < this->getNumRows() - BUpper ) {
    Vec.resize ( _nBands );
    for ( int offs = -BLower; offs <= BUpper; ++offs ) { // this loop will hopefully be unrolled
      Vec[ offs + BLower ].col = RowNum + offs;
      Vec[ offs + BLower ].value = _pData [ dataIndexRO ( RowNum, offs ) ];
    }
  } else {
    Vec.resize ( BLower + ( this->getNumRows() - RowNum ) );
    for ( int j = RowNum - BLower; j < this->getNumRows(); ++j ) {
      Vec[ j - RowNum + BLower ].col = j;
      Vec[ j - RowNum + BLower ].value = _pData [ dataIndex ( RowNum, j ) ];
    }
  }
}

template < typename DataType, int BLower, int BUpper >
void DiagBandMatrix<DataType, BLower, BUpper>::makeSortedRowEntries ( vector<typename Row<DataType>::RowEntry > &Vec, const int RowNum ) const {
  makeRowEntries ( Vec, RowNum );
}


template class DiagBandMatrix < float,       1, 1 >;
template class DiagBandMatrix < double,      1, 1 >;
template class DiagBandMatrix < long double, 1, 1 >;

template class DiagBandMatrix < float,       1, 2 >;
template class DiagBandMatrix < double,      1, 2 >;
template class DiagBandMatrix < long double, 1, 2 >;

template class DiagBandMatrix < float,       2, 1 >;
template class DiagBandMatrix < double,      2, 1 >;
template class DiagBandMatrix < long double, 2, 1 >;

  template class aol::DiagBandMatrix < double, 10, 10 >; // will be used in self-test

}
