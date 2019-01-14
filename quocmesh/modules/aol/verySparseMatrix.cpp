#include "verySparseMatrix.h"

namespace aol {

template <typename DataType>
VerySparseMatrix<DataType>::VerySparseMatrix ( ) : _rows ( ) {
}


template <typename DataType>
VerySparseMatrix<DataType>::~VerySparseMatrix ( ) {
}


template <typename DataType>
VerySparseMatrix<DataType>::VerySparseMatrix ( const VerySparseMatrix<DataType> &other ) : aol::Matrix< DataType > ( other._numRows, other._numCols ), _rows ( other._rows ) {
}


template <typename DataType>
VerySparseMatrix<DataType>& VerySparseMatrix<DataType>::operator= ( const VerySparseMatrix<DataType> &other ) {
  if ( this == &other ) { // beware of self-assignment
    return ( *this );
  }

  if ( other.getNumRows() != this->getNumRows() || other.getNumCols() != this->getNumCols() ) {
    throw Exception ( "VerySparseMatrix::operator= : dimensions don't match.", __FILE__, __LINE__ ); // and we do not want to change the size
  }

  this->_rows = other._rows;

  return ( *this );
}


template <typename DataType>
VerySparseMatrix<DataType>::VerySparseMatrix ( const int Rows, const int Columns ) : aol::Matrix< DataType > ( Rows, Columns ) {
}


template <typename DataType>
void VerySparseMatrix<DataType>::reallocate ( const int Rows, const int Columns ) {
  setZero();
  this->_numRows = Rows;
  this->_numCols = Columns;
}


template <typename DataType>
DataType VerySparseMatrix<DataType>::get ( int I, int J ) const {
#ifdef BOUNDS_CHECK
  if ( I < 0 || J < 0 || I >= static_cast<int> ( this->getNumRows() ) || J >= static_cast<int> ( this->getNumCols() ) ) {
    cerr << I << " " << J << " is out of bounds: " << this->getNumRows() << " " << this->getNumCols() << endl;
    throw aol::Exception ( "aol::VerySparseMatrix::get: Index out of bounds", __FILE__, __LINE__ );
  }
#endif

  if ( _rows.contains ( I ) ) {
    return ( _rows.getRef ( I ).get ( J ) );
  } else {
    return ( aol::NumberTrait<DataType>::zero );
  }
}


template <typename DataType>
void VerySparseMatrix<DataType>::set ( int I, int J, DataType Val ) {
#ifdef BOUNDS_CHECK
  if ( I < 0 || J < 0 || I >= static_cast<int> ( this->getNumRows() ) || J >= static_cast<int> ( this->getNumCols() ) ) {
    cerr << I << " " << J << " is out of bounds: " << this->getNumRows() << " " << this->getNumCols() << endl;
    throw aol::Exception ( "aol::VerySparseMatrix::get: Index out of bounds", __FILE__, __LINE__ );
  }
#endif

  _rows.getRefOrCreate ( I ).set ( J, Val );
}


template <typename DataType>
void VerySparseMatrix<DataType>::add ( int I, int J, DataType Val ) {
#ifdef BOUNDS_CHECK
  if ( I < 0 || J < 0 || I >= static_cast<int> ( this->getNumRows() ) || J >= static_cast<int> ( this->getNumCols() ) ) {
    cerr << I << " " << J << " is out of bounds: " << this->getNumRows() << " " << this->getNumCols() << endl;
    throw aol::Exception ( "aol::VerySparseMatrix::get: Index out of bounds", __FILE__, __LINE__ );
  }
#endif

  _rows.getRefOrCreate ( I ).add ( J, Val );
}


template <typename DataType>
void VerySparseMatrix<DataType>::setZero ( ) {
  _rows.clear();
}


template <typename DataType>
void VerySparseMatrix<DataType>::apply ( const Vector<DataType> &Arg, Vector<DataType> &Dest ) const {
  if ( this->getNumRows() != Dest.size() || this->getNumCols() != Arg.size() ) {
    throw ( Exception ( "aol::VerySparseMatrix::apply: incorrect sizes.", __FILE__, __LINE__ ) );
  }

  for ( typename aol::LookupMap < int, aol::SparseRow<DataType> >::const_iterator rit = _rows.begin(); rit != _rows.end(); ++rit ) {
    Dest[ rit->first ] = (rit->second).mult ( Arg, rit->first );
  }
}


template <typename DataType>
void VerySparseMatrix<DataType>::applyAdd ( const Vector<DataType> &Arg, Vector<DataType> &Dest ) const {
  if ( this->getNumRows() != Dest.size() || this->getNumCols() != Arg.size() ) {
    throw ( Exception ( "aol::VerySparseMatrix::apply: incorrect sizes.", __FILE__, __LINE__ ) );
  }

  for ( typename aol::LookupMap < int, aol::SparseRow<DataType> >::const_iterator rit = _rows.begin(); rit != _rows.end(); ++rit ) {
    Dest[ rit->first ] += (rit->second).mult ( Arg, rit->first );
  }
}


template <typename DataType>
void VerySparseMatrix<DataType>::makeRowEntries ( vector<typename Row<DataType>::RowEntry > &Vec, const int RowNum ) const {
  if ( _rows.contains ( RowNum ) ) {
    _rows.getRef ( RowNum ).makeRowEntries ( Vec, RowNum );
  } else {
    Vec.clear();
  }
}


template <typename DataType>
void VerySparseMatrix<DataType>::makeSortedRowEntries ( vector<typename Row<DataType>::RowEntry > &Vec, const int RowNum ) const {
  if ( _rows.contains ( RowNum ) ) {
    _rows.getRef ( RowNum ).makeSortedRowEntries ( Vec, RowNum );
  } else {
    Vec.clear();
  }
}


template <class DataType>
ostream& aol::VerySparseMatrix<DataType>::printVerySparse ( ostream& os, const aol::Format DataFormatter ) const {
  this->printHead ( os );
  for ( typename aol::LookupMap < int, aol::SparseRow<DataType> >::const_iterator rit = _rows.begin(); rit != _rows.end(); ++rit ) {
    os << "Row " << aol::intFormat ( rit->first ) << ": ";
    std::vector<typename Row<DataType>::RowEntry > vec;
    rit->second.makeRowEntries ( vec, rit->first );
    for ( typename vector<typename Row<DataType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
      os << "(" << aol::intFormat ( rit->first ) << ", " << aol::intFormat ( it->col ) << "): " << DataFormatter ( it->value ) << "   ";
    }
    os << endl;
  }
  return os;
}


template <class DataType>
VerySparseMatrix<DataType>& VerySparseMatrix<DataType>::operator*= ( const DataType Factor ) {
  for ( typename aol::LookupMap < int, aol::SparseRow<DataType> >::iterator rit = _rows.begin(); rit != _rows.end(); ++rit ) {
    (rit->second).scale ( Factor );
  }
  return (*this);
}


template class VerySparseMatrix<float>;
template class VerySparseMatrix<double>;
template class VerySparseMatrix<long double>;

}

