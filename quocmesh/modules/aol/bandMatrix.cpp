#include <bandMatrix.h>

namespace aol {

template < typename _DataType >
GenBandMatrix<_DataType>::GenBandMatrix ( const int Rows, const int Cols, const aol::Vector<int> &Offsets )
  : aol::Matrix<_DataType> ( Rows, Cols ),
    _startApply ( 0 ), _endApply ( 0 ),  _globalToLocal ( 0 ), _localToGlobal ( 0 ),
    _pData ( NULL ),  _nDiags ( Offsets.size() ), _sizeReserved ( 0 ) {

  aol::Vector<int> temp_offsets ( Offsets.size() );
  for ( int i = 0; i < Offsets.size(); ++i ) {
    temp_offsets[i] = Offsets[i];
  }

  this->reallocate ( this->getNumRows(), this->getNumCols(), temp_offsets );
}


template < typename _DataType >
GenBandMatrix<_DataType>::GenBandMatrix ( const int Rows, const int Cols )
  : aol::Matrix<_DataType> ( Rows, Cols ),
    _startApply ( 0 ), _endApply ( 0 ),  _globalToLocal ( 0 ), _localToGlobal ( 0 ),
    _pData ( NULL ), _nDiags ( 0 ) {

  aol::Vector<int> dummy;
  this->reallocate ( this->getNumRows(), this->getNumCols(), dummy );
}


template < typename _DataType >
GenBandMatrix<_DataType>::GenBandMatrix ( const qc::GridStructure & Grid )
  : aol::Matrix<_DataType> ( Grid.getNumberOfNodes(), Grid.getNumberOfNodes() ),
    _startApply ( 0 ), _endApply ( 0 ),  _globalToLocal ( 0 ), _localToGlobal ( 0 ),
    _pData ( NULL ), _nDiags ( 0 ) {

  aol::Vector<int> dummy;
  this->reallocate ( this->getNumRows(), this->getNumCols(), dummy );
}


template < typename _DataType >
GenBandMatrix<_DataType>::GenBandMatrix ( )
  : Matrix<_DataType> ( ),
    _startApply ( 0 ), _endApply ( 0 ),  _globalToLocal ( 0 ), _localToGlobal ( 0 ),
    _pData ( NULL ), _nDiags ( 0 ) {

  aol::Vector<int> dummy;
  this->reallocate ( this->getNumRows(), this->getNumCols(), dummy );
}


template < typename _DataType >
GenBandMatrix<_DataType>::GenBandMatrix ( const GenBandMatrix<_DataType> &Other )
  : Matrix<_DataType> ( Other ),
    _startApply ( Other._startApply ), _endApply ( Other._endApply ),  _globalToLocal ( Other._globalToLocal ), _localToGlobal ( Other._localToGlobal ),
    _pData ( NULL ), _nDiags ( Other._nDiags ) {

  _pData = new DataType[ this->getNumRows() * _nDiags ];
  memcpy ( _pData, Other._pData, this->getNumRows() * _nDiags * sizeof ( DataType ) );
}


template < typename _DataType >
GenBandMatrix<_DataType>::~GenBandMatrix() {
  if ( _pData ) {
    aol::MemoryManager::deallocate ( _pData, _sizeReserved, sizeof ( DataType ) );
  }
}


template < typename _DataType >
/*virtual*/ void GenBandMatrix<_DataType>::reallocate ( const int Rows, const int Cols, const aol::Vector<int>& Offsets ) {
  this->matrixInit ( Rows, Cols );

  if ( _pData != NULL ) {
    aol::MemoryManager::deallocate ( _pData, _sizeReserved, sizeof ( DataType ) );
  }

  _nDiags = Offsets.size();

  _localToGlobal.resize ( _nDiags );
  _localToGlobal.setZero();

  for ( int i = 0; i < _nDiags; ++i ) {
    _localToGlobal[i] = Offsets[i];

    if ( i > 0 ) {
      if ( _localToGlobal[i] < _localToGlobal[i-1] ) {
        cerr << "Offsets for GenBandMatrix not sorted. MakeSortedRowEntries will not work correctly!" << endl;
      }
    }
  }

  _globalToLocal.resize ( 2*this->getNumCols() - 1 );
  _globalToLocal.setZero();
  for ( int i = 0; i < 2*this->getNumCols() - 1; ++i ) {
    _globalToLocal[i] = UNDEFINED_ENTRY;
  }
  for ( int i = 0; i < _nDiags; ++i ) {
    _globalToLocal[ this->getNumCols() - 1 + _localToGlobal[i] ] = i;
  }

  _startApply.resize ( this->getNumRows() );
  _startApply.setZero();
  _endApply.resize ( this->getNumRows() );
  _endApply.setZero();
  for ( int i = 0; i <  this->getNumRows() ; ++i ) {
    if ( _nDiags != 0 ) {
      int j = 0, k = _nDiags;
      while ( ( i + _localToGlobal[ j   ] ) <   0 ) ++j;
      while ( ( i + _localToGlobal[ k-1 ] ) >= this->getNumCols() ) --k;
      _startApply[i] = j;
      _endApply[i]   = k;
    } else {
      _startApply[i] = this->getNumCols();
      _endApply[i]   = 0;
    }
  }

  if ( _nDiags != 0 ) {
    _sizeReserved = this->getNumRows() * _nDiags; // possible change in allocateAtLeast
    _pData = static_cast<DataType*> ( aol::MemoryManager::allocateAtLeast ( _sizeReserved, sizeof ( DataType ) ) );
  } else {
    _pData = NULL;
  }

  setZero();
}



template < typename _DataType >
GenBandMatrix<_DataType>& GenBandMatrix<_DataType>::operator*= ( const DataType Alpha ) {
  for ( int j = 0; j <  this->getNumRows() ; ++j ) {
    for ( int i = 0; i < _nDiags; ++i ) {
      _pData[map_index ( i,j ) ] *= Alpha;
    } // note that the undefined entries remain zero
  }

  return *this;
}


template < typename _DataType >
bool GenBandMatrix<_DataType>::isApproxEqual ( const GenBandMatrix<DataType> &Other, const DataType Epsilon ) {
  if ( !hasSameStructureAs ( Other ) ) {
    std::cerr << "Matrices do not have the same structure " << std::endl;
    return false;
  }

  for ( int j = 0; j <  this->getNumRows() ; ++j ) {
    for ( int i = 0; i < _nDiags; ++i ) {
      const int index = map_index ( i, j );
      if ( fabs ( _pData[index] - Other._pData[index] ) > Epsilon ) {
        std::cerr << "Data differs by " << fabs ( _pData[index] - Other._pData[index] ) << " at pos ("
                  << i << ", " << j << "): "
                  << _pData[index] << " != " << Other._pData[index] << endl;
        return false;
      }
    }
  }
  return true;
}


template < typename _DataType >
GenBandMatrix<_DataType>& GenBandMatrix<_DataType>::addMultiple ( const GenBandMatrix<_DataType> &Mat, DataType Factor ) {
#ifdef DEBUG
  if ( !hasSameStructureAs ( Mat ) ) {
    throw Exception ( "GenBandMatrix::operator+= : trying to add an incompatible matrix", __FILE__, __LINE__ );
  }
#endif
  for ( int i = 0; i < _nDiags; ++i ) {
    for ( int j = 0; j <  this->getNumRows() ; ++j ) {
      const int index = map_index ( i, j );
      _pData[index] += Factor * Mat._pData[index];
    }
  }
  return *this;
}


template < typename _DataType >
GenBandMatrix<_DataType>& GenBandMatrix<_DataType>::operator= ( const GenBandMatrix<DataType> &Mat ) {
  if ( !hasSameStructureAs ( Mat ) ) {
    throw Exception ( "GenBandMatrix::operator= : trying to assign an incompatible matrix", __FILE__, __LINE__ );
  } else {
    // Beware of self-assignment
    if ( this != &Mat ) {
      memcpy ( _pData, Mat._pData, _nDiags *  this->getNumRows()  * sizeof ( DataType ) );
    }
    // do not copy other members
  }
  return *this;
}


template < typename _DataType >
/*virtual*/ void GenBandMatrix<_DataType>::apply ( const Vector<DataType> &src, Vector<DataType> &dst ) const {
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for ( int i = 0; i <  this->getNumRows() ; ++i ) {
    DataType result = aol::ZOTrait<DataType>::zero;
    for ( int j = _startApply[i]; j < _endApply[i]; ++j )
      // note: first indices have different order!
      result += _pData[map_index ( j, i ) ] * src[i + _localToGlobal[j]];
    dst[i] = result;
  }
}


template < typename _DataType >
void GenBandMatrix<_DataType>::applyMasked ( const Vector<DataType> &Arg, Vector<DataType> &Dest,
                                             const BitVector & Mask, IncludeWriteMode applyMode ) const {
  switch ( applyMode ) {
  case INCLUDE_ALL_WRITE_ALL:
    applyMasked<BitMaskFunctorTrue, BitMaskFunctorTrue> ( Arg, Dest, Mask );
    break;
  case INCLUDE_BD_WRITE_INT:
    applyMasked<BitMaskFunctorNegate, BitMaskFunctorIdentity> ( Arg, Dest, Mask );
    break;
  case INCLUDE_ALL_WRITE_INT:
    applyMasked<BitMaskFunctorTrue, BitMaskFunctorIdentity> ( Arg, Dest, Mask );
    break;
  case INCLUDE_INT_WRITE_ALL:
    applyMasked<BitMaskFunctorIdentity, BitMaskFunctorTrue> ( Arg, Dest, Mask );
    break;
  case INCLUDE_INT_WRITE_INT:
    applyMasked<BitMaskFunctorIdentity, BitMaskFunctorIdentity> ( Arg, Dest, Mask );
    break;
  default:
    throw aol::UnimplementedCodeException ( "GenBandMatrix::applyMasked: unknown IncludeWriteMode", __FILE__, __LINE__ );
  }
}


template < typename _DataType >
/*virtual*/ void GenBandMatrix<_DataType>::applyAdd ( const Vector<DataType> &src, Vector<DataType> &dst ) const {
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for ( int i = 0; i <  this->getNumRows() ; ++i ) {
    DataType result = aol::ZOTrait<DataType>::zero;
    for ( int j = _startApply[i]; j < _endApply[i]; ++j ) {
      result += _pData[ map_index ( j, i ) ] * src[ i + _localToGlobal[j] ];
    }
    dst[i] += result;
  }
}


template < typename _DataType >
/*virtual*/ void GenBandMatrix<_DataType>::applyAddMasked ( const Vector<DataType> &Arg, Vector<DataType> &Dest,
                                                            const BitVector & Mask, IncludeWriteMode applyMode ) const {

  // check dimensions
  if ( this->getNumRows() != Dest.size() || this->getNumCols() != Arg.size() ) {
    char errmsg [ 1024 ];
    sprintf ( errmsg, "aol::GenBandMatrix::apply: Cannot apply %d by %d matrix from vector of size %d to vector of size %d.", this->getNumRows(), this->getNumCols(), Arg.size(), Dest.size() );
    throw ( Exception ( errmsg, __FILE__, __LINE__ ) );
  }

  switch ( applyMode ) {
  case INCLUDE_ALL_WRITE_ALL:
    applyAddMasked<BitMaskFunctorTrue, BitMaskFunctorTrue> ( Arg, Dest, Mask );
    break;
  case INCLUDE_BD_WRITE_INT:
    applyAddMasked<BitMaskFunctorNegate, BitMaskFunctorIdentity> ( Arg, Dest, Mask );
    break;
  case INCLUDE_ALL_WRITE_INT:
    applyAddMasked<BitMaskFunctorTrue, BitMaskFunctorIdentity> ( Arg, Dest, Mask );
    break;
  case INCLUDE_INT_WRITE_ALL:
    applyAddMasked<BitMaskFunctorIdentity, BitMaskFunctorTrue> ( Arg, Dest, Mask );
    break;
  case INCLUDE_INT_WRITE_INT:
    applyAddMasked<BitMaskFunctorIdentity, BitMaskFunctorIdentity> ( Arg, Dest, Mask );
    break;
  default:
    throw aol::UnimplementedCodeException ( "GenBandMatrix::applyAddMasked: unknown IncludeWriteMode", __FILE__, __LINE__ );
  }
}


template < typename _DataType >
bool GenBandMatrix<_DataType>::hasSameStructureAs ( const aol::GenBandMatrix<DataType> &Mat ) const {
  bool result = (  this->getNumRows()  == Mat.getNumRows() ) && ( this->getNumCols() == Mat.getNumCols() ) && ( _nDiags == Mat._nDiags );
  if ( result ) {
    for ( int i = 0; i < _nDiags; ++i ) {
      result &= ( _localToGlobal[i] == Mat._localToGlobal[i] );
    }
  }
  return ( result );
}


template < typename _DataType >
void GenBandMatrix<_DataType>::setRowToZero ( const int i ) {
  for ( int d = 0; d < _nDiags; ++d ) {
    _pData[map_index ( d, i ) ] = aol::NumberTrait<DataType>::zero;
  }
}


template < typename _DataType >
void GenBandMatrix<_DataType>::setColToZero ( const int i ) {
  for ( int d = 0; d < _nDiags; ++d ) {
    const int row = i - _localToGlobal[d];
    if ( row >= 0 && row <  this->getNumRows()  ) {
      _pData[map_index ( d, row ) ] = aol::NumberTrait<DataType>::zero;
    }
  }
}

template< typename RealType > const int GenBandMatrix<RealType>::UNDEFINED_ENTRY = std::numeric_limits<int>::min();

template class GenBandMatrix<float>;
template class GenBandMatrix<double>;
template class GenBandMatrix<long double>;


template < typename DataType >
BandMatrix<DataType>::BandMatrix ( const int M, const int N, const int FirstOffset, const int NBands )
  : GenBandMatrix<DataType> ( M, N ),
    _firstOffset ( FirstOffset ), _nBands ( NBands ) {

  cerr << "aol::BandMatrix has not been tested yet." << endl; // remove this if you are sure this class works.

  _firstComplete = aol::Max ( 0, - FirstOffset );
  _lastComplete  = aol::Min ( this->_numRows, this->_numCols, - ( FirstOffset + NBands ) );

  aol::Vector<int> offsets ( _nBands );
  for ( int i = 0; i < _nBands; ++i ) {
    offsets[i] = _firstOffset + i;
  }

  if ( N < 3 ) {
    throw Exception ( "BandMatrix: sizes and offsets incompatible.", __FILE__, __LINE__ );
  } else {
    this->reallocate ( M, N, offsets );
  }

  cerr << "Setting up band matrix " << this->_numRows << " by " << this->_numCols << ", " << _nBands << " diagonals" << endl;
}


template < typename DataType >
/*virtual*/ void BandMatrix<DataType>::apply ( const Vector<DataType> &src, Vector<DataType> &dst ) const {
  // might want to compare sizes
  DataType result = aol::NumberTrait<DataType>::zero;
  for ( int i = 0; i < _firstComplete; ++i ) {
    result = aol::NumberTrait<DataType>::zero;
    for ( int j = 0; j < ( i + _firstOffset + _nBands + 0 ) ; ++j ) {
      //  k = j - i - _firstOffset
      result += this->_pData[ this->map_index ( j - i - _firstOffset, i ) ] * src[j];
    }
    dst[i] = result;
  }
  for ( int i = _firstComplete; i < ( _lastComplete + 1 ); ++i ) {
    result = aol::NumberTrait<DataType>::zero;
    for ( int j = i + _firstOffset; j < ( i + _firstOffset + _nBands + 0 ) ; ++j ) {
      result += this->_pData[ this->map_index ( j - i - _firstOffset, i ) ] * src[j];
    }
    dst[i] = result;
  }
  for ( int i = ( _lastComplete + 1 ); i < this->_numRows; ++i ) {
    result = aol::NumberTrait<DataType>::zero;
    for ( int j = i + _firstOffset; j < this->_numCols ; ++j ) {
      //  k = j - i - _firstOffset
      result += this->_pData[ this->map_index ( j - i - _firstOffset, i ) ] * src[j];
    }
    dst[i] = result;
  }
}


template class BandMatrix<float>;
template class BandMatrix<double>;
template class BandMatrix<long double>;



template < typename DataType >
/*explicit*/ TriBandMatrix<DataType>::TriBandMatrix ( const int N ) : GenBandMatrix<DataType> ( N, N ) {
  aol::Vector< int > offsets ( 3 );
  offsets[0] = -1;
  offsets[1] =  0;
  offsets[2] =  1;

  if ( N < 2 ) {
    throw Exception ( "TriBandMatrix: size too small. Will not work "
                      "and would not be useful.", __FILE__, __LINE__ );
  } else {
    this->reallocate ( N, N, offsets );
  }
}


template < typename DataType >
TriBandMatrix<DataType>::TriBandMatrix ( const int Nx, const int Ny ) : GenBandMatrix<DataType> ( Nx, Nx ) {
  aol::Vector< int > offsets ( 3 );
  offsets[0] = -1;
  offsets[1] =  0;
  offsets[2] =  1;

  if ( Nx < 2  ) {
    throw Exception ( "TriBandMatrix: size too small. Will not work "
                      "and would not be useful.", __FILE__, __LINE__ );
  } else if ( Nx != Ny ) {
    throw Exception ( "TriBandMatrix: needs to be square.", __FILE__, __LINE__ );
  } else {
    this->reallocate ( Nx, Nx, offsets );
  }
}


template < typename DataType >
/*explicit*/ TriBandMatrix<DataType>::TriBandMatrix ( const qc::GridDefinition & Grid ) : GenBandMatrix<DataType> ( Grid ) {
  int n = Grid.getNumberOfNodes();
  aol::Vector< int > offsets ( 3 );
  offsets[0] = -1;
  offsets[1] =  0;
  offsets[2] =  1;

  if ( n < 2 ) {
    throw Exception ( "TriBandMatrix: size too small. Will not work "
                      "and would not be useful.", __FILE__, __LINE__ );
  } else {
    this->reallocate ( n, n, offsets );
  }
}

template < typename DataType >
/*explicit*/ TriBandMatrix<DataType>::TriBandMatrix ( const qc::GridSize<qc::QC_1D> &Size ) : GenBandMatrix<DataType> ( Size.getNumberOfNodes(), Size.getNumberOfNodes() ) {
  int n = Size.getNumberOfNodes();
  aol::Vector< int > offsets ( 3 );
  offsets[0] = -1;
  offsets[1] =  0;
  offsets[2] =  1;

  if ( n < 2 ) {
    throw Exception ( "TriBandMatrix: size too small. Will not work "
                      "and would not be useful.", __FILE__, __LINE__ );
  } else {
    this->reallocate ( n, n, offsets );
  }
}

template < typename DataType >
/*virtual*/ void TriBandMatrix<DataType>::apply ( const Vector<DataType> &Src, Vector<DataType> &Dst ) const {
  // might want to compare sizes
  const int N = this->getNumRows();
  Dst[0] = this->_pData[this->map_index ( 1,0 ) ] * Src[0] + this->_pData[this->map_index ( 2,0 ) ] * Src[1];
  for ( int i = 1; i < this->getNumRows()  - 1; ++i ) {
    Dst[i] = this->_pData[this->map_index ( 0,i ) ] * Src [i-1] + this->_pData[this->map_index ( 1,i ) ] * Src[i] + this->_pData[this->map_index ( 2,i ) ] * Src[i+1];
  }
  Dst[ N - 1 ]  = this->_pData[this->map_index ( 0, N - 1 ) ] * Src[ N -2] + this->_pData[this->map_index ( 1, N-1 ) ] * Src[ N -1];
}


template < typename DataType >
/*virtual*/ void TriBandMatrix<DataType>::applyAdd ( const Vector<DataType> &Src, Vector<DataType> &Dst ) const {
  // might want to compare sizes
  const int N = this->getNumRows();
  Dst[0] += this->_pData[this->map_index ( 1,0 ) ] * Src[0] + this->_pData[this->map_index ( 2,0 ) ] * Src[1];
  for ( int i = 1; i <  N  - 1; ++i ) {
    Dst[i] += this->_pData[this->map_index ( 0,i ) ] * Src [i-1] + this->_pData[this->map_index ( 1,i ) ] * Src[i] + this->_pData[this->map_index ( 2,i ) ] * Src[i+1];
  }
  Dst[ N -1] += this->_pData[this->map_index ( 0, N -1 ) ] * Src[ N -2] + this->_pData[this->map_index ( 1, N -1 ) ] * Src[ N -1];
}


template class TriBandMatrix<float>;
template class TriBandMatrix<double>;
template class TriBandMatrix<long double>;



template<typename DataType>
/*explicit*/ LQuadBandMatrix<DataType>::LQuadBandMatrix ( const int N ) : GenBandMatrix<DataType> ( N, N ) {
  aol::Vector< int > offsets ( 4 );
  offsets[0] = -2;
  offsets[1] = -1;
  offsets[2] =  0;
  offsets[3] =  1;

  if ( N < 4 ) {
    throw Exception ( "LQuadBandMatrix: size too small. Will not work "
                      "and would not be useful.", __FILE__, __LINE__ );
  } else {
    this->reallocate ( N, N, offsets );
  }
}


template<typename DataType>
/*virtual*/ void LQuadBandMatrix<DataType>::apply ( const Vector<DataType> &Src, Vector<DataType> &Dst ) const {
  // might want to compare sizes
  const int N = this->getNumRows();
  Dst[0] = this->_pData[this->map_index ( 2,0 ) ] * Src[0] + this->_pData[this->map_index ( 3,0 ) ] * Src[1];
  Dst[1] = this->_pData[this->map_index ( 1,1 ) ] * Src[0] + this->_pData[this->map_index ( 2,1 ) ] * Src[1] + this->_pData[this->map_index ( 3,1 ) ] * Src[2];
  for ( int i = 2; i < N - 1; ++i ) {
    Dst[i] = this->_pData[this->map_index ( 0,i ) ] * Src [i-2] + this->_pData[this->map_index ( 1,i ) ] * Src[i-1] + this->_pData[this->map_index ( 2,i ) ] * Src[i] + this->_pData[this->map_index ( 3,i ) ] * Src[i+1];
  }
  Dst[N-1] = this->_pData[this->map_index ( 0,N-1 ) ] * Src[N-3] + this->_pData[this->map_index ( 1,N-1 ) ] * Src[N-2] + this->_pData[this->map_index ( 2,N-1 ) ] * Src[N-1];
}


template<typename DataType>
/*virtual*/ void LQuadBandMatrix<DataType>::applyAdd ( const Vector<DataType> &Src, Vector<DataType> &Dst ) const {
  // might want to compare sizes
  const int N = this->getNumRows();
  Dst[0] += this->_pData[this->map_index ( 2,0 ) ] * Src[0] + this->_pData[this->map_index ( 3,0 ) ] * Src[1];
  Dst[1] += this->_pData[this->map_index ( 1,1 ) ] * Src[0] + this->_pData[this->map_index ( 2,1 ) ] * Src[1] + this->_pData[this->map_index ( 3,1 ) ] * Src[2];
  for ( int i = 2; i < N - 1; ++i ) {
    Dst[i] += this->_pData[this->map_index ( 0,i ) ] * Src [i-2] + this->_pData[this->map_index ( 1,i ) ] * Src[i-1] + this->_pData[this->map_index ( 2,i ) ] * Src[i] + this->_pData[this->map_index ( 3,i ) ] * Src[i+1];
  }
  Dst[N-1] += this->_pData[this->map_index ( 0,N-1 ) ] * Src[N-3] + this->_pData[this->map_index ( 1,N-1 ) ] * Src[N-2] + this->_pData[this->map_index ( 2,N-1 ) ] * Src[N-1];
}


template class LQuadBandMatrix<float>;
template class LQuadBandMatrix<double>;
template class LQuadBandMatrix<long double>;

template<typename DataType>
/*explicit*/ RQuadBandMatrix<DataType>::RQuadBandMatrix ( const int N ) : GenBandMatrix<DataType> ( N, N ) {
  aol::Vector< int > offsets ( 4 );
  offsets[0] = -1;
  offsets[1] =  0;
  offsets[2] =  1;
  offsets[3] =  2;

  if ( N < 4 ) {
    throw Exception ( "RQuadBandMatrix: size too small. Will not work "
                      "and would not be useful.", __FILE__, __LINE__ );
  } else {
    this->reallocate ( N, N, offsets );
  }
}


template<typename DataType>
/*virtual*/ void RQuadBandMatrix<DataType>::apply ( const Vector<DataType> &Src, Vector<DataType> &Dst ) const {
  // might want to compare sizes
  const int N = this->getNumRows();
  Dst[0] = this->_pData[this->map_index ( 1,0 ) ] * Src[0] + this->_pData[this->map_index ( 2,0 ) ] * Src[1] + this->_pData[this->map_index ( 3,0 ) ] * Src[2];
  for ( int i = 1; i < N - 2; ++i ) {
    Dst[i] = this->_pData[this->map_index ( 0,i ) ] * Src [i-1] + this->_pData[this->map_index ( 1,i ) ] * Src[i] + this->_pData[this->map_index ( 2,i ) ] * Src[i+1] + this->_pData[this->map_index ( 3,i ) ] * Src[i+2];
  }
  Dst[N-2] = this->_pData[this->map_index ( 0,N-2 ) ] * Src[N-3] + this->_pData[this->map_index ( 1,N-2 ) ] * Src[N-2] + this->_pData[this->map_index ( 2,N-2 ) ] * Src[N-1];
  Dst[N-1] = this->_pData[this->map_index ( 0,N-1 ) ] * Src[N-2] + this->_pData[this->map_index ( 1,N-1 ) ] * Src[N-1];
}


template<typename DataType>
/*virtual*/ void RQuadBandMatrix<DataType>::applyAdd ( const Vector<DataType> &Src, Vector<DataType> &Dst ) const {
  // might want to compare sizes
  const int N = this->getNumRows();
  Dst[0] += this->_pData[this->map_index ( 1,0 ) ] * Src[0] + this->_pData[this->map_index ( 2,0 ) ] * Src[1] + this->_pData[this->map_index ( 3,0 ) ] * Src[2];
  for ( int i = 1; i < N - 2; ++i ) {
    Dst[i] += this->_pData[this->map_index ( 0,i ) ] * Src [i-1] + this->_pData[this->map_index ( 1,i ) ] * Src[i] + this->_pData[this->map_index ( 2,i ) ] * Src[i+1] + this->_pData[this->map_index ( 3,i ) ] * Src[i+2];
  }
  Dst[N-2] += this->_pData[this->map_index ( 0,N-2 ) ] * Src[N-3] + this->_pData[this->map_index ( 1,N-2 ) ] * Src[N-2] + this->_pData[this->map_index ( 2,N-2 ) ] * Src[N-1];
  Dst[N-1] += this->_pData[this->map_index ( 0,N-1 ) ] * Src[N-2] + this->_pData[this->map_index ( 1,N-1 ) ] * Src[N-1];
}


template class RQuadBandMatrix<float>;
template class RQuadBandMatrix<double>;
template class RQuadBandMatrix<long double>;

}
