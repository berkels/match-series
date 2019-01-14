#include <multiVector.h>
#include <gridBase.h>
#include <bzipiostream.h>

template < typename _DataType >
void aol::MultiVector< _DataType >::load ( const char *FileName ) {
  aol::Bzipifstream in ( FileName );
  int value;
  in >> value;
  reallocate ( value, 0 );
  for ( int c = 0; c < numComponents(); ++c ) {
    in >> value;
    this->operator[]( c ).reallocate ( value );
  }
  in.ignore();

  for ( int c = 0; c < numComponents(); ++c ) {
    this->operator[]( c ).loadRaw ( in, qc::PGM_DOUBLE_BINARY );
  }
}

template < typename _DataType >
void aol::MultiVector< _DataType >::save ( const char *FileName ) const {
  aol::Bzipofstream out ( FileName );
  out << numComponents();
  for ( int c = 0; c < numComponents(); ++c )
    out << " " << this->operator[]( c ).size();
  out << endl;
  for ( int c = 0; c < numComponents(); ++c )
    this->operator[]( c ).saveRaw ( out, qc::PGM_DOUBLE_BINARY, this->operator[]( c ).getMinValue(), this->operator[]( c ).getMaxValue() );
}

template < typename _DataType >
void aol::MultiVector< _DataType >::reallocate ( const qc::GridStructure &grid ) {
  this->reallocate ( grid.getDimOfWorld(), grid.getNumberOfNodes() );
}

template <typename _DataType>
typename aol::MultiVector<_DataType>::RealType aol::MultiVector<_DataType>::lpNorm ( RealType p ) const {
  RealType ret = aol::NumberTrait<RealType>::zero;
  for ( int c = 0; c < this->numComponents(); ++c ) {
    ret += this->operator[]( c ).lpNormPowP ( p );
  }

  if ( p == aol::ZOTrait<RealType>::zero ) {
    return static_cast<RealType> ( this->getTotalSize() );
  } else if ( p == aol::ZOTrait<RealType>::one ) {
    return static_cast<RealType> ( ret );
  } else if ( p == static_cast<RealType> ( 2 ) ) {
    return static_cast<RealType> ( sqrt ( ret ) );
  } else                     {
    return static_cast<RealType> ( pow ( ret, static_cast<RealType> ( 1.0 / p ) ) );
  }
}

template <typename _DataType>
void aol::MultiVector<_DataType>::copySplitFrom ( const Vector<DataType> & singleVector ) {
  QUOC_ASSERT ( singleVector.size() == getTotalSize() );
  Vector<int> sizes;
  getSizes(sizes);
  int index = 0;
  for (int i = 0; i < sizes.size(); ++i)
    for (int j = 0; j < sizes[i]; ++j, ++index)
      (*this)[i][j] = singleVector[index];
}

template <typename _DataType>
void aol::MultiVector<_DataType>::copySplitFrom ( const Vector<DataType> & singleVector, const Vector<int> & sizes ) {
  int k = 0;
  reallocate ( 0, 0 );
  for ( int i = 0; i < sizes.size (); ++i ) {
    vecs.push_back ( vec_entry ( new Vector<DataType> ( sizes [i] ) ) );
    for ( int j = 0; j < sizes [i]; ++j, ++k )
      (*this) [i][j] = singleVector [k];
  }
  int rest = singleVector.size() - k;
  if ( rest > 0 ) {
    Vector<DataType>* temp = new Vector<DataType> ( rest );
    for ( int i = 0; i < rest; ++i, ++k )
      ( *temp ) [i] = singleVector [k];
    appendReference ( *temp );
  }
}

template <typename _DataType>
void aol::MultiVector<_DataType>::copySplitTransposeFrom ( const Vector<DataType> & singleVector ) {
  QUOC_ASSERT ( singleVector.size() == getTotalSize() );
  Vector<int> sizes;
  getSizes(sizes);
  if (sizes.getMinValue() != sizes.getMaxValue())
    throw aol::Exception ( "aol::MultiVector::copySplitTransposeFrom only works for component vectors of equal length.", __FILE__, __LINE__ );
  int index = 0;
  for (int i = 0; i < sizes[0]; ++i)
    for (int j = 0; j < sizes.size(); ++j, ++index)
      (*this)[j][i] = singleVector[index];
}

template <typename _DataType>
void aol::MultiVector<_DataType>::getMedianVecOverComponents ( aol::Vector<RealType> &MedianVec, const RealType FillValue ) {
  if ( MedianVec.size() != getEqualComponentSize() )
    throw aol::Exception ( "aol::MultiVector::getMedianVecOverComponents: Size mismatch.", __FILE__, __LINE__ );

  const int numComps = numComponents();

  for ( int j = 0; j < MedianVec.size(); ++j ) {
    aol::Vector<_DataType> medianTemp ( 0 );
    medianTemp.reserve ( numComps );
    for ( int k = 0; k < numComps; ++k ) {
      if ( (*this)[k][j] != aol::NumberTrait<RealType>::Inf )
        medianTemp.pushBack ( (*this)[k][j] );
    }
    if ( medianTemp.size() > 0 )
      MedianVec[j] = medianTemp.getMedianValue();
    else
      MedianVec[j] = FillValue;
  }
}

template <typename _DataType>
void aol::MultiVector<_DataType>::saveToFile ( const char *filename ) const {
  aol::Bzipofstream out ( filename );
  out << aol::VectorFileMagicChar::MultiVector << FileFormatMagicNumber<DataType>::FFType << endl;
  out << "# This is a QuOcMesh file of type " << aol::VectorFileMagicChar::MultiVector << FileFormatMagicNumber<DataType>::FFType << " storing an aol::MultiVector<" << FileFormatMagicNumber<DataType>::FFContainedName << ">" << endl;
  const int numComp = this->numComponents();
  out << numComp << endl;
  for ( int i = 0; i < numComp; ++i ) {
    out << this->operator[](i).size();
    if ( i < ( numComp - 1 ) )
      out << " ";
  }
  out << endl;
  for ( int i = 0; i < numComp; ++i ) {
    const char* buffer = reinterpret_cast<char*> ( this->operator[](i).getData() );
    out.write ( buffer, this->operator[](i).size() * sizeof ( DataType ) );
  }
}


template <typename _DataType>
void aol::MultiVector<_DataType>::loadFromFile ( const char *filename ) {
  aol::Bzipifstream reader ( filename );
  char M = 0;
  int ident = 0;
  reader >> M;
  reader >> ident;
  if ( ( M != aol::VectorFileMagicChar::MultiVector ) || ( ident != FileFormatMagicNumber<DataType>::FFType ) ) {
    cerr << M << ident << ", should be " << aol::VectorFileMagicChar::MultiVector << FileFormatMagicNumber<DataType>::FFType << " (" << FileFormatMagicNumber<DataType>::FFContainedName << ")" << endl;
    throw aol::Exception ( "Illegal magic number for aol::MultiVector", __FILE__, __LINE__ );
  }
  aol::READ_COMMENTS ( reader );
  int numComp;
  reader >> numComp;
  char buffer[1024];
  reader.getline ( buffer, 1 );
  aol::Vector<int> sizes ( numComp );

  for ( int i = 0; i < numComp; ++i ) {
    reader >> sizes[i];
  }
  reader.getline ( buffer, 1 );

  this->reallocate ( sizes );

  for ( int i = 0; i < numComp; ++i ) {
    reader.read ( reinterpret_cast<char*>( this->operator[](i).getData() ), sizes[i] * sizeof( DataType ) );
  }
}


template class aol::MultiVector<signed char>;
template class aol::MultiVector<unsigned char>;
template class aol::MultiVector<short>;
template class aol::MultiVector<unsigned short>;
template class aol::MultiVector<int>;
template class aol::MultiVector<unsigned int>;
template class aol::MultiVector<int64_t>;
template class aol::MultiVector<uint64_t>;
template class aol::MultiVector<float>;
template class aol::MultiVector<double>;
template class aol::MultiVector<long double>;
