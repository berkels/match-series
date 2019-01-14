#include <aol.h>
#include <vec.h>
#include <qmException.h>
#include <gridBase.h>
#include <multiVector.h>
#include <bzipiostream.h>

namespace aol {
template <class T> const T aol::OverflowTrait<T>::max = static_cast<T> ( 255 );

template <> const signed char aol::OverflowTrait<signed char>::max = 127;

template struct aol::OverflowTrait<float>;
template struct aol::OverflowTrait<double>;
template struct aol::OverflowTrait<long double>;

template struct aol::OverflowTrait<unsigned char>;
template struct aol::OverflowTrait<signed char>;
template struct aol::OverflowTrait<unsigned short>;
template struct aol::OverflowTrait<short>;
template struct aol::OverflowTrait<unsigned int>;
template struct aol::OverflowTrait<int>;
template struct aol::OverflowTrait<uint64_t>; // this is unsigned long int, at least on some platforms
template struct aol::OverflowTrait<int64_t>;
#if ( defined __APPLE__ ) || ( BITNESS == 32 ) // This seems to be the only way to get SuiteSparse_long to work under OS X
template struct aol::OverflowTrait<long>;
#endif
}

namespace {

template <typename DataType>
DataType ScalarProduct ( const DataType *a, const DataType *b, const int N ) {
  DataType ret = 0;
  for ( int i = 0 ; i < N ; ++i ) {
    ret += a[i] * b[i];
  }
  return ret;
}

#ifdef USE_SSE
template <>
float ScalarProduct ( const float *a, const float *b, const int N ) {
  float ret = 0;

  const int nLoop = N / 4;

  float *temp = static_cast< float* > ( aol::aligned_memory_allocation ( 4 * sizeof ( float ), 16 ) );

  __m128 m1;
  __m128* m2 = reinterpret_cast< __m128* > ( temp );
  const __m128* pSrc1 = reinterpret_cast< const __m128* > ( a );
  const __m128* pSrc2 = reinterpret_cast< const __m128* > ( b );
  *m2 = _mm_set_ps1 ( 0.0f ); // set all four values of m2 to 0.
  for ( int i = 0; i < nLoop; ++i ) {
    // This code is easy to understand:
    m1 = _mm_mul_ps ( *pSrc1, *pSrc2 );       // *pDest = *pSrc1 * *pSrc2
    *m2 = _mm_add_ps ( m1, *m2 );

    pSrc1++;
    pSrc2++;
  }
  ret = ( temp[0] + temp[1] + temp[2] + temp[3] );
  for ( int i = 4 * nLoop ; i < N ; ++i ) {
    ret += a[i] * b[i];
  }
  aol::aligned_memory_deallocation ( temp );
  return ret;
}

template <>
double ScalarProduct ( const double *a, const double *b, const int N ) {
  double ret = 0;

  const int nLoop = N / 2;

  double *temp = static_cast< double* > ( aol::aligned_memory_allocation ( 2 * sizeof ( double ), 16 ) );
  __m128d m1;
  __m128d* m2 = reinterpret_cast< __m128d* > ( temp );
  const __m128d* pSrc1 = reinterpret_cast< const __m128d* > ( a );
  const __m128d* pSrc2 = reinterpret_cast< const __m128d* > ( b );
  *m2 = _mm_setzero_pd ();  // set both values of m2 to 0.
  for ( int i = 0; i < nLoop; ++i ) {
    // This code is easy to understand:
    m1 = _mm_mul_pd ( *pSrc1, *pSrc2 );       // *pDest = *pSrc1 * *pSrc2
    *m2 = _mm_add_pd ( m1, *m2 );

    pSrc1++;
    pSrc2++;
  }
  ret = ( temp[0] + temp[1] );
  for ( int i = 2 * nLoop ; i < N ; ++i ) {
    ret += a[i] * b[i];
  }
  aol::aligned_memory_deallocation ( temp );
  return ret;
}
#endif

} // end namespace


template <typename _DataType>
void aol::Vector<_DataType>::reserve ( int Length ) {
  if ( Length > _sizeReserved ) {
    DataType* newdata = static_cast<DataType*> ( aol::MemoryManager::allocateAtLeast ( Length, sizeof ( DataType ) ) );
    memset ( newdata, 0, sizeof ( DataType ) * Length );

    if ( _pData ) {
      memcpy ( newdata, _pData, sizeof ( DataType ) * _size );
      // clear remaining part of newly allocated memory:
      if ( _deleteFlag )
        aol::MemoryManager::deallocate ( _pData, _sizeReserved, sizeof ( DataType ) );
    }

    _sizeReserved = Length;
    _pData = newdata;
    _deleteFlag = true;
  }
}

template <typename _DataType>
void aol::Vector<_DataType>::reserveUninit ( int Length ) {
  if ( Length > _sizeReserved ) {
    DataType* newdata = static_cast<DataType*> ( aol::MemoryManager::allocateAtLeast ( Length, sizeof ( DataType ) ) );

    if ( _pData && _deleteFlag ) {
      aol::MemoryManager::deallocate ( _pData, _sizeReserved, sizeof ( DataType ) );
    }

    _sizeReserved = Length;
    _pData = newdata;
    _deleteFlag = true;
  }
}


template <typename _DataType>
void aol::Vector<_DataType>::resize ( const int Length ) {
  if ( Length > _sizeReserved ) {
    reserve ( Length );
    _sizeReserved = _size = Length;
  } else {
    for ( int i = _size; i < Length; ++i ) {
      _pData[i] = aol::ZTrait<DataType>::zero;
    }
    _size = Length;
  }
}

template <typename _DataType>
void aol::Vector<_DataType>::pushBack ( const _DataType Data ) {
  if ( _size + 1 > _sizeReserved ) {
    if ( _size == 0 ) reserve ( 1 );
    else reserve ( 2*_size );
  }
  _pData[_size] = Data;
  _size += 1;
}

template <typename _DataType>
void aol::Vector<_DataType>::pushBackValues ( const aol::Vector<_DataType> &otherVector ) {
  const int currentSize = this->size();
  this->growBy ( otherVector.size() );
  for ( int i = 0; i < otherVector.size(); ++i ) {
    this->set ( currentSize + i, otherVector.get ( i ) );
  }
}


template <typename _DataType>
void aol::Vector<_DataType>::reallocate ( const int Length ) {
  if ( ( _deleteFlag == false ) && ( _pData != NULL ) )
    throw Exception ( "aol::Vector<DataType>::reallocate may not be called on FLAT_COPY vectors!\n", __FILE__, __LINE__ );

  if ( Length > _sizeReserved ) {
    reserve ( Length );
    _sizeReserved = _size = Length;
  } else {
    _size = Length;
  }
  setZero();
}

template< typename _DataType >
void aol::Vector<_DataType>::reallocate ( const qc::GridStructure &grid ) {
  this->reallocate ( grid.getNumberOfNodes() );
}

template< typename _DataType >
void aol::Vector<_DataType>::reallocateClear ( const int Length ) {
  if ( ( _deleteFlag == false ) && ( _pData != NULL ) )
    throw Exception ( "aol::Vector<DataType>::reallocateClear may not be called on FLAT_COPY vectors!\n", __FILE__, __LINE__ );

  aol::MemoryManager::deallocate ( _pData, _sizeReserved, sizeof ( DataType ) );
  _size = Length;
  _sizeReserved = Length;
  _pData = static_cast<DataType*> ( aol::MemoryManager::allocateAtLeast ( _sizeReserved, sizeof ( DataType ) ) );
  memset ( _pData, 0, _sizeReserved * sizeof ( DataType ) );
}

template< typename _DataType >
bool aol::Vector<_DataType>::createOverflowHandledData ( aol::Vector<unsigned char> &buffer,
                                                         const DataType min,
                                                         const DataType max ) const {

  const DataType range = overflowMax - overflowMin;
  int i;

  DataType  dataRange = max - min;

  switch ( overflowHandling ) {
    case aol::CLIP:
      if ( !quietMode ) {
        cerr << "Clipping data to range [";
        if ( sizeof ( unsigned char ) == sizeof ( DataType ) )
          cerr << static_cast< int > ( overflowMin ) << ", " << static_cast<int> ( overflowMax ) << "]" << endl;
        else
          cerr << overflowMin << ", " << overflowMax << "]" << endl;
      }
      for ( i = 0; i < this->_size; ++i ) {
        if ( this->_pData[i] < overflowMin ) {
          buffer[i] = static_cast< unsigned char > ( overflowMin );
        } else {
          if ( this->_pData[i] > overflowMax ) {
            buffer[i] = static_cast< unsigned char > ( overflowMax );
          } else {
            buffer[i] = static_cast< unsigned char > ( this->_pData[i] );
          }
        }
      }
      break;
    case aol::SCALE:
      if ( !quietMode )
        cerr << "Scaling data from range [" << min << ", " << max << "] to [0,255]\n";
      for ( i = 0; i < this->_size; ++i ) {
        buffer[i] = static_cast< unsigned char > ( ( this->_pData[i] - min ) / ( dataRange * 1.0 ) * 255 );
      }
      break;
    case aol::CLIP_THEN_SCALE: {
      if ( !quietMode )
        cerr << "Clipping data to [" << overflowMin << ", " << overflowMax << "] and then scaling to [0,255]\n";
      for ( i = 0; i < this->_size; ++i ) {
        DataType tmp;
        if ( this->_pData[i] < overflowMin )
          tmp = overflowMin;
        else if ( this->_pData[i] > overflowMax )
          tmp = overflowMax;
        else tmp = this->_pData[i];
        buffer[i] = static_cast< unsigned char > ( ( tmp - overflowMin ) / ( 1.0 * range ) * 255 );
      }
    }
    break;
    case aol::REFLECT: {
      if ( !quietMode ) cerr << "Reflecting data at range [" << overflowMin << ", " << overflowMax << "]\n";
      DataType tmp;
      const DataType overflowMin2 = 2 * overflowMin;
      const DataType overflowMax2 = 2 * overflowMax;
      for ( i = 0; i < this->_size; ++i ) {
        tmp = this->_pData[i];
        while ( tmp < overflowMin ) tmp = overflowMin2 - this->_pData[i];
        while ( tmp > overflowMax ) tmp = overflowMax2 - this->_pData[i];
        buffer[i] = static_cast< unsigned char > ( tmp );
      }
    }
    break;
    default:
      ;
  }
  return false;
}

template< typename _DataType >
void aol::Vector<_DataType>::saveRaw ( const char *filename, const qc::SaveType type, const DataType Minimum, const DataType Maximum ) const {
  if ( !quietMode ) cerr << "Writing to uncompressed stream without header...";

  aol::Bzipofstream *out = new aol::Bzipofstream ( filename );

  if ( !out->good() ) {
    delete out;
    throw aol::Exception ( "aol::Vector<DataType>::saveRaw: Cannot open file for writing", __FILE__, __LINE__ );
  }

  saveRaw ( *out, type, Minimum, Maximum );

  delete out;

  if ( !quietMode ) cerr << "done.\n";
}


template< typename _DataType >
void aol::Vector<_DataType>::saveRaw ( ostream &out, const qc::SaveType type, const DataType Minimum, const DataType Maximum ) const {
  int i = 0;

  DataType maximum = Maximum;

  if ( ( type == qc::PGM_UNSIGNED_CHAR_BINARY || type == qc::PGM_UNSIGNED_SHORT_BINARY ) && Maximum == 0 ) {
    cerr << "Cannot scale image to a maximum of 0, scaling to 1 instead.\n";
    // This is necessary to write complete black images
    maximum = 1;
  }

  switch ( type ) {
    case qc::PGM_UNSIGNED_CHAR_ASCII: {
      aol::Vector<unsigned char> buffer ( this->_size );
      bool flag = this->createOverflowHandledData ( buffer, Minimum, maximum );
      for ( i = 0; i < this->_size; ++i ) {
        out <<  static_cast<int> ( buffer[i] ) << " ";
        if ( ( i % 10 ) == 9 ) out << std::endl;
      }
      if ( flag && !quietMode ) cerr << "warning: Data contains values < 0 that are converted to unsigned char\n";
    }
    break;

    case qc::PGM_UNSIGNED_CHAR_BINARY: {
      aol::Vector<unsigned char> buffer ( this->_size );
      bool flag = this->createOverflowHandledData ( buffer, Minimum, maximum );
      out.write ( reinterpret_cast<char*> ( buffer.getData() ), this->_size );
      if ( flag && !quietMode ) cerr << "warning: Data contains values < 0 that are converted to unsigned char\n";
    }
    break;

    case qc::PGM_FLOAT_ASCII: {
      for ( i = 0; i < this->_size; ++i ) {
        out << this->_pData[i] << " ";
        if ( ( i % 10 ) == 9 ) out << std::endl;
      }
    }
    break;

    case qc::PGM_FLOAT_BINARY: {
      if ( ( sizeof ( float ) == sizeof ( DataType ) ) && !(std::numeric_limits<DataType>::is_integer) ) {
        out.write ( reinterpret_cast<char *> ( this->_pData ), this->_size * sizeof ( float ) );
      } else {
        float * buffer = new float[this->_size];
        for ( i = 0; i < this->_size; ++i ) buffer[i] =  static_cast<float> ( this->_pData[i] );
        out.write ( reinterpret_cast<char *> ( buffer ), this->_size * sizeof ( float ) );
        delete[] buffer;
      }
    }
    break;

    case qc::PGM_DOUBLE_BINARY: {
      if ( ( sizeof ( double ) == sizeof ( DataType ) ) && !(std::numeric_limits<DataType>::is_integer) ) {
        out.write ( reinterpret_cast<char *> ( this->_pData ), this->_size * sizeof ( double ) );
      } else {
        double * buffer = new double[this->_size];
        for ( i = 0; i < this->_size; ++i ) buffer[i] =  static_cast<double> ( this->_pData[i] );
        out.write ( reinterpret_cast<char *> ( buffer ), this->_size * sizeof ( double ) );
        delete[] buffer;
      }
    }
    break;

    case qc::PGM_UNSIGNED_SHORT_BINARY: {
      bool flag = false;
      unsigned short *buffer = new unsigned short[this->_size];
      for ( i = 0; i < this->_size; ++i ) {
        DataType temp = this->_pData[i];
        if ( static_cast<double> ( temp ) < 0 ) {
          flag = true;
          // Workaround for icc (replaced - with abs)
          temp = aol::Abs ( this->_pData[i] );
        }
        buffer[i] =  static_cast<unsigned short> ( static_cast<double> ( temp ) / static_cast<double> ( maximum ) * 65535.0 );
      }
      out.write ( reinterpret_cast<char*> ( buffer ), this->_size * sizeof ( unsigned short ) );
      delete[] buffer;
      if ( flag && !this->quietMode ) cerr << "warning: Data contains values < 0 that are converted to unsigned short\n";
    }
    break;

    default:
      throw aol::TypeException ( "aol::Vector<DataType>::saveRaw: invalid type specified", __FILE__, __LINE__ );
  }


  if ( out.good() ) {
    if ( !quietMode )
      cerr << "Successfully wrote to PGM " << type << " stream\n";
  } else {
    const aol::Bzipofstream* pOut = dynamic_cast<const aol::Bzipofstream*>(&out);
    if( pOut != NULL )
      throw aol::Exception ( "aol::Vector<DataType>::saveRaw: Writing to aol::Bzipofstream failed. Is there enough free memory to hold the whole file?", __FILE__, __LINE__ );
    else
      throw aol::Exception ( "aol::Vector<DataType>::saveRaw: Error writing file", __FILE__, __LINE__ );
  }

}

template <typename _DataType>
void aol::Vector<_DataType>::loadRaw ( istream &in, const int Type ) {

  switch ( Type ) {
    case qc::PGM_UNSIGNED_CHAR_ASCII:
      for ( int i = 0; i < this->_size; ++i ) {
        int value;
        in >> value;
        if ( value < 0 || value > 255 )
          throw aol::Exception ( "aol::Vector<DataType>::loadRaw: Illegal ascii number for unsigned char", __FILE__, __LINE__ );
        this->_pData[i] = static_cast<unsigned char> ( value );
      }
      break;
    case qc::PGM_FLOAT_ASCII:
      for ( int i = 0; i < this->_size; ++i ) {
        in >> this->_pData[i];
      }
      break;
    case qc::PGM_UNSIGNED_CHAR_BINARY: {
      unsigned char *buffer = new unsigned char[this->_size];
      in.read ( reinterpret_cast<char*> ( buffer ),  this->_size );
      for ( int i = 0; i < this->_size; ++i ) this->_pData[i] = static_cast<DataType> ( static_cast<unsigned char> ( buffer[i] ) );
      delete[] buffer;
    }
    break;
    case qc::PGM_FLOAT_BINARY:
      if ( ( sizeof ( float ) == sizeof ( DataType ) ) && !(std::numeric_limits<DataType>::is_integer) ) {
        in.read ( reinterpret_cast<char*> ( this->_pData ), this->_size * sizeof ( float ) );
      } else {
        float *buffer = new float[ this->_size];
        in.read ( reinterpret_cast<char*> ( buffer ),  this->_size * sizeof ( float ) );
        for ( int i = 0; i < this->_size; ++i ) this->_pData[i] =  static_cast<DataType> ( buffer[i] );
        delete[] buffer;
      }
      break;
    case qc::PGM_DOUBLE_BINARY:
      if ( ( sizeof ( double ) == sizeof ( DataType ) ) && !(std::numeric_limits<DataType>::is_integer) ) {
        in.read ( reinterpret_cast<char*> ( this->_pData ),  this->_size * sizeof ( double ) );
      } else {
        double *buffer = new double[ this->_size];
        in.read ( reinterpret_cast<char*> ( buffer ),  this->_size * sizeof ( double ) );
        for ( int i = 0; i < this->_size; ++i ) this->_pData[i] =  static_cast<DataType> ( buffer[i] );
        delete[] buffer;
      }
      break;
    case qc::PGM_UNSIGNED_SHORT_BINARY_BIGENDIAN:
    case qc::PGM_UNSIGNED_SHORT_BINARY: {
      unsigned short *buffer = new unsigned short[ this->_size];
      in.read ( reinterpret_cast<char*> ( buffer ),  this->_size * sizeof ( unsigned short ) );
      if ( Type == qc::PGM_UNSIGNED_SHORT_BINARY_BIGENDIAN ) {
        aol::Vector<unsigned short> bufferVec ( buffer, this->_size, aol::FLAT_COPY );
        bufferVec.swapByteOrder();
      }
      for ( int i = 0; i < this->_size; ++i ) this->_pData[i] = static_cast<DataType> ( static_cast<unsigned short> ( buffer[i] ) );
      delete[] buffer;
      break;
    }
    case qc::PGM_SHORT_BINARY: {
      short *buffer = new short[ this->_size];
      in.read ( reinterpret_cast<char*> ( buffer ),  this->_size * sizeof ( short ) );
      for ( int i = 0; i < this->_size; ++i ) this->_pData[i] = static_cast<DataType> ( buffer[i] );
      delete[] buffer;
      break;
    }
    case qc::PGM_UNSIGNED_INT_BINARY: {
      aol::readBinaryData<uint32_t, DataType> ( in, this->_pData, this->_size );
      break;
    }
    case qc::PGM_SIGNED_INT_BINARY: {
      aol::readBinaryData<int32_t, DataType> ( in, this->_pData, this->_size );
      break;
    }
    default:
      throw aol::TypeException ( "Illegal PGM_TYPE", __FILE__, __LINE__ );
  }
  if ( in.fail() )
    throw aol::FileException ( "aol::Vector<DataType>::loadRaw: Reading from istream failed", __FILE__, __LINE__ );
}


template <typename _DataType>
void aol::Vector<_DataType>::saveToFile ( const char *filename ) const {
  aol::Bzipofstream out ( filename );
  out << aol::VectorFileMagicChar::Vector << FileFormatMagicNumber<DataType>::FFType << endl;
  out << "# This is a QuOcMesh file of type " << aol::VectorFileMagicChar::Vector << FileFormatMagicNumber<DataType>::FFType << " storing an aol::Vector<" << FileFormatMagicNumber<DataType>::FFContainedName << ">" << endl;
  out << this->size() << endl;
  const char* buffer = reinterpret_cast<char*> ( this->getData() );
  out.write ( buffer, this->size() * sizeof ( DataType ) );
}


template <typename _DataType>
void aol::Vector<_DataType>::loadFromFile ( const char *filename ) {
  aol::Bzipifstream reader ( filename );
  char M = 0;
  int ident = 0;
  reader >> M;
  reader >> ident;
  if ( ( M != aol::VectorFileMagicChar::Vector ) || ( ident != FileFormatMagicNumber<DataType>::FFType ) ) {
    cerr << M << ident << ", should be " << aol::VectorFileMagicChar::Vector << FileFormatMagicNumber<DataType>::FFType << " (" << FileFormatMagicNumber<DataType>::FFContainedName << ")" << endl;
    throw aol::Exception ( "Illegal magic number for aol::Vector", __FILE__, __LINE__ );
  }
  aol::READ_COMMENTS ( reader );
  int size;
  reader >> size;
  char buffer[1024];
  reader.getline ( buffer, 1 );

  this->reallocate ( size );
  reader.read ( reinterpret_cast<char*>( this->getData() ), size * sizeof( DataType ) );
}




template <typename _DataType>
aol::Vector<_DataType>::Vector ( _DataType* Data, int n, CopyFlag copyFlag )
    : _size ( n ), _sizeReserved ( n ), _deleteFlag ( false ), _pData ( NULL ),
    overflowHandling ( aol::CLIP ), overflowMin ( 0 ), overflowMax ( OverflowTrait<DataType>::max ) {

  switch ( copyFlag ) {
    case FLAT_COPY:
#ifdef USE_SSE
      if ( ( reinterpret_cast<const uintptr_t> ( Data ) % 16 ) != 0 )
        throw Exception ( "aol::Vector<DataType>::Vector( DataType* Data, int n, CopyFlag copyFlag ): With USE_SSE the data pointer of Vector must be 16 byte aligned.\n", __FILE__, __LINE__ );
#endif
      _deleteFlag = false;
      _pData = Data;
      break;

    case DEEP_COPY:
      _pData = static_cast<DataType*> ( aol::MemoryManager::allocateAtLeast ( _sizeReserved, sizeof ( DataType ) ) );
      _deleteFlag = true;
      memcpy ( _pData, Data, _size * sizeof ( DataType ) );
      break;

    default:
      throw aol::Exception ( "aol::Vector<DataType>::Vector( DataType* Data, int n, CopyFlag copyFlag ): illegal copy flag", __FILE__, __LINE__ );
  }
}


template <typename _DataType>
aol::Vector<_DataType>::Vector ( const aol::MultiVector<_DataType>& mv )
    : _deleteFlag ( true ),
    overflowHandling ( aol::CLIP ), overflowMin ( 0 ), overflowMax ( aol::OverflowTrait<DataType>::max ) {
  int i, j, k;

  _size = _sizeReserved = mv.getTotalSize();
  _pData = static_cast<DataType*> ( aol::MemoryManager::allocateAtLeast ( _sizeReserved, sizeof ( DataType ) ) );

  for ( i = 0, k = 0; i < mv.numComponents (); ++i ) {
    const aol::Vector<DataType>& v = mv [i];
    for ( j = 0; j < v.size (); ++j, ++k )
      _pData [k] = v [j];
  }
}

template <typename _DataType>
aol::Vector<_DataType>::Vector ( const std::vector<DataType> & std_vec )
    : _size ( static_cast<int> ( std_vec.size() ) ), _sizeReserved ( static_cast<int> ( std_vec.size() ) ), _deleteFlag ( true ), _pData ( NULL ),
    overflowHandling ( aol::CLIP ), overflowMin ( 0 ), overflowMax ( aol::OverflowTrait<DataType>::max ) {
  _pData = static_cast<DataType*> ( aol::MemoryManager::allocateAtLeast ( _sizeReserved, sizeof ( DataType ) ) );
  for ( int i = 0; i < _size; ++i )
    _pData[i] = std_vec[static_cast<unsigned int> ( i ) ];
}


template <typename _DataType>
aol::Vector<_DataType>::Vector ( aol::Vector<_DataType> const&Vec, CopyFlag copyFlag )
    : Obj ( Vec ), _size ( Vec._size ), _sizeReserved ( Vec._sizeReserved ),
    overflowHandling ( Vec.overflowHandling ), overflowMin ( Vec.overflowMin ), overflowMax ( Vec.overflowMax ) {
  switch ( copyFlag ) {
    case DEEP_COPY:
      _pData = static_cast<DataType*> ( aol::MemoryManager::allocateAtLeast ( _sizeReserved, sizeof ( DataType ) ) );
      _deleteFlag = true;
      memcpy ( _pData, Vec._pData, _size*sizeof ( DataType ) );
      break;

    case FLAT_COPY:
      _pData = Vec._pData;
      _deleteFlag = false;
      //     _size = Vec._size; // this happens in the initialization list.
      //     _sizeReserved = Vec._sizeReserved;
      break;

    case STRUCT_COPY:
      _pData = static_cast<DataType*> ( aol::MemoryManager::allocateAtLeast ( _sizeReserved, sizeof ( DataType ) ) );
      _deleteFlag = true;
      setZero(); // same as in standard constructor: by convention, vectors are initialized with zero.
      break;

    case STRUCT_COPY_UNINIT:
      _pData = static_cast<DataType*> ( aol::MemoryManager::allocateAtLeast ( _sizeReserved, sizeof ( DataType ) ) );
      _deleteFlag = true;
      // do not set to zero
      break;

    default:
      throw aol::Exception ( "Invalid CopyFlag specified", __FILE__, __LINE__ );
      break;
  };
}


template<typename _DataType>
aol::Vector<_DataType>::Vector ( const qc::GridStructure &Grid ) : _size ( Grid.getNumberOfNodes() ), _sizeReserved ( _size ), _deleteFlag ( true ), overflowHandling ( aol::CLIP ), overflowMin ( 0 ), overflowMax ( aol::OverflowTrait<DataType>::max ) {
  _pData = static_cast<DataType*> ( aol::MemoryManager::allocateAtLeast ( _sizeReserved, sizeof ( DataType ) ) );
  if ( !_pData )
    throw OutOfMemoryException ( "Vector<DataType,Realtype>::Vector: Could not allocate memory for Vector.", __FILE__, __LINE__ );

  setZero ();
}


template <typename _DataType>
void aol::Vector<_DataType>::absDiff ( aol::Vector<_DataType> const& Vec1,
                                       aol::Vector<_DataType> const& Vec2 ) {
  if ( Vec1.size() != _size || Vec2.size() != _size ) {
    throw aol::Exception ( "dimensions don't match", "aol::Vector::absDiff" );
  }
  for ( int i = 0; i < _size; ++i ) {
    _pData[ i ] = static_cast< DataType > ( aol::Abs ( Vec1._pData[i] - Vec2._pData[i] ) );
  }
}


template <typename _DataType>
_DataType aol::Vector<_DataType>::operator* ( const aol::Vector<_DataType> &c ) const {
  if ( c.size() != _size ) {
    throw aol::Exception ( "Vector::operator*: Vectorlengths not equal...", __FILE__, __LINE__ );
  } else {
    return ScalarProduct<DataType> ( this->_pData, c.getData(), _size );
  }
  return 0;
}

template <typename _DataType>
_DataType aol::Vector<_DataType>::dotProduct ( const aol::Vector<_DataType> &Vec ) const {
  return this->operator* ( Vec );
}

template <typename _DataType>
void aol::Vector<_DataType>::dump ( bool interactive ) const {
  cout << "( ";
  for ( int i = 0 ; i < _size ; ++i ) {
    printf ( "%f ", static_cast< float > ( _pData[i] ) );
    if ( interactive ) getchar();
  }
  cout << " )" << endl;
}

template <typename _DataType>
void aol::Vector<_DataType>::copyToBuffer ( _DataType* buffer ) const {
  if ( buffer != _pData ) memcpy ( buffer, _pData, sizeof ( DataType ) * _size );
}

template <typename _DataType>
void aol::Vector<_DataType>::readFromBuffer ( const _DataType* const buffer ) {
  if ( buffer != _pData ) memcpy ( _pData, buffer, sizeof ( DataType ) * _size );
}

template <typename _DataType>
aol::Vector<_DataType>& aol::Vector<_DataType>::operator= ( const aol::Vector<_DataType> &Vec ) {
  if ( Vec.size() == _size ) {
    Vec.copyToBuffer ( _pData );

    overflowHandling = Vec.overflowHandling;
    overflowMin = Vec.overflowMin;
    overflowMax = Vec.overflowMax;
  } else
    throw aol::Exception ( "aol::Vector<DataType>::operator= : "
                           "incompatible vector lengths!", __FILE__, __LINE__ );
  return *this;
}

template <typename _DataType>
aol::Vector<_DataType>& aol::Vector<_DataType>::assignFrom ( const aol::BitVector & bitField ) {
  if ( bitField.size() == _size )
    for ( int i = 0; i < _size; ++i )
      ( *this ) [i] = ( bitField[i] ? ZOTrait<DataType>::one : ZOTrait<DataType>::zero );
  else
    throw aol::Exception ( "aol::Vector<DataType>::operator= (const BitVector &):"
                           " incompatible vector lengths!", __FILE__, __LINE__ );
  return *this;
}

template <typename _DataType>
aol::Vector<_DataType> & aol::Vector<_DataType>::copyUnblockedFrom ( const MultiVector<_DataType> & multiVector ) {
  if ( _size != multiVector.getTotalSize() ) {
    cerr << _size << " " << multiVector.getTotalSize() << endl;
    throw aol::Exception ( "aol::Vector::assignFrom ( MultiVector & ) : dimensions do not match.", __FILE__, __LINE__ );
  }

  setZero();
  int components = multiVector.numComponents();
  int n = 0;
  for ( int i = 0; i < components; ++i ) {
    setBlock ( n, multiVector[i] );
    n += multiVector[i].size();
  }
  return *this;
}

template <typename _DataType>
aol::Vector<_DataType>& aol::Vector<_DataType>::operator+= ( const aol::Vector<_DataType> &Vec ) {
  if ( Vec.size() == _size ) {
    for ( int i = 0; i < _size; ++i ) {
      _pData[i] += Vec.get ( i );
    }
  } else {
    throw aol::Exception ( "aol::Vector::operator+= dimensions don't match", __FILE__, __LINE__ );
  }
  return *this;
}

template <typename _DataType>
aol::Vector<_DataType>& aol::Vector<_DataType>::operator-= ( const aol::Vector<_DataType> &Vec ) {
  if ( Vec.size() == _size ) {
    for ( int i = 0; i < _size; ++i ) {
      _pData[i] -= Vec.get ( i );
    }
  } else {
    throw aol::Exception ( "aol::Vector::operator-= dimensions don't match", __FILE__, __LINE__ );
  }
  return *this;
}

template <typename _DataType>
aol::Vector<_DataType>& aol::Vector<_DataType>::operator*= ( const _DataType Value ) {
  DataType *ptr = _pData;
  for ( int i = 0; i < _size; ++i )
    * ( ptr++ ) *= Value;
  return *this;
}

template <typename _DataType>
aol::Vector<_DataType>& aol::Vector<_DataType>::operator/= ( const _DataType Value ) {
  DataType *ptr = _pData;
  for ( int i = 0; i < _size; ++i )
    * ( ptr++ ) /= Value;
  return *this;
}

template <typename _DataType>
aol::Vector<_DataType> & aol::Vector<_DataType>::multMasked ( const _DataType Value, const BitVector & mask, bool invertMask ) {
  DataType *ptr = _pData;
  // Note: ^ is the C notation for XOR. With a short look at the logic table,
  // you can see that (a XOR b) equals a, if b == true, and not(a) if b == false.
  for ( int i = 0; i < _size; ++i ) {
    if ( mask[i] ^ invertMask )
      * ( ptr ) *= Value;
    ptr++;
  }
  return *this;
}

template <typename _DataType>
aol::Vector<_DataType> & aol::Vector<_DataType>::addMasked ( const aol::Vector<_DataType> & Vec, const aol::BitVector & mask, bool invertMask ) {
  if ( Vec.size() != _size )
    throw aol::Exception ( "aol::Vector::addMasked(): dimensions don't match", __FILE__, __LINE__ );

  DataType *ptrThis = _pData;
  DataType *ptrRhs  = Vec._pData;
  // Note: ^ is the C notation for XOR. With a short look at the logic table,
  // you can see that (a XOR b) equals a, if b == true, and not(a) if b == false.
  for ( int i = 0; i < _size; ++i ) {
    if ( mask[i] ^ invertMask )
      * ( ptrThis ) += * ( ptrRhs );
    ptrThis++;
    ptrRhs++;
  }
  return *this;
}

template <typename _DataType>
aol::Vector<_DataType> & aol::Vector<_DataType>::addMultiple ( const aol::Vector<_DataType>& Vec, _DataType Factor ) {
  if ( Vec.size() != _size )
    throw aol::Exception ( "aol::Vector::addMultiple(): dimensions don't match", __FILE__, __LINE__ );

  DataType *ptr = _pData;
  for ( int i = 0; i < _size; ++i ) {
    *ptr += Vec.get ( i ) * Factor;
    ptr++;
  }
  return *this;
}

template <typename _DataType>
aol::Vector<_DataType> & aol::Vector<_DataType>::addMultipleMasked ( const aol::Vector<_DataType> & Vec, _DataType Factor, const aol::BitVector & mask, bool invertMask ) {
  if ( Vec.size() != _size )
    throw aol::Exception ( "aol::Vector::addMultipleMasked(): dimensions don't match", __FILE__, __LINE__ );

  DataType *ptrThis = _pData;
  DataType *ptrRhs  = Vec._pData;
  // Note: ^ is the C notation for XOR. With a short look at the logic table,
  // you can see that (a XOR b) equals a, if b == true, and not(a) if b == false.
  for ( int i = 0; i < _size; ++i ) {
    if ( mask[i] ^ invertMask )
      * ( ptrThis ) += * ( ptrRhs ) * Factor;
    ptrThis++;
    ptrRhs++;
  }
  return *this;
}

template <typename _DataType>
void aol::Vector<_DataType>::setSum ( const aol::Vector<_DataType>& Vec1,
                                      const aol::Vector<_DataType>& Vec2, _DataType Factor ) {
  if ( Vec1.size() != _size || Vec2.size() != _size ) {
    cerr << "ERROR in add: incompatible vectorlengths...\n";
    return;
  }
  DataType *ptr = _pData;
  for ( int i = 0; i < _size; ++i ) {
    *ptr = Vec1.get ( i ) + Vec2.get ( i ) * Factor;
    ptr++;
  }
}


template <typename _DataType>
typename aol::Vector<_DataType>::RealType aol::Vector<_DataType>::lpNormPowP ( RealType p ) const {
  RealType ret = 0.0;
  DataType *ptr = _pData;
  if ( p == static_cast<RealType> ( 0 ) ) {
    return static_cast<RealType> ( _size );
  } else if ( p == static_cast<RealType> ( 1 ) ) {
    for ( int i = 0; i < _size; ++i ) ret += static_cast< RealType > ( aol::Abs ( *ptr++ ) );
  } else if ( p == static_cast<RealType> ( 2 ) ) {
    for ( int i = 0; i < _size; ++i ) {
      ret += aol::Sqr ( *ptr );
      ptr++;
    }
  } else                      {
    for ( int i = 0; i < _size; ++i ) ret += static_cast< RealType > (  pow ( static_cast< RealType > ( aol::Abs ( *ptr++ ) ), static_cast< RealType > ( p ) ) );
  }
  return static_cast<RealType> ( ret );
}

template <typename _DataType>
typename aol::Vector<_DataType>::RealType aol::Vector<_DataType>::lpNorm ( RealType p ) const {
  RealType ret = lpNormPowP ( p );
  if ( p == aol::ZOTrait<RealType>::zero ) {
    return static_cast<RealType> ( _size );
  } else if ( p == aol::ZOTrait<RealType>::one ) {
    return static_cast<RealType> ( ret );
  } else if ( p == static_cast<RealType> ( 2 ) ) {
    return static_cast<RealType> ( sqrt ( ret ) );
  } else                     {
    return static_cast<RealType> ( pow ( ret, static_cast<RealType> ( 1.0 / p ) ) );
  }
}

template <typename _DataType>
_DataType aol::Vector<_DataType>::sum ( int p ) const {
  int i;
  DataType ret = static_cast< DataType > ( 0.0 );
  DataType *ptr = _pData;

  switch ( p ) {
    case 0:
      return static_cast< DataType > ( _size );
    case 1:
      for ( i = 0; i < _size; ++i )
        ret += static_cast< DataType > ( *ptr++ );
      break;
    case 2:
      for ( i = 0; i < _size; ++i )
        ret += aol::Sqr ( *ptr++ );
      break;
    default:
      throw aol::Exception ( "aol::Vector<DataType>::sum( int power ) not implemented for power != (0, 1, 2)", __FILE__, __LINE__ );
      //  for(i=0; i<N; ++i)
      //    ret += (DataType)pow((RealType)*ptr++, (RealType)p);
      //  break;
  }
  return ret;
}

template <typename _DataType>
_DataType aol::Vector<_DataType>::sumWeighted ( aol::Vector<_DataType> const& weight ) const {
#ifdef BOUNDS_CHECK
  if ( weight.size() < _size ) {
    char error[1024];
    sprintf ( error, "Vector<DataType>::sumWeighted: weight vector not long enough (%d), should be %d at least.\n", weight.size(), _size );
    throw OutOfBoundsException ( error, __FILE__, __LINE__ );
  }
#endif
  _DataType ret = 0.;
  _DataType *ptr = _pData;
  for ( int i = 0; i < _size; ++i )
    ret += *ptr++ * weight[i];
  return ret;
}

template <typename _DataType>
void aol::Vector<_DataType>::getBlock ( int start, aol::Vector<_DataType>& block ) const {
  int siz = block.size ();

  if ( start + siz > size () )
    throw aol::Exception ( "aol::Vector<T,R>::getBlock: Parameter mismatch: Block does not fit", __FILE__, __LINE__ );

  for ( int i = 0; i < siz; ++i )
    block.set ( i, get ( start + i ) );
}

template <typename _DataType>
void aol::Vector<_DataType>::setBlock ( int start, const aol::Vector<_DataType>& block ) {
  int siz = block.size ();

  if ( start + siz > size () )
    throw aol::Exception ( "aol::Vector<T,R>::getBlock: Parameter mismatch: Block does not fit", __FILE__, __LINE__ );

  for ( int i = 0; i < siz; ++i )
    set ( start + i, block.get ( i ) );
}

template <typename _DataType>
void aol::Vector<_DataType>::addBlock ( int start, const aol::Vector<_DataType>& block ) {
  int siz = block.size ();

  if ( start + siz > size () )
    throw aol::Exception ( "aol::Vector<T,R>::addBlock: Parameter mismatch: Block does not fit", __FILE__, __LINE__ );

  for ( int i = 0; i < siz; ++i )
    add ( start + i, block.get ( i ) );
}

template <typename _DataType>
typename aol::Vector<_DataType>::RealType aol::Vector<_DataType>::getWeightedMeanValue ( const aol::Vector<RealType> &Weights ) const {
  const int numVals = this->size();

  RealType mean = 0;
  for ( int i = 0; i < numVals; ++i )
    mean += Weights[i] * get(i);
  mean /= Weights.sum();

  return mean;
}

template <typename _DataType>
typename aol::Vector<_DataType>::RealType aol::Vector<_DataType>::getMedianValue ( ) const {
  const int sz = this->size();
  std::vector< DataType > stdVector ( sz );
  for ( int i = 0; i < sz; ++i ) {
    stdVector[i] = (*this)[i];
  }
  std::sort ( stdVector.begin(), stdVector.end() );
  RealType ret;
  if ( sz == 0 ) {
    throw aol::Exception ( "Cannot compute median of empty vector", __FILE__, __LINE__ );
  } else if ( sz % 2 == 1 ) { // for odd size, take central value
    ret = static_cast<RealType> ( stdVector[ sz/2 ] );
  } else { // for even size, take arithmetic mean of two central values
    ret = static_cast<RealType> ( stdVector[ sz/2 - 1 ] + stdVector [ sz/2 ] ) / 2;
  }
  return ( ret );
}

template <typename _DataType>
typename aol::Vector<_DataType>::RealType aol::Vector<_DataType>:: getWeightedMedianValue( const aol::Vector<RealType> &Weights ) const {
  const int numVals = this->size();

  if ( numVals != Weights.size() )
    throw Exception ( "aol::Vector<DataType>::getWeightedMedianValue: numVals != Weights.size() !\n", __FILE__, __LINE__ );

  if ( Weights.getMinValue() <= 0 )
    throw Exception ( "aol::Vector<DataType>::getWeightedMedianValue: Nonpositive weights are not supported!\n", __FILE__, __LINE__ );

  std::vector<std::pair<_DataType, RealType> > valuesAndWeights;
  valuesAndWeights.reserve ( numVals );
  for ( int i = 0; i < numVals; ++i )
    valuesAndWeights.push_back ( std::pair<_DataType, RealType> ( this->get ( i ), Weights[i] ) );

  std::sort( valuesAndWeights.begin(), valuesAndWeights.end() );
  const RealType halfOfTotalWeight = 0.5 * Weights.sum();

  RealType weight = 0;
  for ( int i = 0; i < numVals; ++i ) {
    weight += valuesAndWeights[i].second;
    if ( ( weight > halfOfTotalWeight ) || ( i == ( numVals - 1 ) ) )
      return valuesAndWeights[i].first;
    else if ( aol::appeqAbsolute ( weight, halfOfTotalWeight ) ) {
      return ( valuesAndWeights[i].first * valuesAndWeights[i].second + valuesAndWeights[i+1].first * valuesAndWeights[i+1].second ) / ( valuesAndWeights[i].second +  valuesAndWeights[i+1].second );
    }
  }
  return valuesAndWeights[numVals-1].first;
}

template <typename _DataType>
void aol::Vector<_DataType>::sortValues ( ) {
  _DataType* begin = getData();
  _DataType* end = getData() + size();
  std::make_heap ( begin, end );
  std::sort_heap ( begin, end );
}

template <typename _DataType>
int aol::Vector<_DataType>::indexOfFirstOccurence ( const DataType Value ) const {
  for ( int i = 0; i < this->size(); ++i ) {
    if ( Value == this->get ( i ) ) {
      return ( i );
    }
  }
  return ( -1 ); // not found
}

template <typename _DataType>
void aol::Vector<_DataType>::erase ( const int Index ) {
  for ( int i = Index + 1; i < this->size(); ++i ) {
    this->set ( i - 1, this->get ( i ) );
  }
  this->growBy ( -1 );
}

template <typename _DataType>
void aol::Vector<_DataType>::insert ( const int Index, _DataType value ) {
  this->growBy ( 1 );
  for ( int i = this->size () - 1; i > Index; --i ) {
    this->set ( i, this->get ( i - 1 ) );
  }

  this->set ( Index, value );
}

template <typename _DataType>
void aol::Vector<_DataType>::eraseFirstOccurence ( const DataType Value ) {
  const int pos = indexOfFirstOccurence ( Value );
  if ( pos != -1 ) {
    this->erase ( pos );
  } // else do nothing (value not found)
}


template <typename _DataType>
unsigned int aol::Vector<_DataType>::crc32OfData ( ) const {
  return ( aol::crc32 ( this->getData(), this->size() * sizeof ( DataType ) ) );
}

template <typename _DataType>
void aol::Vector<_DataType>::createHistogram ( aol::Vector<int> &/*Histo*/ ) const {
  throw Exception ( "aol::Vector<DataType>::createHistogram may only be called if DataType is unsigned char!\n", __FILE__, __LINE__ );
}

// GCC wants this specialization to be put into the namespace.
namespace aol {

template <>
void Vector<unsigned char>::createHistogram ( aol::Vector<int> &Histo ) const {
  const int numPixels = this->size();
  Histo.reallocate ( 256 );
  for ( int i = 0; i < numPixels; ++i )
    ++Histo [ this->get ( i ) ];

}

}

template <typename DataType>
bool aol::Vector<DataType>::quietMode =
#ifdef VERBOSE
  false;
#else
  true;
#endif

template <typename _DataType> const aol::Format* aol::Vector<_DataType>::_pFormat = &mixedFormat;
template <typename _DataType> bool aol::Vector<_DataType>::PrettyFormat = true;

// template class aol::Vector<bool>; // aol::Vector<bool> should not be used ( += doesn't make sense with bool ), use aol::BitVector instead.
// template class aol::Vector<char>; // DO NOT USE Vector<char>: char may be signed or unsigned, depending on your platform.
template class aol::Vector<signed char>;
template class aol::Vector<unsigned char>;
template class aol::Vector<short>;
template class aol::Vector<unsigned short>;
template class aol::Vector<int>;
template class aol::Vector<unsigned int>;
template class aol::Vector<int64_t>;
template class aol::Vector<uint64_t>;
#if ( defined __APPLE__ ) || ( BITNESS == 32 ) // This seems to be the only way to get SuiteSparse_long to work under OS X
template class aol::Vector<long>;
#endif
template class aol::Vector<float>;
template class aol::Vector<double>;
template class aol::Vector<long double>;
