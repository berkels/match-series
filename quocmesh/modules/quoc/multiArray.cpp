#include <scalarArray.h>
#include <multiArray.h>
#include <pngInterface.h>

namespace qc {


template <typename DataType, int imagedim>
#ifdef USE_LIB_PNG
void loadPNGToMArray ( qc::MultiArray<DataType, 2, imagedim> &MArg, const char *FileName, const bool StripAlpha ) {
  PNGLoadInterface pngInterface( FileName, false, StripAlpha );
  MArg.reallocate ( pngInterface.getWidth(), pngInterface.getHeight() );
  pngInterface.writeDataToMultiArray<DataType>( MArg );
#else
void loadPNGToMArray ( qc::MultiArray<DataType, 2, imagedim> &, const char *, const bool ) {
  throw aol::Exception ( "Reading PNG files requires libpng. Compile with -DUSE_LIB_PNG.", __FILE__, __LINE__ );
#endif //USE_LIB_PNG
}

template <typename DataType, int imagedim>
#ifdef USE_LIB_PNG
void savePNGFromMArray ( const qc::MultiArray<DataType, 2, imagedim> &MArg, const char *FileName ) {
  PNGSaveInterface pngInterface( FileName );
  pngInterface.loadDataFromMultiArray<DataType>( MArg );
#else
void savePNGFromMArray ( const qc::MultiArray<DataType, 2, imagedim> &, const char * ) {
  throw aol::Exception ( "Writing PNG files requires libpng. Compile with -DUSE_LIB_PNG.", __FILE__, __LINE__ );
#endif //USE_LIB_PNG
}

/**
 * This struct + static function construction is a workaround to C++ limitations.
 *
 * \author Berkels
 */
template <typename DataType, int rangedim, int imagedim>
struct doPNGLoadSave {
  static void loadPNG ( qc::MultiArray<DataType, rangedim, imagedim> &/*MArg*/, const char * /*FileName*/ ) {
    throw aol::Exception( "Only the template specializations of doPNGLoadSave should be used.", __FILE__, __LINE__);
  }
};

template <typename DataType>
struct doPNGLoadSave<DataType, 2, 2> {
  static void loadPNG ( qc::MultiArray<DataType, 2, 2> &MArg, const char *FileName ) {
    loadPNGToMArray ( MArg, FileName, false );
  }

  static void savePNG ( const qc::MultiArray<DataType, 2, 2> &/*MArg*/, const char * /*FileName*/ ) {
    throw aol::UnimplementedCodeException( "doPNGLoadSave<DataType, 2, 2> PNGSaveInterface::savePNG is not implemented.", __FILE__, __LINE__);
  }
};

template <typename DataType>
struct doPNGLoadSave<DataType, 2, 3> {
  static void loadPNG ( qc::MultiArray<DataType, 2, 3> &MArg, const char *FileName ) {
    loadPNGToMArray ( MArg, FileName, true );
  }

  static void savePNG ( const qc::MultiArray<DataType, 2, 3> &MArg, const char *FileName ) {
    savePNGFromMArray ( MArg, FileName );
  }
};

template <typename DataType>
struct doPNGLoadSave<DataType, 2, 4> {
  static void loadPNG ( qc::MultiArray<DataType, 2, 4> &MArg, const char *FileName ) {
    loadPNGToMArray ( MArg, FileName, false );
  }

  static void savePNG ( const qc::MultiArray<DataType, 2, 4> &MArg, const char *FileName ) {
    savePNGFromMArray ( MArg, FileName );
  }
};

template <typename DataType>
struct doPNGLoadSave<DataType, 2, 5> {
  static void loadPNG ( qc::MultiArray<DataType, 2, 5> &MArg, const char *FileName ) {
    loadPNGToMArray ( MArg, FileName, false );
  }

  static void savePNG ( const qc::MultiArray<DataType, 2, 5> &MArg, const char *FileName ) {
    savePNGFromMArray ( MArg, FileName );
  }
};

template <typename DataType>
struct doPNGLoadSave<DataType, 2, 1> {
  static void loadPNG ( qc::MultiArray<DataType, 2, 1> &MArg, const char *FileName ) {
    MArg[0].loadPNG ( FileName );
  }

  static void savePNG ( const qc::MultiArray<DataType, 2, 1> &MArg, const char *FileName ) {
    MArg[0].savePNG ( FileName );
  }
};

template <class DataType, int rangedim, int imagedim>
void qc::MultiArray<DataType, rangedim, imagedim>::loadPNG ( const char *FileName ) {
  doPNGLoadSave<DataType, rangedim, imagedim>::loadPNG ( *this, FileName );
}

template <class DataType, int rangedim, int imagedim>
void qc::MultiArray<DataType, rangedim, imagedim>::savePNG ( const char *FileName ) const {
  doPNGLoadSave<DataType, rangedim, imagedim>::savePNG ( *this, FileName );
}

template <class DataType, int rangedim, int imagedim>
void qc::MultiArray<DataType, rangedim, imagedim>::loadPPM ( const char *FileName ) {
  if ( ( rangedim != 2 ) || ( imagedim != 3 ) )
    throw aol::UnimplementedCodeException ( "qc::MultiArray<DataType>::loadPPM: Only implemented for the case \"2, 3\"", __FILE__, __LINE__ );

  std::ifstream in ( FileName, std::ios::binary );

  char tmp[4];
  int inWidth, inHeight, inMaxColor;

  in.get ( tmp, 3 );
  if ( tmp[0] != 'P' )
    throw aol::TypeException ( "qc::MultiArray<DataType>::loadPPM: wrong file format", __FILE__, __LINE__ );

  const char type = tmp[1];

  switch ( tmp[1] ) {
    case '3':
    case '6':
      break;
    default:
      throw aol::TypeException ( "qc::MultiArray<DataType>::loadPPM: wrong magic number", __FILE__, __LINE__ );
  }

  aol::READ_COMMENTS ( in );
  in >> inWidth;
  if ( in.fail() )
    throw aol::FileFormatException ( "qc::MultiArray<DataType>::loadPPM: error reading width", __FILE__, __LINE__ );

  aol::READ_COMMENTS ( in );
  in >> inHeight;
  if ( in.fail() )
    throw aol::FileFormatException ( "qc::MultiArray<DataType>::loadPPM: error reading height", __FILE__, __LINE__ );

  aol::READ_COMMENTS ( in );
  in >> inMaxColor;
  if ( in.fail() )
    throw aol::FileFormatException ( "qc::MultiArray<DataType>::loadPPM: error reading max color", __FILE__, __LINE__ );

  in.ignore();

  this->reallocate ( inWidth, inHeight );

  for ( int j = 0; j < inHeight; ++j ) {
    for ( int i = 0; i < inWidth; ++i ) {
      for ( int c = 0; c < 3; ++c ) {
        if ( type == '3' )
          in >> this->comp ( c ).getReference ( i, j );
        else
          this->comp ( c ).set ( i, j, aol::readBinaryData<unsigned char, DataType> ( in ) );
      }
    }
  }
}

} // namespace qc

template void qc::MultiArray<double, 1, 3>::loadPNG ( const char *FileName );

template class qc::MultiArray<double, 2, 1>;
template class qc::MultiArray<double, 2, 2>;
template class qc::MultiArray<float, 2, 3>;
template class qc::MultiArray<double, 2, 3>;
template class qc::MultiArray<long double, 2, 3>;
template class qc::MultiArray<unsigned char, 2, 3>;
template class qc::MultiArray<float, 2, 4>;
template class qc::MultiArray<double, 2, 4>;
template class qc::MultiArray<long double, 2, 4>;

template class qc::MultiArray<double, 2, 5>;
