#include <pngInterface.h>

#ifdef USE_LIB_PNG

qc::PNGLoadInterface::PNGLoadInterface ( const char *FileName, const bool QuietMode, const bool StripAlpha )
  : _quietMode ( QuietMode ) {
  _filePointer = fopen ( FileName, "rb" );
  if ( !_filePointer ) {
    string errorMessage = "Could not open file ";
    errorMessage += FileName;
    errorMessage += " for reading.\n";
    throw aol::Exception ( errorMessage, __FILE__, __LINE__ );
  }
  const int numberOfBytesToCheck = 4;
  png_byte header[numberOfBytesToCheck];
  const int numberOfElementsRead = static_cast<int> ( fread ( header, 1, numberOfBytesToCheck, _filePointer ) );
  if ( numberOfElementsRead != numberOfBytesToCheck )
    throw aol::IOException ( "Error reading header from file", __FILE__, __LINE__ );

  const bool is_png = !png_sig_cmp ( header, static_cast< png_size_t > ( 0 ), numberOfBytesToCheck );
  if ( !is_png ) {
    throw aol::Exception ( "File is not a png", __FILE__, __LINE__ );
  }

  _pngPointer = png_create_read_struct
    ( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );

  if ( !_pngPointer )
    throw aol::Exception ( "Executing png_create_read_struct failed", __FILE__, __LINE__ );

  _infoPointer = png_create_info_struct ( _pngPointer );

  if ( !_infoPointer ) {
    png_destroy_read_struct ( &_pngPointer, static_cast< png_infopp > ( NULL ), static_cast< png_infopp > ( NULL ) );
    throw aol::Exception ( "Executing png_create_info_struct failed", __FILE__, __LINE__ );
  }

  _endInfoPointer = png_create_info_struct ( _pngPointer );
  if ( !_endInfoPointer ) {
    png_destroy_read_struct ( &_pngPointer, &_infoPointer, static_cast< png_infopp > ( NULL ) );
    throw aol::Exception ( "Executing png_create_info_struct failed", __FILE__, __LINE__ );
  }

  if ( setjmp ( png_jmpbuf ( _pngPointer ) ) ) {
    png_destroy_read_struct ( &_pngPointer, &_infoPointer, &_endInfoPointer );
    fclose ( _filePointer );
    throw aol::Exception ( "An error occured", __FILE__, __LINE__ );
  }

  png_init_io ( _pngPointer, _filePointer );

  png_set_sig_bytes ( _pngPointer, numberOfBytesToCheck );

  // PNG_TRANSFORM_EXPAND is necessary to handle PNGs that use a palette.
  png_read_png ( _pngPointer, _infoPointer, ( StripAlpha ? PNG_TRANSFORM_STRIP_ALPHA : 0 )|PNG_TRANSFORM_EXPAND, NULL );

  _width = png_get_image_width ( _pngPointer, _infoPointer );
  _height = png_get_image_height ( _pngPointer, _infoPointer );
  _bit_depth = png_get_bit_depth ( _pngPointer, _infoPointer );
  _channels = png_get_channels ( _pngPointer, _infoPointer );
  // channels - number of channels of info for the color type
  // (valid values are 1 (GRAY, PALETTE), 2 (GRAY_ALPHA), 3 (RGB), 4 (RGB_ALPHA or RGB + filler byte))
  if ( !_quietMode ) {
    cerr << FileName << " read.\n";
    cerr << "width = " << _width << ", height = " << _height << ", bit_depth = " << _bit_depth << ", channels = " << _channels << endl;
  }

  if ( _bit_depth == 16 )
    throw aol::UnimplementedCodeException ( "qc::PNGLoadInterface::PNGLoadInterface: Reading 16-bit per channel PNG data not implemented", __FILE__, __LINE__ );
}

template< typename DataType >
void qc::PNGLoadInterface::writeDataToScalarArray2d ( qc::ScalarArray<DataType, qc::QC_2D> &Array ) {
  png_bytepp row_pointers = png_get_rows ( _pngPointer, _infoPointer );

  for ( int x = 0; x < _width; x++ ) {
    for ( int y = 0; y < _height; y++ ) {
      double temp = 0.;
      for ( int channel = 0; channel < _channels; channel++ )
        temp += static_cast<double> ( row_pointers[y][_channels*x+channel] );
      temp /= static_cast<double> ( _channels );
      Array.set ( x, y, static_cast<DataType> ( temp ) );
    }
  }
}

template void qc::PNGLoadInterface::writeDataToScalarArray2d ( qc::ScalarArray<signed char, qc::QC_2D> &Array );
template void qc::PNGLoadInterface::writeDataToScalarArray2d ( qc::ScalarArray<unsigned char, qc::QC_2D> &Array );
template void qc::PNGLoadInterface::writeDataToScalarArray2d ( qc::ScalarArray<short, qc::QC_2D> &Array );
template void qc::PNGLoadInterface::writeDataToScalarArray2d ( qc::ScalarArray<unsigned short, qc::QC_2D> &Array );
template void qc::PNGLoadInterface::writeDataToScalarArray2d ( qc::ScalarArray<int, qc::QC_2D> &Array );
template void qc::PNGLoadInterface::writeDataToScalarArray2d ( qc::ScalarArray<unsigned int, qc::QC_2D> &Array );
template void qc::PNGLoadInterface::writeDataToScalarArray2d ( qc::ScalarArray<float, qc::QC_2D> &Array );
template void qc::PNGLoadInterface::writeDataToScalarArray2d ( qc::ScalarArray<double, qc::QC_2D> &Array );
template void qc::PNGLoadInterface::writeDataToScalarArray2d ( qc::ScalarArray<long double, qc::QC_2D> &Array );

qc::PNGLoadInterface::~PNGLoadInterface() {
  png_destroy_read_struct ( &_pngPointer, &_infoPointer, &_endInfoPointer );
  fclose ( _filePointer );
}

qc::PNGSaveInterface::PNGSaveInterface ( const char *FileName ) {
  _filePointer = fopen ( FileName, "wb" );
  if ( !_filePointer ) {
    string errorMessage = "Could not open file ";
    errorMessage += FileName;
    errorMessage += " for writing.\n";
    throw aol::Exception ( errorMessage, __FILE__, __LINE__ );
  }
  _pngPointer = png_create_write_struct
    ( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );
  if ( !_pngPointer )
    throw aol::Exception ( "Executing png_create_write_struct failed", __FILE__, __LINE__ );

  _infoPointer = png_create_info_struct ( _pngPointer );
  if ( !_infoPointer ) {
    png_destroy_write_struct ( &_pngPointer, static_cast< png_infopp > ( NULL ) );
    throw aol::Exception ( "Executing png_create_info_struct failed", __FILE__, __LINE__ );
  }

  if ( setjmp ( png_jmpbuf ( _pngPointer ) ) ) {
    png_destroy_write_struct ( &_pngPointer, &_infoPointer );
    fclose ( _filePointer );
    throw aol::Exception ( "An error occured", __FILE__, __LINE__ );
  }

  png_init_io ( _pngPointer, _filePointer );
}

qc::PNGSaveInterface::~PNGSaveInterface() {
  png_destroy_write_struct ( &_pngPointer, &_infoPointer );
  fclose ( _filePointer );
}

template< typename DataType >
void qc::PNGSaveInterface::loadDataFromScalarArray2d ( const qc::ScalarArray<DataType, qc::QC_2D> &Array ) {
  if ( ( Array.getNumX() == 0 ) || ( Array.getNumY() == 0 ) ) {
    cerr << Array.getNumX() << " " << Array.getNumY() << endl;
    throw aol::Exception ( "Cannot save ScalarArray of size 0 as PNG.", __FILE__, __LINE__ );
  }

  png_set_IHDR ( _pngPointer, _infoPointer, Array.getNumX(), Array.getNumY(), 8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
    PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT );

  png_bytepp row_pointers = new png_bytep[ Array.getNumY() ];

  aol::Vector<unsigned char> vecDummy ( Array.getNumX() * Array.getNumY() );
  Array.createOverflowHandledData ( vecDummy );
  unsigned char *dummy = vecDummy.getData();

  for ( int i = 0; i < Array.getNumY(); i++ )
    row_pointers[i] = dummy + sizeof ( unsigned char ) * i * Array.getNumX();

  png_set_rows ( _pngPointer, _infoPointer, row_pointers );

  png_write_png ( _pngPointer, _infoPointer, PNG_TRANSFORM_IDENTITY, NULL );

  delete[] row_pointers;
}

template void qc::PNGSaveInterface::loadDataFromScalarArray2d ( const qc::ScalarArray<signed char, qc::QC_2D> &Array );
template void qc::PNGSaveInterface::loadDataFromScalarArray2d ( const qc::ScalarArray<unsigned char, qc::QC_2D> &Array );
template void qc::PNGSaveInterface::loadDataFromScalarArray2d ( const qc::ScalarArray<short, qc::QC_2D> &Array );
template void qc::PNGSaveInterface::loadDataFromScalarArray2d ( const qc::ScalarArray<unsigned short, qc::QC_2D> &Array );
template void qc::PNGSaveInterface::loadDataFromScalarArray2d ( const qc::ScalarArray<int, qc::QC_2D> &Array );
template void qc::PNGSaveInterface::loadDataFromScalarArray2d ( const qc::ScalarArray<unsigned int, qc::QC_2D> &Array );
template void qc::PNGSaveInterface::loadDataFromScalarArray2d ( const qc::ScalarArray<float, qc::QC_2D> &Array );
template void qc::PNGSaveInterface::loadDataFromScalarArray2d ( const qc::ScalarArray<double, qc::QC_2D> &Array );
template void qc::PNGSaveInterface::loadDataFromScalarArray2d ( const qc::ScalarArray<long double, qc::QC_2D> &Array );

#endif // USE_LIB_PNG
