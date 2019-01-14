#ifndef __PNGINTERFACE_H
#define __PNGINTERFACE_H

#ifdef USE_LIB_PNG
#include <png.h>
#include <multiArray.h>

namespace qc {

/**
 * Interface class to load PNG files into quocmesh data structures.
 *
 * Needs libpng.
 *
 * \author Berkels
 */
class PNGLoadInterface {
  FILE *_filePointer;
  png_structp _pngPointer;
  png_infop _infoPointer;
  png_infop _endInfoPointer;
  int _width;
  int _height;
  int _bit_depth;
  int _channels;
  bool _quietMode;
public:
  PNGLoadInterface ( const char *FileName, const bool QuietMode = false, const bool StripAlpha = true );

  ~PNGLoadInterface();

  template< typename DataType >
  void writeDataToScalarArray2d ( qc::ScalarArray<DataType, qc::QC_2D> &Array );

  template< typename DataType, int MArrayImageDim >
  void writeDataToMultiArray ( qc::MultiArray<DataType, 2, MArrayImageDim> &VArray ) {
    if ( ( MArrayImageDim != 3 ) && ( MArrayImageDim != 4 ) && !( ( MArrayImageDim == 2 ) && ( _channels == 2 ) ) )
      throw aol::UnimplementedCodeException( "PNGLoadInterface::writeDataToMultiArray is not implemented for MArrayImageDim != 3, 4", __FILE__, __LINE__);

    png_bytepp row_pointers = png_get_rows ( _pngPointer, _infoPointer );

    for ( int x = 0; x < _width; x++ ) {
      for ( int y = 0; y < _height; y++ ) {
        double temp = 0.;
        for ( int channel = 0; channel < aol::Min ( _channels, MArrayImageDim ); channel++ ) {
          temp = static_cast<double> ( row_pointers[y][_channels*x+channel] );
          VArray.comp ( channel ).set ( x, y, static_cast<DataType> ( temp ) );
        }
      }
    }

    // If the input png is grayscale with alpha channel and we want RGBA as output,
    // we copy the alpha channel to the forth component of the MultiArray.
    if ( ( _channels == 2 ) && ( MArrayImageDim == 4 ) ) {
      VArray[3] = VArray[1];
    }

    // If the input png is grayscale, we copy the gray channel to all
    // components of the MultiArray. In case of ( _channels == 2 ) && ( MArrayImageDim == 3 )
    // this purposely removes the alpha channel.
    if ( ( _channels <= 2 ) && ( ( MArrayImageDim == 3 ) || ( MArrayImageDim == 4 ) ) ) {
      VArray[1] = VArray[0];
      VArray[2] = VArray[0];
    }

    // If the input image doesn't have an alpha channel, but the output is supposed to have one,
    // assume that everything is opaque.
    if ( ( _channels != 4 ) && ( _channels != 2 ) && ( MArrayImageDim == 4 ) )
      VArray[3].setAll ( 255 );
  }

  int getWidth() {
    return _width;
  }
  int getHeight() {
    return _height;
  }
};

/**
 * Interface class to save PNG files from quocmesh data structures.
 *
 * Needs libpng.
 *
 * \author Berkels
 */
class PNGSaveInterface {
  FILE *_filePointer;
  png_structp _pngPointer;
  png_infop _infoPointer;
public:
  PNGSaveInterface ( const char *FileName );

  ~PNGSaveInterface();

  template< typename DataType >
  void loadDataFromScalarArray2d ( const qc::ScalarArray<DataType, qc::QC_2D> &Array );

  template< typename DataType, int MArrayImageDim >
  void loadDataFromMultiArray ( const qc::MultiArray<DataType, 2, MArrayImageDim> &VArray ) {
    if ( ( MArrayImageDim != 3 ) && ( MArrayImageDim != 4 ) )
      throw aol::UnimplementedCodeException( "PNGSaveInterface::loadDataFromMultiArray is not implemented for MArrayImageDim != 3, 4", __FILE__, __LINE__);

    const int numX = VArray[0].getNumX();
    const int numY = VArray[0].getNumY();
    for( int i = 1; i < VArray.numComponents(); i++ ){
      if ( (numX != VArray[i].getNumX()) || (numY != VArray[i].getNumY()) )
        throw aol::Exception ( "All arrays have to be of the same dimensions.", __FILE__, __LINE__ );
    }
    if ( ( numX == 0 ) || ( numY == 0 ) ) {
      cerr << numX << " " << numY << endl;
      throw aol::Exception ( "Cannot save ScalarArray of size 0 as PNG.", __FILE__, __LINE__ );
    }
    png_set_IHDR ( _pngPointer, _infoPointer, numX, numY, 8, ( ( MArrayImageDim == 3 ) ? PNG_COLOR_TYPE_RGB : PNG_COLOR_TYPE_RGB_ALPHA ), PNG_INTERLACE_NONE,
      PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT );

    png_bytepp row_pointers = new png_bytep[ MArrayImageDim*numY ];
    unsigned char *dummy = new unsigned char[ MArrayImageDim*numX*numY ];
    aol::Vector<unsigned char> vecDummyChannelData ( numX*numY );
    unsigned char *dummyChannelData = vecDummyChannelData.getData();

    // This is not the most efficient way to implement this, but using
    // createOverflowHandeledData greatly reduces code duplication.
    for( int channel = 0; channel < MArrayImageDim; channel++ ){
      VArray[channel].createOverflowHandledData ( vecDummyChannelData );
      for ( int i = 0; i < numX*numY; i++ )
        dummy[MArrayImageDim*i+channel] = dummyChannelData[i];
    }

    for ( int i = 0; i < numY; i++ )
      row_pointers[i] = dummy + sizeof ( unsigned char ) * MArrayImageDim * i * numX;

    png_set_rows ( _pngPointer, _infoPointer, row_pointers );

    png_write_png ( _pngPointer, _infoPointer, PNG_TRANSFORM_IDENTITY, NULL );

    delete[] row_pointers;
    delete[] dummy;
  }
};


} // namespace qc
#endif // USE_LIB_PNG

#endif  // __PNGINTERFACE_H
