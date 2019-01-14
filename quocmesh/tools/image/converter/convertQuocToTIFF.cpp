/**
 * \file
 * \brief Converts one (or multiple) 2D Quoc array(s) to a single TIFF with 16 bit precision (unsigned short).
 *
 * Usage: convertQuocToTIFF InputFile [numSlices]
 *
 * If used with a single argument, it assumes that the argument is a 2D array file.
 * If used with two arguments, it assumes that the first argument is a filename mask for
 * several 2D array files and the second argument is the number of slices.
 *
 * \author Berkels
 */

#include <aol.h>

#ifdef USE_EXTERNAL_CIMG

#include <scalarArray.h>
#include <cimgIncludes.h>

typedef double RType;

int main ( int argc, char **argv ) {
  try {
    if ( ( argc != 2 ) && ( argc != 3 ) ) {
      cerr << "USAGE: " << argv[0] << "  <InputFile> [numSlices]" << endl;
      return EXIT_FAILURE;
    }

    const string inputFileName = argv[1];
    const int numSlices = ( argc >= 3 ) ? atoi ( argv[2] ) : 0;

    if ( numSlices <= 1 ) {
      qc::ScalarArray<unsigned short, qc::QC_2D> imageArray ( inputFileName.c_str() );

      cimg_library::CImg<unsigned short> cimg (imageArray.getData(), imageArray.getNumX(), imageArray.getNumY() );
      cimg.save_tiff ( ( aol::getBaseFileName ( inputFileName ) + ".tif" ).c_str() );
    }
    else {
      const aol::Vec3<int> sliceSize = qc::getSizeFromArrayFile ( aol::strprintf ( inputFileName.c_str(), 0 ) );
      qc::ScalarArray<unsigned short, qc::QC_3D> imageArray ( sliceSize[0], sliceSize[1], numSlices );
      imageArray.loadSlices ( inputFileName.c_str(), qc::QC_Z, 0, numSlices - 1 );

      cimg_library::CImg<unsigned short> cimg ( imageArray.getData(), imageArray.getNumX(), imageArray.getNumY(), imageArray.getNumZ() );
      // Remove the "%...d" mask for the slice index from the output filename.
      const int posOfPercent = inputFileName.find ( '%' );
      const int posOfD = inputFileName.find ( 'd', posOfPercent );
      const string outputBaseFileName = aol::getBaseFileName ( inputFileName.substr ( 0, posOfPercent ) +  inputFileName.substr ( posOfD + 1, inputFileName.size() ) );
      cimg.save_tiff ( ( aol::strprintf ( outputBaseFileName.c_str(), 0 ) + ".tif" ).c_str() );
    }
  }
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}

#else // USE_EXTERNAL_CIMG

int main ( int, char** ) {
  cerr << "Needs to be compiled with the external cimg.\n";
  return 0;
}

#endif // USE_EXTERNAL_CIMG
