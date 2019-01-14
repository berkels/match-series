/**
 * \file
 * \brief Converts a tiff file (including those with 16 bit precision) to the quoc format.
 *        If the tiff file contains multiple slices, the slices are written as separate 2D
 *        quoc arrays.

 * Usage: convertTIFFToQuoc InputFile
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
    if ( argc != 2 ) {
      cerr << "USAGE: " << argv[0] << "  <InputFile>" << endl;
      return EXIT_FAILURE;
    }

    const string inputFileName = argv[1];
    const string outputBaseFileName = aol::getBaseFileName ( inputFileName );

    cerr << "Reading TIFF from " << inputFileName << endl;

    cimg_library::CImg<RType> cimg;
    cimg.load_tiff ( inputFileName.c_str() );

    if ( cimg.depth () > 1 ) {
      qc::ScalarArray<RType, qc::QC_3D> a ( cimg.width(), cimg.height(), cimg.depth (), cimg.data(), aol::FLAT_COPY );
      a.saveSlices ( ( outputBaseFileName + "_%03d" + qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str(), qc::QC_Z, qc::PGM_DOUBLE_BINARY );
    }
    else {
      qc::ScalarArray<RType, qc::QC_2D> a ( cimg.width(), cimg.height() , cimg.data(), aol::FLAT_COPY );
      a.save ( ( outputBaseFileName + qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str(), qc::PGM_DOUBLE_BINARY );
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
