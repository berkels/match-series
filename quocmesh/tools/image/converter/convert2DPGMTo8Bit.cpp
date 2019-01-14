/**
 * \file
 * \brief Loads 2D quoc dataset with floating point values and saves it as two 8 bit grey value pgm images
 *        (one of them mapping [0, 1] to [0, 255], the other mapping the full range to [0, 255]).
 *
 * Usage: convert2DPGMTo8Bit inputfile
 *
 * \author Berkels
 */

#include <aol.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <quocTimestepSaver.h>

typedef float RType;
int main ( int argc, char **argv ) {

  try {
    char inputFileName[1024];
    char outputFileName[1024];

    if ( ( argc < 2 ) || ( argc > 6 ) ){
      cerr << "USAGE: " << argv[0] << " <input_file> [divideInputBy255] [onlySaveScaled] [removeInf] [saturatedEntryPercentage]" << endl;
      return EXIT_FAILURE;
    }
    if ( argc >= 2 ) {
      sprintf ( inputFileName, "%s",  argv[1] );
    }
    const bool divideInputBy255 = ( argc >= 3 ) ? ( atoi(argv[2]) != 0 ) : false;
    const bool onlySaveScaled = ( argc >= 4 ) ? ( atoi(argv[3]) != 0 ) : false;
    const bool removeInf = ( argc >= 5 ) ? ( atoi(argv[4]) != 0 ) : false;
    const RType saturatedEntryPercentage = ( argc >= 6 ) ? ( atof(argv[5]) != 0 ) : 0;
    cerr << "Reading 2D PGM from " << inputFileName << endl;

    qc::ScalarArray<RType, qc::QC_2D> a ( inputFileName );
    if ( divideInputBy255 ) {
      cerr << "Dividing input by 255\n";
      a /= 255;
    }
    if ( removeInf ) {
      const int numEntries = a.size();
      aol::BitVector infMask ( numEntries );
      for ( int i = 0; i < numEntries; ++i )
        infMask.set ( i, aol::isInf ( a[i] ) );

      const RType minValue = a.getMinValueMasked ( infMask, true );
      const RType maxValue = a.getMaxValueMasked ( infMask, true );

      a.setAllMasked ( ( minValue + maxValue ) / 2, infMask );
    }
    sprintf ( outputFileName, "%s_8.png", inputFileName );
    a.setQuietMode ( true );
    a.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0., 1. );
    if ( onlySaveScaled == false )
      a.savePNG ( outputFileName );

    qc::DefaultArraySaver<RType, qc::QC_2D> saver ( true, true );
    if ( saturatedEntryPercentage > 0 )
      saver.setEnhanceContrastSaturationPercentage ( saturatedEntryPercentage );
    saver.saveStep( a, -1, inputFileName );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
