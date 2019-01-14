/**
 * \file
 * \brief Deforms the 2D quoc input array as specified by the input displacement field and saves the result with double precision.
 *
 * \author Berkels
 */

#include <aol.h>
#include <op.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <multiArray.h>
#include <levelSetDrawer.h>
#include <timestepSaver.h>
#include <configurators.h>
#include <deformations.h>
#include <dm3Import.h>

typedef double RType;
typedef qc::RectangularGridConfigurator<RType, qc::QC_2D, aol::GaussQuadrature<RType, qc::QC_2D, 3> > ConfType;

int main ( int argc, char **argv ) {

  try {
    if ( ( argc < 5 ) || ( argc > 10 ) ) {
      cerr << "USAGE: " << argv[0] << " <input_png> <input_file def_x> <input_file def_y> <out_file> [invertDev] [NNinterp] [keepDMX] [removeNaNs] [xDerivNormThreshold]" << endl;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_FAILURE;
    }

    const string inputImageFileName = argv[1];
    const string inputFileNameDefX = argv[2];
    const string inputFileNameDefY = argv[3];
    const string outputFileName = argv[4];

    const bool useInverseDeformation = ( argc > 5 ) ? ( atoi ( argv[5] ) != 0 ) : false;
    const bool NearestNeighborInterpolation = ( argc > 6 ) ? ( atoi ( argv[6] ) != 0 ) : false;
    const bool keepDMX = ( argc > 7 ) ? ( atoi ( argv[7] ) != 0 ) : false;
    const bool removeNaNs = ( argc > 8 ) ? ( atoi ( argv[8] ) != 0 ) : false;
    const RType xDerivNormThreshold = ( argc > 9 ) ? atof ( argv[9] ) : 0;

    qc::ScalarArray<RType, qc::QC_2D> dataArray ( inputImageFileName.c_str() );
    if ( removeNaNs ) {
      for ( int i = 0; i < dataArray.size(); ++i ) {
        if ( aol::isNaN( dataArray[i] ) )
          dataArray[i] = 0;
      }
    }
    qc::MultiArray<RType, qc::QC_2D> deformation ( inputFileNameDefX, inputFileNameDefY );
    ConfType::InitType grid ( dataArray.getSize() );
    qc::ScalarArray<RType, qc::QC_2D> deformedArray ( grid );

    if ( useInverseDeformation == false ) {
      qc::DeformImage<ConfType>( dataArray, grid, deformedArray, deformation, true, 0, NearestNeighborInterpolation );
      if ( xDerivNormThreshold > 0 )
        qc::applyStretchMute<RType, qc::QC_2D> ( deformedArray, deformation, xDerivNormThreshold, 0 );
    }
    else
      qc::InvDeformImage<ConfType>( dataArray, grid, deformedArray, deformation );

    if ( keepDMX
        && ( aol::fileNameEndsWith( inputImageFileName.c_str(), "dm3" ) || aol::fileNameEndsWith( inputImageFileName.c_str(), "dm4" ) ) ) {
      qc::DM3Reader dmreader( inputImageFileName );
      dmreader.saveQuocDataInDM3Container ( deformedArray, outputFileName.c_str() );
    }
    else
      deformedArray.save ( outputFileName.c_str(), qc::PGM_DOUBLE_BINARY );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
