/**
 * \file
 * \brief Computes a*image + b (works both in 2d and 3d). In 2D, an integer coordinate offset is supported.
 *
 * Usage: multShiftImg inputImg outputImg a b [xoffset2D yoffset2D]
 *
 * \warning The output-image will only be saved with double values if the suffix is '.bz2'!
 *
 * \author Nemitz
 */

#include <quoc.h>
#include <scalarArray.h>
#include <scalarArray.h>

#include <qmException.h>
#include <auxiliary.h>

typedef double RealType;
using namespace aol::color;

int main ( int argc, char **argv ) {
  try {
    if ( argc < 5 ) {
      cerr << "Computes a*image + b (either in 2d or in 3d).\n";
      cerr << "Attention: If the suffix of the output is '.bz2', it will be saved with double values, otherwise not!\n\n";
      cerr << aol::color::red << "usage: " << argv[0] << " input output a b [xoffset2D yoffset2D]\n" << aol::color::reset;
      return EXIT_FAILURE;
    }

    const double a = atof ( argv[3] );
    const double b = atof ( argv[4] );

    const qc::Dimension dim = qc::getDimensionFromArrayFile ( argv[1] );


    // ------------------------- --------------------------------------------
    // ------------------------- 3D-Daten -----------------------------------
    // ----------------------------------------------------------------------
    if ( dim == qc::QC_3D ) {
      const qc::ScalarArray<RealType, qc::QC_3D> inputimg ( argv[1] );
      qc::ScalarArray<RealType, qc::QC_3D> outputimg ( inputimg );


      cerr << blue << "Multiplying 3d-volume with " << red << a << blue << " and adding " << red << b << blue << " to it...";

      outputimg *= a;
      outputimg.addToAll ( b );

      cerr << "done! Saving...\n" << reset;

      outputimg.save ( argv[2], qc::SaveTypeTrait<RealType>::BinarySaveType );
    }


    // ------------------------- --------------------------------------------
    // ------------------------- 2D-Daten -----------------------------------
    // ----------------------------------------------------------------------
    if ( dim == qc::QC_2D ) {
      const qc::ScalarArray<RealType, qc::QC_2D> inputimg ( argv[1] );
      qc::ScalarArray<RealType, qc::QC_2D> outputimg ( inputimg );

      if ( argc >= 7 ) {
        const qc::CoordType offset ( atoi ( argv[5] ), atoi ( argv[6] ) );
        cerr << "Shifting coordinates by offset " << offset << endl;
        outputimg.shiftByOffsetFrom ( offset, inputimg );
      }

      cerr << blue << "Multiplying 2d-image with " << red << a << blue << " and adding " << red << b << blue << " to it...";

      outputimg *= a;
      outputimg.addToAll ( b );

      cerr << "done! Saving...\n" << reset;

      qc::recognizeEndingAndSave2d ( outputimg, argv[2] );
    }

    cerr << blue << "done! Thanx for using multShiftImg.cpp.\n";
    cerr << reset;
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
}
