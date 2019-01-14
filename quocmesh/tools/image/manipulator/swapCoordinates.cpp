/** \file
 *  \brief Swap coordinates
 *
 *  Swap coordinates within an image (2D / 3D), i.e. rotate 90 degrees while changing size of the dimensions accordingly. Does NOT resample. \\
 *  Usage: swapCoordinates inputPath outputPath coord1 coord2 \\
 *
 *  \author Niklas Mevenkamp
 */

#include <aol.h>
#include <quoc.h>
#include <scalarArray.h>

typedef double RealType;

template <qc::Dimension Dim, typename PictureType>
void swapCoordinates ( const PictureType &Input, PictureType &Output, const int Coord1, const int Coord2 ) {
  if ( Dim == qc::QC_2D || Dim == qc::QC_3D ) {
    if ( Coord1 < 0 || Coord1 > 2 || Coord2 < 0 || Coord2 > 2 ) throw aol::Exception ( "Coordinates invalid!", __FILE__, __LINE__ );
    if ( Dim == qc::QC_2D && ( Coord1 > 1 || Coord2 > 1 ) ) throw aol::Exception ( "Coordinates invalid for 2D!", __FILE__, __LINE__ );
    
    if ( Coord1 != Coord2 ) {
      int coord1 = aol::Min<int> ( Coord1, Coord2 );
      int coord2 = aol::Max<int> ( Coord1, Coord2 );
      
      aol::Vec3<int> size = Input.getSize ( );
      const int tmp = size[Coord1];
      size[Coord1] = size[Coord2];
      size[Coord2] = tmp;
      
      qc::GridSize<Dim> gridSize ( size );
      Output.reallocate ( gridSize );
      
      for ( int z=0; z<Output.getNumZ ( ) ; ++z ) {
        for ( int y=0; y<Output.getNumY ( ) ; ++y ) {
          for ( int x=0; x<Output.getNumX ( ) ; ++x ) {
            qc::CoordType xOut ( x, y, z );
            qc::CoordType xIn ( xOut );
            xIn[coord1] = xOut[coord2];
            xIn[coord2] = xOut[coord1];
            Output.set ( xOut, Input.get ( xIn ) );
          }
        }
      }
    }
  }
}

// get the first character of the header
char getFirstHeaderChar ( const char* fileName ) {
  aol::ipfstream file ( fileName );
  qc::ArrayHeader header;
  qc::ReadArrayHeader ( file, header );
  file.close();
  
  return header.magic[1];
}

int main ( int argc, char **argv ) {
  try {
    if ( argc < 6 ) {
      cerr << "usage: " << argv[0] << " <inputPath> <outputPath> <dimension> <coord1> <coord2>\n";
      return EXIT_FAILURE;
    }
    
    const int dim = atoi ( argv[3] );

    if ( dim != 2 && dim != 3 ) cerr << "Dimension needs to be either 2 or 3!\n";
    else {
      const int coord1 = atoi ( argv[4] );
      const int coord2 = atoi ( argv[5] );
      
      // ------------------------- 3D data -----------------------------------
      if ( dim == 3 ) {
        qc::ScalarArray<RealType, qc::QC_3D> inputImg ( argv[1] ), outputImg;
        swapCoordinates<qc::QC_3D> ( inputImg, outputImg, coord1, coord2 );
        outputImg.save ( argv[2], qc::SaveTypeTrait<RealType>::BinarySaveType );
      }

      // ------------------------- 2D-Data -----------------------------------
      if ( dim == 2 ) {
        qc::ScalarArray<RealType, qc::QC_2D> inputImg ( argv[1] ), outputImg;
        swapCoordinates<qc::QC_3D> ( inputImg, outputImg, coord1, coord2 );
        outputImg.save ( argv[2], qc::SaveTypeTrait<RealType>::BinarySaveType );
      }
    }
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
}
