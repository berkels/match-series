#ifndef __IMAGETOOLS_H
#define __IMAGETOOLS_H

#include <array.h>
#include <matrix.h>
#include <kernel3d.h>
#include <auxiliary.h>
#include <bzipiostream.h>
#include <multiArray.h>
#include <rectangularGrid.h>
#include <gnuplotter.h>
#include <parameterParser.h>
#include <convolution.h>
#include <rgbColorMap.h>
#include <connectedComponents.h>
#include <generator.h>

namespace qc {

enum ConcatType { LEFT_RIGHT, TOP_BOTTOM };

/** Class that allows the concatenation of ScalarArray<QC_2D>-objects for the purpose of output as a PGM image.
 *  It the ScalarArrays are of different sizes the images will be aligned top and right in the output image.
 *  Any background will be filled with the specified color.
 *  Please note that the recursive use of ConcatImage is possible. I.e. you can concatenate ConcatImages with
 *  this class as well.
 *  Please note further that this class is meaningful only for output as PGM image. In fact to save memory the
 *  implementation does some nasty things in terms of memory management.
 */

template<typename DataType, typename RealType = float>
class ConcatImage : public ScalarArray<DataType, qc::QC_2D> {
protected:
  const ScalarArray<DataType, qc::QC_2D> &_image1, &_image2;
  const ConcatType          _type;
  const int                 _spacing;
  const unsigned char       _bgColor;
public:
  ConcatImage ( ScalarArray<DataType, qc::QC_2D> &image1,
                ScalarArray<DataType, qc::QC_2D> &image2,
                const ConcatType type = LEFT_RIGHT,
                const int spacing = 1, const unsigned char bgColor = 255 ) :
      ScalarArray<DataType, qc::QC_2D> ( 1, 1 ),
      _image1 ( image1 ), _image2 ( image2 ), _type ( type ), _spacing ( spacing ), _bgColor ( bgColor ) {
    switch ( type ) {
    case LEFT_RIGHT :
      this->numX = image1.getNumX() + spacing + image2.getNumX();
      this->numY = aol::Max ( image1.getNumY(), image2.getNumY() );
      break;
    case TOP_BOTTOM:
      this->numX = aol::Max ( image1.getNumX(), image2.getNumX() );
      this->numY = image1.getNumY() + spacing + image2.getNumY();
      break;
    }
  }

  virtual bool createOverflowHandeledData ( unsigned char *tmp, DataType * /*min*/ = NULL, DataType * /*max*/ = NULL ) const {
    const int w1 = _image1.getNumX();
    const int w2 = _image2.getNumX();
    const int h1 = _image1.getNumY();
    const int h2 = _image2.getNumY();

    const int w = this->getNumX();
    const int h = this->getNumY();

    typedef unsigned char UCHAR;

    UCHAR *img1 = new UCHAR[w1*h1];
    UCHAR *img2 = new UCHAR[w2*h2];

    _image1.createOverflowHandeledData ( img1, NULL, NULL );
    _image2.createOverflowHandeledData ( img2, NULL, NULL );

    memset ( tmp, _bgColor, w*h );

    int x, y;

    switch ( _type ) {
    case LEFT_RIGHT: {
      UCHAR *ptr = tmp, *ptr1 = img1, *ptr2 = img2;
      const int c = aol::Min ( h1, h2 );
      const int d = aol::Max ( h1, h2 );
      for ( y = 0; y < c; ++y ) {
        for ( x = 0; x < w1; ++x ) *ptr++ = *ptr1++;
        ptr += _spacing;
        for ( x = 0; x < w2; ++x ) *ptr++ = *ptr2++;
      }
      if ( c == d ) break;
      if ( c == h1 ) {
        for ( y = c; y < d; ++y ) {
          ptr += ( w1 + _spacing );
          for ( x = 0; x < w2; ++x ) *ptr++ = *ptr2++;
        }
      } else {
        for ( y = c; y < d; ++y ) {
          for ( x = 0; x < w1; ++x ) *ptr++ = *ptr1++;
          ptr += ( w2 + _spacing );
        }
      }
    }
    break;
    case TOP_BOTTOM: {
      UCHAR *ptr = tmp, *ptr1 = img1, *ptr2 = img2;
      const int c = w - w1;
      const int d = w - w2;
      for ( y = 0; y < h1; ++y ) {
        for ( x = 0; x < w1; ++x ) *ptr++ = *ptr1++;
        ptr += c;
      }
      ptr += ( _spacing * w );
      for ( y = 0; y < h2; ++y ) {
        for ( x = 0; x < w2; ++x ) *ptr++ = *ptr2++;
        ptr += d;
      }
    }
    break;
    };

    delete[] img1;
    delete[] img2;

    return false;
  }

};

template <typename RealType, typename GridType>
void write1dImage ( const GridType &Grid, const aol::Vector< RealType > &Image, const char* baseFileName, const bool /*WriteSlices*/ = 0, const Comp /*Direction*/ = qc::QC_X ) {
  if ( Grid.getDimOfWorld() == 1 ) {
    char fileName[1024];
    sprintf ( fileName, "%s%s", baseFileName, qc::getDefaultArraySuffix( qc::QC_1D ) );
    ScalarArray<RealType, qc::QC_1D> ImageArray ( Image, Grid.getNumX() );
    ImageArray.save ( fileName, PGM_DOUBLE_BINARY );

    aol::PlotDataFileHandler<RealType> plotHandler;
    plotHandler.generateFunctionPlot ( ImageArray );
    aol::Plotter<RealType> plotter;
    plotter.addPlotCommandsFromHandler ( plotHandler );
    plotter.set_outfile_base_name ( baseFileName );
    plotter.genPNG();
  }
  else
    throw aol::Exception ( "You can't execute write1dImage on a non-1D grid!\n" );
  }

template <typename RealType, typename GridType>
void write2dImage ( const GridType &Grid, const aol::Vector< RealType > &Image, const char* baseFileName, const bool /*WriteSlices*/ = 0, const Comp /*Direction*/ = qc::QC_X, const bool SaveAsDouble = false ) {
  if ( Grid.getDimOfWorld() == 2 ) {
    ScalarArray<RealType, qc::QC_2D> ImageArray ( Image, Grid.getNumX(), Grid.getNumY() );
    if ( SaveAsDouble )
      ImageArray.save ( ( string ( baseFileName ) + qc::getDefaultArraySuffix( qc::QC_2D ) ).c_str(), qc::PGM_DOUBLE_BINARY );
    else {
      ImageArray.setQuietMode ( true );
      ImageArray.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0., 1. );
      ImageArray.savePNG ( ( string ( baseFileName ) + ".png" ).c_str() );
    }
  }
  else
    throw aol::Exception ( "You can't execute write2dImage on a 3D grid!\n" );
}

template <typename RealType, typename GridType>
void write3dImage ( const GridType &Grid, const aol::Vector< RealType > &Image, const char* baseFileName, const bool WriteSlices = 0, const Comp Direction = qc::QC_X ) {
  if ( Grid.getDimOfWorld() == 3 ) {
    char fileName[1024];
    qc::ScalarArray<RealType, qc::QC_3D> *ImageArray = NULL;
    RealType minValue = Image.getMinValue();
    if ( ( minValue < 0 ) && WriteSlices ) {
      ImageArray = new qc::ScalarArray<RealType, qc::QC_3D> ( Image, Grid.getNumX(), Grid.getNumY(), Grid.getNumZ(), aol::DEEP_COPY );
      ImageArray->addToAll ( - minValue );
    } else {
      ImageArray = new qc::ScalarArray<RealType, qc::QC_3D> ( Image, Grid.getNumX(), Grid.getNumY(), Grid.getNumZ() );
    }
    ImageArray->setQuietMode ( true );
    if ( WriteSlices == 0 ) {
      sprintf ( fileName, "%s.dat.bz2", baseFileName );
      ImageArray->save ( fileName, PGM_DOUBLE_BINARY );
    } else {
      sprintf ( fileName, "%s_%%03d.pgm", baseFileName );
      ImageArray->saveSlices ( fileName, Direction, qc::PGM_UNSIGNED_CHAR_BINARY, NULL, aol::CLIP_THEN_SCALE, 0, 1 );
    }
    delete ImageArray;
  }
  else
    throw aol::Exception ( "You can't execute write3dImage on a 2D grid!\n" );
}

template <typename RealType, typename GridType>
void writeImage ( const GridType &Grid, const aol::Vector< RealType > &Image, const char* baseFileName, const bool WriteSlices = 0, const Comp Direction = qc::QC_X, const bool SaveAsDouble = false ) {
  if ( ( Grid.getDimOfWorld() != 2 ) && ( SaveAsDouble ) )
    throw aol::Exception ( "qc::writeImage: SaveAsDouble is currently only implemented in 2D", __FILE__, __LINE__ );

  if ( Grid.getDimOfWorld() == 1 ) {
    write1dImage<RealType, GridType> ( Grid, Image, baseFileName, WriteSlices, Direction );
  }
  else if ( Grid.getDimOfWorld() == 2 ) {
    write2dImage<RealType, GridType> ( Grid, Image, baseFileName, WriteSlices, Direction, SaveAsDouble );
  }
  else if ( Grid.getDimOfWorld() == 3 ) {
    write3dImage<RealType, GridType> ( Grid, Image, baseFileName, WriteSlices, Direction );
  }
  else
    throw aol::Exception ( "qc::writeImage: illegal dimension", __FILE__, __LINE__ );
}

template <typename RealType, typename GridType>
void writeImage ( const GridType &Grid, const aol::MultiVector< RealType > &MImage, const char* baseFileName, const bool WriteSlices = 0, const Comp Direction = qc::QC_X, const bool SaveAsDouble = false ) {
  char filename[1024];
  for ( int i = 0; i < MImage.numComponents(); i++ ) {
    sprintf ( filename, "%s_%d", baseFileName, i );
    qc::writeImage<RealType, GridType> ( Grid, MImage[i], filename, WriteSlices, Direction, SaveAsDouble );
  }
}

template <typename RealType, typename GridType>
void writeImage ( const GridType &Grid, const qc::MultiArray< RealType, 2, 1 > &MImage, const char* baseFileName, const bool WriteSlices = 0, const Comp Direction = qc::QC_X, const bool SaveAsDouble = false ) {
  qc::writeImage<RealType, GridType> ( Grid, MImage[0], baseFileName, WriteSlices, Direction, SaveAsDouble );
}

template <typename RealType, typename GridType>
void writeImage ( const GridType &/*Grid*/, const qc::MultiArray< RealType, 2, 3 > &MImage, const char* baseFileName, const bool /*WriteSlices*/ = 0, const Comp /*Direction*/ = qc::QC_X, const bool /*SaveAsDouble*/ = false ) {
  qc::MultiArray< RealType, 2, 3 > image ( MImage, aol::FLAT_COPY );
  image.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0., 1. );
  image.savePNG ( aol::strprintf ( "%s.png", baseFileName ).c_str() );
}

/**
 * \author Berkels
 */
template <typename RealType>
void plotLineOf2DArray ( const ScalarArray<RealType, qc::QC_2D> &Array, const int LineNumber, const char *BaseName, const aol::PlotOutFileType OutType ){
  const int numX = Array.getNumX();
  aol::Vector<RealType> slice ( numX );
  for ( int i = 0; i < numX; ++i )
    slice[i] = Array.get ( i, LineNumber );

  aol::Plotter<RealType> plotter;
  plotter.set_outfile_base_name( BaseName );
  aol::PlotDataFileHandler<RealType> plotHandler;
  plotHandler.generateFunctionPlot( slice );
  plotter.addPlotCommandsFromHandler( plotHandler );
  plotter.genPlot( OutType );
}

/**
 * Maps the scalar values of a ScalarArray from [0,1] to color values in a MultiArray
 * using the standard HSV hue scale.
 *
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
void convertScalarToColor ( const qc::ScalarArray<RealType, Dim> &Scalar, qc::MultiArray<unsigned char, Dim, 3> &BufArray ) {

  RealType rgb[3] = {0, 0, 0}, hsv[3];

  for ( int i = 0; i < Scalar.size(); ++i ) {
    RealType b = Scalar[i];

    while ( b > 1. ) b -= 1.;
    while ( b < 0. ) b += 1.;

    // A weighting to alter the brightness could be added here.
    const RealType w = 1;

    hsv[0] = b;
    hsv[1] = w;
    hsv[2] = w;

    aol::RGBColorMap<RealType>::hsv2rgb ( hsv, rgb );

    BufArray[0][i] = static_cast< unsigned char> ( rgb[0] * 255 );
    BufArray[1][i] = static_cast< unsigned char> ( rgb[1] * 255 );
    BufArray[2][i] = static_cast< unsigned char> ( rgb[2] * 255 );
  }
}
  
/**
  * Maps the scalar values of a ScalarArray to color values in a MultiArray
  * using the specified color transformation
  *
  * \author Tatano
  */
template <typename RealType, qc::Dimension Dim>
void convertScalarToColorUsingcolorTrans ( const qc::ScalarArray<RealType, Dim> &Scalar, qc::MultiArray<unsigned char, Dim, 3> &BufArray, const aol::colorTrans &colorTransformation, const RealType Min = aol::NumberTrait<RealType>::getNaN(), const RealType Max = aol::NumberTrait<RealType>::getNaN() ) {
  const RealType maxAbsScalar = Scalar.getMaxAbsValue();
  const RealType min = aol::isNaN ( Min ) ? -maxAbsScalar : Min;
  const RealType max = aol::isNaN ( Max ) ? maxAbsScalar : Max;
  const aol::RGBColorMap<RealType> hsvMap ( min, max, colorTransformation );
    
  for ( int i = 0; i < Scalar.size(); ++i ) {
    aol::Vec3< RealType > color;
    hsvMap.scalarToColor ( Scalar[i], color );
    for ( int j = 0; j < 3; ++j )
      BufArray[j][i] = static_cast< unsigned char> ( color[j] * 255 );
  }
}

/**
 * Algorithm for bicubic interpolation.
 *
 * \author Paul Breeuwsma
 *
 * Converted from Java by Berkels.
 * Original Java source available at http://www.paulinternet.nl/?page=bicubic
 */
template <typename RealType>
class CachedBicubicInterpolator {
  aol::Mat<4,4, RealType> _a;

public:
  void updateCoefficients ( const aol::Mat<4,4, RealType> &P ) {
    _a[0][0] = P[1][1];
    _a[0][1] = -.5*P[1][0] + .5*P[1][2];
    _a[0][2] = P[1][0] - 2.5*P[1][1] + 2*P[1][2] - .5*P[1][3];
    _a[0][3] = -.5*P[1][0] + 1.5*P[1][1] - 1.5*P[1][2] + .5*P[1][3];
    _a[1][0] = -.5*P[0][1] + .5*P[2][1];
    _a[1][1] = .25*P[0][0] - .25*P[0][2] - .25*P[2][0] + .25*P[2][2];
    _a[1][2] = -.5*P[0][0] + 1.25*P[0][1] - P[0][2] + .25*P[0][3] + .5*P[2][0] - 1.25*P[2][1] + P[2][2] - .25*P[2][3];
    _a[1][3] = .25*P[0][0] - .75*P[0][1] + .75*P[0][2] - .25*P[0][3] - .25*P[2][0] + .75*P[2][1] - .75*P[2][2] + .25*P[2][3];
    _a[2][0] = P[0][1] - 2.5*P[1][1] + 2*P[2][1] - .5*P[3][1];
    _a[2][1] = -.5*P[0][0] + .5*P[0][2] + 1.25*P[1][0] - 1.25*P[1][2] - P[2][0] + P[2][2] + .25*P[3][0] - .25*P[3][2];
    _a[2][2] = P[0][0] - 2.5*P[0][1] + 2*P[0][2] - .5*P[0][3] - 2.5*P[1][0] + 6.25*P[1][1] - 5*P[1][2] + 1.25*P[1][3] + 2*P[2][0] - 5*P[2][1] + 4*P[2][2] - P[2][3] - .5*P[3][0] + 1.25*P[3][1] - P[3][2] + .25*P[3][3];
    _a[2][3] = -.5*P[0][0] + 1.5*P[0][1] - 1.5*P[0][2] + .5*P[0][3] + 1.25*P[1][0] - 3.75*P[1][1] + 3.75*P[1][2] - 1.25*P[1][3] - P[2][0] + 3*P[2][1] - 3*P[2][2] + P[2][3] + .25*P[3][0] - .75*P[3][1] + .75*P[3][2] - .25*P[3][3];
    _a[3][0] = -.5*P[0][1] + 1.5*P[1][1] - 1.5*P[2][1] + .5*P[3][1];
    _a[3][1] = .25*P[0][0] - .25*P[0][2] - .75*P[1][0] + .75*P[1][2] + .75*P[2][0] - .75*P[2][2] - .25*P[3][0] + .25*P[3][2];
    _a[3][2] = -.5*P[0][0] + 1.25*P[0][1] - P[0][2] + .25*P[0][3] + 1.5*P[1][0] - 3.75*P[1][1] + 3*P[1][2] - .75*P[1][3] - 1.5*P[2][0] + 3.75*P[2][1] - 3*P[2][2] + .75*P[2][3] + .5*P[3][0] - 1.25*P[3][1] + P[3][2] - .25*P[3][3];
    _a[3][3] = .25*P[0][0] - .75*P[0][1] + .75*P[0][2] - .25*P[0][3] - .75*P[1][0] + 2.25*P[1][1] - 2.25*P[1][2] + .75*P[1][3] + .75*P[2][0] - 2.25*P[2][1] + 2.25*P[2][2] - .75*P[2][3] - .25*P[3][0] + .75*P[3][1] - .75*P[3][2] + .25*P[3][3];
  }

  RealType getValue ( const RealType X, const RealType Y ) const {
    RealType x2 = aol::Sqr( X );
    RealType x3 = x2 * X;
    RealType y2 = aol::Sqr( Y );
    RealType y3 = y2 * Y;

    return ( _a[0][0] + _a[0][1] * Y + _a[0][2] * y2 + _a[0][3] * y3 ) +
           ( _a[1][0] + _a[1][1] * Y + _a[1][2] * y2 + _a[1][3] * y3 ) * X +
           ( _a[2][0] + _a[2][1] * Y + _a[2][2] * y2 + _a[2][3] * y3 ) * x2 +
           ( _a[3][0] + _a[3][1] * Y + _a[3][2] * y2 + _a[3][3] * y3 ) * x3;
  }
};

/**
 * \author Berkels
 */
template <typename RealType>
RealType naiveNearestNeighborSearchDistance ( const std::set<aol::Vec<2, short> > &PointSet, const aol::Vec<2, short> &Point ) {
  RealType distance = aol::NumberTrait<RealType>::Inf;
  for ( std::set<aol::Vec<2, short> >::const_iterator it = PointSet.begin(); it != PointSet.end(); ++it ) {
    // short is too small to store the distance.
    const aol::Vec2<RealType> distanceVec ( (*it)[0] - Point[0], (*it)[1] - Point[1] );
    distance = aol::Min ( distance, distanceVec.normSqr() );
  }
  return sqrt ( distance );
}

/**
 * \author Berkels
 */
template <typename RealType>
void findMaxima ( const qc::ScalarArray<RealType, qc::QC_2D> &Image, const int PatchSize, aol::RandomAccessContainer<aol::Vec2<short> > &Maxima ) {
  const int offset = PatchSize/2;
  const qc::CoordType lower ( PatchSize, PatchSize );
  const qc::CoordType upper ( Image.getNumX() - 2*PatchSize + 1, Image.getNumY() - 2*PatchSize + 1 );
  qc::ScalarArray<RealType, qc::QC_2D> temp ( PatchSize, PatchSize );
  for ( qc::RectangularIterator<qc::QC_2D> it ( lower, upper ); it.notAtEnd(); ++it ) {
    Image.copyBlockTo ( *it, temp );
    if ( temp.get ( offset, offset ) >= temp.getMaxValue() )
      Maxima.pushBack ( aol::Vec2<short> ( (*it)[0] + offset, (*it)[1] + offset ) );
  }
}

/**
 * \author Berkels
 */
void cleanMask ( qc::BitArray<qc::QC_2D> &Mask, const int DropComponentsSmallerThan, const bool DropBoundaryComponents );

/**
  * \author Tatano
  */
void cleanMask3D ( qc::BitArray<qc::QC_3D> &Mask, const int DropComponentsSmallerThan, const bool DropBoundaryComponents );
  
/**
 * \author Tatano
 */
void getBiggestComponent ( qc::BitArray<qc::QC_2D> &Mask );

/**
 * \author Berkels
 */
void convertMaskToSet ( const qc::BitArray<qc::QC_2D> &Mask, std::set<aol::Vec<2, short> > &MaskSet );

/**
 * Very simple volume rendering using orthogonal projection and a basic emission / absorption model.
 *
 * \author Berkels
 */
template <typename RealType>
void renderVolume ( const qc::ScalarArray<RealType, qc::QC_3D> &Volume, const qc::Comp Direction, qc::ScalarArray<RealType, qc::QC_2D> &Projection ){
  int tmpNumX = 0, tmpNumY = 0, numOfSlices = 0;

  switch ( Direction ) {
    case qc::QC_X:
      tmpNumX = Volume.getNumY();
      tmpNumY = Volume.getNumZ();
      numOfSlices = Volume.getNumX();
      break;
    case qc::QC_Y:
      tmpNumX = Volume.getNumX();
      tmpNumY = Volume.getNumZ();
      numOfSlices = Volume.getNumY();
      break;
    case qc::QC_Z:
      tmpNumX = Volume.getNumX();
      tmpNumY = Volume.getNumY();
      numOfSlices = Volume.getNumZ();
      break;
    default:
      throw aol::UnimplementedCodeException ( "Unknown qc::Comp", __FILE__, __LINE__ );
  }

  Projection.reallocate ( tmpNumX, tmpNumY );

  // The following may look like an awful example of code duplication, but
  // the order of the loops is important for the execution speed. Furthermore,
  // it can't hurt not to have the switch in the inner most loop.
  switch ( Direction ) {
    case qc::QC_X:
      for ( int y = 0; y < tmpNumY; ++y ) {
        for ( int x = 0; x < tmpNumX; ++x ) {
          for ( int z = 0; z <numOfSlices; ++z ) {
            const RealType value = Volume.get ( z, x, y );
            const RealType opacity = value;
            Projection.set ( x, y, opacity * value + ( 1 - opacity ) * Projection.get ( x, y ) );
          }
        }
      }
      break;
    case qc::QC_Y:
      for ( int y = 0; y < tmpNumY; ++y ) {
        for ( int z = 0; z <numOfSlices; ++z ) {
          for ( int x = 0; x < tmpNumX; ++x ) {
            const RealType value = Volume.get ( x, z, y );
            const RealType opacity = value;
            Projection.set ( x, y, opacity * value + ( 1 - opacity ) * Projection.get ( x, y ) );
          }
        }
      }
      break;
    case qc::QC_Z:
      for ( int z = 0; z <numOfSlices; ++z ) {
        for ( int y = 0; y < tmpNumY; ++y ) {
          for ( int x = 0; x < tmpNumX; ++x ) {
            const RealType value = Volume.get ( x, y, z );
            const RealType opacity = value;
            Projection.set ( x, y, opacity * value + ( 1 - opacity ) * Projection.get ( x, y ) );
          }
        }
      }
      break;
    default:
      throw aol::UnimplementedCodeException ( "Unknown qc::Comp", __FILE__, __LINE__ );
  }
}

/**
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
void resampleArrayToQuadraticArrayKeepingAspectRatio ( const qc::ScalarArray<RealType, Dim> &InputArray, qc::ScalarArray<RealType, Dim> &OutputArray, const RealType FillValue = 0 ) {
  // OutputArray needs to be quadatic, otherwise the code below won't keep the aspect ratio.
  GridSize<Dim> sizeChecker ( OutputArray );
  sizeChecker.quadraticOrDie ();

  GridSize<Dim> paddedSize ( static_cast<short> ( aol::Max ( InputArray.getNumX(), InputArray.getNumY(), InputArray.getNumZ() ) ) );
  qc::ScalarArray<RealType, Dim> paddedArray ( paddedSize );
  paddedArray.padFrom ( InputArray, FillValue );
  OutputArray.resampleFrom ( paddedArray );
}

/**
 * \author Berkels
 */
template <typename RealType>
void joinTwo2DArraysVertically( const qc::ScalarArray<RealType, qc::QC_2D> &Arg1, const qc::ScalarArray<RealType, qc::QC_2D> &Arg2, qc::ScalarArray<RealType, qc::QC_2D> &Dest, const int Spacing = 0 ) {
  if( Arg1.getNumX() != Arg2.getNumX() )
    throw aol::Exception( "Arg1.getNumX() != Arg2.getNumX() !", __FILE__, __LINE__);

  Dest.reallocate( Arg1.getNumX(), Arg1.getNumY() + Spacing + Arg2.getNumY() );
  Dest.pasteFrom ( Arg1, 0, 0 );
  Dest.pasteFrom ( Arg2, 0, Arg1.getNumY() + Spacing );

  for( int j = 0; j < Spacing; j++ ) {
    for( int i = 0; i < Arg1.getNumX(); i++ ) {
      Dest.set( i, j + Arg1.getNumY(), 0.5 );
    }
  }
}

/**
 * \author Berkels
 */
template <typename RealType>
void shrinkAndPad2DArray ( qc::ScalarArray<RealType, qc::QC_2D> &Array, const RealType ScaleFactor, const RealType FillValue = 0 ) {
  qc::ScalarArray<RealType, qc::QC_2D> arrayScaled ( static_cast<int> ( Array.getNumX() * ScaleFactor ), static_cast<int> ( Array.getNumY() * ScaleFactor ) );
  arrayScaled.resampleFrom ( Array );
  Array.padFrom ( arrayScaled, FillValue );
}

/**
 * \author Berkels
 */
template <typename RealType>
void crop2DImage ( qc::ScalarArray<RealType, qc::QC_2D> &Image, const aol::ParameterParser &Parser ) {
  if ( Parser.checkAndGetBool ( "cropInput" ) ) {
    qc::ScalarArray<RealType, qc::QC_2D> imageCropped ( Parser.getInt ( "cropSizeX" ), Parser.getInt ( "cropSizeY" ) );
    Image.copyBlockTo ( Parser.getInt ( "cropStartX" ), Parser.getInt ( "cropStartY" ) , imageCropped );
    Image.reallocate ( imageCropped );
    Image = imageCropped;
  }
}
  
/**
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
void applyStretchMute ( qc::ScalarArray<RealType, qc::QC_2D> &DeformedImage,
                        const qc::MultiArray<RealType, Dim> &Phi,
                        const RealType XDerivNormThreshold,
                        const RealType ExtensionConstant = aol::NumberTrait<RealType>::Inf ) {
  qc::GridSize<Dim> gridSize = qc::GridSize<Dim>::createFrom ( DeformedImage );
  qc::ScalarArray<RealType, Dim> xDerivNorm ( gridSize );
  qc::ScalarArrayHelper<qc::ScalarArray<RealType, Dim> >::getPointWiseXDerivNorm ( Phi, xDerivNorm );
  typename qc::BitArray<Dim> maskedPositions ( gridSize );
  maskedPositions.thresholdFrom ( xDerivNorm, XDerivNormThreshold );
  maskedPositions.erodeByOne();
  maskedPositions.dilateByOne();
  for ( int k = 0; k < xDerivNorm.size(); ++k )
    if ( maskedPositions[k] )
      DeformedImage[k] = ExtensionConstant;
}

/**
 * \author Berkels
 */
int computeMean ( const aol::Vector<int> &Histo, int StartIndex, int EndIndex );

/**
 * \author Berkels
 */
int computeFirstMoment ( const aol::Vector<int> &Histo, int StartIndex, int EndIndex );

/**
 * \author Berkels
 */
template <typename RealType>
RealType computeEpsilonEq ( const RealType XdX, const RealType XdY, const RealType YdX, const RealType YdY ) {
  const RealType epsilonXX = XdX + 0.5 * ( aol::Sqr ( XdX ) + aol::Sqr ( YdX ) );
  const RealType epsilonYY = YdY + 0.5 * ( aol::Sqr ( XdY ) + aol::Sqr ( YdY ) );
  const RealType epsilonXY = 0.5 * ( XdY + YdX + XdX*XdY + YdX*YdY );
  return sqrt ( aol::Sqr ( epsilonXX ) + 2*aol::Sqr ( epsilonXY ) + aol::Sqr ( epsilonYY ) );
}

/**
 * \author Berkels
 */
template <typename RealType>
void computeEpsilonEq ( const qc::MultiArray<RealType, qc::QC_2D> &GradU,
                        const qc::MultiArray<RealType, qc::QC_2D> &GradV,
                        qc::ScalarArray<RealType, qc::QC_2D> &EpsilonEq ) {
  for ( int i = 0; i < GradU[0].size(); ++i )
    EpsilonEq[i] = computeEpsilonEq ( GradU[0][i], GradU[1][i], GradV[0][i], GradV[1][i] );
}
  
/**
 *\brief Finds the curve points of a black-white image of the curve and stores them in a MultiVector
 * \author Tatano
 */
template <typename ConfiguratorType>
  void getCurveFromImage(const typename ConfiguratorType::InitType &grid, const qc::ScalarArray<typename ConfiguratorType::RealType, qc::QC_2D> &image, aol::MultiVector<typename ConfiguratorType::RealType> &points){
  qc::ScalarArray<typename ConfiguratorType::RealType, qc::QC_2D> temp(image, aol::DEEP_COPY);
  temp.scaleValuesTo01();
  points.reallocate(2, 0);
  const ConfiguratorType config(grid);
  typename ConfiguratorType::InitType::OldAllNodeIterator fnit;
  for( fnit = grid._nBeginIt; fnit != grid._nEndIt; ++fnit){
    if ( temp.get( ( *fnit ) [0], ( *fnit ) [1] ) == 1.0 ){
      points[0].pushBack( ( *fnit ) [0] * grid.H() );
      points[1].pushBack( ( *fnit ) [1] * grid.H() );
    }
  }
}
  
/**
 *\brief Sorts the points of a curve C given as MultiVector of points c=(x,y)
 * \author Tatano
 */
template <typename RealType>
void sortCurvePoints(const aol::MultiVector<RealType> &points, aol::MultiVector<RealType> &sortedPoints){
  sortedPoints.reallocate(points);
  aol::MultiVector<RealType> temp(points);
  aol::Vector<RealType> dist(points[0].size()-1);
  aol::Vec2<RealType> pos;
  pos[0] = sortedPoints[0][0] = points[0][0];
  pos[1] = sortedPoints[1][0] = points[1][0];
  int l=1;
  temp[0].erase(0);
  temp[1].erase(0);
  int minInd = 0;
  while (temp[0].size()>0) {
    for (int i=0; i<temp[0].size(); i++) {
      dist[i] = sqrt(aol::Sqr(pos[0] - temp[0][i]) + aol::Sqr(pos[1] - temp[1][i]));
    }
    minInd = dist.getMinIndexAndValue().first;
    pos[0] = sortedPoints[0][l] = temp[0][minInd];
    pos[1] = sortedPoints[1][l] = temp[1][minInd];
    l++;
    temp[0].erase(minInd);
    temp[1].erase(minInd);
    dist.reallocate(temp[0].size());
  }
 /* dist.reallocate(points[0].size()-1);
  RealType distLP = sqrt(aol::Sqr(sortedPoints[0][points[0].size()-1] - sortedPoints[0][0]) + aol::Sqr(sortedPoints[1][points[0].size()-1] - sortedPoints[1][0]) );
  int iter = 0;
  while (distLP > 0.1 && iter<points[0].size()-1 ) {
    pos[0] = sortedPoints[0][points[0].size()-1];
    pos[1] = sortedPoints[1][points[0].size()-1];
    for (int i=0; i<points[0].size()-1; i++) {
      dist[i] = sqrt(aol::Sqr(pos[0] - sortedPoints[0][i]) + aol::Sqr(pos[1] - sortedPoints[1][i]));
    }
    cerr << dist << endl;
    minInd = dist.getMinIndexAndValue().first;
    cerr << minInd << endl;
    int posInt = points[0].size()-1;
    sortedPoints[0].insert(minInd, pos[0]);
    sortedPoints[1].insert(minInd, pos[1]);
    sortedPoints[0].erase(posInt);
    sortedPoints[1].erase(posInt);
    distLP = sqrt(aol::Sqr(sortedPoints[0][points[0].size()-1] - sortedPoints[0][0]) + aol::Sqr(sortedPoints[1][points[0].size()-1] - sortedPoints[1][0]) );
    iter++;
    dist.reallocate(points[0].size()-1);
  }*/
  /*dist.reallocate(temp[0].size());
  for (int i=0; i<dist.size(); i++) {
    dist[i] = sqrt(aol::Sqr(pos[0] - temp[0][i]) + aol::Sqr(pos[1] - temp[1][i]));
  }
  minInd = dist.getMinIndexAndValue().first;
  sortedPoints[0][1] = temp[0][minInd];
  sortedPoints[1][1] = temp[1][minInd];
  temp[0].erase(minInd);
  temp[1].erase(minInd);
  l=2;
  dist.reallocate( points[0].size() - temp[0].size());
  aol::Vector<RealType> minDist(temp[0].size());
  aol::Vector<int> minIndDist(temp[0].size());
  while (temp[0].size()>0) {
    for (int j=0; j<temp[0].size(); j++) {
      for (int i=0; i<dist.size(); i++)
        dist[i] = sqrt(aol::Sqr(temp[0][j] - sortedPoints[0][i]) + aol::Sqr(temp[1][j] - sortedPoints[1][i]));
      minDist[j] = dist.getMinValue();
      minIndDist[j] = dist.getMinIndexAndValue().first;
    }
    minInd = minDist.getMinIndexAndValue().first;
    l = minIndDist[minInd];
    RealType distPrec = sqrt(aol::Sqr(temp[0][minInd] - sortedPoints[0][l-1]) + aol::Sqr(temp[1][minInd] - sortedPoints[1][l-1]));
    RealType distSuc = sqrt(aol::Sqr(temp[0][minInd] - sortedPoints[0][l]) + aol::Sqr(temp[1][minInd] - sortedPoints[1][l]));
    if (distPrec < distSuc) {
      sortedPoints[0].insert(l, temp[0][minInd]);
      sortedPoints[1].insert(l, temp[1][minInd]);
    }else{
      sortedPoints[0].insert(l+1, temp[0][minInd]);
      sortedPoints[1].insert(l+1, temp[1][minInd]);
    }
    sortedPoints[0].erase(points[0].size());
    sortedPoints[1].erase(points[1].size());
    temp[0].erase(minInd);
    temp[1].erase(minInd);
    dist.reallocate(sortedPoints[0].size() - temp[0].size());
    minDist.reallocate(temp[0].size());
    minIndDist.reallocate(temp[0].size());
  }*/
}

} // end of namespace qc.

#endif
