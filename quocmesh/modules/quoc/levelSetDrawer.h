#ifndef __LEVELSETDRAWER_H
#define __LEVELSETDRAWER_H

#include <gridBase.h>
#include <isolineIterator2d.h>
#include <scalarArray.h>
#include <multiArray.h>
#include <geom.h>

namespace qc {

// this is just the first prototype, it should implement behaviour
// as it is known from the colorbar (firstline, step, min, max) later
template <typename ConfiguratorType>
class LevelSetDrawer {
  const typename ConfiguratorType::InitType &_grid;
  typedef typename ConfiguratorType::RealType RealType;
public:
  LevelSetDrawer ( const typename ConfiguratorType::InitType &Grid )
      : _grid ( Grid ) {}

  void draw ( const ScalarArray<RealType, qc::QC_2D> &LevelSetFunction,
              ScalarArray<RealType, qc::QC_2D> &Image,
              RealType IsoValue = 0. ) const {
    IsoLineManager2d<ConfiguratorType> isoManager ( _grid, LevelSetFunction );

    isoManager.setIsoValue ( IsoValue );

    ScalarArray<RealType, qc::QC_2D> dist ( _grid );

    dist.setAll ( 10. ); // just must be larger than diam of reference element

    for ( IsoLineIterator2d<ConfiguratorType> isoIt ( isoManager.begin(), false ); isoIt != isoManager.end(); isoIt++ ) {
      aol::LineSegment<RealType, qc::QC_3D> segment ( isoIt->points[0], isoIt->points[1] );

      for ( int i = 0; i < 2; i++ ) {
        for ( int j = 0; j < 2; j++ ) {
          short x = isoIt->intersect[i][j][0];
          short y = isoIt->intersect[i][j][1];
          aol::Vec3<RealType> pt ( x, y, 0. );

          RealType d = segment.dist ( pt ) / sqrt ( 2. );

          if ( dist.get ( x, y ) > d ) {
            dist.set ( x, y, d );
            // _tmp[0].set( x, y, (1.-d) + _u_data.get(x,y)*( 2.*d-1.) ); // invert
            Image.set ( x, y, ( 1. - d ) + Image.get ( x, y ) * d );
          }
        }
      }
    }
  }

  /*
   * Draws a levelline of a levelset function onto a color image. ColorChannel
   * specifies the channel in which the levelline is drawn, e.g. if
   * ColorChannel == 0, the levelline is drawn red.
   *
   * \author Berkels
   */
  void draw ( const ScalarArray<RealType, qc::QC_2D> &LevelSetFunction,
              MultiArray<RealType, 2, 3> &VImage,
              RealType IsoValue,
              const int ColorChannel ) const {
    const int numX = LevelSetFunction.getNumX();
    const int numY = LevelSetFunction.getNumY();
    const RealType vImageMaxValue = VImage.getMaxValue();

    if ( vImageMaxValue < 0. )
      throw aol::Exception ( "LevelSetDrawer::draw should not be used on images with only negative values.", __FILE__, __LINE__ );

    qc::ScalarArray<RealType, qc::QC_2D> levellineImage ( numX, numY );

    draw ( LevelSetFunction, levellineImage, IsoValue );

    for ( int y = 0; y < numY; y++ ) {
      for ( int x = 0; x < numX; x++ ) {
        RealType levellineImageValue = aol::Clamp( levellineImage.get ( x, y ), aol::ZOTrait<RealType>::zero, aol::ZOTrait<RealType>::one );
        // If levellineImageValue is zero, we don't have to alter VImage at this position.
        if( levellineImageValue != 0. ){
          // If levellineImageValue exceeds a certain threshold, we want the line to be drawn
          // in the pure color chosen by ColorChannel, therefore we have to set the values
          // in the other channels to zero.
          if ( levellineImageValue > 0.7 ){
            VImage[(ColorChannel+2)%3][y*numX+x] = 0.;
            VImage[(ColorChannel+1)%3][y*numX+x] = 0.;
          }
          // Add the levellineValue into the chosen ColorChannel of VImage and clamp the
          // result into [0,vImageMaxValue]. levellineImageValue lives in [0.,1.], so we
          // recale it to be in [0,vImageMaxValue].
          VImage[ColorChannel%3][y*numX+x]
            = aol::Clamp(vImageMaxValue*levellineImageValue + VImage[ColorChannel%3][y*numX+x], aol::ZOTrait<RealType>::zero, vImageMaxValue );
        }
      }
    }
  }

  void drawToGnuplot ( const ScalarArray<RealType, qc::QC_2D> &LevelSetFunction,
                       ofstream &OutFile,
                       RealType IsoValue = 0.,
                       const bool ScaleToUnitDomain = false,
                       const bool MirrorAtYAxis = false ) const {
    RealType h_x = 1.0;
    RealType h_y = 1.0 ;
    if ( ScaleToUnitDomain ){
      h_x /= static_cast<RealType> ( LevelSetFunction.getNumX() - 1 );
      h_y /= static_cast<RealType> ( LevelSetFunction.getNumY() - 1 );
    }

    std::vector<qc::IsoFragmentInf<RealType, 2> > segments;

    for ( int i = 0; i <= 0; i++ ) {
      IsoLineManager2d<ConfiguratorType> isoManager ( _grid, LevelSetFunction );
      isoManager.setIsoValue ( IsoValue + static_cast< RealType > ( i ) * 0.5 / LevelSetFunction.getNumX() );

      for ( IsoLineIterator2d<ConfiguratorType> isoIt = isoManager.begin(); isoIt != isoManager.end(); isoIt++ )
        segments.push_back ( *isoIt );
    }

    aol::RandomAccessContainer< aol::RandomAccessContainer<aol::Vec2<RealType> > > curves;
    while ( segments.size() > 0 ) {
      curves.pushBack ( aol::RandomAccessContainer<aol::Vec2<RealType> > () );
      extractCurve ( segments, curves[curves.size()-1] );
    }

    for ( int k = 0; k < curves.size(); ++k ) {
      aol::RandomAccessContainer<aol::Vec2<RealType> > &curve = curves[k];
      for ( int i = 0; i < curve.size(); ++i ) {
        if ( MirrorAtYAxis )
          curve[i][1] = LevelSetFunction.getNumY()-1-curve[i][1];
        curve[i][0] = h_x*curve[i][0];
        curve[i][1] = h_y*curve[i][1];
        OutFile << curve[i][0] << " " << curve[i][1] << endl;
      }
        OutFile << endl;
    }
  }
  void drawToGnuplot ( const ScalarArray<RealType, qc::QC_2D> &LevelSetFunction,
                       const char *FileName,
                       RealType IsoValue = 0.,
                       const bool ScaleToUnitDomain = false,
                       const bool MirrorAtYAxis = false ) const {
    ofstream out ( FileName );
    drawToGnuplot ( LevelSetFunction, out, IsoValue, ScaleToUnitDomain, MirrorAtYAxis );
    out.close();
  }

private:
  void extractCurve ( std::vector<qc::IsoFragmentInf<RealType, 2> > &Segments, aol::RandomAccessContainer<aol::Vec2<RealType> > &Curve ) const {
    for ( int i = 0; i < 2; ++i )
      Curve.pushBack ( aol::Vec2<RealType> ( Segments[0].points[i][0], Segments[0].points[i][1] ) );

    Segments.erase ( Segments.begin() );

    while ( Curve[0] != Curve[Curve.size()-1] ) {
      bool nextPointFound = false;
      for ( typename std::vector<qc::IsoFragmentInf<RealType, 2> >::iterator it = Segments.begin(); it != Segments.end(); ++it ) {
        // Invalid segment, remove it and restart search for the next segment.
        if ( it->points[0] == it->points[1] ) {
          Segments.erase ( it );
          nextPointFound = true;
          break;
        }
        aol::Vec2<RealType> startPoint( it->points[0][0], it->points[0][1] );
        aol::Vec2<RealType> endPoint( it->points[1][0], it->points[1][1] );
        if ( Curve[Curve.size()-1] == startPoint ) {
          Curve.pushBack ( endPoint );
          Segments.erase ( it );
          nextPointFound = true;
          break;
        }
        else if ( Curve[Curve.size()-1] == endPoint ) {
          Curve.pushBack ( startPoint );
          Segments.erase ( it );
          nextPointFound = true;
          break;
        }
      }
      if ( nextPointFound == false )
        break;
    }
  }
};

}

#endif
