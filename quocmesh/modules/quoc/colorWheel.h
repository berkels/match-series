#ifndef __COLORWHEEL_H
#define __COLORWHEEL_H

#include <scalarArray.h>
#include <multiArray.h>
#include <rgbColorMap.h>

namespace qc {

template <typename RealType>
class ColorWheel {
  const int _numX, _numY;
  qc::ScalarArray<RealType, qc::QC_2D> _scaledInputArray;
  const RealType _lowerAngle, _upperAngle;
  const aol::RGBColorMap<RealType> _colorMap;
public:
  ColorWheel ( const qc::ScalarArray<RealType, qc::QC_2D> &AngleArray, const RealType LowerAngle, const RealType UpperAngle, const bool Rescale )
      : _numX ( AngleArray.getNumX() ),
      _numY ( AngleArray.getNumY() ),
      _scaledInputArray ( AngleArray ),
      _lowerAngle ( ( Rescale ) ? _scaledInputArray.getMinValue() : LowerAngle ),
      _upperAngle ( ( Rescale ) ? _scaledInputArray.getMaxValue() : UpperAngle ),
      _colorMap ( 0., 1., aol::HSV_MAGENTA_TO_RED ) {

    if ( Rescale ) {
      _scaledInputArray.addToAll ( - _lowerAngle );
    }

    _scaledInputArray /= _scaledInputArray.getMaxValue();
  }
  void saveColoredWheel ( const char *Filename ) {
    qc::MultiArray<unsigned char, 2, 3> bufArray ( _numX, _numY );

    convertScalarToColor ( _scaledInputArray, bufArray, _lowerAngle, _upperAngle );

    colwheel ( bufArray, 40*_numX / 513, _numY - 80*_numY / 513, 60*_numX / 513, -aol::NumberTrait<long double>::pi / 4, aol::NumberTrait<long double>::pi / 4, _numX, _numY );
    FILE* fp = fopen ( Filename, "wb" );
    bufArray.save ( fp );
    fclose ( fp );
  }
  void saveGrayWheel ( const char *Filename ) {
    qc::ScalarArray<RealType, qc::QC_2D> gray ( _scaledInputArray );
    colwheel_gray ( gray, 40*_numX / 513, _numY - 80*_numY / 513, 60*_numX / 513, _lowerAngle, _upperAngle, _numX, _numY );
    gray.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0., 1. );
    gray.save ( Filename, qc::PGM_UNSIGNED_CHAR_BINARY );
  }
private:
  void colwheel ( qc::MultiArray<unsigned char, 2, 3> &BufArray, int cx, int cy, int r, RealType a1, RealType a2, int numX, int numY );
  void colwheel_gray ( qc::ScalarArray<RealType, qc::QC_2D> &out, int cx, int cy, int r, RealType a1, RealType a2, int numX, int numY );
  void convertScalarToColor ( const qc::ScalarArray<RealType, qc::QC_2D> &Scalar, qc::MultiArray<unsigned char, 2, 3> &BufArray, const RealType a1, const RealType a2 ) {
    const int numX = Scalar.getNumX();
    const int numY = Scalar.getNumY();
    for ( int y = 0; y < numY; y++ ) {
      for ( int x = 0; x < numX; x++ ) {
        RealType b = Scalar.get ( x, y );


        b = ( a2 - a1 ) * b + a1;

        if ( b > aol::NumberTrait<long double>::pi / 4 ) {
          b = aol::NumberTrait<long double>::pi / 4;
        }

        if ( b < -aol::NumberTrait<long double>::pi / 4 ) {
          b = -aol::NumberTrait<long double>::pi / 4;
        }

        b = -b / ( aol::NumberTrait<long double>::pi / 2 ) + 0.5;


        aol::Vec3< RealType > color;
        _colorMap.scalarToColor ( b, color );
        for ( int i = 0; i < 3; ++i )
          BufArray[i][y*numX+x] = static_cast<unsigned char> ( color[ i ] * 255 );
      }
    }
  }
};

template <typename RealType>
void ColorWheel<RealType>::colwheel ( qc::MultiArray<unsigned char, 2, 3> &BufArray, int cx, int cy, int r, RealType a1, RealType a2, int numX, int /*numY*/ ) {

  aol::Vec3< RealType > rgb;
  int i = 0;

  a1 /= 2 * aol::NumberTrait<long double>::pi;
  a2 /= 2 * aol::NumberTrait<long double>::pi;

  for ( int x = cx - r;x <= cx + r;x++ ) {
    for ( int y = cy - r;y <= cy + r;y++, i++ ) {
      RealType vx = x - cx, vy = y - cy;
      RealType  n = sqrt ( vx * vx + vy * vy );
      if ( n > r )
        continue;
      if ( n > 1.0e-6 ) {
        vx /= n;
        vy /= n;
        RealType angle = - ( ( vy >= 0 ) ? acosf ( vx ) : aol::NumberTrait<long double>::pi + acosf ( -vx ) ) / ( 2 * aol::NumberTrait<long double>::pi ); //angle of v with y-axis, in [0..1]
        if ( angle < 0. )
          angle += 1.;
        RealType b = angle;

        if ( b > a2 )
          b -= 1.;
        if ( b < a1 || b > a2 ) {
          rgb[0] = rgb[1] = rgb[2] = 0;
          continue;
        } else {
          RealType a = ( b - a1 ) / ( a2 - a1 );
          _colorMap.scalarToColor ( a, rgb );
        }
      } else {
        rgb[0] = rgb[1] = rgb[2] = 0;
        continue;
      }

      int pos = ( y * numX + x );
      BufArray[0][pos]   = static_cast< unsigned char > ( rgb[0] * 255 );
      BufArray[1][pos]   = static_cast< unsigned char > ( rgb[1] * 255 );
      BufArray[2][pos]   = static_cast< unsigned char > ( rgb[2] * 255 );
    }
  }
}


template <typename RealType>
void ColorWheel<RealType>::colwheel_gray ( qc::ScalarArray<RealType, qc::QC_2D> &out, int cx, int cy, int r, RealType a1, RealType a2, int /*numX*/, int /*numY*/ ) {
  int i = 0;
  a1 /= 2 * aol::NumberTrait<long double>::pi;
  a2 /= 2 * aol::NumberTrait<long double>::pi;

  for ( int x = cx - r;x <= cx + r;x++ ) {
    for ( int y = cy - r;y <= cy + r;y++, i++ ) {
      RealType vx = x - cx, vy = y - cy;
      RealType  n = sqrt ( vx * vx + vy * vy );
      RealType a = 0;
      if ( n > r )
        continue;
      if ( n > 1.0e-6 ) {
        vx /= n;
        vy /= n;
        RealType angle = ( ( vy >= 0 ) ? acosf ( vx ) : aol::NumberTrait<long double>::pi + acosf ( -vx ) ) / ( 2 * aol::NumberTrait<long double>::pi ); //angle of v with y-axis, in [0..1]
        if ( angle < 0. )
          angle += 1.;
        RealType b = angle;



        if ( b > a2 )
          b -= 1.;
        if ( b < a1 || b > a2 ) {
          continue;
        } else {
          a = ( b - a1 ) / ( a2 - a1 );

        }
      } else {
        continue;
      }

      out.set ( x, y, a );
    }
  }
}

/**
 * \author Berkels, Droske
 */
template <typename RealType>
class FullColorWheel {
  const int _numX, _numY;
  qc::ScalarArray<RealType, qc::QC_2D> _scaledInputArray;
  const qc::ScalarArray<RealType, qc::QC_2D> *_weight;
  const RealType _lowerAngle, _upperAngle;
  bool _fadeWheel;
  bool _fadeToBlack;
public:
  FullColorWheel ( const qc::ScalarArray<RealType, qc::QC_2D> &AngleArray, const RealType LowerAngle, const RealType UpperAngle )
    : _numX ( AngleArray.getNumX() ),
      _numY ( AngleArray.getNumY() ),
      _scaledInputArray ( AngleArray ), _weight ( NULL ),
      _lowerAngle ( LowerAngle ),
      _upperAngle ( UpperAngle ),
      _fadeWheel ( false ),
      _fadeToBlack ( true ) {

    //! The values in [_lowerAngle,_upperAngle] are scaled to [0,1]. AngleArray should not contain values outside [_lowerAngle,_upperAngle].
    _scaledInputArray.addToAll( -_lowerAngle );
    _scaledInputArray /= (_upperAngle-_lowerAngle);
  }

  void setWeight ( const qc::ScalarArray<RealType, qc::QC_2D> &weight ) {
    _weight = &weight;
  }
  
  void setFadeWheel ( const bool FadeWheel ) {
    _fadeWheel = FadeWheel;
  }
  
  void setFadeToBlack ( const bool FadeToBlack ) {
    _fadeToBlack = FadeToBlack;
  }

  void saveColoredWheel ( const std::string filename, const bool ShowColorWheel = true ) {
    qc::MultiArray<unsigned char, 2, 3> bufArray ( _numX, _numY );

    convertScalarToColor ( _scaledInputArray, bufArray );

    const int minNumXY = aol::Min( _numX, _numY );
    if ( ShowColorWheel )
      colwheel ( bufArray, 80 * minNumXY / 513, _numY - 80 * minNumXY / 513, 60 * minNumXY / 513, _lowerAngle, _upperAngle, _numX, _numY );
    if ( aol::fileNameEndsWith ( filename.c_str(), ".png" ) )
      bufArray.savePNG ( filename.c_str() );
    else {
      FILE* fp = fopen ( filename.c_str(), "wb" );
      bufArray.save ( fp );
      fclose ( fp );
    }
  }

private:
  void colwheel ( qc::MultiArray<unsigned char, 2, 3> &BufArray, int cx, int cy, int r, RealType a1, RealType a2, int numX, int numY );

  void convertScalarToColor ( const qc::ScalarArray<RealType, qc::QC_2D> &Scalar, qc::MultiArray<unsigned char, 2, 3> &BufArray ) {

    const int numX = Scalar.getNumX();
    const int numY = Scalar.getNumY();

    RealType rgb[3] = {0, 0, 0}
                      , hsv[3];

    for ( int y = 0; y < numY; y++ ) {
      for ( int x = 0; x < numX; x++ ) {
        RealType b = Scalar.get ( x, y );

        while ( b > 1. ) b -= 1.;
        while ( b < 0. ) b += 1.;

        const RealType w = _weight ? aol::Clamp( _weight->get ( x, y ), aol::ZOTrait<RealType>::zero, aol::ZOTrait<RealType>::one ) : 1.;

        hsv[0] = b;
        hsv[1] = w;
        hsv[2] = _fadeToBlack ? w : 1;

        aol::RGBColorMap<RealType>::hsv2rgb ( hsv, rgb );

        BufArray[0][y*numX+x] = static_cast< unsigned char> ( rgb[0] * 255 );
        BufArray[1][y*numX+x] = static_cast< unsigned char> ( rgb[1] * 255 );
        BufArray[2][y*numX+x] = static_cast< unsigned char> ( rgb[2] * 255 );
      }
    }
  }

};

template <typename RealType>
void FullColorWheel<RealType>::colwheel ( qc::MultiArray<unsigned char, 2, 3> &BufArray, int cx, int cy, int r, RealType a1, RealType a2, int numX, int ) {

  RealType rgb[3] = {0, 0, 0}, hsv[3];
  int i = 0;

  a1 /= 2 * aol::NumberTrait<long double>::pi;
  a2 /= 2 * aol::NumberTrait<long double>::pi;

  for ( int x = cx - r;x <= cx + r;x++ ) {
    for ( int y = cy - r;y <= cy + r;y++, i++ ) {
      RealType vx = x - cx, vy = y - cy;
      RealType  n = sqrt ( vx * vx + vy * vy );
      if ( n > r )
        continue;
      if ( n > 1.0e-6 ) {
        vx /= n;
        vy /= n;
        RealType angle = - ( ( vy >= 0 ) ? acos ( vx ) : aol::NumberTrait<long double>::pi + acos ( -vx ) ) / ( 2 * aol::NumberTrait<long double>::pi ); //angle of v with y-axis, in [0..1]

        if ( angle < 0. )
          angle += 1.;
        RealType b = angle;

        if ( b > a2 )
          b -= 1.;

        if ( b < a1 || b > a2 ) {
          rgb[0] = rgb[1] = rgb[2] = 0;
          continue;
        } else {
          RealType a = ( b - a1 ) / ( a2 - a1 );

          hsv[0] = a;     //0=red,0.25=vio,0.5=green,0.75=blue,1=red again
          hsv[1] = _fadeWheel ? n/r : 1;     //0=white, 1=colored
          hsv[2] = ( _fadeWheel && _fadeToBlack ) ? n/r : 1;     //luminance 0=black,1=colored
          aol::RGBColorMap<RealType>::hsv2rgb ( hsv, rgb );
        }
      } else {
        rgb[0] = rgb[1] = rgb[2] = _fadeToBlack ? 0 : 1;
        if ( !_fadeWheel )
          continue;
      }

      int pos = ( y * numX + x );
      if( pos >= BufArray[0].size() ){
        pos = BufArray[0].size() - 1;
        cerr << "Warning: Out of bounds access in FullColorWheel::colwheel prevented\n";
      }
      BufArray[0][pos]   = static_cast< unsigned char> ( rgb[0] * 255 );
      BufArray[1][pos]   = static_cast< unsigned char> ( rgb[1] * 255 );
      BufArray[2][pos]   = static_cast< unsigned char> ( rgb[2] * 255 );
    }
  }
}

/**
 * Writes a 2D vector field color coded based on the angle the vectors have with respect
 * to the x-axis.
 *
 * A negative factor value will set the factor to max |v(x)|
 *
 * Does a mirroing at the x-axis to reflect the mirroring the ScalarArray2D load and save
 * functions do. Needs to be consistent with qc::WriteVectorFieldAsGnuplotFile
 *
 * \author Berkels, Droske
 */
template <typename RealType>
void writeColorField ( const qc::Array<RealType> &vx, const qc::Array<RealType> &vy, const string &filename, RealType factor = -1, const bool ShowColorWheel = true,
                       const bool FadeWheel = false,
                       const bool FadeToBlack = true ) {

  qc::ScalarArray<RealType, qc::QC_2D> cw    ( vx, aol::STRUCT_COPY );
  qc::ScalarArray<RealType, qc::QC_2D> weight ( vx, aol::STRUCT_COPY );

  RealType max = 0.;

  aol::Vec2<int> size ( vx.getNumX(), vy.getNumY() );

  for ( int x = 0; x < size[0]; x++ ) {
    for ( int y = 0; y < size[1]; y++ ) {
      aol::Vec2<RealType> v ( vx.get ( x, y ), vy.get ( x, y ) );
      max = aol::Max ( max, v.norm() );
    }
  }

#ifdef VERBOSE
  cerr << "factor = " << factor << " max = " << max << endl;
#endif
  if ( factor < 0 ) factor = max;
#ifdef VERBOSE
  cerr << "factor = " << factor << " max = " << max << endl;
#endif


  for ( int x = 0; x < size[0]; x++ ) {
    for ( int y = 0; y < size[1]; y++ ) {
      aol::Vec2<RealType> v ( vx.get ( x, y ), vy.get ( x, y ) );

      RealType w = v.norm();
      if ( w > factor ) {
        w = factor;
      }

      w = w / factor;

      weight.set ( x, y, w );

      RealType alpha = 0;

      if ( v.normSqr() > 1e-8 ) {
        RealType n = v.norm();
        v[0] /= n;
        v[1] /= n;

        // -acos is in [-pi,0], shift this to [pi,2pi].
        alpha = - acos ( v[0] ) + 2 * aol::NumberTrait<long double>::pi;
        if ( v[1] < 0. ) {
          alpha = 2 * aol::NumberTrait<long double>::pi - alpha;
        }
      }
      cw.set ( x, y, alpha );
    }
  }

  qc::FullColorWheel<RealType> wheel ( cw, 0, 2*aol::NumberTrait<long double>::pi );
  wheel.setWeight ( weight );
  wheel.setFadeWheel ( FadeWheel );
  wheel.setFadeToBlack ( FadeToBlack );
  wheel.saveColoredWheel ( filename.c_str(), ShowColorWheel );
}

} // end namespace qc
#endif
