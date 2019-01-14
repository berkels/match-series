#ifndef __RGBCOLORMAP_H
#define __RGBCOLORMAP_H

#include <aol.h>
#include <multiVector.h>
#include <smallVec.h>

namespace aol {

enum colorTrans { BLACK_TO_WHITE,
                  WHITE_TO_BLACK,
                  BLUE_TO_RED,
                  HSV_BLUE_TO_RED,
                  HSV_BLUE_TO_MAGENTA,
                  HSV_RED_TO_BLUE,
                  TOPOGRAPHIC,
                  WHITE_TO_GREEN_TO_BLUE_TO_RED,
                  WHITE_TO_BLUE,
                  HSV_MAGENTA_TO_RED,
                  UNDEFINED };

/** Routines associated with color mapping
 *
 * <ul><li>Mapping from scalar values in range [Min, Max] to three-channel color with predefined color maps</li>
 *     <li>rgb to hsv</li>
 *     <li>hsv to rgb</li>
 *     <li>rgb color string</li></ul>
 *
 *  \author Schwen
 */
template< typename RealType >
class RGBColorMap {
protected:
  RealType _min, _max;
  aol::MultiVector<RealType> _colorValues;

public:
  RGBColorMap ( ) : _min ( aol::NumberTrait<RealType>::NaN ), _max ( aol::NumberTrait<RealType>::NaN ), _colorValues() {
  }

  // copy constructor and assignment operator should do the correct thing

  RGBColorMap ( const RealType Min, const RealType Max, const colorTrans ColorTrans = BLACK_TO_WHITE );

  //! map scalar value to color triple, returning color value (slow due to temporary object)
  inline aol::Vec3<RealType> scalarToColor ( const RealType val ) const {
    aol::Vec3<RealType> col;
    scalarToColor ( val, col );
    return ( col );
  }

  //! map scalar value to color triple
  inline void scalarToColor ( const RealType val, aol::Vec3<RealType> &color ) const {
    for ( short i = 0; i < 3; ++i )
      color[i] = _colorValues[i].interpolateInRange ( val, _min, _max );
  }

  //! Simple routine which converts from a HSV triple to an RGB triple. Both have to be in the range 0..1,0..1,0..1
  inline static void hsv2rgb ( const RealType* hsv, RealType* rgb )  {
    const int hueCase = static_cast< int > ( hsv[0] * 6 );
    const RealType frac = 6 * hsv[0] - hueCase;
    const RealType lx  = hsv[2] * ( 1.0 - hsv[1] );
    const RealType ly  = hsv[2] * ( 1.0 - hsv[1] * frac );
    const RealType lz  = hsv[2] * ( 1.0 - hsv[1] * ( 1.0 - frac ) );

    switch ( hueCase ) {
    case 0:
    case 6:
      rgb[0] = hsv[2];
      rgb[1] = lz;
      rgb[2] = lx;
      break;  /* 0<hue<1/6   */
    case 1:
      rgb[0] = ly;
      rgb[1] = hsv[2];
      rgb[2] = lx;
      break;  /* 1/6<hue<2/6 */
    case 2:
      rgb[0] = lx;
      rgb[1] = hsv[2];
      rgb[2] = lz;
      break;  /* 2/6<hue<3/6 */
    case 3:
      rgb[0] = lx;
      rgb[1] = ly;
      rgb[2] = hsv[2];
      break;  /* 3/6<hue/4/6 */
    case 4:
      rgb[0] = lz;
      rgb[1] = lx;
      rgb[2] = hsv[2];
      break;  /* 4/6<hue<5/6 */
    case 5:
      rgb[0] = hsv[2];
      rgb[1] = lx;
      rgb[2] = ly;
      break;  /* 5/6<hue<1   */
    default:
      throw aol::UnimplementedCodeException ( "hsv2rgb: unhandled switch case!", __FILE__, __LINE__ );
    }
  }

  //! Simple routine which converts from an RGB triple to an HSV triple. Both have to be in the range 0..1,0..1,0..1
  inline static void rgb2hsv ( const RealType* rgb, RealType* hsv )  {
    const RealType M = aol::Max<RealType> ( rgb[0], rgb[1], rgb[2] );
    const RealType m = aol::Min<RealType> ( rgb[0], rgb[1], rgb[2] );
    const RealType d = M - m;

    hsv[2] = M;    //value == max(r,g,b)
    hsv[1] = ( M > 0.00001 ) ? d / M : 0; //saturation

    if ( hsv[1] == 0 )
      hsv[0] = 0;  //achromatic case, hue is 0 by convention
    else     //chromatic case
    {
      if ( rgb[0] == M )
        hsv[0] = ( rgb[1] - rgb[2] ) / d;
      else if ( rgb[1] == M )
        hsv[0] = 2 + ( rgb[2] - rgb[0] ) / d;
      else
        hsv[0] = 4 + ( rgb[0] - rgb[1] ) / d;

      hsv[0] /= 6;
      if ( hsv[0] < 0 )
        hsv[0] += 1;
    }
  }

  /*
   *\author Tatano
   * Note: The transformation is taken from Sridhar, Colour Image Processing, Oxford University Press, 2011 (pag. 352)
   */
  //! Simple routine which converts from an RGB triple to an HSI triple. Both have to be in the range 0..1,0..1,0..1
  inline static void rgb2hsi ( const RealType* rgb, RealType* hsi )  {
    const RealType m = aol::Min<RealType> ( rgb[0], rgb[1], rgb[2] );
    const RealType sum = rgb[0] + rgb[1] + rgb[2];
    const RealType theta = acos(0.5 * (rgb[0]-rgb[1]+rgb[0]-rgb[2])/(sqrt(aol::Sqr(rgb[0]-rgb[1])+(rgb[0]-rgb[2])*(rgb[1]-rgb[2]))+std::numeric_limits<RealType>::epsilon()) );
    
    if (rgb[2] > rgb[1]) {
      hsi[0] = (360. - theta)/360.;
    }else
      hsi[0] = theta/360.;
    if (sum == 0.) {
      hsi[1] = 0.;
    }else
      hsi[1] = 1. - 3. * m/sum;
    hsi[2] = 1./3.0 * sum;
  }

  /*
   *\author Tatano
   * Note: The transformation is taken from Sridhar, Colour Image Processing, Oxford University Press, 2011 (pag. 353)
   */
  //! Simple routine which converts from an HSI triple to an RGB triple. Both have to be in the range 0..1,0..1,0..1
  inline static void hsi2rgb ( RealType* hsi, RealType* rgb )  {
    hsi[0] *= 360;
    if ( hsi[0] >= 0. && hsi[0] < 120. ) {
      rgb[0] = hsi[2] * ( 1. + hsi[1] * cos(hsi[0])/cos(60. - hsi[0]) );
      rgb[2] = hsi[2] * ( 1. - hsi[1] );
      rgb[1] = 3. * hsi[2] - ( rgb[0] + rgb[2] );
    }
    else if ( hsi[0] >= 120. && hsi[0] < 240. ) {
      hsi[0] -= 120.;
      rgb[1] = hsi[2] * ( 1. + hsi[1] * cos(hsi[0])/cos(60. - hsi[0]) );
      rgb[0] = hsi[2] * ( 1. - hsi[1] );
      rgb[2] = 3. * hsi[2] - ( rgb[0] + rgb[1] );
    }
    else if ( hsi[0] >= 240. && hsi[0] < 360. ) {
      hsi[0] -= 240.;
      rgb[2] = hsi[2] * ( 1. + hsi[1] * cos(hsi[0])/cos(60. - hsi[0]) );
      rgb[1] = hsi[2] * ( 1. - hsi[1] );
      rgb[0] = 3. * hsi[2] - ( rgb[2] + rgb[1] );
    }
  }
  
  /*
   *\author Tatano
   * Note: The transformation is taken from Sridhar, Colour Image Processing, Oxford University Press, 2011 (pag. 346)
   */
  //! Simple routine which converts from an RGB triple to an gray value. Both have to be in the range 0..1,0..1,0..1
  inline static void rgb2gray ( const RealType* rgb, RealType &gray )  {
    gray = 0.2125 * rgb[0] + 0.7154 * rgb[1] + 0.072 * rgb[2];
  }

  /*
   *\author Tatano
   * Note: The transformation is taken from Sridhar, Colour Image Processing, Oxford University Press, 2011 (pag. 359)
   */
  //! Simple routine which converts from an RGB triple to an YCbCr triple. Both have to be in the range 0..1,0..1,0..1
  inline static void rgb2YCbCr ( const RealType* rgb, RealType* YCbCr )  {
    YCbCr[0] =   0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2];
    YCbCr[1] = - 0.169 * rgb[0] - 0.331 * rgb[1] + 0.500 * rgb[2];
    YCbCr[2] =   0.500 * rgb[0] - 0.419 * rgb[1] - 0.081 * rgb[2];
  }

  /*
   *\author Tatano
   * Note: The transformation is taken from Sridhar, Colour Image Processing, Oxford University Press, 2011 (pag. 360)
   */
  //! Simple routine which converts from an YCbCr triple to an RGB triple. Both have to be in the range 0..1,0..1,0..1
  inline static void YCbCr2rgb ( const RealType* YCbCr, RealType* rgb )  {
    rgb[0] =  1.0 * YCbCr[0] + 0.000 * YCbCr[1] + 1.403 * YCbCr[2];
    rgb[1] =  1.0 * YCbCr[0] - 0.344 * YCbCr[1] - 0.714 * YCbCr[2];
    rgb[2] =  1.0 * YCbCr[0] + 1.773 * YCbCr[1] + 0.000 * YCbCr[2];
  }
  
  /*
   *\author Mevenkamp
   * Note: transformation taken from: http://www.easyrgb.com/index.php?X=MATH&H=02#text2
   */
  //! Simple routine which converts from an RGB triple (in range [0,255]) to an XYZ triple.
  inline static void rgb2xyz ( const RealType* rgb, RealType* xyz ) {
    RealType var_rgb[3] = {static_cast<RealType>(rgb[0] / 255.0), static_cast<RealType>(rgb[1] / 255.0), static_cast<RealType>(rgb[2] / 255.0)};
    
    for ( int c=0; c<3 ; ++c ) {
      if ( var_rgb[c] > 0.04045 ) var_rgb[c] = pow ( ( ( rgb[c] + 0.055f ) / 1.055f ), static_cast<RealType>(2.4) );
      else var_rgb[c] = var_rgb[c] / 12.92;
      var_rgb[c] *= 100.0;
    }
    
    //Observer = 2°, Illuminant = D65
    xyz[0] = var_rgb[0] * 0.4124 + var_rgb[1] * 0.3576 + var_rgb[2] * 0.1805;
    xyz[1] = var_rgb[0] * 0.2126 + var_rgb[1] * 0.7152 + var_rgb[2] * 0.0722;
    xyz[2] = var_rgb[0] * 0.0193 + var_rgb[1] * 0.1192 + var_rgb[2] * 0.9505;
  }
  
  /*
   *\author Mevenkamp
   * Note: transformation is a concatenation of rgb2xyz above and the transformation from XYZ to CIE L*ab taken from: http://www.easyrgb.com/index.php?X=MATH&H=07#text7
   */
  //! Simple routine which converts from an RGB triple (in range [0,255]) to a CIE L*ab triple.
  inline static void rgb2CIELab ( const RealType* rgb, RealType* CIELab ) {
    RealType xyz[3];
    rgb2xyz ( rgb, xyz );
    
    // Observer = 2°, Illuminant = D65
    xyz[0] /= 95.047;
    xyz[1] /= 100.00;
    xyz[2] /= 108.883;
    
    for ( int c=0; c<3 ; ++c ) {
      if ( xyz[c] > 0.008856 ) xyz[c] = pow ( xyz[c], static_cast<RealType>( 1.0 / 3.0 ) );
      else xyz[c] = 7.787 * xyz[c] + 16.0 / 116.0;
    }
    
    CIELab[0] = ( 116 * xyz[1] ) - 16;
    CIELab[1] = 500 * ( xyz[0] - xyz[1] );
    CIELab[2] = 200 * ( xyz[1] - xyz[2] );
  }
  
  //! \brief Returns the color values.
  inline const aol::MultiVector<RealType>& getColorValues () const {
    return _colorValues;
  }

  inline const RealType getMin () const {
    return _min;
  }

  inline const RealType getMax () const {
    return _max;
  }

  //! Generate rgb hex notation used in html, gnuplot, ...
  static string gnuplotColorString( const RealType* rgb );
  static string gnuplotColorString( const aol::MultiVector < RealType >& rgb, int i );

};

}

#endif
