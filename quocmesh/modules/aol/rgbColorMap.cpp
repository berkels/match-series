#include <rgbColorMap.h>

namespace aol {

template< typename RealType >
RGBColorMap<RealType>::RGBColorMap ( const RealType Min, const RealType Max, const colorTrans ColorTrans /* = BLACK_TO_WHITE */ ) : _min ( Min ), _max ( Max ), _colorValues ( 3, 0 ) {
  switch ( ColorTrans ) {
  case BLACK_TO_WHITE: {
    _colorValues.reallocate ( 3, 2 );
    _colorValues[0][0] = 0.0;  _colorValues[1][0] = 0.0;  _colorValues[2][0] = 0.0; // black
    _colorValues[0][1] = 1.0;  _colorValues[1][1] = 1.0;  _colorValues[2][1] = 1.0; // white
    break;
  }

  case WHITE_TO_BLACK: {
    _colorValues.reallocate ( 3, 2 );
    _colorValues[0][0] = 1.0;  _colorValues[1][0] = 1.0;  _colorValues[2][0] = 1.0; // white
    _colorValues[0][1] = 0.0;  _colorValues[1][1] = 0.0;  _colorValues[2][1] = 0.0; // black
    break;
  }

  case BLUE_TO_RED: {
    _colorValues.reallocate ( 3, 4 );
    _colorValues[0][0] = 0.0;  _colorValues[1][0] = 0.0;  _colorValues[2][0] = 1.0; // blue
    _colorValues[0][1] = 0.4;  _colorValues[1][1] = 0.0;  _colorValues[2][1] = 0.6;
    _colorValues[0][2] = 0.6;  _colorValues[1][2] = 0.0;  _colorValues[2][2] = 0.4;
    _colorValues[0][3] = 1.0;  _colorValues[1][3] = 0.0;  _colorValues[2][3] = 0.0; // red
    break;
  }

  case HSV_BLUE_TO_RED: {
    _colorValues.reallocate ( 3, 5 );
    _colorValues[0][0] = 0.0;  _colorValues[1][0] = 0.0;  _colorValues[2][0] = 1.0; // blue
    _colorValues[0][1] = 0.0;  _colorValues[1][1] = 1.0;  _colorValues[2][1] = 1.0; // cyan
    _colorValues[0][2] = 0.0;  _colorValues[1][2] = 1.0;  _colorValues[2][2] = 0.0; // green
    _colorValues[0][3] = 1.0;  _colorValues[1][3] = 1.0;  _colorValues[2][3] = 0.0; // yellow
    _colorValues[0][4] = 1.0;  _colorValues[1][4] = 0.0;  _colorValues[2][4] = 0.0; // red
    break;
  }

  case HSV_BLUE_TO_MAGENTA: {
    _colorValues.reallocate ( 3, 6 );
    _colorValues[0][0] = 0.0;  _colorValues[1][0] = 0.0;  _colorValues[2][0] = 1.0; // blue
    _colorValues[0][1] = 0.0;  _colorValues[1][1] = 1.0;  _colorValues[2][1] = 1.0; // cyan
    _colorValues[0][2] = 0.0;  _colorValues[1][2] = 1.0;  _colorValues[2][2] = 0.0; // green
    _colorValues[0][3] = 1.0;  _colorValues[1][3] = 1.0;  _colorValues[2][3] = 0.0; // yellow
    _colorValues[0][4] = 1.0;  _colorValues[1][4] = 0.0;  _colorValues[2][4] = 0.0; // red
    _colorValues[0][5] = 1.0;  _colorValues[1][5] = 0.0;  _colorValues[2][5] = 1.0; // magenta
    break;
  }

  case HSV_RED_TO_BLUE: {
    _colorValues.reallocate ( 3, 5 );
    _colorValues[0][0] = 1.0;  _colorValues[1][0] = 0.0;  _colorValues[2][0] = 0.0; // red
    _colorValues[0][1] = 1.0;  _colorValues[1][1] = 1.0;  _colorValues[2][1] = 0.0; // yellow
    _colorValues[0][2] = 0.0;  _colorValues[1][2] = 1.0;  _colorValues[2][2] = 0.0; // green
    _colorValues[0][3] = 0.0;  _colorValues[1][3] = 1.0;  _colorValues[2][3] = 1.0; // cyan
    _colorValues[0][4] = 0.0;  _colorValues[1][4] = 0.0;  _colorValues[2][4] = 1.0; // blue
    break;
  }

  case HSV_MAGENTA_TO_RED: {
    _colorValues.reallocate ( 3, 6 );
    _colorValues[0][0] = 1.0;  _colorValues[1][0] = 0.0;  _colorValues[2][0] = 1.0; // magenta
    _colorValues[0][1] = 0.0;  _colorValues[1][1] = 1.0;  _colorValues[2][1] = 0.0; // green
    _colorValues[0][2] = 1.0;  _colorValues[1][2] = 1.0;  _colorValues[2][2] = 0.0; // yellow
    _colorValues[0][3] = 0.0;  _colorValues[1][3] = 0.0;  _colorValues[2][3] = 1.0; // blue
    _colorValues[0][4] = 0.0;  _colorValues[1][4] = 1.0;  _colorValues[2][4] = 1.0; // cyan
    _colorValues[0][5] = 1.0;  _colorValues[1][5] = 0.0;  _colorValues[2][5] = 0.0; // red
    break;
  }

  case TOPOGRAPHIC: {
    _colorValues.reallocate ( 3, 12 );
    _colorValues[0][ 0] = 0.0;  _colorValues[1][ 0] = 0.0 ;  _colorValues[2][ 0] = 0.5 ;
    _colorValues[0][ 1] = 0.0;  _colorValues[1][ 1] = 0.0 ;  _colorValues[2][ 1] = 1.0 ;
    _colorValues[0][ 2] = 0.4;  _colorValues[1][ 2] = 0.4 ;  _colorValues[2][ 2] = 1.0 ;
    _colorValues[0][ 3] = 0.6;  _colorValues[1][ 3] = 0.6 ;  _colorValues[2][ 3] = 1.0 ;
    _colorValues[0][ 4] = 0.3;  _colorValues[1][ 4] = 1.0 ;  _colorValues[2][ 4] = 0.3 ;
    _colorValues[0][ 5] = 0.0;  _colorValues[1][ 5] = 1.0 ;  _colorValues[2][ 5] = 0.0 ;
    _colorValues[0][ 6] = 0.0;  _colorValues[1][ 6] = 0.5 ;  _colorValues[2][ 6] = 0.0 ;
    _colorValues[0][ 7] = 0.3;  _colorValues[1][ 7] = 0.2 ;  _colorValues[2][ 7] = 0.0 ;
    _colorValues[0][ 8] = 0.5;  _colorValues[1][ 8] = 0.3 ;  _colorValues[2][ 8] = 0.0 ;
    _colorValues[0][ 9] = 0.8;  _colorValues[1][ 9] = 0.5 ;  _colorValues[2][ 9] = 0.0 ;
    _colorValues[0][10] = 1.0;  _colorValues[1][10] = 0.7 ;  _colorValues[2][10] = 0.5 ;
    _colorValues[0][11] = 1.0;  _colorValues[1][11] = 0.85;  _colorValues[2][11] = 0.75;
    break;
  }

  case WHITE_TO_GREEN_TO_BLUE_TO_RED: {
    _colorValues.reallocate ( 3, 4 );
    _colorValues[0][0] = 1.0;  _colorValues[1][0] = 1.0;  _colorValues[2][0] = 1.0; // white
    _colorValues[0][1] = 0.0;  _colorValues[1][1] = 1.0;  _colorValues[2][1] = 0.0; // green
    _colorValues[0][2] = 0.0;  _colorValues[1][2] = 0.0;  _colorValues[2][2] = 1.0; // blue
    _colorValues[0][3] = 1.0;  _colorValues[1][3] = 0.0;  _colorValues[2][3] = 0.0; // red
    break;
  }

  case WHITE_TO_BLUE: {
    _colorValues.reallocate ( 3, 2 );
    _colorValues[0][0] = 1.0;  _colorValues[1][0] = 1.0;  _colorValues[2][0] = 1.0; // white
    _colorValues[0][1] = 0.0;  _colorValues[1][1] = 0.0;  _colorValues[2][1] = 1.0; // blue
    break;
  }

  default:
    throw aol::Exception ( "aol::RGBColorMap: illegal colorTrans", __FILE__, __LINE__ );
  }
}


template< typename RealType>
/*static*/ string RGBColorMap<RealType>::gnuplotColorString( const RealType* rgb ) {
  stringstream s;
  s << "\"#";
  s << hex << std::setfill ( '0' ) << std::setw ( 6 ) << ( static_cast<int> ( 255. * rgb[0] ) << 16 )
    + ( static_cast<int> ( 255. * rgb[1] ) << 8  )
    + ( static_cast<int> ( 255. * rgb[2] )       );
  s << "\"";
  return s.str();
}

//! \brief Generate RGB hex notation from color values and index.
//! \param[in] rgb RGB color values as given by getColorValues.
//! \param[in] i   The index of the color to convert.
template< typename RealType>
string RGBColorMap<RealType>::gnuplotColorString( const aol::MultiVector < RealType >& rgb, const int i ) {
  stringstream s;
  s << "\"#";
  s << hex << std::setfill ( '0' ) << std::setw ( 6 ) << ( static_cast<int> ( 255. * rgb[0][i] ) << 16 )
    + ( static_cast<int> ( 255. * rgb[1][i] ) << 8  )
    + ( static_cast<int> ( 255. * rgb[2][i] )       );
  s << "\"";
  return s.str();
}


template class RGBColorMap<float>;
template class RGBColorMap<double>;
template class RGBColorMap<long double>;

}
