/* *******************************************************************
 * Implementation of EpsWriter
 * ******************************************************************* */

#include "epswriter.h"

namespace aol {

EpsWriter::EpsWriter ( const char* filename, colorTrans colTrans, unsigned int numColors )
  : _colTrans( colTrans ), _numColors( numColors ) {
  // cerr << "Writing eps to file " << filename << endl;
  output = fopen ( filename, "w" );
  if ( output ) {
    writePreamble1();
    writeColors ( colTrans, numColors );
    writePreamble2();
  } else {
    cerr << "Error - could not open eps file " << filename << endl;
  }
}

void EpsWriter::close() {
  writePostamble();
  fclose ( output );
}

void EpsWriter::writeLine ( int x1, int y1, int x2, int y2, int thickness, int color ) {
  fprintf ( output, "%d slw n %d %d m %d %d l gs col%d s gr\n", thickness, x1, y1, x2, y2, color );
}

void EpsWriter::writeFilledTrapezoid ( int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4, int color ) {
  fprintf ( output, "0 slw n %d %d m %d %d l %d %d l %d %d l cp gs col%d 1.0 shd ef gr\n", x1, y1, x2, y2, x3, y3, x4, y4, color );
}

void EpsWriter::writePreamble1(unsigned int sizeX, unsigned int sizeY) {
  // this preamble is taken from xfig eps export
  fprintf ( output, "%%!PS-Adobe-2.0 EPSF-2.0\n" );
  //   fprintf( output, "%%%%Title: eps file \n");
  fprintf ( output, "%%%%BoundingBox: 0 0 %u %u\n", sizeX, sizeY );
  fprintf ( output, "%%%%Magnification: 1.0000\n" );
  fprintf ( output, "%%%%EndComments\n" );

  fprintf ( output, "/$F2psDict 200 dict def\n" );
  fprintf ( output, "$F2psDict begin\n" );
  fprintf ( output, "$F2psDict /mtrx matrix put\n" );
  fprintf ( output, "/col-1 {0 setgray} bind def\n" );
}

void EpsWriter::writePreamble2() {
  // other preliminary stuff
  fprintf ( output, "end\n" );
  fprintf ( output, "save\n" );
  fprintf ( output, "newpath 0 1025 moveto 0 0 lineto 1025 0 lineto 1025 1025 lineto closepath clip newpath\n" );
  fprintf ( output, " 0 1025 translate\n" );
  fprintf ( output, "1 -1 scale\n" );
  fprintf ( output, "\n" );
  fprintf ( output, "/cp {closepath} bind def\n" );
  fprintf ( output, "/ef {eofill} bind def\n" );
  fprintf ( output, "/gr {grestore} bind def\n" );
  fprintf ( output, "/gs {gsave} bind def\n" );
  fprintf ( output, "/sa {save} bind def\n" );
  fprintf ( output, "/rs {restore} bind def\n" );
  fprintf ( output, "/l {lineto} bind def\n" );
  fprintf ( output, "/m {moveto} bind def\n" );
  fprintf ( output, "/rm {rmoveto} bind def\n" );
  fprintf ( output, "/n {newpath} bind def\n" );
  fprintf ( output, "/s {stroke} bind def\n" );
  fprintf ( output, "/sh {show} bind def\n" );
  fprintf ( output, "/slc {setlinecap} bind def\n" );
  fprintf ( output, "/slj {setlinejoin} bind def\n" );
  fprintf ( output, "/slw {setlinewidth} bind def\n" );
  fprintf ( output, "/srgb {setrgbcolor} bind def\n" );
  fprintf ( output, "/rot {rotate} bind def\n" );
  fprintf ( output, "/sc {scale} bind def\n" );
  fprintf ( output, "/sd {setdash} bind def\n" );
  fprintf ( output, "/ff {findfont} bind def\n" );
  fprintf ( output, "/sf {setfont} bind def\n" );
  fprintf ( output, "/scf {scalefont} bind def\n" );
  fprintf ( output, "/sw {stringwidth} bind def\n" );
  fprintf ( output, "/tr {translate} bind def\n" );
  fprintf ( output, "/tnt {dup dup currentrgbcolor\n" );
  fprintf ( output, "  4 -2 roll dup 1 exch sub 3 -1 roll mul add\n" );
  fprintf ( output, "  4 -2 roll dup 1 exch sub 3 -1 roll mul add\n" );
  fprintf ( output, "  4 -2 roll dup 1 exch sub 3 -1 roll mul add srgb}\n" );
  fprintf ( output, "  bind def\n" );
  fprintf ( output, "/shd {dup dup currentrgbcolor 4 -2 roll mul 4 -2 roll mul\n" );
  fprintf ( output, "  4 -2 roll mul srgb} bind def\n" );
  fprintf ( output, "/$F2psBegin {$F2psDict begin /$F2psEnteredState save def} def\n" );
  fprintf ( output, "/$F2psEnd {$F2psEnteredState restore end} def\n" );
  fprintf ( output, "\n" );
  fprintf ( output, "$F2psBegin\n" );
  fprintf ( output, "10 setmiterlimit\n" );
  fprintf ( output, "0 slj 0 slc\n" );
  fprintf ( output, " 0.001024 0.001024 sc\n" );
}

void EpsWriter::writeColors ( aol::Vector<float>* red, aol::Vector<float>* green, aol::Vector<float>* blue ) {
  if ( red && green && blue ) {
    if ( ( red->size() == green->size() ) && ( green->size() == blue->size() ) ) {
      for ( int i = 0; i < red->size(); i++ ) {
        fprintf ( output, "/col%d {%f %f %f srgb} bind def\n", i, red->get ( i ), green->get ( i ), blue->get ( i ) );
      }
    } else {
      cerr << "Sizes of r, g, b color values incompatible" << endl;
    }
  } else { // use some standard color palette
    fprintf ( output, "/col0 {0.0 0.0 0.0 srgb} bind def\n" );
    fprintf ( output, "/col1 {0.0 0.0 1.0 srgb} bind def\n" );
    fprintf ( output, "/col2 {0.0 1.0 0.0 srgb} bind def\n" );
    fprintf ( output, "/col3 {0.0 1.0 1.0 srgb} bind def\n" );
    fprintf ( output, "/col4 {1.0 0.0 0.0 srgb} bind def\n" );
    fprintf ( output, "/col5 {1.0 0.0 1.0 srgb} bind def\n" );
    fprintf ( output, "/col6 {1.0 1.0 0.0 srgb} bind def\n" );
    fprintf ( output, "/col7 {1.0 1.0 1.0 srgb} bind def\n" );
  }
}

  void EpsWriter::writeColors ( colorTrans colTrans, unsigned int numColors  ) {
  if ( colTrans != UNDEFINED ) {
    aol::RGBColorMap<RealType> colorMap(0.0, 1.0, colTrans);
    for ( unsigned int i = 0; i < numColors; i++ ) {
      fprintf ( output, "/col%d {%f %f %f srgb} bind def\n", i, colorMap.scalarToColor(static_cast<RealType>(i) / static_cast<RealType>(numColors - 1))[0],
                                                    colorMap.scalarToColor(static_cast<RealType>(i) / static_cast<RealType>(numColors - 1))[1],
                                                    colorMap.scalarToColor(static_cast<RealType>(i) / static_cast<RealType>(numColors - 1))[2] );
    }
  } else { // use some standard color palette
    fprintf ( output, "/col0 {0.0 0.0 0.0 srgb} bind def\n" );
    fprintf ( output, "/col1 {0.0 0.0 1.0 srgb} bind def\n" );
    fprintf ( output, "/col2 {0.0 1.0 0.0 srgb} bind def\n" );
    fprintf ( output, "/col3 {0.0 1.0 1.0 srgb} bind def\n" );
    fprintf ( output, "/col4 {1.0 0.0 0.0 srgb} bind def\n" );
    fprintf ( output, "/col5 {1.0 0.0 1.0 srgb} bind def\n" );
    fprintf ( output, "/col6 {1.0 1.0 0.0 srgb} bind def\n" );
    fprintf ( output, "/col7 {1.0 1.0 1.0 srgb} bind def\n" );
  }
}

void EpsWriter::writePostamble() {
  fprintf ( output, "$F2psEnd\n" );
  fprintf ( output, "rs\n" );
  fprintf ( output, "showpage\n" );
}

// TODO: Merge the functions below with aol::RGBColorMap
void EpsWriter::setTopographical256 ( aol::Vector<float>* red, aol::Vector<float>* green, aol::Vector<float>* blue ) {
  aol::Vector<float> redv ( 12 ), greenv ( 12 ), bluev ( 12 );
  redv[ 0] = 0.0;  greenv[ 0] = 0.0;  bluev [ 0] = 0.5;
  redv[ 1] = 0.0;  greenv[ 1] = 0.0;  bluev [ 1] = 1.0;
  redv[ 2] = 0.4;  greenv[ 2] = 0.4;  bluev [ 2] = 1.0;
  redv[ 3] = 0.6;  greenv[ 3] = 0.6;  bluev [ 3] = 1.0;
  redv[ 4] = 0.3;  greenv[ 4] = 1.0;  bluev [ 4] = 0.3;
  redv[ 5] = 0.0;  greenv[ 5] = 1.0;  bluev [ 5] = 0.0;
  redv[ 6] = 0.0;  greenv[ 6] = 0.5;  bluev [ 6] = 0.0;
  redv[ 7] = 0.3;  greenv[ 7] = 0.2;  bluev [ 7] = 0.0;
  redv[ 8] = 0.5;  greenv[ 8] = 0.3;  bluev [ 8] = 0.0;
  redv[ 9] = 0.8;  greenv[ 9] = 0.5;  bluev [ 9] = 0.0;
  redv[10] = 1.0;  greenv[10] = 0.7;  bluev [10] = 0.5;
  redv[11] = 1.0;  greenv[11] = 0.85; bluev [11] = 0.75;

  red->resize ( 256 );
  green->resize ( 256 );
  blue->resize ( 256 );

  for ( int i = 0; i < 256; i++ ) {
    red->set (   i, redv.interpolateInRange   ( i / 255.0f, -1.0f, 1.0f ) );
    green->set ( i, greenv.interpolateInRange ( i / 255.0f, -1.0f, 1.0f ) );
    blue->set (  i, bluev.interpolateInRange  ( i / 255.0f, -1.0f, 1.0f ) );
  }
}

void EpsWriter::setBlueToRed256 ( aol::Vector<float>* red, aol::Vector<float>* green, aol::Vector<float>* blue ) {
  aol::Vector<float> redv ( 4 ), greenv ( 4 ), bluev ( 4 );

  redv[ 0] = 0.0; greenv[ 0] = 0.0; bluev [ 0] = 1.0;
  redv[ 1] = 0.4; greenv[ 1] = 0.0; bluev [ 1] = 0.6;
  redv[ 2] = 0.6; greenv[ 2] = 0.0; bluev [ 2] = 0.4;
  redv[ 3] = 1.0; greenv[ 3] = 0.0; bluev [ 3] = 0.0;

  red->resize ( 256 );
  green->resize ( 256 );
  blue->resize ( 256 );

  for ( int i = 0; i < 256; i++ ) {
    red->set (   i, redv.interpolateInRange   ( i / 255.0f, -1.0f, 1.0f ) );
    green->set ( i, greenv.interpolateInRange ( i / 255.0f, -1.0f, 1.0f ) );
    blue->set (  i, bluev.interpolateInRange  ( i / 255.0f, -1.0f, 1.0f ) );
  }
}


//equal to HSV_BLUE_TO_RED
void EpsWriter::setColortrans256 ( aol::Vector<float>* red, aol::Vector<float>* green, aol::Vector<float>* blue ) {
  aol::Vector<float> redv ( 5 ), greenv ( 5 ), bluev ( 5 );

  // hsv blue to red
  redv[ 0] = 0.0; greenv[ 0] = 0.0; bluev [ 0] = 1.0; // blue
  redv[ 1] = 0.0; greenv[ 1] = 1.0; bluev [ 1] = 1.0; // cyan
  redv[ 2] = 0.0; greenv[ 2] = 1.0; bluev [ 2] = 0.0; // green
  redv[ 3] = 1.0; greenv[ 3] = 1.0; bluev [ 3] = 0.0; // yellow
  redv[ 4] = 1.0; greenv[ 4] = 0.0; bluev [ 4] = 0.0; // red

  red->resize ( 256 );
  green->resize ( 256 );
  blue->resize ( 256 );

  for ( int i = 0; i < 256; i++ ) {
    red->set (   i, redv.interpolateInRange   ( i / 255.0f, -1.0f, 1.0f ) );
    green->set ( i, greenv.interpolateInRange ( i / 255.0f, -1.0f, 1.0f ) );
    blue->set (  i, bluev.interpolateInRange  ( i / 255.0f, -1.0f, 1.0f ) );
  }
}

void EpsWriter::setWhiteToBlack256 ( aol::Vector<float>* red, aol::Vector<float>* green, aol::Vector<float>* blue ) {
  aol::Vector<float> redv ( 2 ), greenv ( 2 ), bluev ( 2 );

  // greyscale
  redv[ 0] = 1.0; greenv[ 0] = 1.0; bluev [ 0] = 1.0; // white
  redv[ 1] = 0.0; greenv[ 1] = 0.0; bluev [ 1] = 0.0; // black

  red->resize ( 256 );
  green->resize ( 256 );
  blue->resize ( 256 );

  for ( int i = 0; i < 256; i++ ) {
    red->set (   i, redv.interpolateInRange   ( i / 255.0f, 0.0f, 1.0f ) );
    green->set ( i, greenv.interpolateInRange ( i / 255.0f, 0.0f, 1.0f ) );
    blue->set (  i, bluev.interpolateInRange  ( i / 255.0f, 0.0f, 1.0f ) );
  }
}

void EpsWriter::writeScalePPM ( colorTrans colTrans, const char* filename ) {
  FILE * outdat;
  const unsigned int numColors = 256;
  outdat = fopen ( filename, "w" );
  aol::RGBColorMap<RealType> colorMap(0.0, 1.0, colTrans);
  fprintf ( outdat, "P3\n%d %d\n255\n", 1, numColors );
  for ( int i = numColors ; i >= 0 ; i-- ) { // want zero at bottom of image
    fprintf ( outdat, "%d %d %d\n", static_cast<int> ( colorMap.scalarToColor(static_cast<RealType>(i) / static_cast<RealType>(numColors))[0] ),
                              static_cast<int> ( colorMap.scalarToColor(static_cast<RealType>(i) / static_cast<RealType>(numColors))[1] ),
                                    static_cast<int> ( colorMap.scalarToColor(static_cast<RealType>(i) / static_cast<RealType>(numColors))[2] ) );
  }
  fclose ( outdat );
}


void EpsWriter::setRotationMatrix ( aol::FullMatrix<double> &rot, double xang_deg, double yang_deg, double zang_deg ) {
  // rot is set to the rotation matrix specified by the three angles (in degrees)
  aol::FullMatrix<double> xrot ( 3, 3 ), yrot ( 3, 3 ), zrot ( 3, 3 );

  const double
    xang = xang_deg * 0.0174329252, // pi / 180
    yang = yang_deg * 0.0174329252,
    zang = zang_deg * 0.0174329252;

  xrot.set ( 0, 0,        1.0 );   xrot.set ( 0, 1,        0.0 );   xrot.set ( 0, 2,        0.0 );
  xrot.set ( 1, 0,        0.0 );   xrot.set ( 1, 1,  cos ( xang ) );   xrot.set ( 1, 2, -sin ( xang ) );
  xrot.set ( 2, 0,        0.0 );   xrot.set ( 2, 1,  sin ( xang ) );   xrot.set ( 2, 2,  cos ( xang ) );

  yrot.set ( 0, 0,  cos ( yang ) );   yrot.set ( 0, 1,        0.0 );   yrot.set ( 0, 2, -sin ( yang ) );
  yrot.set ( 1, 0,        0.0 );   yrot.set ( 1, 1,        1.0 );   yrot.set ( 1, 2,        0.0 );
  yrot.set ( 2, 0,  sin ( yang ) );   yrot.set ( 2, 1,        0.0 );   yrot.set ( 2, 2,  cos ( yang ) );

  zrot.set ( 0, 0,  cos ( zang ) );   zrot.set ( 0, 1, -sin ( zang ) );   zrot.set ( 0, 2,        0.0 );
  zrot.set ( 1, 0,  sin ( zang ) );   zrot.set ( 1, 1,  cos ( zang ) );   zrot.set ( 1, 2,        0.0 );
  zrot.set ( 2, 0,        0.0 );   zrot.set ( 2, 1,        0.0 );   zrot.set ( 2, 2,        1.0 );

  cerr << "Please ignore the following warnings. The matrices are only 3x3 ... " << endl;
  rot  = xrot;
  rot *= yrot;
  rot *= zrot;
  cerr << "Thank you for your patience :-)" << endl;
}

}
