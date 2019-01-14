#ifndef __EPSWRITER_H
#define __EPSWRITER_H

#include<aol.h>
#include<vec.h>
#include<matrix.h>
#include <rgbColorMap.h>

namespace aol {

/**
 *  Class for writing color eps files.
 *  The square [0,1] is discretized as [0, 1e6]
 *  The eps dialect used here is due to reverse engineering of xfig exports.
 *  \todo redesign class avoiding Vector pointers, maybe move parts of code elsewhere
 *  \author schwen
 **/
class EpsWriter {
public:
  typedef double RealType;

  EpsWriter ( const char* filename, colorTrans colorTrans = UNDEFINED, unsigned int numColors = 256 );

  ~EpsWriter() {
    close();
  }

  void writeLine ( int x1, int y1, int x2, int y2, int thickness, int color );

  void writeFilledTrapezoid ( int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4, int color );


  static void setTopographical256 ( aol::Vector<float>* red, aol::Vector<float>* green, aol::Vector<float>* blue );

  static void setBlueToRed256 ( aol::Vector<float>* red, aol::Vector<float>* green, aol::Vector<float>* blue );

  static void setColortrans256 ( aol::Vector<float>* red, aol::Vector<float>* green, aol::Vector<float>* blue );

  static void setWhiteToBlack256 ( aol::Vector<float>* red, aol::Vector<float>* green, aol::Vector<float>* blue );

  static void setRotationMatrix ( aol::FullMatrix<double> &rot, double xang_deg, double yang_deg, double zang_deg );

  static void writeScalePPM ( colorTrans colTrans, const char* filename );

 protected:
  colorTrans _colTrans;
  unsigned int _numColors;

 private:
  FILE * output;

  EpsWriter();

  void writePreamble1(unsigned int sizeX = 1025, unsigned int sizeY = 1025);
  void writePreamble2();
  void writeColors ( aol::Vector<float>* red, aol::Vector<float>* green, aol::Vector<float>* blue );
  void writeColors ( colorTrans colorTrans, unsigned int numColors  );
  void writePostamble();

  void close(); // may want to allow explicite closing, but then need to make sure closed only once!

};

}

#endif
