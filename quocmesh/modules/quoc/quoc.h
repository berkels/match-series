#ifndef __QUOC_H
#define __QUOC_H

#include <aol.h>
#include <smallVec.h>

namespace qc {

struct ArrayHeader {
  char magic[4];
  int numX, numY, numZ;
  double max;
};

/** \brief Magic numbers for pgm and pgm style files used for saving
 *  The type conventions for 1D, 2D and 3D PGM are:
 *  <table>
 *  <tr><td>type 2</td><td>ASCII unsigned char data PGM</td></tr>
 *  <tr><td>type 5</td><td>RAW unsigned char data PGM</td></tr>
 *  <tr><td>type 7</td><td>ASCII float data PGM</td></tr>
 *  <tr><td>type 8</td><td>RAW float data PGM</td></tr>
 *  <tr><td>type 9</td><td>RAW double data PGM </td></tr>
 *  <tr><td>type 10</td><td>RAW unsigned short data PGM </td></tr>
 *  <tr><td>type 11</td><td>RAW signed short data PGM </td></tr>
 *  <tr><td>type 33</td><td>RAW signed int data PGM </td></tr>
 *  <tr><td>type 34</td><td>RAW unsigned int data PGM </td></tr>
 *  <tr><td>type 211</td><td>PNG </td></tr>
 *  </table>
 *  Note that there is an enum available that should be used instead of those integers.
 *  Types 2, 5 and 211 can be directly viewed using e.g. xv or gimp for 2D data (i.e. on Array2D).
 *  Type 211 cannot be used in 3D and does not allow for comments.
 *  The floating point types (7, 8, 9) don't use the overflow handling. Instead the values are
 *  saved directly, just casted to the selected floating point type.
 */
enum SaveType {
  PGM_UNSIGNED_CHAR_ASCII    =  2,                     //!< can be opened in image processing programs (xv, gimp, ...)
  PGM_UNSIGNED_CHAR_BINARY   = aol::FF_UNSIGNED_CHAR,  //!< can be opened in image processing programs (xv, gimp, ...)
  PGM_FLOAT_ASCII            =  7,
  PGM_FLOAT_BINARY           = aol::FF_FLOAT,
  PGM_DOUBLE_BINARY          = aol::FF_DOUBLE,
  PGM_UNSIGNED_SHORT_BINARY  = aol::FF_UNSIGNED_SHORT,
  PGM_SHORT_BINARY           = aol::FF_SIGNED_SHORT,
  PGM_UNSIGNED_INT_BINARY    = aol::FF_UNSIGNED_INT,
  PGM_SIGNED_INT_BINARY      = aol::FF_SIGNED_INT,
  PNG_2D                     = 211,
  PGM_UNSIGNED_SHORT_BINARY_BIGENDIAN = 212,
  NETCDF                     = 213
};

/**
 * Given a qc::SaveType, return "sizeof" of the underlying data type.
 *
 * \author Berkels
 */
int getSizeOfSaveType ( const qc::SaveType Type );

/** Type of traversal
 */
enum tMode {
  PREFIX,   /**< callback is called before traversal of children */
  POSTFIX,  /**< callback is called after traversal of children */
  INFIX     /**< callback is called only on leaf elements */
};

/** Possible values for dimension of world
 */
enum Dimension {
  QC_1D = 1,    /**< dimension 1 */
  QC_2D = 2,    /**< dimension 2 */
  QC_3D = 3     /**< dimension 3 */
};

/** possible components that can be used
 *  to identify a component in position vectors.
 *  The hard coded values ensure that one can use
 *  the enum values as index (do not modify!).
 */
enum Comp {
  QC_X = 0,
  QC_Y = 1,
  QC_Z = 2
};

/**
 * Given a qc::SaveType, return the suffix usually used in names of files with this type,
 * e.g. returns .pgm for PGM_UNSIGNED_CHAR_ASCII. The optional argument Dim has to be
 * specified to get dimension dependent suffixes like .q2bz.
 *
 * \author Berkels
 */
const char *getDefaulSuffixOfSaveType ( const qc::SaveType Type, const qc::Dimension Dim = qc::QC_3D );

typedef aol::Vec3<short> CoordType;


/** Trait for default binary and ASCII save types for 2d/3d arrays and different data types
 *  To be specialized only for those types where a corresponding SaveType exists.
 * \author Schwen
 */
template< typename DataType >
class SaveTypeTrait {
public:
  static const qc::SaveType AsciiSaveType;
  static const qc::SaveType BinarySaveType;
};

bool getDefaultArrayCompressionState( const qc::Dimension Dim );
void setDefaultArrayCompressionState( const qc::Dimension Dim, bool State );
void resetDefaultArrayCompressionState( const qc::Dimension Dim );
  
const char *getDefaultArraySuffix ( const qc::Dimension Dim );

//! computes largest p such that 2^p <= Num
inline unsigned logBaseTwo ( int Num ) {
  if ( Num < 1 ) {
    cerr << "ERROR: logBaseTwo: invalid parameter: " << Num << endl;
    return 0;
  }

  int p = 0;
  while ( ( 1 << p ) < Num ) {
    p++;
  }
  if ( 1 << p != Num ) {
    p--;
  }
  return p;
}

//! computes smallest p such that 2^p >= Num
inline unsigned logBaseTwoCeil ( int Num ) {
  if ( Num < 1 ) {
    cerr << "ERROR: logBaseTwoCeil: invalid parameter: " << Num << endl;
    return 0;
  }

  int p = 0;
  while ( ( 1 << p ) < Num ) {
    p++;
  }
  return p;
}

//! inverse-lexicographical combination of three integers
inline int ILexCombine3 ( const int x, const int y, const int z, const int numX, const int numY ) {
  return (  ( z * numY +  y ) * numX + x );
}

//! inverse-lexicographical combination of two integers
inline int ILexCombine2 ( const int x, const int y, const int numX ) {
  return ( y * numX + x );
}

//! @enum Doc...
enum DiffVarType {
  DIFF_X,
  DIFF_Y,
  DIFF_Z,
  DIFF_XX,
  DIFF_YY,
  DIFF_ZZ,
  DIFF_YZ,
  DIFF_XZ,
  DIFF_XY
};

}

#endif
