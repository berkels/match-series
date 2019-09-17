#include <quoc.h>

namespace qc {

int getSizeOfSaveType ( const qc::SaveType Type ) {
  switch ( Type ) {
  case PGM_UNSIGNED_CHAR_ASCII:
  case PGM_UNSIGNED_CHAR_BINARY:
    return 1;
    break;
  case PGM_UNSIGNED_SHORT_BINARY:
  case PGM_SHORT_BINARY:
    return 2;
    break;
  case PGM_FLOAT_ASCII:
  case PGM_FLOAT_BINARY:
  case PGM_UNSIGNED_INT_BINARY:
  case PGM_SIGNED_INT_BINARY:
    return 4;
    break;
  case PGM_DOUBLE_BINARY:
    return 8;
    break;
  default:
    throw aol::Exception ( "qc::getSizeOfSaveType: Unsupported SaveType", __FILE__, __LINE__ );
    break;
  }
}

const char *getDefaulSuffixOfSaveType ( const qc::SaveType Type, const qc::Dimension Dim ) {
  switch ( Type ) {
    case PGM_UNSIGNED_CHAR_ASCII:
    case PGM_UNSIGNED_CHAR_BINARY:
      return ".pgm";
      break;
    case PGM_UNSIGNED_SHORT_BINARY:
    case PGM_SHORT_BINARY:
    case PGM_FLOAT_ASCII:
    case PGM_FLOAT_BINARY:
    case PGM_UNSIGNED_INT_BINARY:
    case PGM_SIGNED_INT_BINARY:
    case PGM_DOUBLE_BINARY:
      return getDefaultArraySuffix ( Dim );
      break;
    case PNG_2D:
      return ".png";
      break;
    case NETCDF:
      return ".nc";
      break;
    default:
      throw aol::Exception ( "qc::getDefaulSuffixOfSaveType: Unsupported SaveType", __FILE__, __LINE__ );
      break;
  }
}

static bool defaultArrayCompression_2D_enabled = true;
static bool defaultArrayCompression_3D_enabled = true;
  
bool getDefaultArrayCompressionState( const qc::Dimension Dim ) {
  if ( Dim == qc::QC_1D ) return false;
  else if ( Dim == qc::QC_2D ) return defaultArrayCompression_2D_enabled;
  else if ( Dim == qc::QC_3D ) return defaultArrayCompression_3D_enabled;
  else throw aol::Exception ( "qc::getDefaultArrayCompressionState: Unsupported Dimension", __FILE__, __LINE__ );
}
  
void setDefaultArrayCompressionState( const qc::Dimension Dim, bool State ) {
  if ( Dim == qc::QC_2D ) defaultArrayCompression_2D_enabled = State;
  else if ( Dim == qc::QC_3D ) defaultArrayCompression_3D_enabled = State;
  else throw aol::Exception ( "qc::getDefaultArrayCompressionState: Unsupported Dimension", __FILE__, __LINE__ );
}
  
void resetDefaultArrayCompressionState( const qc::Dimension Dim ) {
  if ( Dim == qc::QC_2D ) defaultArrayCompression_2D_enabled = true;
  else if ( Dim == qc::QC_3D ) defaultArrayCompression_3D_enabled = true;
  else throw aol::Exception ( "qc::getDefaultArrayCompressionState: Unsupported Dimension", __FILE__, __LINE__ );
}
  
const char *getDefaultArraySuffix ( const qc::Dimension Dim ) {
  if ( Dim == qc::QC_3D ) {
    if ( defaultArrayCompression_3D_enabled )
      return ".q3bz";
    else
      return ".q3";
  } else if ( Dim == qc::QC_2D ) {
    if ( defaultArrayCompression_2D_enabled )
      return ".q2bz";
    else
      return ".q2";
  } else if ( Dim == qc::QC_1D )
    return ".q1";
  else
    return ".dat.bz2";
}

template<> const qc::SaveType qc::SaveTypeTrait< unsigned char  >::AsciiSaveType = qc::PGM_UNSIGNED_CHAR_ASCII;
template<> const qc::SaveType qc::SaveTypeTrait< float          >::AsciiSaveType = qc::PGM_FLOAT_ASCII;

template<> const qc::SaveType qc::SaveTypeTrait< unsigned char  >::BinarySaveType = qc::PGM_UNSIGNED_CHAR_BINARY;
template<> const qc::SaveType qc::SaveTypeTrait< short          >::BinarySaveType = qc::PGM_SHORT_BINARY;
template<> const qc::SaveType qc::SaveTypeTrait< unsigned short >::BinarySaveType = qc::PGM_UNSIGNED_SHORT_BINARY;
template<> const qc::SaveType qc::SaveTypeTrait< float          >::BinarySaveType = qc::PGM_FLOAT_BINARY;
template<> const qc::SaveType qc::SaveTypeTrait< double         >::BinarySaveType = qc::PGM_DOUBLE_BINARY;

}
