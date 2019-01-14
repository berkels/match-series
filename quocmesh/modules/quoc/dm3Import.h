#ifndef __DM3IMPORT_H
#define __DM3IMPORT_H

#include <scalarArray.h>

namespace qc {

/**
 * Class to parse version 3 and 4 of the Gatan Digital Micrograph file format.
 *
 * \note The DM3/DM4 format parsing should work with most DM3/DM4 files, but the actual data extraction is
 *       rather hackish and only works with a small subset of files compliant with the DM3/DM4 format.
 *
 * \author Berkels
 */
class DM3Reader {
  struct DM3ImageData {
    int num;
    int offset;
    int dataSize;
    int dataType;
    int numX;
    int numY;
    int numZ;
    int declaredType;
    string name;
    DM3ImageData ( )
      : num ( -1 ),
        offset ( -1 ),
        dataSize ( -1 ),
        dataType ( -1 ),
        numX ( -1 ),
        numY ( -1 ),
        numZ ( -1 ),
        declaredType ( -1 ) { }

    qc::SaveType getSaveType ( ) const {
      if ( dataType == 3 )
        return qc::PGM_SIGNED_INT_BINARY;
      else if ( dataType == 2 )
        return qc::PGM_SHORT_BINARY;
      else if ( dataType == 4 )
        return qc::PGM_UNSIGNED_SHORT_BINARY;
      else if ( dataType == 5 )
        return qc::PGM_UNSIGNED_INT_BINARY;
      else if ( dataType == 6 )
        return qc::PGM_FLOAT_BINARY;
      else if ( dataType == 10 )
        return qc::PGM_UNSIGNED_CHAR_BINARY;
      else
        throw aol::UnimplementedCodeException ( aol::strprintf ( "Data type %d not implemented ", dataType ).c_str(), __FILE__, __LINE__ );
    }
  };

  std::vector<std::string> _currentTag;
  std::vector<DM3ImageData> _dataEntries;
  ifstream _in;
  int _version;

  int64_t readInterger ( ifstream &In ) {
    if ( _version == 3 )
      return aol::swapByteOrder<int32_t> ( aol::readBinaryData<int32_t, int32_t> ( In ) );
    else
      return aol::swapByteOrder<int64_t> ( aol::readBinaryData<int64_t, int64_t> ( In ) );
  }

public:
  DM3Reader ( const string &InputFileName )
   : _in ( InputFileName.c_str(), ios::binary ), _version ( 3 ) {

    if ( _in.good() == false )
      throw aol::FileException ( aol::strprintf ( "Cannot open file %s for reading", InputFileName.c_str() ).c_str(), __FILE__, __LINE__ );

    _version = aol::swapByteOrder<int32_t> ( aol::readBinaryData<int32_t, int32_t> ( _in ) );
    if ( ( _version != 3 ) && ( _version != 4 ) )
      throw aol::IOException ( aol::strprintf ( "Unexpected version (got %d), expected 3 or 4", _version ).c_str(), __FILE__, __LINE__ );

     /*const int fileLength =*/ readInterger ( _in );

    const int byteOrder = aol::swapByteOrder<int32_t> ( aol::readBinaryData<int32_t, int32_t> ( _in ) );
    if ( byteOrder != 1 )
      throw aol::UnimplementedCodeException ( "Reading of detected byte order not implemented", __FILE__, __LINE__ );

    readSortedClosed ( _in );

    /*const int numTags =*/ readInterger ( _in );

    do {
      if ( readTagOrTagDir ( _in ) == 0 )
        break;

      _currentTag.clear();
    } while ( true );
  }

  const DM3ImageData& getDataEntry ( ) const {
    if ( _dataEntries.size() > 0 )
      return _dataEntries.back();

    throw aol::Exception ( "No data entry found.", __FILE__, __LINE__ );
  }


  template <typename DataType>
  void exportDataToScalarArray ( qc::ScalarArray<DataType, qc::QC_2D> &Array ) {
    const DM3ImageData& dataEntry = getDataEntry();
    _in.seekg ( dataEntry.offset );
    Array.loadRaw ( _in, dataEntry.getSaveType(), dataEntry.numX, dataEntry.numY );
  }

  template <typename DataType>
  void exportDataToScalarArray ( qc::ScalarArray<DataType, qc::QC_3D> &Array ) {
    const DM3ImageData& dataEntry = getDataEntry();
    _in.seekg ( dataEntry.offset );
    Array.loadRaw ( _in, dataEntry.getSaveType(), dataEntry.numX, dataEntry.numY, ( dataEntry.numZ > 1 ) ? dataEntry.numZ : 1 );
  }

  void saveDataAsScalarArray ( const string &OutBaseName, const bool WriteSlices = false ) {
    const DM3ImageData& dataEntry = getDataEntry();
    if ( dataEntry.numZ > 1 ) {
      qc::ScalarArray<double, qc::QC_3D> a;
      exportDataToScalarArray ( a );
      if ( WriteSlices )
        a.saveSlices( ( string ( OutBaseName ) + "_%03d" + qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str(), qc::QC_Z, qc::PGM_DOUBLE_BINARY );
      else
        a.save ( ( string ( OutBaseName ) + qc::getDefaultArraySuffix ( qc::QC_3D ) ).c_str(), qc::PGM_DOUBLE_BINARY );
    }
    else {
      qc::ScalarArray<double, qc::QC_2D> a;
      exportDataToScalarArray ( a );
      a.save ( ( string ( OutBaseName ) + qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str(), qc::PGM_DOUBLE_BINARY );
    }
  }

  template <typename RealType>
  void saveQuocDataInDM3Container ( const qc::ScalarArray<RealType, qc::QC_2D> &QuocData, const string &OutBaseName ) {
    const int offset = getDataEntry().offset;

    // get length of file:
    _in.seekg ( 0, _in.end );
    const int length = _in.tellg();
    _in.seekg ( 0, _in.beg );
    const string outFileName = ( OutBaseName + ( ( _version == 3 ) ? ".dm3" : ".dm4" ) );
    std::ofstream out ( outFileName.c_str(), ios::binary );

    if ( out.good() == false )
      throw aol::FileException ( aol::strprintf ( "qc::DM3Reader::saveQuocDataInDM3Container: Cannot open \"%s\" for writing", outFileName.c_str() ).c_str(), __FILE__, __LINE__ );

    for ( int i = 0; i < offset; ++i )
      out.put ( _in.get() );

    qc::ScalarArray<RealType, qc::QC_2D> paddedArray ( getDataEntry().numX, getDataEntry().numY );
    paddedArray.padFrom ( QuocData );
    for ( int i = 0; i < getDataEntry().dataSize; ++i )
      writeValue ( out, getDataEntry().dataType, paddedArray[i] );

    _in.seekg ( out.tellp(), _in.beg );
    for ( int i = out.tellp(); i < length; ++i )
      out.put ( _in.get() );
  }
  
  template <typename RealType>
  void saveQuocDataInDM3Container ( const qc::ScalarArray<RealType, qc::QC_3D> &QuocData, const string &OutBaseName ) {
    const int offset = getDataEntry().offset;
    
    // get length of file:
    _in.seekg ( 0, _in.end );
    const int length = _in.tellg();
    _in.seekg ( 0, _in.beg );
    const string outFileName = ( OutBaseName + ( ( _version == 3 ) ? ".dm3" : ".dm4" ) );
    std::ofstream out ( outFileName.c_str(), ios::binary );
    
    if ( out.good() == false )
      throw aol::FileException ( aol::strprintf ( "qc::DM3Reader::saveQuocDataInDM3Container: Cannot open \"%s\" for writing", outFileName.c_str() ).c_str(), __FILE__, __LINE__ );
    
    for ( int i = 0; i < offset; ++i )
      out.put ( _in.get() );
    
    qc::ScalarArray<RealType, qc::QC_3D> paddedArray ( getDataEntry().numX, getDataEntry().numY, getDataEntry().numZ );
    paddedArray.padFrom ( QuocData );
    for ( int i = 0; i < getDataEntry().dataSize; ++i )
      writeValue ( out, getDataEntry().dataType, paddedArray[i] );
    
    _in.seekg ( out.tellp(), _in.beg );
    for ( int i = out.tellp(); i < length; ++i )
      out.put ( _in.get() );
  }
protected:
  void readTagName ( ifstream &in, const int DirNum ) {
    const int tagNameLength = aol::swapByteOrder<int16_t> ( aol::readBinaryData<int16_t, int16_t> ( in ) );
    if ( tagNameLength < 0 )
      throw aol::IOException ( "negative tag name length.", __FILE__, __LINE__ );

    char *tagName = new char[tagNameLength+1];
    aol::readBinaryData<char, char> ( in, tagName, tagNameLength );
    tagName[tagNameLength] = 0;
    string tagNameString = tagName;
    if ( tagNameLength > 0 )
      _currentTag.push_back ( tagName );
    else if ( DirNum > -1 )
      _currentTag.push_back ( aol::strprintf ( "%d", DirNum ) );
    else
      _currentTag.push_back ( "nameless" );
    delete[] tagName;

    if ( ( _currentTag.size() == 3 )
      && ( _currentTag[0] == "ImageList" )
      && ( _currentTag[2] == "ImageData" ) ) {
      _dataEntries.push_back( DM3ImageData() );
      _dataEntries.back().num = atoi ( _currentTag[1].c_str() );
    }
  }

  string getCurrentTag ( ) const {
    string tag;
    for ( unsigned int i = 0; i < _currentTag.size(); ++i ) {
      tag += _currentTag[i];
      if ( i < _currentTag.size() - 1 )
        tag += ".";
    }
    return tag;
  }

  void readSortedClosed ( ifstream &in ) {
    /*const int sorted =*/ aol::readBinaryData<unsigned char, int> ( in );
    /*const int closed =*/ aol::readBinaryData<unsigned char, int> ( in );
  }

  template <typename DataType>
  DataType parseValue ( ifstream &in ) {
    return aol::readBinaryData<DataType, DataType> ( in );
  }

  void readValue ( ifstream &in, const int Type ) {
    switch( Type ) {
    case 2:
      parseValue<int16_t> ( in );
      break;
    case 3:
      parseValue<int32_t> ( in );
      break;
    case 4:
      parseValue<uint16_t> ( in );
      break;
    case 5:
      parseValue<uint32_t> ( in );
      break;
    case 6:
      parseValue<float> ( in );
      break;
    case 7:
      parseValue<double> ( in );
      break;
    case 8:
    case 10:
      parseValue<signed char> ( in );
      break;
    case 11:
    case 12:
      parseValue<int64_t> ( in );
      break;
    default:
      throw aol::UnimplementedCodeException ( aol::strprintf ( "type %d not implemented ", Type ).c_str(), __FILE__, __LINE__ );
    }
  }

  template<typename RealType>
  void writeValue ( ofstream &out, const int Type, const RealType Value ) {
    switch( Type ) {
      case 3:
        aol::writebinary<int32_t> ( out, Value );
        break;
      case 4:
        aol::writebinary<uint16_t> ( out, Value );
        break;
      case 5:
        aol::writebinary<uint32_t> ( out, Value );
        break;
      case 6:
        aol::writebinary<float> ( out, Value );
        break;
      case 10:
        aol::writebinary<uint8_t> ( out, aol::Clamp<RealType> ( Value, 0, 255 ) );
        break;
      default:
        throw aol::UnimplementedCodeException ( aol::strprintf ( "type %d not implemented ", Type ).c_str(), __FILE__, __LINE__ );
    }
  }

  void readTag ( ifstream &in ) {
    // Read "%%%%"
    for ( int i = 0; i < 4; ++i )
      if ( aol::readBinaryData<unsigned char, unsigned char> ( in ) != '%' )
        throw aol::IOException ( "Unexpected character read", __FILE__, __LINE__ );

    const int64_t infoSize = readInterger( in );

    aol::Vector<int64_t> info ( infoSize );
    for ( int i = 0; i < infoSize; ++i )
      info[i] = readInterger ( in );

    if ( infoSize < 1 )
      throw aol::IOException ( "Info size must be positive", __FILE__, __LINE__ );
    else if ( infoSize == 1 ) {
      if ( ( _dataEntries.size() > 0 ) && ( _currentTag.size() == 5 )
        && ( _currentTag[0] == "ImageList" )
        && ( _dataEntries.back().num == atoi ( _currentTag[1].c_str() ) )
        && ( _currentTag[2] == "ImageData" )
        && ( _currentTag[3] == "Dimensions" )
        && ( info[0] == 5 ) ) {
        const int dir = atoi ( _currentTag[4].c_str() );
        const int numXYZ = aol::readBinaryData<uint32_t, int> ( in );
        if ( dir == 0 )
          _dataEntries.back().numX = numXYZ;
        else if ( dir == 1 )
          _dataEntries.back().numY = numXYZ;
        else if ( dir == 2 )
          _dataEntries.back().numZ = numXYZ;
      }
      else if ( ( _dataEntries.size() > 0 ) && ( _currentTag.size() == 4 )
        && ( _currentTag[0] == "ImageList" )
        && ( _dataEntries.back().num == atoi ( _currentTag[1].c_str() ) )
        && ( _currentTag[2] == "ImageData" )
        && ( _currentTag[3] == "DataType" )
        && ( info[0] == 5 ) ) {
        _dataEntries.back().declaredType = aol::readBinaryData<uint32_t, int> ( in );
      }
      else
        readValue ( in, info[0] );
    }
    else {
      const int numberType = info[0];

      switch( numberType ) {
      // 15 = group of data
      case 15:
        readGroup ( in, info );
        break;
      // 20 = array
      case 20:
        if ( info[1] == 15 )
          readArrayOfGroups ( in, info );
        else {
          if ( ( _dataEntries.size() > 0 ) && ( _currentTag.size() == 4 )
            && ( _currentTag[0] == "ImageList" )
            && ( _dataEntries.back().num == atoi ( _currentTag[1].c_str() ) )
            && ( _currentTag[2] == "ImageData" )
            && ( _currentTag[3] == "Data" ) ) {
              _dataEntries.back().offset = in.tellg();
              _dataEntries.back().dataType = info[1];
              _dataEntries.back().dataSize = info[2];
          }

          for ( int i = 0; i < info[2]; ++i )
            readValue ( in, info[1] );
        }
        break;
      default:
        throw aol::UnimplementedCodeException ( aol::strprintf ( "numberType %d not implemented", numberType ).c_str(), __FILE__, __LINE__ );
      }
    }
    _currentTag.pop_back();
  }

  void readTagDirectory ( ifstream &in ) {
    readSortedClosed ( in );
    const int64_t numTags = readInterger ( in );
    for ( int i = 0; i < numTags; ++i )
      readTagOrTagDir ( in, i );
    _currentTag.pop_back();
  }

  int readTagOrTagDir ( ifstream &in, const int DirNum = -1 ) {
    const int tagIdentifier = aol::readBinaryData<unsigned char, int> ( in );

    if ( tagIdentifier == 0 )
      return 0;

    readTagName ( in, DirNum );
    if ( _version == 4 )
      aol::swapByteOrder<int64_t> ( aol::readBinaryData<int64_t, int64_t> ( in ) );

    if ( tagIdentifier == 21 )
      readTag ( in );
    else if  ( tagIdentifier == 20 )
      readTagDirectory ( in );
    else
      throw aol::IOException ( aol::strprintf ( "Invalid tag identifier %d", tagIdentifier ), __FILE__, __LINE__ );

    return tagIdentifier;
  }

  void readGroup ( ifstream &in, const aol::Vector<int64_t> &info ) {
    if ( info[1] != 0 )
      throw aol::IOException ( "length of groupname has to be zero", __FILE__, __LINE__ );

    const int numGroupEntries = info[2];

    if ( ( info.size() - 3 ) != 2*numGroupEntries )
      throw aol::IOException ( "inconsistent sizes", __FILE__, __LINE__ );

    aol::Vector<int> groupEntryTypes ( numGroupEntries );
    for ( int i = 0; i < numGroupEntries; ++i ) {
      if ( info[3+2*i] != 0 )
        throw aol::IOException ( "length of fieldname has to be zero", __FILE__, __LINE__ );
      groupEntryTypes[i] = info[4+2*i];
    }

    for ( int i = 0; i < numGroupEntries; ++i )
      readValue ( in, groupEntryTypes[i] );
  }

  void readArrayOfGroups ( ifstream &in, const aol::Vector<int64_t> &info ) {
    if ( info[2] != 0 )
      throw aol::IOException ( "length of groupname has to be zero", __FILE__, __LINE__ );

    const int numGroupEntries = info[3];

    if ( ( info.size() - 5 ) != 2*numGroupEntries )
      throw aol::IOException ( "inconsistent sizes", __FILE__, __LINE__ );

    aol::Vector<int> groupEntryTypes ( numGroupEntries );
    for ( int i = 0; i < numGroupEntries; ++i ) {
      if ( info[4+2*i] != 0 )
        throw aol::IOException ( "length of fieldname has to be zero", __FILE__, __LINE__ );
      groupEntryTypes[i] = info[5+2*i];
    }

    for ( int j = 0; j < info[4+2*numGroupEntries]; ++j )
      for ( int i = 0; i < numGroupEntries; ++i )
        readValue ( in, groupEntryTypes[i] );
  }

};

} // namespace qc

#endif // __DM3IMPORT_H
