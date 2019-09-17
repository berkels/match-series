/**
 * \file
 * \brief Converts a FEI TIA SER file to 2D ScalarArrays with double precision.
 *
 * C++ port of http://ncem.lbl.gov/software_exchange/serReader.m with some adjustments (original Matlab code by Peter Ercius)
 * and an update ported from https://bitbucket.org/ercius/openncem
 *
 * Usage: convertSERToQuoc InputFile
 *
 * \author Berkels
 */

#include <dm3Import.h>

qc::SaveType getSerSaveType ( const int DataType ) {
  switch ( DataType ) {
    case 1:
      return qc::PGM_UNSIGNED_CHAR_BINARY;
    case 2:
      return qc::PGM_UNSIGNED_SHORT_BINARY;
    case 3:
      return qc::PGM_UNSIGNED_INT_BINARY;
    case 4:
      throw aol::UnimplementedCodeException ( "Reading signed char data from SER files is not implemented yet.", __FILE__, __LINE__ );
    case 5:
      return qc::PGM_SHORT_BINARY;
    case 6:
      return qc::PGM_SIGNED_INT_BINARY;
    case 7:
      return qc::PGM_FLOAT_BINARY;
    case 8:
      return qc::PGM_DOUBLE_BINARY;
    default:
      throw aol::Exception ( "Unsupported data type", __FILE__, __LINE__ );
  }
}

int main ( int argc, char **argv ) {

  try {
    if ( argc < 2 ) {
      cerr << "USAGE: " << argv[0] << "  <InputFile>" << endl;
      return EXIT_FAILURE;
    }

    const string inFileName = argv[1];

    std::ifstream in ( inFileName.c_str(), std::ios::binary );

    if ( in.good() == false )
      throw aol::IOException ( aol::strprintf ( "Error opening file \"%s\"", inFileName.c_str() ).c_str(), __FILE__, __LINE__ );

    const int byteOrder = aol::readBinaryData<int16_t, int> ( in );
    if ( byteOrder != 0x4949 )
      throw aol::Exception ( "Unsupported byte order.", __FILE__, __LINE__ );

    /*const int seriesID =*/ aol::readBinaryData<int16_t, int> ( in );
    const int seriesVersion = aol::readBinaryData<int16_t, int> ( in );
    const int dataTypeID = aol::readBinaryData<int32_t, int> ( in );

    const int tagTypeID = aol::readBinaryData<int32_t, int> ( in );
    if ( tagTypeID == 0x4152 )
      cerr << "Tag is time only\n";
    else if ( tagTypeID == 0x4142 )
      cerr << "Tag is 2D with time\n";
    else
      cerr << "Unkown tag\n";

    // the number of data elements in the original data set
    const int totalNumberElements = aol::readBinaryData<int32_t, int> ( in );
    // the number of data elements written to the file
    const int validNumberElements = aol::readBinaryData<int32_t, int> ( in );
    // the offset (in bytes) to the beginning of the data offset array
    const int offsetArrayOffset = aol::readBinaryData<int32_t, int> ( in );
    // the number of dimensions of the indices (not the data)
    const int numberDimensions = aol::readBinaryData<int32_t, int> ( in );
    cerr << " numberDimensions " << numberDimensions << endl;

    // Dimension arrays (starts at byte offset 30)
    for ( int i = 0; i < numberDimensions; ++i ) {
      cerr << "dimension " << i << ":\n";
      cerr << "dimension size  = " << aol::readBinaryData<int32_t, int> ( in ) << endl;

      // calibration value at element calibrationElement
      /*const double calibrationOffset =*/ aol::readBinaryData<double, double> ( in );
      // calibration delta between elements of the series
      /*const double calibrationDelta =*/ aol::readBinaryData<double, double> ( in );
      // the element in the series with a calibraion value of calibrationOffset
      /*const int calibrationElement =*/ aol::readBinaryData<int32_t, int> ( in );

      {
        // length of the description string
        const int descriptionLength = aol::readBinaryData<int32_t, int> ( in );
        // describes the dimension
        char *description = new char[descriptionLength+1];
        description[descriptionLength] = 0;
        aol::readBinaryData<char, char> ( in, description, descriptionLength );
        cerr << "description = " << description << endl;
        delete[] description;
      }

      {
        // length of the units string
        const int unitsLength = aol::readBinaryData<int32_t, int> ( in );
        // name of the units in this dimension
        char *units = new char[unitsLength+1];
        units[unitsLength] = 0;
        aol::readBinaryData<char, char> ( in, units, unitsLength );
        cerr << "units = " << units << endl;
        delete[] units;
        
      }
    }

    // Get arrays containing byte offsets for data and tags
    // seek in the file to the offset arrays
    in.seekg ( offsetArrayOffset );
    // Data offset array (the byte offsets of the individual data elements)
    aol::Vector<int> dataOffsetArray ( totalNumberElements );
    if ( seriesVersion <= 528 )
      aol::readBinaryData<int32_t, int> ( in, dataOffsetArray.getData(), totalNumberElements );
    else
      aol::readBinaryData<int64_t, int> ( in, dataOffsetArray.getData(), totalNumberElements );

    // Tag Offset Array (byte offsets of the individual data tags)
    aol::Vector<int> tagOffsetArray ( totalNumberElements );
    if ( seriesVersion <= 528 )
      aol::readBinaryData<int32_t, int> ( in, tagOffsetArray.getData(), totalNumberElements );
    else
      aol::readBinaryData<int64_t, int> ( in, tagOffsetArray.getData(), totalNumberElements );
    switch ( dataTypeID ) {
      // 1D data
      case 0x4120:
        {
          aol::MultiVector<double> data ( validNumberElements, 0 );
          for ( int i = 0; i < validNumberElements; ++i ) {
            in.seekg ( dataOffsetArray[i] );
            // calibration value at element calibrationElement
            /*const double calibrationOffset =*/ aol::readBinaryData<double, double> ( in );
            // calibration delta between elements of the array
            /*const double calibrationDelta =*/ aol::readBinaryData<double, double> ( in );
            // element in the array with calibration value of calibrationOffset
            /*const int calibrationElement =*/ aol::readBinaryData<int32_t, int> ( in );

            const int dataType = aol::readBinaryData<int16_t, int> ( in );

            // number of elements in the array
            const int num = aol::readBinaryData<int32_t, int> ( in );
            data[i].resize( num );
            data[i].loadRaw ( in, getSerSaveType ( dataType ) );
          }
          
          if ( data.allDimsEqual() ) {
            qc::ScalarArray<double, qc::QC_2D> a ( data[0].size(), validNumberElements);
            a.copyUnblockedFrom ( data );
            a.save ( aol::strprintf ( "%s%s", aol::getBaseFileName( inFileName ).c_str(), qc::getDefaultArraySuffix( qc::QC_2D ) ).c_str(), qc::PGM_DOUBLE_BINARY );
          }
          else {
            for ( int i = 0; i < validNumberElements; ++i ) {
              qc::ScalarArray<double, qc::QC_1D> a ( data[i], aol::FLAT_COPY );
              a.save ( aol::strprintf ( "%s%04i%s", aol::getBaseFileName( inFileName ).c_str(), i, qc::getDefaultArraySuffix( qc::QC_1D ) ).c_str(), qc::PGM_DOUBLE_BINARY );
            }
            
          }
        }
        break;
      // Get data from 2D elements
      case 0x4122:
        {
          for ( int i = 0; i < validNumberElements; ++i ) {
            in.seekg ( dataOffsetArray[i] );
            // calibration at element calibrationElement along x
            /*const double calibrationOffsetX =*/ aol::readBinaryData<double, double> ( in );
            /*const double calibrationDeltaX =*/ aol::readBinaryData<double, double> ( in );
            // element in the array along x with calibration value of calibrationOffset
            /*const int calibrationElementX =*/ aol::readBinaryData<int32_t, int> ( in );
            /*const double calibrationOffsetY =*/ aol::readBinaryData<double, double> ( in );
            /*const double calibrationDeltaY =*/ aol::readBinaryData<double, double> ( in );
            /*const int calibrationElementY =*/ aol::readBinaryData<int32_t, int> ( in );
            const int dataType = aol::readBinaryData<int16_t, int> ( in );

            const int numX = aol::readBinaryData<int32_t, int> ( in );
            const int numY = aol::readBinaryData<int32_t, int> ( in );
            qc::ScalarArray<double, qc::QC_2D> a ( numX, numY );
            // SER files seemed to be vertically flipped compared to standard pixel images like PGM.
            qc::ScalarArray<double, qc::QC_2D> aFlipped ( numX, numY );
            aFlipped.loadRaw ( in, getSerSaveType ( dataType ), numX, numY );
            a.flipFrom ( aFlipped, qc::QC_Y );
            a.save ( aol::strprintf ( "%s%03i%s", aol::getBaseFileName( inFileName ).c_str(), i, qc::getDefaultArraySuffix( qc::QC_2D ) ).c_str(), qc::PGM_DOUBLE_BINARY );

            switch ( tagTypeID ) {
              case 0x4142:
                throw aol::UnimplementedCodeException ( "Reading Time-and-Position tags not implemented yet.", __FILE__, __LINE__ );
                break;
                // Get data from "Time-only" tags
              case 0x4152:
              {
                in.seekg ( tagOffsetArray[i] );
                // type of the tag (should be 0x4152)
                const int tagTypeID2 = aol::readBinaryData<int32_t, int> ( in );
                if ( tagTypeID2 != 0x4152 )
                  throw aol::Exception ( "Time tag mismatch.", __FILE__, __LINE__ );

                // number of seconds since Jan. 1, 1970
                const time_t timeTag = aol::readBinaryData<int32_t, time_t> ( in );
                char timeString[256];
                strftime ( timeString, 256, "%H:%M on %A, %d %B %Y", localtime ( &timeTag ) );
                cerr << "Element " << i << " created " << timeString << endl;
              }
                break;
              default:
                throw aol::Exception ( "Unexpected time tag type.", __FILE__, __LINE__ );
            }
          }
        }
        break;
      default:
        throw aol::Exception ( "Unexpected dimension.", __FILE__, __LINE__ );
    }
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
