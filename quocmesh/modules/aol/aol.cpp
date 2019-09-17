#include <aol.h>
#include <qmException.h>
#include <polynomial.h>
#include <indexMapper.h>


namespace aol {

template<class T> T returnArgument ( const T a ) {
  return a;
}

template float returnArgument<float> ( const float a );
template double returnArgument<double> ( const double a );
template long double returnArgument<long double> ( const long double a );
template short returnArgument<short> ( const short a );
template unsigned short returnArgument<unsigned short> ( const unsigned short a );
template int returnArgument<int> ( const int a );
template signed char returnArgument<signed char> ( const signed char a );
template unsigned char returnArgument<unsigned char> ( const unsigned char a );
template unsigned int returnArgument<unsigned int> ( const unsigned int a );
template int64_t returnArgument<int64_t> ( const int64_t a );
template uint64_t returnArgument<uint64_t> ( const uint64_t a );
#if ( defined __APPLE__ ) || ( BITNESS == 32 ) // This seems to be the only way to get SuiteSparse_long to work under OS X
template long returnArgument<long> ( const long a );
#endif

void assertAlert(string assertion, string file, int line, string message) {
  if (message == "")
    throw Exception(string ("Assertion failed: ") + assertion, file, line);
  else
    throw Exception(message + " (" + assertion + ")", file, line);
}

template <> const bool ZTrait<bool>::zero = false;
template <> const bool ZOTrait<bool>::one = true;

template <> const char ZTrait<char>::zero = 0;
template <> const char ZOTrait<char>::one = 1;

template <> const unsigned char ZTrait<unsigned char>::zero = 0;
template <> const unsigned char ZOTrait<unsigned char>::one = 1;

template <> const signed char ZTrait<signed char>::zero = 0;
template <> const signed char ZOTrait<signed char>::one = 1;

template <> const short ZTrait<short>::zero = 0;
template <> const short ZOTrait<short>::one = 1;

template <> const int ZTrait<int>::zero = 0;
template <> const int ZOTrait<int>::one = 1;

template <> const unsigned int ZTrait<unsigned int>::zero = 0;
template <> const unsigned int ZOTrait<unsigned int>::one = 1;

template <> const unsigned short ZTrait<unsigned short>::zero = 0;
template <> const unsigned short ZOTrait<unsigned short>::one = 1;

template <> const int64_t ZTrait<int64_t>::zero = 0;
template <> const int64_t ZOTrait<int64_t>::one = 1;

template <> const uint64_t ZTrait<uint64_t>::zero = 0;
template <> const uint64_t ZOTrait<uint64_t>::one = 1;

#if ( defined __APPLE__ ) || ( BITNESS == 32 ) // This seems to be the only way to get SuiteSparse_long to work under OS X
template <> const long ZTrait<long>::zero = 0;
template <> const long ZOTrait<long>::one = 1;
#endif

template <> const float ZTrait<float>::zero = 0.f;
template <> const float ZOTrait<float>::one = 1.f;

template <> const double ZTrait<double>::zero = 0.0;
template <> const double ZOTrait<double>::one = 1.0;

template <> const long double ZTrait<long double>::zero = 0.0l;
template <> const long double ZOTrait<long double>::one = 1.0l;

template <> const std::complex<double> ZTrait<std::complex<double> >::zero = 0.0;
template <> const std::complex<double> ZOTrait<std::complex<double> >::one = 1.0;


template <> float NumberTrait<float>::getNaN ( )             { return ( std::numeric_limits<float>::quiet_NaN() ); }
template <> double NumberTrait<double>::getNaN ( )           { return ( std::numeric_limits<double>::quiet_NaN() ); }
template <> long double NumberTrait<long double>::getNaN ( ) { return ( std::numeric_limits<long double>::quiet_NaN() ); }

template <> float NumberTrait<float>::getInf ( )             { return ( std::numeric_limits<float>::infinity() ); }
template <> double NumberTrait<double>::getInf ( )           { return ( std::numeric_limits<double>::infinity() ); }
template <> long double NumberTrait<long double>::getInf ( ) { return ( std::numeric_limits<long double>::infinity() ); }

template <> float NumberTrait<float>::getPi ( )              { return ( static_cast<float>( 3.1415926535897932384626433832795029L ) ); }
template <> double NumberTrait<double>::getPi ( )            { return ( static_cast<double>( 3.1415926535897932384626433832795029L ) ); }
template <> long double NumberTrait<long double>::getPi ( )  { return ( static_cast<long double> ( 3.1415926535897932384626433832795029L ) ); }

template <> const float NumberTrait<float>::NaN = NumberTrait<float>::getNaN();
template <> const float NumberTrait<float>::Inf = NumberTrait<float>::getInf();
template <> const float NumberTrait<float>::pi  = NumberTrait<float>::getPi();

template <> const double NumberTrait<double>::NaN = NumberTrait<double>::getNaN();
template <> const double NumberTrait<double>::Inf = NumberTrait<double>::getInf();
template <> const double NumberTrait<double>::pi  = NumberTrait<double>::getPi();

template <> const long double NumberTrait<long double>::NaN = NumberTrait<long double>::getNaN();
template <> const long double NumberTrait<long double>::Inf = NumberTrait<long double>::getInf();
template <> const long double NumberTrait<long double>::pi  = NumberTrait<long double>::getPi();

template <> const char MaxInitializerTrait<char>::MaxInitializer = numeric_limits<char>::min();
template <> const unsigned char MaxInitializerTrait<unsigned char>::MaxInitializer = numeric_limits<unsigned char>::min();
template <> const signed char MaxInitializerTrait<signed char>::MaxInitializer = numeric_limits<signed char>::min();
template <> const short MaxInitializerTrait<short>::MaxInitializer = numeric_limits<short>::min();
template <> const unsigned short MaxInitializerTrait<unsigned short>::MaxInitializer = numeric_limits<unsigned short>::min();
template <> const int MaxInitializerTrait<int>::MaxInitializer = numeric_limits<int>::min();
template <> const unsigned int MaxInitializerTrait<unsigned int>::MaxInitializer = numeric_limits<unsigned int>::min();
template <> const int64_t MaxInitializerTrait<int64_t>::MaxInitializer = numeric_limits<int64_t>::min();
template <> const uint64_t MaxInitializerTrait<uint64_t>::MaxInitializer = numeric_limits<uint64_t>::min();
#if ( defined __APPLE__ ) || ( BITNESS == 32 ) // This seems to be the only way to get SuiteSparse_long to work under OS X
template <> const long MaxInitializerTrait<long>::MaxInitializer = numeric_limits<long>::min();
#endif
template <> const float MaxInitializerTrait<float>::MaxInitializer = -NumberTrait<float>::getInf();
template <> const double MaxInitializerTrait<double>::MaxInitializer = -NumberTrait<double>::getInf();
template <> const long double MaxInitializerTrait<long double>::MaxInitializer = -NumberTrait<long double>::getInf();

template <> const char* TemplateToStringTrait<int>::GetString() { return "int"; }
template <> const char* TemplateToStringTrait<float>::GetString() { return "float"; }
template <> const char* TemplateToStringTrait<double>::GetString() { return "double"; }

const aol::SimpleFormat scientificFormat ( 17, 10, ios::scientific | ios::left | ios::showpos );
const aol::SimpleFormat longScientificFormat ( 25, 17, ios::scientific | ios::left | ios::showpos );
const aol::SimpleFormat shortFormat ( 6, 2, ios::fixed | ios::left | ios::showpos );
const aol::SimpleFormat intFormat ( 4, 0, ios::fixed | ios::right );
const aol::SimpleFormat longIntFormat ( 8, 0, ios::fixed | ios::right );
const aol::SimpleFormat fillFormat ( 6, 0, ios::fixed | ios::right, '0' );
const aol::MixedFormat mixedFormat ( 3, 6 );
const aol::MixedFormat detailedFormat ( 8, 12 );



string SimpleFormat::operator () ( const long double val ) const {
  ostringstream oss;
  oss.fill ( fillchar );
  oss << setw ( width ) << setprecision ( precision ) << setiosflags ( format );
  if ( val != 0.0 )
#if ( defined (__MINGW32_VERSION) || defined(__MINGW64__) ) && ( GCC_VERSION == 40601 )
    // Long double output is broken in MinGW 4.6.1. As workaround use double instead.
    oss << ( static_cast<double> ( val ) );
#else
    oss << val;
#endif
  else oss << resetiosflags ( ios::fixed | ios::scientific ) << 0;
  return oss.str ();
}

string SimpleFormat::operator () ( const int64_t val ) const {
  ostringstream oss;
  oss.fill ( fillchar );
  oss << setw ( width ) << setprecision ( precision ) << setiosflags ( format ) << resetiosflags ( ios::fixed | ios::scientific | ios::showpos ) << val;
  return oss.str ();
}

string SimpleFormat::operator () ( const uint64_t val ) const {
  ostringstream oss;
  oss.fill ( fillchar );
  oss << setw ( width ) << setprecision ( precision ) << setiosflags ( format ) << resetiosflags ( ios::fixed | ios::scientific | ios::showpos ) << val;
  return oss.str ();
}

string MixedFormat::operator () ( const long double val ) const {
  ostringstream oss;
  long double max = pow ( 10.0l, static_cast<long double> ( ( before - 1 ) ) );
  long double eps = pow ( 10.0l, static_cast<long double> ( -after ) );
  long double rval = eps * rint ( val / eps );
  long double l, r;
  r = fabs ( aol::Modf ( rval, &l ) );
  if ( val == std::numeric_limits<long double>::infinity () )
    oss << string ( before - 2, ' ' ) << setw ( 2 + 1 + after ) << setiosflags ( ios::left ) << " INF";
  else if ( val == - std::numeric_limits<long double>::infinity () )
    oss << string ( before - 2, ' ' ) << setw ( 2 + 1 + after ) << setiosflags ( ios::left ) << "-INF";
  else if ( !aol::isFinite ( val ) )
    oss << string ( before - 1, ' ' ) << setw ( 1 + 1 + after ) << setiosflags ( ios::left ) << "NAN";
  else if ( val == 0.0 )
    oss << setw ( before ) << setiosflags ( ios::right ) << 0 << string ( 1 + after, ' ' );
  else if ( ( fabs ( l ) >= max ) || ( l == 0.0 && r < sqrt ( eps ) ) ) {
// ios::scientific under MinGW / VC++ produces different output than with GCC under Linux.
#if defined (__MINGW32_VERSION) || defined(__MINGW64__) || defined (_MSC_VER)
    //! \todo Improve the rounding.
    const int exponent = static_cast<int> ( floor ( log10( aol::Abs ( val ) ) ) );
    const double coefficient = val * pow ( 10., static_cast<double> ( -exponent ) );
    if ( after > 4 )
      oss << aol::strprintf ( "%*.*f", before+after-3, after - 4, coefficient );
    else if ( ( after == 4 ) || ( after == 3 ) ) {
      oss << aol::strprintf ( "%*d", before, static_cast<int> ( floor ( coefficient ) ) );
      if ( after == 4 )
        oss << ".";
    }
    else
      throw ( aol::Exception ( "Scientific notation needs \"after >= 3\"!", __FILE__, __LINE__ ) );
    oss << aol::strprintf ( "e%s%02d", ( exponent >= 0 ) ? "+" : "-", aol::Abs ( exponent ) );
#else
    oss << setprecision ( aol::Max ( after - 4, 0 ) ) << setiosflags ( ios::left | ios::scientific )
    << string ( before - 2, ' ' ) << ( val >= 0.0 ? " " : "" ) << val;
    if ( after == 4 ) {
      // Fix because scientific with precision "0" does not print a "."
      string str = oss.str ();
      return string ( str, 0, before ) + "." + string ( str, before, string::npos );
    }
#endif
  }
  else {
#if !defined (__MINGW32_VERSION) && !defined(__MINGW64__) // MinGW doesn't seem to support long double precision with printf.
    oss << aol::strprintf ( "%*.*Lf", before+after+1, after, val );
#else
    oss << aol::strprintf ( "%*.*f", before+after+1, after, static_cast<double>(val) );
#endif
  }
  return oss.str ();
}

string MixedFormat::operator () ( const int64_t val ) const {
  ostringstream oss;
  oss << setw ( before ) << setiosflags ( ios::right ) << val << string ( 1 + after, ' ' );
  return oss.str ();
}

string MixedFormat::operator () ( const uint64_t val ) const {
  ostringstream oss;
  oss << setw ( before ) << setiosflags ( ios::right ) << val << string ( 1 + after, ' ' );
  return oss.str ();
}

string generateCurrentTimeAndDateString ( ) {
  char   today[256];
  time_t now = time ( NULL );
  strftime ( today, 256, "%H:%M on %A, %d %B %Y", localtime ( &now ) );
  return string ( today );
}

const char* getFileName( const char* FileNameWithPath ){
  const char* temp = strrchr(FileNameWithPath, '/');
  // No occurrence of '/', i.e. FileNameWithPath doesn't contain a Unix path.
  if ( temp == NULL ) {
#if defined ( _WIN32 ) || defined ( _WIN64 )
    const char* tempTwo = strrchr(FileNameWithPath, '\\');
    // No occurrence of '\\', i.e. FileNameWithPath also doesn't contain a Windows path.
    if ( tempTwo == NULL )
      return FileNameWithPath;
    else
      return (tempTwo + 1);
#else
    return FileNameWithPath;
#endif
  }
  else
    return (temp + 1);
}

string getBaseFileName ( const string &FileNameWithPath ) {
  string outputFileNameBase = aol::getFileName ( FileNameWithPath.c_str() );

  size_t found = outputFileNameBase.rfind ( '.' );
  if ( found != string::npos )
    outputFileNameBase = outputFileNameBase.substr ( 0, found );

  return outputFileNameBase;
}

string getPath ( const string &FileNameWithPath ) {
  string path = FileNameWithPath;
  size_t pos = path.rfind ( aol::getFileName ( FileNameWithPath.c_str() ) );
  if ( pos != std::string::npos )
    path.erase ( pos );
  return path;
}
  
string expandTildeOrEnvVar ( const string &InputString ) {
  if ( InputString.length() == 0 )
    return InputString;

  if ( InputString[0] == '~' ) {
    const char *homeVar = getenv ( "HOME" );
    string tempString = string ( homeVar ) + ( InputString.c_str() + 1 );
    return tempString;
  }
  else if ( ( InputString[0] == '$' ) && ( InputString.length() > 1 ) && ( InputString[1] == '{' ) ){
    string tempString = InputString;
    std::size_t pos = tempString.find('}');
    if ( pos != std::string::npos ) {
      const string varName = tempString.substr ( 2, pos - 2 );
      const char *varValue = getenv( varName.c_str() );
      if ( varValue )
        return varValue + tempString.substr ( pos + 1 );
      else
        throw Exception ( aol::strprintf ( "Undefined environment variable '%s'", varName.c_str() ).c_str(), __FILE__, __LINE__ );
    }
    else
      throw Exception ( aol::strprintf ( "Invalid environment variable format, full string is '%s'", InputString.c_str() ).c_str(), __FILE__, __LINE__ );
  }
  else
    return InputString;
}
  
bool fileNameEndsWith ( const char* FileName, const char* Ending ){
  // Trim trailing whitespace to ignore it in the comparison.
  string s = FileName;
  string whitespaces (" \t\f\v\n\r");
  size_t found;
  found = s.find_last_not_of(whitespaces);
  if ( found != string::npos )
    s.erase(found+1);
  else
    s.clear();

  if ( strlen ( s.c_str() ) <= strlen ( Ending ) )
    return false;

  return ( strcasecmp ( s.c_str() + strlen ( s.c_str() ) - strlen ( Ending ), Ending ) == 0 );
}

bool hasBzipSuffix ( const char *FileName ) {
  return ( fileNameEndsWith ( FileName, ".bz2" ) || fileNameEndsWith ( FileName, ".q2bz" ) || fileNameEndsWith ( FileName, ".q3bz" ) );
}

bool hasZlibSuffix ( const char *FileName ) {
  return fileNameEndsWith ( FileName, ".zz" );
}

bool hasGzipSuffix ( const char *FileName ) {
  return ( fileNameEndsWith ( FileName, ".gz" ) || fileNameEndsWith ( FileName, ".svgz" ) );
}

string getOutFileName ( const string &inFile, bool fromEnd ) {

  int firstPt, secPt;

  if ( fromEnd ) {
    secPt   = inFile.rfind ( '.' );
    firstPt = inFile.rfind ( '.', secPt - 1 );
  } else {
    firstPt = inFile.find ( '.' );
    secPt   = inFile.find ( '.', firstPt + 1 );
  }

  // if inFileName doesn't contain a dot:
  if ( fromEnd ? ( secPt == -1 ) : ( firstPt == -1 ) ) {
    return string ( inFile ) + ".0001";
  }

  // if inFileName contains only one dot, put "0001" in front of ending:
  if ( fromEnd ? ( firstPt == -1 ) : ( secPt == -1 ) ) {
    string before ( inFile, 0, ( fromEnd ? secPt : firstPt ) + 1 );
    string after ( inFile, ( fromEnd ? secPt : firstPt ), inFile.size() );
    return before + "0001" + after;
  }

  // if inFileName contains at least two dots,
  // increment number between the first two dots:
  if ( firstPt != -1 &&
       secPt   != -1 ) {
    string before ( inFile, 0, firstPt + 1 );
    string between ( inFile, firstPt + 1, secPt - ( firstPt + 1 ) );
    string after ( inFile, secPt, inFile.size() );
    int origBetweenSize = between.size();

    int num = atoi ( between.c_str() );
    num++;
    between = to_string ( num );
    int newBetweenSize = between.size();

    for ( int i = 0; i < origBetweenSize - newBetweenSize; ++i ) {
      between.insert ( 0, "0" );
    }
    string outFile = before + between + after;
    return outFile;
  }

  throw ( aol::Exception ( "Case shouldn't happen.", __FILE__, __LINE__ ) );
  return ( "" );
}

string strprintf(const char * format, ...) {
  // declare variable argument list
  va_list az;

  // copy from my input into this list
  // (second argument is nothing really
  // used, but the last named argument
  // of myself)
  va_start(az, format);

  // give this argument list variable
  // to vscprintf instead of my own
  // arguments (that is in what functions like vsprintf
  // differ from functions like sprintf)
  const int sizeNeeded = aol::vscprintf(format, az) + 1;

  // restore stack into clean state:
  va_end(az);

  if ( sizeNeeded < 0 )
    throw aol::Exception ( "strprintf encountered an error when determining the necessary buffer size\n", __FILE__, __LINE__ );

  char *buffer = new char[sizeNeeded];

  va_start(az, format);
  vsprintf (buffer, format, az);
  va_end(az);

  // automatic return type conversion into string:
  string ret = buffer;
  delete[] buffer;
  return ret;
}

void checkNextLineOrString ( istream &In, const char *ExpectedString, const bool Line ) {
  string temp;
  if ( Line )
    getline ( In, temp );
  else
    In >> temp;

  if ( temp.compare ( ExpectedString ) ) {
    string errorMessage = aol::strprintf ( "Unknown header: Expected \n\"%s\"\n, got \n\"%s\".\n", ExpectedString, temp.c_str() );
    throw aol::Exception ( errorMessage, __FILE__, __LINE__ );
  }
}

int countDigitsOfNumber ( const int N ) {
  return ( N != 0 ) ? static_cast<int>(floor(log10(static_cast<double>(aol::Abs(N)))))+1 : 1;
}

unsigned int crc32 ( const void* const ptr, const unsigned int size ) {
  const unsigned char* const data = static_cast< const unsigned char* >( ptr );
  const unsigned int crc32Poly = 0xEDB88320;
  unsigned int dwarf = 0xffffffff;
  for ( uint64_t i = 0; i < 8 * size; ++i ) {
    const unsigned int entry = i / ( 8 * sizeof( unsigned char ) );
    const unsigned char shift = ( i % ( sizeof ( unsigned char ) * 8 ) );
    const bool currentBit = ( data[ entry ] >> shift ) & 1;
    if ( static_cast<bool> ( dwarf & 1 ) == currentBit ) {
      dwarf = dwarf >> 1;
    } else {
      dwarf = ( dwarf >> 1 ) ^ crc32Poly;
    }
  }
  return ( dwarf );
}


void StopWatch::start() {
  delta = 0.0;
  getWallClockTime( started, _startedMiliseconds );
  cont();
}

void StopWatch::cont() {
  t1 = getRuntimeSoFar();
}

void StopWatch::stop() {
  t2 = getRuntimeSoFar();
  getWallClockTime( finished, _finishedMiliseconds );
  delta += ( t2 - t1 );
}

double StopWatch::total_runtime() {
  return ( getRuntimeSoFar() );
}

const char* StopWatch::startedAt() {
  strftime ( _startedString, 512, "%H:%M on %A, %d %B %Y", localtime ( &started ) );
  return ( _startedString );
}

const char* StopWatch::stoppedAt() {
  strftime ( _stoppedString, 512, "%H:%M on %A, %d %B %Y", localtime ( &finished ) );
  return ( _stoppedString );
}

bool StopWatch::_suppressOpenMPWarning = false;


void printSelfTestSuccessMessage ( const string &message ) {
  cerr << aol::color::ok << endl
       << "--------------------------------------------------------------------------------" << endl
       << message << endl
       << "--------------------------------------------------------------------------------" << endl
       << aol::color::reset;
}

//! for consistent output of failed selfTests (message should be 80 characters wide)
void printSelfTestFailureMessage ( const string &message ) {
  cerr << aol::color::error << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << message << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << aol::color::reset;
}

bool checkForBenchmarkArguments ( const int argn, char **argv, string &ResultFilename ) {
  if ( ( argn == 4 ) && ( string ( argv [1] ) == string ( "bench" ) ) && ( string (argv [2]) == string ( "file" ) ) ) {
    ResultFilename = argv [3];
    return true;
  }
  else
    return false;
}

//! read and ignore comments in saved files
void READ_COMMENTS ( istream &in ) {
  if ( in.fail() )
    throw aol::Exception ( "aol::READ_COMMENTS: called on stream with in.fail()", __FILE__, __LINE__ );

  do {
    char _tmp={' '}, _temp[256]={' '};
    int  _done=0;
    while (_tmp == ' ' || _tmp == '\n' || _tmp == '\t' || _tmp == '\r')
      _tmp = static_cast<char> (in.get());
    in.unget();
    _done=0;
    while ( !_done ) {
      _tmp = static_cast<char> ( in.get() );
      if ( _tmp=='#' ) {
        in.getline(_temp, 255);
      } else {
        in.unget();
        _done=1;
      }
    }
  } while( 0 );

  if ( in.fail() )
    throw aol::Exception ( "aol::READ_COMMENTS: error reading comment", __FILE__, __LINE__ );
}

void copyFile ( const char* source, const char* dest ) {
  ifstream in ( source, ios::binary );
  ofstream out ( dest, ios::binary );

  out << in.rdbuf ();

  out.close ();
  in.close ();
}

void printHGChangesetID ( ostream& os ) {
#ifdef HG_CHANGESET_ID
  os << "HG changeset: " << HG_CHANGESET_ID << std::endl;
#else
  os << "Could not determine HG changeset ID. Is USE_MERCURIAL enabled?" << std::endl;
#endif
}



template<> const aol::FileFormatType aol::FileFormatMagicNumber< bool           >::FFType = aol::FF_BOOL;
template<> const aol::FileFormatType aol::FileFormatMagicNumber< signed char    >::FFType = aol::FF_SIGNED_CHAR;
template<> const aol::FileFormatType aol::FileFormatMagicNumber< unsigned char  >::FFType = aol::FF_UNSIGNED_CHAR;
template<> const aol::FileFormatType aol::FileFormatMagicNumber< signed short   >::FFType = aol::FF_SIGNED_SHORT;
template<> const aol::FileFormatType aol::FileFormatMagicNumber< unsigned short >::FFType = aol::FF_UNSIGNED_SHORT;
template<> const aol::FileFormatType aol::FileFormatMagicNumber< signed int     >::FFType = aol::FF_SIGNED_INT;
template<> const aol::FileFormatType aol::FileFormatMagicNumber< unsigned int   >::FFType = aol::FF_UNSIGNED_INT;
template<> const aol::FileFormatType aol::FileFormatMagicNumber< int64_t        >::FFType = aol::FF_SIGNED_64INT;
template<> const aol::FileFormatType aol::FileFormatMagicNumber< uint64_t       >::FFType = aol::FF_UNSIGNED_64INT;
#if ( defined __APPLE__ ) || ( BITNESS == 32 ) // This seems to be the only way to get SuiteSparse_long to work under OS X
template<> const aol::FileFormatType aol::FileFormatMagicNumber< long           >::FFType = aol::FF_SIGNED_64INT;
#endif
template<> const aol::FileFormatType aol::FileFormatMagicNumber< float          >::FFType = aol::FF_FLOAT;
template<> const aol::FileFormatType aol::FileFormatMagicNumber< double         >::FFType = aol::FF_DOUBLE;
template<> const aol::FileFormatType aol::FileFormatMagicNumber< long double    >::FFType = aol::FF_LONG_DOUBLE;


template<> const string aol::FileFormatMagicNumber< bool           >::FFContainedName = "bool";
template<> const string aol::FileFormatMagicNumber< signed char    >::FFContainedName = "signed char";
template<> const string aol::FileFormatMagicNumber< unsigned char  >::FFContainedName = "unsigned char";
template<> const string aol::FileFormatMagicNumber< signed short   >::FFContainedName = "signed short";
template<> const string aol::FileFormatMagicNumber< unsigned short >::FFContainedName = "unsigned short";
template<> const string aol::FileFormatMagicNumber< signed int     >::FFContainedName = "signed int";
template<> const string aol::FileFormatMagicNumber< unsigned int   >::FFContainedName = "unsigned int";
template<> const string aol::FileFormatMagicNumber< int64_t        >::FFContainedName = "int64_t";
template<> const string aol::FileFormatMagicNumber< uint64_t       >::FFContainedName = "uint64_t";
#if ( defined __APPLE__ ) || ( BITNESS == 32 ) // This seems to be the only way to get SuiteSparse_long to work under OS X
template<> const string aol::FileFormatMagicNumber< long           >::FFContainedName = "long";
#endif
template<> const string aol::FileFormatMagicNumber< float          >::FFContainedName = "float";
template<> const string aol::FileFormatMagicNumber< double         >::FFContainedName = "double";
template<> const string aol::FileFormatMagicNumber< long double    >::FFContainedName = "long double";

} // namespace aol
