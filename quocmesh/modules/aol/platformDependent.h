#ifndef __PLATFORMDEPENDENT_H
#define __PLATFORMDEPENDENT_H

/**
 * \file
 *
 * This file is intended to collect small pieces of platform dependent
 * code. Long platform dependent code pieces don't need to be put in here,
 * e.g. aol::HashSet.
 */

#include <vector>
// Necessary for std::ifstream/ofstream (aol.h may not be included!).
#include <fstream>
#include <qmException.h>
// do not include <aol.h> to avoid problems with circular dependency.

#ifdef __GNUC__
// With this define, GCC_VERSION is 40301 for GCC 4.3.1 and can be conveniently checked (e.g. __WARNING_ON)
#define GCC_VERSION ( __GNUC__ * 10000 \
                      + __GNUC_MINOR__ * 100 \
                      + __GNUC_PATCHLEVEL__ )
#if (__GNUC__ >= 3) && !defined(_LIBCPP_VERSION)
#define HAVE_STDIO_FILEBUF 1
#endif
#endif

#ifdef HAVE_STDIO_FILEBUF
#include <ext/stdio_filebuf.h>
using  namespace __gnu_cxx;
#endif

#if !defined ( _WIN32 ) && !defined ( _WIN64 )
#include<sys/times.h>
#include<sys/resource.h>
#else
#include<time.h>
#endif

#ifndef _MSC_VER
#include <stdint.h>
#endif

#ifdef __MINGW32_VERSION
// Necessary for va_list (see aol::vscprintf).
#include <cstdarg>
#endif

// Mechanism to disable and enable warnings
// Usage example: WARNING_OFF ( old-style-cast ) and WARNING_ON ( old-style-cast )
//
// Convention: Disable only warnings that are enabled in the QuocMeshes
//  by default, and re-enable them afterwards
//
// Note: Support for diagnostic pragmas was added in GCC 4.2.0.
#if ( defined ( __GNUC__ ) ) && ( GCC_VERSION >= 40200 ) && !( defined ( __INTEL_COMPILER ) ) && !( defined ( __NVCC__ ) )
#define __PRAGMA(P) _Pragma(#P)
#define __WARNING_OFF(WARN) __PRAGMA(GCC diagnostic ignored #WARN)
#define __WARNING_ON(WARN) __PRAGMA(GCC diagnostic warning #WARN)
#define WARNING_OFF(WARN) __WARNING_OFF(-W ## WARN)
#define WARNING_ON(WARN) __WARNING_ON(-W ## WARN)
#else
#define WARNING_OFF(WARN)
#define WARNING_ON(WARN)
#endif

/*
 * windows specific settings
 */

#ifdef __CYGWIN__

#define rint(a)   double(int(a+0.5)); //round to nearest integer

#ifdef _B
#undef _B
#endif
#ifdef _N
#undef _N
#endif
#ifdef _U
#undef _U
#endif

#endif // __CYGWIN__

#if defined(_MSC_VER)
// Warning 4661 is disabled, because the implementation of the class SurfMesh is split over several cpp files.
#pragma warning ( disable : 4661 )
// Warning 4355 is disabled, because the implementation of the class aol::DofMask uses 'this' in the base member initializer list
#pragma warning ( disable : 4355 )

// This allows to include hash_set with recent VC++ versions.
#if ( _MSC_VER >= 1900 )
#define _SILENCE_STDEXT_HASH_DEPRECATION_WARNINGS 1
#endif

#define snprintf _snprintf
#define strcasecmp _stricmp

#include <float.h>
#include <limits>

typedef signed __int8 int8_t;
typedef unsigned __int8 uint8_t;
typedef signed __int16 int16_t;
typedef unsigned __int16 uint16_t;
typedef signed __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef signed __int64 int64_t;
typedef unsigned __int64 uint64_t;

#pragma warning( disable : 4305 )
#define popen(s1,s2) _popen(s1,s2)
#define pclose(pf)  _pclose(pf)

#define rint(a)   double(int(a+0.5)); //round to nearest integer
#define copysign(a,b) _copysign(a,b)

inline double cbrt ( double x ) {
  if ( fabs(x) < std::numeric_limits<double>::epsilon() )
    return 0.;

  // pow is only guaranteed for work positive x.
  if (x > 0.0)
    return pow(x, 1./3.);
  else
    return -pow(-x, 1./3.);
}

#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2      1.570796326794896619
#endif

#ifndef M_LN2
#define M_LN2       0.69314718055994530942
#endif

// check if stlport STL is used:
#ifdef _STLP_CALL

// The compiler has problems with overloaded pow-function, if stlport's STL is used,
// so we need the following two template functions.

//! pow function for input arguments of same type
template<class T> T pow ( T a, T b ) {
  return ( T ) pow ( double ( a ), double ( b ) );
}

//! pow function for input arguments of different type
template<class T1, class T2> double pow ( T1 a, T2 b ) {
  return pow ( double ( a ), double ( b ) );
}

#endif //_STLP_CALL

#endif // _MSC_VER

namespace aol{

//! Function for platform indepedent aligned memory allocation
//! \author Berkels
void* aligned_memory_allocation ( const size_t MemorySize, const size_t Alignment );

//! Function to deallocate the memory allocated by aligned_memory_allocation
//! \author Berkels
void aligned_memory_deallocation ( void* Pointer );

void makeDirectory ( const char *DirectoryName, bool verbose = true );

void recursivelyMakeDirectory ( const char *DirectoryName, bool verbose = false );

//! Checks, if directory exists
//! \author Toelkes
bool directoryExists ( const std::string& directoryName );

//! Generates a unique directory name using dirBaseName as a base
//! \author Toelkes
bool generateUniqueDirectoryName ( std::string &dirBaseName );

//! Creates and opens a temporary file. TempFileName has to contain a template filename ending with `XXXXXX'.
//! (This does not apply to VC++, here the original contents of TempFileName are ignored.
//! The temporary filename is written into TempFileName.
//! The file has to be deleted manually, after TempFileStream is closed.
void generateTemporaryFile ( char *TempFileName, std::ofstream &TempFileStream );

//! Returns the current working directory with slashes as seperators between subdirectories
void getCurrentDirectoryUnixStyle ( char *CurrentDirectoryName, const int MaxNameLength = 128);

void setCurrentDirectory ( std::string dir );

//! Creates a vector with the names of all files (and possibly subdirectories) of the directory \arg Dir.
//! \author Berkels
void createDirectoryListing ( const char *Dir, std::vector<std::string> &DirList, const bool IncludeDirectories = false );

//! Test for file existence
bool fileExists ( std::string filename );

//! Returns the size of the file and -1 if the file does not exist.
//! \author Berkels
int getSizeOfFile ( const std::string &Filename );

//! Appends the specified directory name to the PATH environment variable.
//! \author Berkels
void appendSearchPath ( const char *DirectoryName );

//! Calls system( "PAUSE" ), if called on a platform where the output window is automatically closed on exit (MinGW for example)
void callSystemPauseIfNecessaryOnPlatform ();

//! Writes the current wall clock time in seconds to Seconds and the miliseconds part to Miliseconds.
void getWallClockTime( time_t &Seconds, unsigned short &Miliseconds );

//! Returns the current runtime in seconds as double.
double getRuntimeSoFar();

//! Find out whether we are currently using elecric fence malloc debugger
bool areWeUsingEFence ();

//! Choose which memory is counted by memusage
enum MemoryUsageFlags { VIRTUAL_MEMORY = 1,                 //!< Size of virtual memory space of this process, including resident, swap and pages not currently mapped to either
                        RESIDENT_MEMORY = 2,                //!< Size of memory pages in RAM
                        SWAP_MEMORY = 4,                    //!< Size of memory pages on swap devices
                        DATA_MEMORY = 8,                    //!< Size of data segment
                        STACK_MEMORY = 16,                  //!< Size of program stack
                        MAX_VIRTUAL_MEMORY = 32,            //!< Maximum of virtual memory
                        MAX_RESIDENT_MEMORY = 64,           //!< Maximum of resident memory
                        MEMORY_MANAGER_MEMORY = 128,        //!< Memory used by vector manager for buffering deleted Vectors
                        MALLOCED_MEMORY = 256,              //!< Memory allocated by malloc()
                        MMAPPED_MEMORY = 512                //!< Memory allocated by mmap()
};

//! Returns number of bytes used by program, by default the size of the heap is returned
//! Works only on Linux systems. MALLOCED_MEMORY + MMAPPED_MEMORY gives the most exact figure, but flows over at 4 GiB.
int64_t memusage ( const MemoryUsageFlags flags = static_cast<MemoryUsageFlags> ( RESIDENT_MEMORY + SWAP_MEMORY ) );

//! Returns the number of characters (not including the trailing '\0') which would
//! have been written to the final string by vsprintf.
//! \author Berkels
int vscprintf(const char *format, va_list ap);

//! pfstream is a GNU extension and therefore unknown to the intel compiler
//! We supply a substitution here
//! @ingroup Input
class ipfstream :
#ifdef __GNUC__
#if __GNUC__ >= 3
      public std::istream {
#else
      public std::ifstream {
#endif
#else // OTHER compiler
      public std::ifstream {
#endif

private:
  bool isPipe;
  FILE *f;
#ifdef _WIN32
  bool isTemp;
  char* TempFileName;
#endif

  FILE *openFile ( const char *FileName );
#ifdef HAVE_STDIO_FILEBUF
  stdio_filebuf<char>* buf;
#endif

public:
  ipfstream ( const char *name );

  ~ipfstream() {
    close();
  }

  void close() {
#ifdef HAVE_STDIO_FILEBUF
    // Detach before closing file/pipe
    rdbuf ( NULL );
    if ( buf ) {
      delete buf;
      buf = NULL;
    }
#else
#if defined(__GNUC__) && !(__GNUC__ >= 3)
    ifstream::close();
#endif
    if ( f ) {
#ifdef _WIN32
      fclose ( f );
#else
      if ( isPipe ) pclose ( f );
      else fclose ( f );
#endif
      f = NULL;
#ifdef _WIN32
      if ( isTemp ) {
        remove ( TempFileName );
        delete[] TempFileName;
      }
#endif // _WIN32
    }
#endif
  }
};

//! pfstream is a GNU extension and therefore unknown to the intel compiler
//! We supply a substitution here
//! @ingroup Output
class opfstream :
#ifdef __GNUC__
#if __GNUC__ >= 3
      public std::ostream {
#else
      public std::ofstream {
#endif
#else // OTHER compiler
      public std::ofstream {
#endif
private:
#if defined(__MINGW32_VERSION) || defined(__CYGWIN__) // new gcc under Windows
  bool isTemp;
  char* TempFileName;
#endif
  bool isPipe;
  FILE *f;

#ifdef HAVE_STDIO_FILEBUF
  stdio_filebuf<char>* buf;
#endif

public:
#ifdef __GNUC__
#if __GNUC__ < 3  // old gcc
  opfstream ( const char *name ) : std::ofstream() {
    if ( ( isPipe = ( name[0] == '|' ) ) ) f = popen ( name + 1, "w" );
    else f = fopen ( name, "w" );
    if ( f ) attach ( fileno ( f ) );
  }
#else // new gcc
#if defined(__MINGW32_VERSION) || defined(__CYGWIN__) // new gcc under Windows
  opfstream ( const char *name ) :
      std::ostream ( NULL ), isTemp ( name[0] == '|' ), isPipe ( false )
#else
  opfstream ( const char *name ) :
      std::ostream ( NULL ), isPipe ( name[0] == '|' ), f ( isPipe ? popen ( name + 1, "w" ) : fopen ( name, "w" ) )
#ifdef HAVE_STDIO_FILEBUF
          , buf ( f ? new stdio_filebuf<char> ( f, std::ios::out ) : NULL )
#endif
#endif
  {
#if defined(__MINGW32_VERSION) || defined(__CYGWIN__) // new gcc under Windows
    if ( isTemp ) {
      TempFileName = new char[256];
      sprintf ( TempFileName, "%s.pgm", name + 12 );
    }

    f = ( isTemp ? fopen ( TempFileName, "wb" ) : fopen ( name, "wb" ) );
    buf = ( f ? new stdio_filebuf<char> ( f, std::ios::out ) : NULL );
#endif
    if ( !f ) throw Exception ( std::string ( "Could not open file " ) + name, __FILE__, __LINE__ );
#ifdef HAVE_STDIO_FILEBUF
    rdbuf ( buf );
#endif
  }
#endif

#else // other compiler
#ifdef _WIN32
  opfstream ( const char *name ) : std::ofstream() {
    isPipe = ( name[0] == '|' );
    f = ( isPipe ? popen ( name + 1, "w" ) : fopen ( name, "w" ) );
// Hack the make this class somewhat work with VC++ 2005 SP1
// Will be removed once Bzipofstream is implemented.
#if defined(_MSC_VER)
    open ( name, std::ios::binary );
#endif
  }
#else
  opfstream ( const char *name ) : std::ofstream ( ( isPipe = ( name[0] == '|' ) ) ? ( f = popen ( name + 1, "w" ) ) : ( f = fopen ( name, "w" ) ) ) {}
#endif //_WIN32
#endif //compiler

  ~opfstream();

  void close() {
#ifdef __GNUC__
#if __GNUC__ < 3
    std::ofstream::close();
#endif
#endif
  }
};

bool isErfTestOK ( );

//! table of ansi color codes
namespace color {
// example: cerr << color::error << "error message noticed by everybody :-)" << color::reset;
extern const std::string reset;
extern const std::string invert;
extern const std::string black;
extern const std::string red;
extern const std::string green;
extern const std::string brown;  // aka yellow
extern const std::string blue;
extern const std::string purple; // aka magenta
extern const std::string cyan;
extern const std::string light_grey;
extern const std::string dark_grey;
extern const std::string light_red;
extern const std::string light_green;
extern const std::string yellow; // aka light yellow
extern const std::string light_blue;
extern const std::string pink;   // aka light magenta
extern const std::string light_cyan;
extern const std::string white;
extern const std::string beep;
extern const std::string error;
extern const std::string ok;
extern const std::string residuum;
}

} // namespace aol

#endif
