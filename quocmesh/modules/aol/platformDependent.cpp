#ifdef USE_LIB_BZ2
#if (defined _WIN32) || defined(__CYGWIN__)
#include <bzlib.h>
// Unfortunately bzlib.h defines a lot of crap under some platforms.
// We need to kill the defines here.
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#endif
#endif

#ifdef GNU
#include <sys/types.h>
#include <sys/stat.h>
#if !defined ( _WIN32 ) && !defined ( _WIN64 )
#include <sys/user.h>
#else
#ifndef __MINGW64__
#include <ddk/ntddk.h>
#endif
#endif
#endif

#ifdef __APPLE__
#include <mach/mach.h>
#include <unistd.h>
#endif // MACOSX

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <aol.h>
#include <vec.h>

#if defined(__MINGW32_VERSION) || defined(__MINGW64__) || defined(__CYGWIN__)
#include <io.h>
#include <fcntl.h>
#endif
#include <sys/stat.h>
#if defined(_MSC_VER)
#include <sys/types.h>
#include <sys/stat.h>
#include <direct.h>
#include <windows.h>
#endif
#include <sys/timeb.h>

#ifdef _WIN32
#include <psapi.h>
#endif

#include <dirent.h>

namespace aol{

#if !defined ( _WIN32 ) && !defined ( _WIN64 )
#include <dlfcn.h>
#endif

#ifndef __APPLE__
#include <malloc.h>
#else
#include <sys/malloc.h>
#endif

#ifdef USE_SSE
#include <xmmintrin.h>
void* aligned_memory_allocation ( const size_t MemorySize, const size_t Alignment ) {
#if defined (__MINGW32_VERSION) || defined(__MINGW64__)
  return __mingw_aligned_malloc ( MemorySize, Alignment );
#else
#ifdef __CYGWIN__
  return memalign ( Alignment, MemorySize );
#else
  void* p = NULL;
  if ( posix_memalign ( &p, Alignment, MemorySize ) != 0 )
    throw aol::OutOfMemoryException ( "aol::aligned_memory_allocation: Allocate not successful.", __FILE__, __LINE__ );

  return p;
#endif // __CYGWIN__
#endif // __MINGW32_VERSION
}

void aligned_memory_deallocation ( void* Pointer ) {
#if defined (__MINGW32_VERSION) || defined(__MINGW64__)
  __mingw_aligned_free ( Pointer );
#else
  free ( Pointer );
#endif // __MINGW32_VERSION
}
#endif

#if defined(_MSC_VER)
// Visual C++ uses UNICODE and therefore wide chars instead of chars for certain functions.
// In the following are some functions to conveniently convert between wchar and char (and vice versa).
// wide_to_narrow if taken from http://www.codeguru.com/forum/archive/index.php/t-336106.html
CHAR wide_to_narrow(WCHAR w)
{
  // simple typecast
  // works because UNICODE incorporates ASCII into itself
  return CHAR(w);
}

WCHAR narrow_to_wide(CHAR w)
{
  return WCHAR(w);
}

void convertCharToWChar( const CHAR *Arg, WCHAR *Dest ){
  std::transform(Arg, Arg + strlen( Arg ) + 1, Dest, narrow_to_wide);
}

void convertWCharToChar( const WCHAR *Arg, CHAR *Dest ){
  std::transform(Arg, Arg + wcslen( Arg ) + 1 , Dest, wide_to_narrow);
}

wstring wstringFromCharPointer( const CHAR *Arg ){
  WCHAR dest[1024];
  convertCharToWChar( Arg, dest );
  return wstring( dest );
}

#endif

void convertFilenameToUnixStyle ( char *CurrentDirectoryName ) {
  int i = 0;
  while( CurrentDirectoryName[i] != '\0' ){
    if( CurrentDirectoryName[i] == '\\' )
      CurrentDirectoryName[i] = '/';
    i++;
  }
}

void makeDirectory ( const char *DirectoryName, bool verbose ) {
#if defined(_WIN32) && !defined(__MINGW32_VERSION) && !defined(__CYGWIN__)
  string systemCommand = "md \"";
  systemCommand += DirectoryName;
  systemCommand += "\"";
  // VC++ demands subdirectories to be separated by '\\' instead of '/'.
  std::replace( systemCommand.begin(), systemCommand.end(), '/', '\\' );
  system ( systemCommand.c_str() );
  if ( verbose )
    cerr << "Created directory " << DirectoryName << endl;
#else
  struct stat directory;
  const int statReturn = stat ( DirectoryName, &directory );

  if ( ( statReturn == - 1 ) || !S_ISDIR ( directory.st_mode ) ) {
#if defined(__MINGW32_VERSION) || defined(__MINGW64__)
    mkdir ( DirectoryName );
#else
    string systemCommand = "mkdir \"";
    systemCommand += DirectoryName;
    systemCommand += "\"";
    if ( system ( systemCommand.c_str() ) != EXIT_SUCCESS )
      cerr << "aol::makeDirectory: Calling '" << systemCommand << "' returned an error.\n";
#endif
    if ( verbose )
      cerr << "Created directory " << DirectoryName << endl;
  } else {
    if ( verbose )
      cerr << "Directory " << DirectoryName << " already exists\n";
  }
#endif
}

  void recursivelyMakeDirectory ( const char *DirectoryName, bool verbose ) {
#if defined(_WIN32) && !defined(__MINGW32_VERSION) && !defined(__MINGW64__) && !defined(__CYGWIN__)
  std::string dir = DirectoryName;
  // Use the Windows path separator
  std::replace( dir.begin(), dir.end(), '/', '\\' );
  if ( dir.size() == 0 )
    return;

  // Make sure the path ends with a separator
  if ( dir[dir.size()-1] != '\\' )
    dir += '\\';

  char folder[MAX_PATH];
  const char *end;
  ZeroMemory(folder, MAX_PATH * sizeof(char));
  end = strchr ( dir.c_str(), '\\');

  while ( end != NULL ) {
    strncpy ( folder, dir.c_str(), end - dir.c_str() + 1 );
    if ( !CreateDirectoryA(folder, NULL ))
    {
	  if ( ( GetLastError() == ERROR_ALREADY_EXISTS ) && verbose )
        cerr << folder << " already exists" << endl;
	}
	else if ( verbose )
		cerr << folder << "created" << endl;
    end = strchr(++end, '\\');
  }
  return;
#else
  struct stat directory;
  const int statReturn = stat ( DirectoryName, &directory );

  if ( ( statReturn == - 1 ) || !S_ISDIR ( directory.st_mode ) ) {
#if defined(__MINGW32_VERSION) || defined(__MINGW64__)
  throw aol::Exception ( "Not implemented for mingw", __FILE__, __LINE__, __FUNCTION__ );
#else
    string systemCommand = "mkdir -p \"";
    systemCommand += DirectoryName;
    systemCommand += "\"";
    if ( system ( systemCommand.c_str() ) != EXIT_SUCCESS )
      cerr << "aol::makeDirectory: Calling mkdir returned an error.\n";
#endif
    if ( verbose )
      cerr << "Created directory " << DirectoryName << endl;
  } else {
    if ( verbose )
      cerr << "Directory " << DirectoryName << " already exists\n";
  }
#endif
}

  bool directoryExists ( const std::string& directoryName ) {
    struct stat directory;
    const int statReturn = stat ( directoryName.c_str (), &directory );

    return !( ( statReturn == -1 ) || !S_ISDIR ( directory.st_mode ) );
  }

  bool generateUniqueDirectoryName ( std::string &dirBaseName ) {
    std::string uniqueDirName = dirBaseName;
    unsigned int counter = 1;

    while ( directoryExists ( uniqueDirName ) ) {
      uniqueDirName = ( dirBaseName + '_' + aol::to_string ( counter++ ) );
    }

    if ( dirBaseName != uniqueDirName ) {
      dirBaseName = uniqueDirName;
      return true;
    }

    return false;
  }

void generateTemporaryFile ( char *TempFileName, ofstream &TempFileStream ) {
#if defined(_MSC_VER)
  // based on http://msdn2.microsoft.com/en-us/library/aa363875.aspx
  DWORD dwRetVal;
  DWORD dwBufSize=1024;
  UINT uRetVal;
  WCHAR szTempName[1024];
  WCHAR lpPathBuffer[1024];

  // Get the temp path.
  dwRetVal = GetTempPath(dwBufSize,     // length of the buffer
                         lpPathBuffer); // buffer for path
  if (dwRetVal > dwBufSize || (dwRetVal == 0))
  {
      printf ("GetTempPath failed with error %lu.\n", GetLastError());
      throw Exception ( "Could not generate a temporary file.", __FILE__, __LINE__ );
  }

  // Create a temporary file.
  uRetVal = GetTempFileName( lpPathBuffer,                          // directory for tmp files
                             wstringFromCharPointer("NEW").c_str(), // temp file name prefix
                             0,                                     // create unique name
                             szTempName);                           // buffer for name

  if (uRetVal == 0)
  {
    printf ("GetTempFileName failed with error %lu.\n", GetLastError());
    throw Exception ( "Could not generate a temporary file.", __FILE__, __LINE__ );
  }
  convertWCharToChar( szTempName, TempFileName );
  convertFilenameToUnixStyle( TempFileName );
  TempFileStream.open ( TempFileName );
#else
  if ( TempFileName[0] == 0 )
    throw Exception ( "Can't generate a temporary file with empty file name mask.", __FILE__, __LINE__ );
  int fileDescriptor = 0;
#if defined(__MINGW32_VERSION) || defined(__MINGW64__)
  if ( mktemp ( TempFileName ) ) {
    do
      fileDescriptor = open ( TempFileName, O_CREAT | O_EXCL, S_IREAD | S_IWRITE );
    while ( ! ( fileDescriptor == -1 && errno == EEXIST ) && mktemp ( TempFileName ) );
  } else
    fileDescriptor = -1;
#else
  fileDescriptor = mkstemp ( TempFileName );
#endif
  if ( fileDescriptor == -1 )
    throw Exception ( "Could not generate a temporary file.", __FILE__, __LINE__ );

  TempFileStream.open ( TempFileName );
  close ( fileDescriptor );
#endif
}

void getCurrentDirectoryUnixStyle ( char *CurrentDirectoryName, const int MaxNameLength ) {
  if (
#if defined(_MSC_VER)
    _getcwd
#else
    getcwd
#endif
      (CurrentDirectoryName, MaxNameLength) == NULL )
    cerr << "aol::getCurrentDirectoryUnixStyle: Failed to get the directory.\n";
#if defined(_WIN32)
  convertFilenameToUnixStyle( CurrentDirectoryName );
#endif
}

void setCurrentDirectory ( string dir ) {
#if (defined(_MSC_VER))
  _chdir ( dir.c_str() );
#elif ( defined ( GNU ) ) || ( defined (__clang__) )
  if ( chdir ( dir.c_str() ) != 0 )
    cerr << "aol::setCurrentDirectory: Failed to set the directory.\n";
#else
  throw UnimplementedCodeException("setCurrentDirectory", __FILE__, __LINE__);
#endif
}


#if !defined ( __MINGW32_VERSION ) && !defined (__MINGW64__)
/**
 * Helper function for aol::createDirectoryListing.
 *
 * \author Berkels
 */
bool isDirEntDirectory ( dirent *Entry ) {
#if defined (  _MSC_VER )
  if ( Entry->d_type & DT_DIR )
#else
  if ( Entry->d_type == 4 )
#endif
    return true;

  return false;
}
#endif

void createDirectoryListing ( const char *Dir, std::vector<std::string> &DirList, const bool IncludeDirectories ) {
  DIR *dirp;
  struct dirent *entry;

  if ( ( dirp = opendir( Dir ) ) != NULL ) {
    while ( ( entry = readdir ( dirp ) ) != NULL ) {
      if ( IncludeDirectories ||
#if defined ( __MINGW32_VERSION ) || defined (__MINGW64__)
      ( fileExists ( ( string ( Dir ) + entry->d_name ) ) == true )
#else
      ( isDirEntDirectory ( entry ) == false )
#endif
       )
        DirList.push_back ( entry->d_name );
    }
    std::sort ( DirList.begin(), DirList.end() );
    closedir ( dirp );
  }
}

bool fileExists ( string filename ) {
#if defined ( GNU ) || defined ( __clang__ )
  struct stat buf;
  return !stat ( filename.c_str (), &buf ) && S_ISREG ( buf.st_mode );
#else
#if (defined(_MSC_VER))
  struct _stat buf;
  return !_stat( filename.c_str (), &buf ) && (buf.st_mode & _S_IFREG);
#else
  throw UnimplementedCodeException("fileExists", __FILE__, __LINE__);
#endif
#endif
}

int getSizeOfFile ( const std::string &Filename ) {
  if ( aol::fileExists ( Filename ) ) {
    struct stat filestatus;
    stat( Filename.c_str(), &filestatus );
    return filestatus.st_size;
  }
  else
    return -1;
}

void appendSearchPath ( const char *DirectoryName ) {
  std::string path;
  const char *pathVar = getenv ( "PATH" );
  if ( pathVar ) {
    path = pathVar;
    path += ":";
  }
  path += DirectoryName;
#if defined(__MINGW32_VERSION) || defined(__MINGW64__) || defined(_MSC_VER)
  throw UnimplementedCodeException ( "appendSearchPath not implemented under VC++", __FILE__, __LINE__ );
#else
  setenv ( "PATH", path.c_str(), 1 );
#endif
}

void callSystemPauseIfNecessaryOnPlatform (){
#if defined(__MINGW32_VERSION) || defined(__MINGW64__) || defined(_MSC_VER)
  const char *var = getenv ( "QUOC_NO_SYSTEM_PAUSE" );
  if ( ( var == NULL ) || ( strcmp ( var, "ON" ) ) ) {
    cerr << "Press enter to continue. . .\n";
    cin.get();
  }
#endif
}

void getWallClockTime( time_t &Seconds, unsigned short &Miliseconds ){
#if (defined(_WIN32) && !defined(__CYGWIN__))
  _timeb t;
  _ftime( &t );
#else
  timeb t;
  ftime( &t );
#endif
  Seconds = t.time;
  Miliseconds = t.millitm;
}

// MinGW defines CLOCKS_PER_SEC with an old style.
#ifdef _WIN32
WARNING_OFF(old-style-cast)
#endif
double getRuntimeSoFar() {
  double time_sec = 0.0;

#if defined ( _WIN32 ) || defined ( _WIN64 )
  time_sec = static_cast<double> ( clock() ) / CLOCKS_PER_SEC;
#else
  struct rusage ru_self, ru_children;
  getrusage ( RUSAGE_SELF, &ru_self );
  getrusage ( RUSAGE_CHILDREN, &ru_children );

  // the following is probably the right thing to do: sum up user and system time of both program and child processes called
  time_sec += 1.0 * ( ru_self.ru_utime.tv_sec ) +  1.0e-6 * ( ru_self.ru_utime.tv_usec );          // self user time
  time_sec += 1.0 * ( ru_self.ru_stime.tv_sec ) +  1.0e-6 * ( ru_self.ru_stime.tv_usec );          // self system time
  time_sec += 1.0 * ( ru_children.ru_utime.tv_sec ) +  1.0e-6 * ( ru_children.ru_utime.tv_usec );  // children user time
  time_sec += 1.0 * ( ru_children.ru_stime.tv_sec ) +  1.0e-6 * ( ru_children.ru_stime.tv_usec );  // children system time
#endif

  return ( time_sec );
}
#ifdef _WIN32
WARNING_ON(old-style-cast)
#endif

bool areWeUsingEFence () {
#if ( ( defined(GNU) || defined(__clang__) ) && !defined(_WIN32) && !defined(_WIN64))
  return dlsym (dlopen (NULL, RTLD_LAZY), "EF_ALLOW_MALLOC_0") != NULL;
#else
  throw UnimplementedCodeException ( "areWeUsingEFence", __FILE__, __LINE__ );
#endif
}

int64_t memusage ( const MemoryUsageFlags flags ) {
#if ( ( defined(GNU) || defined ( __INTEL_COMPILER ) || defined ( __clang__ ) ) && !defined(_WIN32) && !defined(_WIN64) && !defined(__APPLE__))
  int64_t mem = 0;
  char cline [1024];

  std::ifstream procfs ( "/proc/self/status" );

  while (procfs.getline (cline, 1000)) {
    std::string line (cline);

    if (flags &      VIRTUAL_MEMORY) { if (line.find ("VmSize:") == 0) { int p = line.find ("kB"); mem += convert<int> (line.substr (7,p-7)); } }
    if (flags &     RESIDENT_MEMORY) { if (line.find ("VmRSS:")  == 0) { int p = line.find ("kB"); mem += convert<int> (line.substr (6,p-6)); } }
    if (flags &         SWAP_MEMORY) { if (line.find ("VmSwap:") == 0) { int p = line.find ("kB"); mem += convert<int> (line.substr (7,p-7)); } }
    if (flags &         DATA_MEMORY) { if (line.find ("VmData:") == 0) { int p = line.find ("kB"); mem += convert<int> (line.substr (7,p-7)); } }
    if (flags &        STACK_MEMORY) { if (line.find ("VmStk:")  == 0) { int p = line.find ("kB"); mem += convert<int> (line.substr (6,p-6)); } }
    if (flags &  MAX_VIRTUAL_MEMORY) { if (line.find ("VmPeak:") == 0) { int p = line.find ("kB"); mem += convert<int> (line.substr (7,p-7)); } }
    if (flags & MAX_RESIDENT_MEMORY) { if (line.find ("VmHWM:")  == 0) { int p = line.find ("kB"); mem += convert<int> (line.substr (6,p-6)); } }

  }

  mem *= 1024;

  if (flags & MEMORY_MANAGER_MEMORY) {
#ifndef DO_NOT_USE_MEMORYMANAGER
    mem += aol::MemoryManager::memoryManagerMemoryUsage ();
#endif
  }

  char buffer[4001] = {0};

  {
    // malloc_stats prints to stderr, redirect this to a pipe

    // mallinfo uses unsigned ints which can overflow on 64bit systems
    // malloc_info is 64bit safe, but does not count mmapped memory
    // malloc_get_state is utterly system dependent
    // malloc_stats works internally also with unsigned ints, but (in contrast to mallinfo) this could theoretically be improved
    int out_pipe[2];
    int saved_stderr;
    saved_stderr = dup (STDERR_FILENO);
    if ( pipe (out_pipe) != 0 )
      cerr << "aol::memusage: Error opening pipe.\n";
    dup2 (out_pipe[1], STDERR_FILENO);
    close (out_pipe[1]);
    malloc_stats ();
    // Standard error seems to be unavailable here, so use cout instead of cerr.
    if ( read (out_pipe[0], buffer, 4000) <= 0 )
      cout << "aol::memusage: Error reading pipe output.\n";
    dup2 (saved_stderr, STDERR_FILENO);
  }

  std::string buf ( buffer );
  std::string::size_type p = 0, pe;
  int64_t mema = 0;
  // Attention: With openmp there may be multiple arenas!
  while ( (p = buf.find ("Arena", p)) != string::npos) {
    p = buf.find("in use bytes", p);
    p = buf.find ("=", p);
    pe = buf.find ("\n", p);
    mema += convert<int64_t> (buf.substr (p+1, pe-p-1));
  }
  p = buf.find ("Total (incl. mmap):");
  p = buf.find("in use bytes", p);
  p = buf.find ("=", p);
  pe = buf.find ("\n", p);
  int64_t memtot = convert<int64_t> (buf.substr (p+1, pe-p-1));

  if (flags & MALLOCED_MEMORY) {
    mem += mema;
    cerr << aol::color::red << "aol::memusage(MALLOCED_MEMORY) flows over at 4 GiB" << aol::color::reset << endl;
  }
  if (flags & MMAPPED_MEMORY) {
    mem += memtot - mema;
    cerr << aol::color::red << "aol::memusage(MMAPPED_MEMORY) flows over at 4 GiB"  << aol::color::reset << endl;
  }

  return mem;
#else
#ifdef _WIN32
  HANDLE process = OpenProcess ( PROCESS_QUERY_INFORMATION|PROCESS_VM_READ, FALSE, GetCurrentProcessId() );
  if (NULL == process)
    return 0;

  uint64_t mem = 0;

  PROCESS_MEMORY_COUNTERS pmc;
  if ( GetProcessMemoryInfo( process, &pmc, sizeof(pmc)) )
    mem = pmc.WorkingSetSize;

  CloseHandle( process );

  return mem;
#elif defined ( __APPLE__ )
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

  if ( KERN_SUCCESS != task_info(mach_task_self(), TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&t_info), &t_info_count)) {
    return -1;
  }
  return t_info.resident_size;
#else
  throw UnimplementedCodeException ( "memusage", __FILE__, __LINE__ );
#endif
  aol::doNothingWithArgumentToPreventUnusedParameterWarning ( flags );
  return -1;
#endif
}

int vscprintf(const char *format, va_list ap) {
#if defined(_MSC_VER)
  return _vscprintf( format, ap );
#else
  // C99 compliant code (not supported by VC++ 2008)
  return vsnprintf ( NULL, 0, format, ap );
#endif
}

FILE* ipfstream::openFile ( const char *FileName ) {
#ifndef _WIN32
  char fullCommand[4096];
#endif
  if ( fileNameEndsWith ( FileName, ".bz2" ) ) {
#ifdef _WIN32
#ifdef USE_LIB_BZ2
    FILE *tempf;
    BZFILE *tempb;
    TempFileName = new char[256];
    unsigned char buf[200000];
    sprintf ( TempFileName, "%s.pgm", FileName );
    int bzerror;
    tempf = fopen ( FileName, "rb" );
    if ( !tempf ){
      string errorMessage = "Could not open bz2 file ";
      errorMessage += FileName;
      errorMessage += " for reading.\n";
      throw Exception ( errorMessage, __FILE__, __LINE__ );
    }
    ofstream tempfile ( TempFileName, ofstream::binary );
    if ( !tempfile ) {
      fclose ( tempf );
      throw Exception ( "Could not open tempfile for decompressing.", __FILE__, __LINE__ );
    }

    tempb = BZ2_bzReadOpen ( &bzerror, tempf, 0, 0, NULL, 0 );
    while ( bzerror == BZ_OK ) {
      const int bytesRead = BZ2_bzRead ( &bzerror, tempb, buf, sizeof ( buf ) );
      tempfile.write ( reinterpret_cast<char*> ( buf ), bytesRead );
    }
    tempfile.close();
    BZ2_bzReadClose ( &bzerror, tempb );
    fclose ( tempf );
    f = fopen ( TempFileName, "rb" );
    isPipe = false;
    isTemp = true;
#else
    throw Exception ( "Reading bz2 compressed files without using bzlib under windows is impossible", __FILE__, __LINE__ );
#endif // USE_LIB_BZ2
#else
    sprintf ( fullCommand, "bunzip2 -c %s 2>/dev/null", FileName );
    f = popen ( fullCommand, "r" );
    isPipe = true;
#endif // _WIN32
  } else if ( fileNameEndsWith ( FileName, ".gz" ) ) {
#ifdef _WIN32
    throw Exception ( "Reading gz compressed files under windows is not implemented", __FILE__, __LINE__ );
#else
    sprintf ( fullCommand, "gunzip -c %s 2>/dev/null", FileName );
    f = popen ( fullCommand, "r" );
    isPipe = true;
#endif // _WIN32
  } else {
#ifdef _WIN32
    f = fopen ( FileName, "rb" );
    isTemp = false;
#else
    f = fopen ( FileName, "r" );
#endif
    isPipe = false;
  }
  return f;
}

#ifdef __GNUC__
#if __GNUC__ < 3  // old gcc
ipfstream::ipfstream ( const char *name ) : ifstream() {
  openFile ( name );
  if ( f ) attach ( fileno ( f ) );
  else throw Exception ( string ( "Could not open file " ) + name, __FILE__, __LINE__ );
}
#else // new gcc
ipfstream::ipfstream ( const char *name ) :
    istream ( NULL ), isPipe ( false ), f ( openFile ( name ) )
#ifdef HAVE_STDIO_FILEBUF
        , buf ( f ? new stdio_filebuf<char> ( f, ios::in ) : NULL )
#endif
    // isPipe will be set by openFile
{
  if ( !f ) throw Exception ( string ( "Could not open file " ) + name, __FILE__, __LINE__ );
#ifdef HAVE_STDIO_FILEBUF
  rdbuf ( buf );
#endif
}
#endif

#else // other compiler
#ifdef _WIN32
ipfstream::ipfstream ( const char *name ) : ifstream() {
  openFile ( name );
  if ( !f ) throw Exception ( string ( "Could not open file " ) + name, __FILE__, __LINE__ );
  if( isTemp )
    open( TempFileName, ios::binary );
  else
    open( name, ios::binary );
}
#else
ipfstream::ipfstream ( const char *name ) : ifstream ( f = openFile ( name ) ) {
  if ( !f ) throw Exception ( string ( "Could not open file " ) + name, __FILE__, __LINE__ );
}
#endif //_WIN32
#endif //compiler

opfstream::~opfstream() {
  if ( f ) {
    flush();
#ifdef HAVE_STDIO_FILEBUF
    // Detach before closing file/pipe
    rdbuf ( NULL );
    if ( buf ) {
      delete buf;
      buf = NULL;
    }
#endif
#if defined(__MINGW32_VERSION) || defined(__CYGWIN__)
    fclose ( f );
#else
    if ( isPipe ) pclose ( f );
    else fclose ( f );
#endif
    f = NULL;
#if defined(__MINGW32_VERSION) || defined(__CYGWIN__) // new gcc under Windows
    if ( isTemp ) {
#ifdef USE_LIB_BZ2
      FILE *outfile;
      BZFILE *b;
      // copy TempFileName to outfileName without the temporary suffix ".pgm"
      string tempString ( TempFileName );
      string tempString2 ( tempString.begin(), tempString.end() - 4 );
      const char* outfileName = tempString2.c_str();

      unsigned char buf[200000];
      int bzerror;
      FILE *infile;
      infile = fopen ( TempFileName, "rb" );
      outfile = fopen ( outfileName, "wb" );
      b = BZ2_bzWriteOpen ( &bzerror, outfile, 9, 0, 0 );
      while ( !feof ( infile ) ) {
        const int bytesRead = fread ( reinterpret_cast<char*> ( buf ), 1, sizeof ( buf ), infile );
        BZ2_bzWrite ( &bzerror, b, buf, bytesRead );
      }
      unsigned int temp;
      fclose ( infile );
      BZ2_bzWriteClose ( &bzerror, b, 0, &temp, &temp );
      fclose ( outfile );
      remove ( TempFileName );
      delete[] TempFileName;
#else
      throw Exception ( "Writing bz2 compressed files without using bzlib under windows is impossible", __FILE__, __LINE__ );
#endif // USE_LIB_BZ2
    }
#endif //__MINGW32_VERSION
  }
}

bool isErfTestOK ( ) {
#if defined GNU
  const int N = 100;
  double maxDiff = -1.0;
  for ( int i = 0; i < N; ++i ) {
    const double
      xval = -2.0 + i / ( 5.0 * (N-1) ),
      difference = fabs ( erf ( xval ) - aol::Erf ( xval ) );
    if ( difference > maxDiff ) {
      maxDiff = difference;
    }
  }
  return ( maxDiff < 1.5e-7 );
#else
  // erf not provided, so no test possible
  return ( true );
#endif
}

namespace color {
#if !defined(__MINGW32_VERSION) && !defined(_MSC_VER)
const string reset       = "\033[0;0m";
const string invert      = "\033[0;7m";
const string black       = "\033[0;30m";
const string red         = "\033[0;31m";
const string green       = "\033[0;32m";
const string brown       = "\033[0;33m";
const string blue        = "\033[0;34m";
const string purple      = "\033[0;35m";
const string cyan        = "\033[0;36m";
const string light_grey  = "\033[0;37m";
const string dark_grey   = "\033[1;30m";
const string light_red   = "\033[1;31m";
const string light_green = "\033[1;32m";
const string yellow      = "\033[1;33m";
const string light_blue  = "\033[1;34m";
const string pink        = "\033[1;35m";
const string light_cyan  = "\033[1;36m";
const string white       = "\033[1;37m";
const string beep        = "\007";
#else // A standard windows console doesnt support this type of color setting
const string reset       = "";
const string invert      = "";
const string black       = "";
const string red         = "";
const string green       = "";
const string brown       = "";
const string blue        = "";
const string purple      = "";
const string cyan        = "";
const string light_grey  = "";
const string dark_grey   = "";
const string light_red   = "";
const string light_green = "";
const string yellow      = "";
const string light_blue  = "";
const string pink        = "";
const string light_cyan  = "";
const string white       = "";
const string beep        = "";
#endif // __MINGW32_VERSION
const string error       = beep + red;
const string ok          = green;
const string residuum    = blue;
}

} // namespace aol
