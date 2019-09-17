#! \file
#! \brief Main cmake configuration file

CMAKE_MINIMUM_REQUIRED ( VERSION 2.8 )

PROJECT ( quocmesh )

ENABLE_TESTING()

# cmake 2.8.11.2 introduced a new policy regarding automatically linking Qt executables
# to qtmaintarget on Windows. We have to explicitly select the new behavior to
# prevent warnings about this policy change.
IF ( POLICY CMP0020 )
  CMAKE_POLICY ( SET CMP0020 NEW )
ENDIF ( )

# cmake 3.0 introduced a new policy that enables MACOSX_RPATH by default.
# The new behavior should be safe for us to use, so enable it to prevent warnings.
IF ( POLICY CMP0042 )
  CMAKE_POLICY ( SET CMP0042 NEW )
ENDIF ( )

# ADD_COMPILE_DEFINITIONS needs cmake 3.12. Supply a replacement for older versions.
IF ( ${CMAKE_VERSION} VERSION_LESS "3.12.0" )
  MACRO ( ADD_COMPILE_DEFINITIONS )
    SET ( Definitions "${ARGN}" )
    FOREACH ( Definition IN LISTS Definitions )
      ADD_COMPILE_OPTIONS ( "-D${Definition}" )
    ENDFOREACH ( )
  ENDMACRO ( ADD_COMPILE_DEFINITIONS )
ENDIF ( )

SET ( CMAKE_MODULE_PATH )
IF ( WIN32 )
  LIST ( APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/windows/" )
ENDIF ( )
LIST ( APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/" )
# Also allow internal cmake modules.
LIST ( APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/internal/cmake/ )

# Generate a string that characterizes the current platform.
SITE_NAME( QUOC_HOSTNAME )
SET ( QUOC_PLATFORM_INFO "${QUOC_HOSTNAME}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_VERSION}-${CMAKE_CXX_COMPILER_ID}" )

# If the generator supports it, allow duplicate custom targets. This allows
# to create short aliases for the executable targets in STDPROJECTS.
IF ( "${CMAKE_GENERATOR}" MATCHES Make )
  SET_PROPERTY ( GLOBAL PROPERTY ALLOW_DUPLICATE_CUSTOM_TARGETS TRUE )
  SET ( CREATE_SHORT_EXECUTABLE_TARGET_ALIASES 1 )
ELSE ( )
  SET ( CREATE_SHORT_EXECUTABLE_TARGET_ALIASES 0 )
ENDIF ( )

# When using an IDE, turn on folder support to allow the grouping of targets in folders.
IF ( MSVC OR XCODE )
  SET_PROPERTY ( GLOBAL PROPERTY USE_FOLDERS ON )
ENDIF ( )

# If the user didn't create a selection file yet, create one by copying the default one.
IF ( NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/cmake.selection )
  FILE ( COPY cmake.selection.default DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  FILE ( RENAME ${CMAKE_CURRENT_BINARY_DIR}/cmake.selection.default ${CMAKE_CURRENT_BINARY_DIR}/cmake.selection )
ENDIF ( )

# Look for a Python interpreter.
FIND_PACKAGE ( PythonInterp QUIET )

# Allow the user to deactivate the vector manager. Since some options (like USE_OPENMP) may need to turn off
# the vector manager, we define the option here but only turn off the vector manager after the other options
# had a chance to overwrite DO_NOT_USE_MEMORYMANAGER.
#! \cmakeoption{Deactivate memory manager,OFF (i.e. memory manager <u>is</u> used)}
OPTION ( DO_NOT_USE_MEMORYMANAGER "" OFF )

#! \cmakeoption{Generate header inclusion test,OFF}
OPTION ( GENERATE_INCLUDE_TEST "In development cmake replacement of util/testIndividualInclusion.sh" OFF )

#! \cmakeoption{Use OpenMP,OFF}
#! \note for gcc and msvc compiler only
OPTION ( USE_OPENMP "" OFF )

# First configure some compiler specific things.
IF ( MSVC )
  ADD_COMPILE_DEFINITIONS ( _UNICODE _CRT_SECURE_NO_WARNINGS _SCL_SECURE_NO_WARNINGS )
  # Turn off the "conversion from 'A' to 'B', possible loss of data" warnings. Since our gcc configuration
  # doesn't show them, huge amounts of them have accumulated in the code and it would be a lot of work
  # (with little gain) to fix them
  ADD_COMPILE_OPTIONS ( /wd4244 )
  SET ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:libcmt.lib" )
  SET ( CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /NODEFAULTLIB:msvcrt.lib" )

  IF ( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    # With x64, the object file for tpCFEGrid.cpp gets so big that we need to enable large object files.
    SET ( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /bigobj" )
  ENDIF ( )

  IF ( USE_OPENMP )
    IF ( CMAKE_CXX_COMPILER_ID STREQUAL "Clang" )
      ADD_COMPILE_OPTIONS ( -Xclang -fopenmp )
      FIND_LIBRARY ( OMP_LIBRARY NAMES omp libomp )
      IF ( NOT OMP_LIBRARY )
        MESSAGE ( FATAL_ERROR "Could not find libomp." )
      ENDIF ()
      LIST ( APPEND SYSTEM_LIBRARIES ${OMP_LIBRARY} )
    ELSE ( )
      ADD_COMPILE_OPTIONS ( /openmp )
    ENDIF ( )
    SET ( DO_NOT_USE_MEMORYMANAGER 1 )
  ENDIF ( USE_OPENMP )

  # VC++ doesn't have dirent.h, so try to find it. For this we supply a VC++ compatible version in external/dirent.
  FIND_PATH ( MSVC_DIRENT_INCLUDE_DIR dirent.h external/dirent )
  IF ( MSVC_DIRENT_INCLUDE_DIR )
    INCLUDE_DIRECTORIES ( SYSTEM ${MSVC_DIRENT_INCLUDE_DIR} )
  ELSE ( )
    MESSAGE ( FATAL_ERROR "Could not find dirent.h." )
  ENDIF ( )
ENDIF ( MSVC )

IF ( CMAKE_COMPILER_IS_GNUCXX )
  #! \cmakeoption{Parse gcc error output,ON}
  #! \note for gcc compiler only
  OPTION ( PARSE_GCC_ERRORS "Make cmake use our GCC error parsing script like colorize.pl" ON )
  IF ( PARSE_GCC_ERRORS )
    #  To let the script know where the source and thus the error parse scripts are sneak in the path as argument.
    SET_PROPERTY ( GLOBAL PROPERTY RULE_LAUNCH_COMPILE "bash ${CMAKE_CURRENT_SOURCE_DIR}/util/cmakeParseError.sh ${CMAKE_CURRENT_SOURCE_DIR}" )
  ENDIF ( )

  # Figure out which GCC version we are using.
  EXECUTE_PROCESS ( COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE QUOC_GCC_VERSION )

  ADD_COMPILE_OPTIONS ( -fmessage-length=0 -Wcast-qual -Wall -Wextra -Wno-conversion -Winit-self -Wcast-align -Wundef -Wswitch-default -Wredundant-decls -pedantic )
  ADD_COMPILE_DEFINITIONS ( GNU )
  # The following flags are only valid for C++.
  SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wold-style-cast -Woverloaded-virtual -Wnon-virtual-dtor" )
  IF ( MINGW )
    IF ( NOT ( ${QUOC_GCC_VERSION} VERSION_LESS 5.1 ) AND ( CMAKE_SIZEOF_VOID_P EQUAL 8 ) )
      # With x64, the object file for tpCFEGrid.cpp gets so big that we need to enable large object files.
      SET ( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wa,-mbig-obj" )
    ENDIF ( )
  ELSE ( )
    ADD_COMPILE_OPTIONS ( -fPIC -ansi )

    #! \cmakeoption{Use efence when in debug mode,ON}
    #! \note for gcc compiler only
    OPTION ( USE_EFENCE_IN_DEBUG_MODE "In debug configurations link efence if enabled." OFF )
    IF ( USE_EFENCE_IN_DEBUG_MODE )
      SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} debug -lefence )
    ENDIF ()
  ENDIF ( )
  # G++ compiler flags only used in release mode.
  SET ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wdisabled-optimization -funroll-loops" )

  # aol::areWeUsingEFence needs this lib (not available under Windows).
  IF ( NOT WIN32 )
    SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} -ldl )
  ENDIF ( )

  # Support for "march=native" was introduced in GCC version 4.2.0.
  # The GCC 4.2.1 version that comes with OS X doesn't seem to support this.
  # Newer GCC versions do support this under OS X if clang is used as assembler, but this leads to strange problems.
  # For instance, the tpcfe selfTest fails, so don't use this as default there.
  IF ( NOT ( ${QUOC_GCC_VERSION} VERSION_LESS 4.2 ) AND ( NOT APPLE ) )
    #! \cmakeoption{gcc -march flag,native for gcc version 4.2 or greater}
    SET ( GCC_SYSTEM_TYPE "native" CACHE STRING "select system type for gcc -march=xxx" )
  ELSE ( )
    SET ( GCC_SYSTEM_TYPE "" CACHE STRING "select system type for gcc -march=xxx" )
  ENDIF ( )

  IF ( NOT ( GCC_SYSTEM_TYPE STREQUAL "" ) )
    # To use "march=native" with a CPU that supports AVX under OS X, clang has to be used as assembler.
    # The GNU assembler for OS X is too old to support these extensions.
    IF ( APPLE AND ( GCC_SYSTEM_TYPE STREQUAL "native" ) )
      SET ( CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Wa,-q" )
      SET ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wa,-q" )
    ENDIF ( )
    SET ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=${GCC_SYSTEM_TYPE}" )
    SET ( CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=${GCC_SYSTEM_TYPE}" )
  ENDIF ( )

  # MinGW GCC 4.6.1 warns about some things that are obviously correct, so turn off the corresponding warning.
  # The same is true for GCC 4.6.3 under Ubuntu, so we turn if off for GCC 4.6.0 to 4.6.3.
  IF ( ( NOT ( ${QUOC_GCC_VERSION} VERSION_LESS 4.6 ) ) AND ( ${QUOC_GCC_VERSION} VERSION_LESS 4.6.4 ) )
    ADD_COMPILE_OPTIONS ( -Wno-array-bounds )
  ENDIF ( )

  #! \cmakeoption{Activate SSE,ON}
  #! \note for gcc compiler only
  OPTION ( USE_SSE "" ON )
  IF ( USE_SSE )
    SET ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -msse2 -DUSE_SSE" )
  ENDIF ( )

  #! \cmakeoption{Use MPI,OFF}
  #! \note for gcc compiler only
  OPTION ( USE_MPI "" OFF )
  IF ( USE_MPI)
     ADD_COMPILE_DEFINITIONS ( MPI__)
     FIND_PACKAGE ( MPI REQUIRED )
     INCLUDE_DIRECTORIES ( SYSTEM ${MPI_INCLUDE_PATH} )
     LIST ( APPEND SYSTEM_LIBRARIES ${MPI_LIBRARIES} )
  ENDIF ( USE_MPI )

  IF ( USE_OPENMP )
    ADD_COMPILE_OPTIONS ( -fopenmp )
    SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} -lgomp )
    # With OpenMP the vector manager does more harm then good, so turn it off
    # (overriding the user setting but keeping the cache value).
    SET ( DO_NOT_USE_MEMORYMANAGER 1 )
  ENDIF ( USE_OPENMP )

  #! \cmakeoption{Use gprof profiler,OFF}
  #! \note for gcc compiler only
  OPTION ( USE_PROFILING "" OFF )
  IF ( USE_PROFILING )
    ADD_COMPILE_OPTIONS ( -pg )
    SET ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg" )
    IF ( APPLE )
      MESSAGE ( FATAL_ERROR "The -pg option seems to be broken under OS X. Thus, profiling with GCC is not supported here. Instruments.app is a potential profiling alternative." )
    ENDIF ( )
  ENDIF ( USE_PROFILING )

  #! \cmakeoption{Use libstdc++ profile mode,OFF}
  #! \note for gcc compiler only
  OPTION ( USE_STD_PROFILING "Profiling and automatic suggestions for improvement when using STL data structures." OFF )
  IF ( USE_STD_PROFILING )
    SET ( CMAKE_BUILD_TYPE RelWithDebInfo )
    STRING ( REPLACE "-O2" "-O0" CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
    STRING ( REPLACE "-O2" "-O0" CMAKE_C_FLAGS_RELWITHDEBINFO ${CMAKE_C_FLAGS_RELWITHDEBINFO})
    ADD_COMPILE_OPTIONS ( -fno-inline )
    ADD_COMPILE_DEFINITIONS ( _GLIBCXX_PROFILE )
  ENDIF ( USE_STD_PROFILING )

  IF ( NOT ( ${QUOC_GCC_VERSION} VERSION_LESS 4.6 ) )
    #! \cmakeoption{Enable loop optimizer,OFF}
    #! \note for gcc compiler version 4.6 or greater only
    OPTION ( USE_OPT_GRAPHITE "Enable the graphite loop optimizer" OFF )
    IF ( USE_OPT_GRAPHITE )
      SET ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -floop-interchange -floop-strip-mine -floop-block" )
    ENDIF ( USE_OPT_GRAPHITE )
  ENDIF ( )

  IF ( NOT ( ${QUOC_GCC_VERSION} VERSION_LESS 4.5 ) )
    #! \cmakeoption{Use link time optimization,OFF}
    #! \note for gcc compiler version 4.5 or greater only
    OPTION ( USE_OPT_LTO "Use link time optimization (LTO)" OFF )
    IF ( USE_OPT_LTO )
      SET ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto" )
      SET ( CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -flto" )
      SET ( CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} -flto" )
      SET ( CMAKE_MODULE_LINKER_FLAGS_RELEASE "${CMAKE_MODULE_LINKER_FLAGS_RELEASE} -flto" )
    ENDIF ( USE_OPT_LTO )
  ENDIF ( )

  IF ( NOT ( ${QUOC_GCC_VERSION} VERSION_LESS 4.7 ) )
    #! \cmakeoption{Use agressive optimization flags,OFF}
    #! \note for gcc compiler version 4.7 or greater only
    OPTION ( USE_OPT_AGGRESSIVE "Enable more aggressive optimization flags, at the moment -Ofast which enables -ffast-math among other flags." OFF )
    IF ( USE_OPT_AGGRESSIVE )
      SET ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast" )
    ENDIF ( USE_OPT_AGGRESSIVE )
  ENDIF ( )
ENDIF ( CMAKE_COMPILER_IS_GNUCXX )

#! \cmakeoption{Build and use DUMA,OFF}
OPTION ( BUILD_AND_USE_DUMA "Build and use DUMA. If the sources are not available try to download them." OFF )
IF ( BUILD_AND_USE_DUMA )
  # If the duma source is not the the external directory, try to download it.
  IF ( NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/duma-source/duma.c )
    SET ( DUMA_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/duma-source )
    SET ( DUMA_ARCHIVE ${DUMA_SOURCE_DIR}/duma.tar.gz )
    FILE ( DOWNLOAD http://sourceforge.net/projects/duma/files/latest/download ${DUMA_ARCHIVE} SHOW_PROGRESS )
    EXECUTE_PROCESS ( COMMAND ${CMAKE_COMMAND} -E tar xzf ${DUMA_ARCHIVE} WORKING_DIRECTORY ${DUMA_SOURCE_DIR} )

    # Find out in which subdirectory the source was extracted
    FILE ( GLOB DUMA_SUBDIRS RELATIVE ${DUMA_SOURCE_DIR} ${DUMA_SOURCE_DIR}/*)
    FOREACH ( DUMA_SUBDIR ${DUMA_SUBDIRS})
      IF ( IS_DIRECTORY ${DUMA_SOURCE_DIR}/${DUMA_SUBDIR} )
        SET ( DUMA_EXTRACTED_SOURCE_DIR ${DUMA_SOURCE_DIR}/${DUMA_SUBDIR} )
        BREAK ( )
      ENDIF ( )
    ENDFOREACH ( )

    # Copy all files from this subdirectory into the directory we want to have the source in.
    FILE ( GLOB DUMA_FILES ${DUMA_EXTRACTED_SOURCE_DIR}/*.c ${DUMA_EXTRACTED_SOURCE_DIR}/*.cpp ${DUMA_EXTRACTED_SOURCE_DIR}/*.h )
    FOREACH ( DUMA_FILE ${DUMA_FILES} )
      FILE ( COPY ${DUMA_FILE} DESTINATION ${DUMA_SOURCE_DIR} )
    ENDFOREACH()

    # Remove the extracted subdirectory and the downloaded archie, we have extracted all the files we need.
    EXECUTE_PROCESS ( COMMAND ${CMAKE_COMMAND} -E remove_directory ${DUMA_EXTRACTED_SOURCE_DIR} WORKING_DIRECTORY ${DUMA_SOURCE_DIR} )
    FILE ( REMOVE ${DUMA_ARCHIVE} )
  ENDIF ( )

  ADD_SUBDIRECTORY ( external/duma-source )
  INCLUDE_DIRECTORIES ( SYSTEM external/duma-source ${CMAKE_CURRENT_BINARY_DIR}/external/duma-source )
  LIST ( APPEND SYSTEM_LIBRARIES ${EXTERNAL_DUMA_LIBRARIES} )
  ADD_COMPILE_DEFINITIONS ( USE_DUMA )
ENDIF ( )

#! \cmakeoption{Build and use TetGen,OFF}
OPTION ( BUILD_AND_USE_TETGEN "" OFF )
IF ( BUILD_AND_USE_TETGEN )
  ADD_SUBDIRECTORY ( external/tetgen )
  INCLUDE_DIRECTORIES ( SYSTEM external/tetgen ${CMAKE_CURRENT_BINARY_DIR}/external/tetgen )
  LIST ( APPEND SYSTEM_LIBRARIES ${EXTERNAL_TET_LIBRARIES} )
  ADD_COMPILE_DEFINITIONS ( USE_TETGEN )
  ADD_COMPILE_DEFINITIONS ( TETLIBRARY )
ENDIF ( )

IF ( DO_NOT_USE_MEMORYMANAGER )
  ADD_COMPILE_DEFINITIONS ( DO_NOT_USE_MEMORYMANAGER )
ENDIF ( )

# Determine bitness
IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  ADD_COMPILE_DEFINITIONS( BITNESS=64 )
ELSEIF ( CMAKE_SIZEOF_VOID_P EQUAL 4 )
  ADD_COMPILE_DEFINITIONS( BITNESS=32 )
ELSE ( )
  MESSAGE ( WARNING "Could not determine bitness from size of void pointer being ${CMAKE_SIZEOF_VOID_P}." )
ENDIF ( )

# aol::memusage needs this lib under Windows.
IF ( WIN32 )
  SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} psapi )
ENDIF ( )

# Addtional flags for debugging
SET ( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DBOUNDS_CHECK -DDEBUG -DDO_NOT_USE_MEMORYMANAGER" )

#! \cmakeoption{Use STL debug mode,OFF}
OPTION ( USE_STLDEBUG "Bounds checking etc. for STL" OFF )
IF ( USE_STLDEBUG )
  SET ( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG" )
ENDIF ( USE_STLDEBUG )

IF ( CMAKE_COMPILER_IS_GNUCXX )
  IF ( NOT ( ${QUOC_GCC_VERSION} VERSION_LESS 4.6 ) )
    #! \cmakeoption{Use C++11,OFF}
    #! \note for gcc compiler version 4.7 or greater or Clang only
    OPTION ( USE_C++11 "Make the new c++ standard available" OFF )
    IF ( USE_C++11 )
      MESSAGE ( WARNING "Option USE_C++11 is deprecated, please set 'C++_STANDARD=C++11' instead" )
    ENDIF ( )
  ENDIF ( )
ELSE ( )
  IF ( CMAKE_CXX_COMPILER_ID STREQUAL "Clang" )
    IF ( APPLE OR WIN32 )
      # Turn off some warnings that were introduced with recent Clang versions. Note that the Clang variant
      # that comes with XCode uses a completely different version number scheme.
      IF ( CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 8.1 OR ( CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 3.9 AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.0 ) )
        ADD_COMPILE_OPTIONS ( -Wno-undefined-var-template -Wno-null-dereference )
      ENDIF ( )
    ENDIF ( )
    # additional libraries only required for linux?
    IF ( ( NOT WIN32 ) AND ( NOT APPLE ) )
      LIST ( APPEND SYSTEM_LIBRARIES -lstdc++ -lm -ldl )
    ENDIF ( )
    IF ( MSVC )
      ADD_COMPILE_DEFINITIONS ( UNICODE )
      SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fms-compatibility-version=19" )
    ENDIF ( )
    OPTION ( USE_C++11 "" OFF )
    IF ( USE_C++11 )
      # Unlike GCC, Clang doesn't allow the option -std=c++11 to be used with C code.
      IF ( NOT MSVC )
        SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
      ENDIF ( )
      ADD_COMPILE_DEFINITIONS ( USE_CPP11 )
      IF ( APPLE )
        # Furthermore, to use unordered_set with clang it is also necessary to use libc++.
        SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++" )
        # Moreover, at least under OS X 10.8.5, we need to explicitly link libc++.
        LIST ( APPEND SYSTEM_LIBRARIES -lc++ )
      ENDIF ( )
    ENDIF ( )
  ENDIF ( )
ENDIF ( )


#! \cmakeoption{C++ standard,default}
#! meant to supersede USE_C++11 option
SET ( C++_STANDARD "default" CACHE STRING "C++ standard to use (default means compiler's default will be used)")
# add possible values for C++_STANDARD option, always default and c++98
SET_PROPERTY ( CACHE C++_STANDARD PROPERTY STRINGS default C++98 )

IF ( CMAKE_COMPILER_IS_GNUCXX )
  # add further values depending on version
  # C++11 was available earlier than gcc 4.6 but we tested 4.6
  IF ( NOT ( ${QUOC_GCC_VERSION} VERSION_LESS 4.6 ) )
    SET_PROPERTY ( CACHE C++_STANDARD APPEND PROPERTY STRINGS C++11 )
    # Similarly C++14 was certainly available earlier
    IF ( NOT( ${QUOC_GCC_VERSION} VERSION_LESS 5.1 ) )
      SET_PROPERTY ( CACHE C++_STANDARD APPEND PROPERTY STRINGS C++14 )
    ENDIF()
  ENDIF()

  # set flags and defines
  # check default
  IF ( C++_STANDARD STREQUAL "default" )
    # do not add a compiler flag
    REMOVE_DEFINITIONS ( -DUSE_CPP11 -DUSE_CPP14 )
    # gcc 5 defaults to -std=gnu11 only for C code
  ENDIF()

  # check C++98
  IF ( C++_STANDARD STREQUAL "C++98" )
    ADD_COMPILE_OPTIONS ( -std=c++98 )
    REMOVE_DEFINITIONS ( -DUSE_CPP11 -DUSE_CPP14 )
  ENDIF()

  # check C++11
  IF ( ( C++_STANDARD STREQUAL "C++11" ) OR ( USE_C++11 ) )
    IF ( ${QUOC_GCC_VERSION} VERSION_LESS 4.7 )
      ADD_COMPILE_OPTIONS ( -std=c++0x )
    ELSE ( )
      ADD_COMPILE_OPTIONS ( -std=c++11 )
    ENDIF ( )
    ADD_COMPILE_DEFINITIONS ( USE_CPP11 )
    REMOVE_DEFINITIONS ( -DUSE_CPP14 )
  ENDIF()

  # check C++14
  IF ( C++_STANDARD STREQUAL "C++14" )
    ADD_COMPILE_OPTIONS ( -std=c++14 )
    ADD_COMPILE_DEFINITIONS ( USE_CPP11 USE_CPP14 )
  ENDIF()
ENDIF()


# Look for libbz2. If we don't find it, build it from the source in external/bz2.
IF ( NOT MINGW )
  FIND_PACKAGE ( BZip2 QUIET )
ENDIF ( )
IF ( BZIP2_FOUND )
  INCLUDE_DIRECTORIES ( ${BZIP2_INCLUDE_DIR} )
  SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${BZIP2_LIBRARIES} )
ELSE ( BZIP2_FOUND )
  ADD_SUBDIRECTORY( external/bz2 )
  INCLUDE_DIRECTORIES ( ${EXTERNAL_BZIP2_INCLUDE_DIR} )
  SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${EXTERNAL_BZIP2_LIBRARIES} )
ENDIF ( BZIP2_FOUND )

ADD_COMPILE_DEFINITIONS ( USE_LIB_BZ2 )

# Look for zlib. If it is not found it could be build from sources, see BUILD_ZLIB_AND_LIBPNG.
FIND_PACKAGE ( ZLIB QUIET )
IF ( ZLIB_FOUND )
  INCLUDE_DIRECTORIES ( ${ZLIB_INCLUDE_DIRS} )
  SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${ZLIB_LIBRARIES} )
  ADD_COMPILE_DEFINITIONS ( USE_ZLIB )
ENDIF ( )

#! \cmakeoption{Use png library,ON}
OPTION ( USE_PNG "" OFF )
IF ( USE_PNG )
  #! \cmakeoption{Build libz and libpng,OFF}
  OPTION ( BUILD_ZLIB_AND_LIBPNG "Build libpng and zlib (which is required for libpng) from sources." OFF )
  IF ( NOT BUILD_ZLIB_AND_LIBPNG )
    FIND_PACKAGE ( PNG REQUIRED )
    IF ( PNG_FOUND )
      INCLUDE_DIRECTORIES ( ${PNG_INCLUDE_DIR} )
      SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${PNG_LIBRARIES} )
      ADD_COMPILE_DEFINITIONS ( USE_LIB_PNG )
    ENDIF ( PNG_FOUND )
  ELSE ( )
    # If the zlib source is not the the external directory, try to download it.
    IF ( NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/zlib-source/zlib.h )
      SET ( ZLIB_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/zlib-source )
      SET ( ZLIB_ARCHIVE ${ZLIB_SOURCE_DIR}/zlib.tar.gz )
      FILE ( DOWNLOAD https://zlib.net/zlib-1.2.11.tar.gz ${ZLIB_ARCHIVE} SHOW_PROGRESS )
      EXECUTE_PROCESS ( COMMAND ${CMAKE_COMMAND} -E tar xzf ${ZLIB_ARCHIVE} WORKING_DIRECTORY ${ZLIB_SOURCE_DIR} )

      SET ( ZLIB_EXTRACTED_SOURCE_DIR ${ZLIB_SOURCE_DIR}/zlib-1.2.11 )

      # Copy all files from this subdirectory into the directory we want to have the source in.
      FILE ( GLOB ZLIB_FILES ${ZLIB_EXTRACTED_SOURCE_DIR}/*.c ${ZLIB_EXTRACTED_SOURCE_DIR}/*.h )
      FOREACH ( ZLIB_FILE ${ZLIB_FILES} )
       FILE ( COPY ${ZLIB_FILE} DESTINATION ${ZLIB_SOURCE_DIR} )
      ENDFOREACH()

      # Remove the extracted subdirectory and the downloaded archive, we have extracted all the files we need.
      EXECUTE_PROCESS ( COMMAND ${CMAKE_COMMAND} -E remove_directory ${ZLIB_EXTRACTED_SOURCE_DIR} WORKING_DIRECTORY ${ZLIB_SOURCE_DIR} )
      FILE ( REMOVE ${ZLIB_ARCHIVE} )
    ENDIF ( )

    ADD_SUBDIRECTORY( external/zlib-source )
    INCLUDE_DIRECTORIES ( ${EXTERNAL_ZLIB_INCLUDE_DIR} )
    SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${EXTERNAL_ZLIB_LIBRARIES} )

    # If the libpng source is not the the external directory, try to download it.
    IF ( NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/png-source/png.h )
      SET ( LIBPNG_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/png-source )
      SET ( LIBPNG_ARCHIVE ${LIBPNG_SOURCE_DIR}/libpng.tar.gz )
      FILE ( DOWNLOAD https://ftp-osl.osuosl.org/pub/libpng/src/libpng14/libpng-1.4.22.tar.gz ${LIBPNG_ARCHIVE} SHOW_PROGRESS )
      EXECUTE_PROCESS ( COMMAND ${CMAKE_COMMAND} -E tar xzf ${LIBPNG_ARCHIVE} WORKING_DIRECTORY ${LIBPNG_SOURCE_DIR} )

      SET ( LIBPNG_EXTRACTED_SOURCE_DIR ${LIBPNG_SOURCE_DIR}/libpng-1.4.22 )

      # Copy all files from this subdirectory into the directory we want to have the source in.
      FILE ( GLOB LIBPNG_FILES ${LIBPNG_EXTRACTED_SOURCE_DIR}/*.c ${LIBPNG_EXTRACTED_SOURCE_DIR}/*.h )
      FOREACH ( LIBPNG_FILE ${LIBPNG_FILES} )
       FILE ( COPY ${LIBPNG_FILE} DESTINATION ${LIBPNG_SOURCE_DIR} )
      ENDFOREACH()

      # Remove the extracted subdirectory and the downloaded archive, we have extracted all the files we need.
      EXECUTE_PROCESS ( COMMAND ${CMAKE_COMMAND} -E remove_directory ${LIBPNG_EXTRACTED_SOURCE_DIR} WORKING_DIRECTORY ${LIBPNG_SOURCE_DIR} )
      FILE ( REMOVE ${LIBPNG_ARCHIVE} )
    ENDIF ( )

    ADD_SUBDIRECTORY( external/png-source )
    INCLUDE_DIRECTORIES ( SYSTEM ${EXTERNAL_PNG_INCLUDE_DIR} )
    SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${EXTERNAL_PNG_LIBRARIES} )
    ADD_COMPILE_DEFINITIONS ( USE_LIB_PNG )
  ENDIF ( )
ENDIF ( )

#! \cmakeoption{Build and use wavelet library,OFF}
OPTION ( BUILD_AND_USE_WAVELET "Build and use wavelet1d. If the sources are not available try to download them." OFF )
IF ( BUILD_AND_USE_WAVELET )
  # If the wavelet1d source is not in the external directory, try to download it and apply patch.
  IF ( NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/wavelet1d-source/wavelet.cpp )
  	SET ( WAVELET_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/wavelet1d-source )
  	SET ( WAVELET_ARCHIVE ${WAVELET_SOURCE_DIR}/wavelet1d.rar )
  	FILE ( DOWNLOAD http://wavelet1d.googlecode.com/files/wavelet-03.rar ${WAVELET_ARCHIVE} SHOW_PROGRESS )
  	EXECUTE_PROCESS ( COMMAND ${CMAKE_COMMAND} -E tar xzf ${WAVELET_ARCHIVE} WORKING_DIRECTORY ${WAVELET_SOURCE_DIR} )

  	# Copy all files into the directory we want to have the source in
  	FILE ( COPY ${WAVELET_SOURCE_DIR}/wavelet-03/src/wavelet.cpp DESTINATION ${WAVELET_SOURCE_DIR} )
	FILE ( COPY ${WAVELET_SOURCE_DIR}/wavelet-03/src/wavelet.h DESTINATION ${WAVELET_SOURCE_DIR} )

	# Remove the downloaded and extracted archive, we have extracted all the files we need.
    EXECUTE_PROCESS ( COMMAND ${CMAKE_COMMAND} -E remove_directory ${WAVELET_SOURCE_DIR}/wavelet-03 WORKING_DIRECTORY ${WAVELET_SOURCE_DIR} )
    FILE ( REMOVE ${WAVELET_ARCHIVE} )

    #Apply patch
    IF ( NOT WIN32 )
      EXECUTE_PROCESS ( COMMAND patch -p0 -i ${WAVELET_SOURCE_DIR}/wavelet.patch ${WAVELET_SOURCE_DIR}/wavelet.cpp )
    ENDIF ( )
  ENDIF ( )

  ADD_SUBDIRECTORY( external/wavelet1d-source )
  INCLUDE_DIRECTORIES ( SYSTEM external/wavelet1d-source )
  LIST ( APPEND SYSTEM_LIBRARIES ${EXTERNAL_WAVELET_LIBRARIES} )
  ADD_COMPILE_DEFINITIONS ( USE_LIB_WAVELET )
ENDIF ( )

#! \cmakeoption{Build and use wavelet library,OFF}
OPTION ( BUILD_AND_USE_HUNGARIAN "Build and use hungarian. If the sources are not available try to download them." OFF )
IF ( BUILD_AND_USE_HUNGARIAN )
  # If the hungarian source is not in the external directory, try to download it and apply patch.
  IF ( NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/hungarian-source/Hungarian.cpp )
  	SET ( HUNGARIAN_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/hungarian-source )
  	SET ( HUNGARIAN_ARCHIVE ${HUNGARIAN_SOURCE_DIR}/hungarian-v2.0.zip )
  	FILE ( DOWNLOAD http://robotics.usc.edu/~lantao/codes/hungarian-v2.0.zip ${HUNGARIAN_ARCHIVE} SHOW_PROGRESS )
  	EXECUTE_PROCESS ( COMMAND ${CMAKE_COMMAND} -E tar xzf ${HUNGARIAN_ARCHIVE} WORKING_DIRECTORY ${HUNGARIAN_SOURCE_DIR} )

  	# Copy all files into the directory we want to have the source in
  	FILE ( COPY ${HUNGARIAN_SOURCE_DIR}/hungarian-v2.0/Hungarian.h DESTINATION ${HUNGARIAN_SOURCE_DIR} )
	FILE ( COPY ${HUNGARIAN_SOURCE_DIR}/hungarian-v2.0/Hungarian.cpp DESTINATION ${HUNGARIAN_SOURCE_DIR} )
	FILE ( COPY ${HUNGARIAN_SOURCE_DIR}/hungarian-v2.0/Assignment.h DESTINATION ${HUNGARIAN_SOURCE_DIR} )
	FILE ( COPY ${HUNGARIAN_SOURCE_DIR}/hungarian-v2.0/Assignment.cpp DESTINATION ${HUNGARIAN_SOURCE_DIR} )
	FILE ( COPY ${HUNGARIAN_SOURCE_DIR}/hungarian-v2.0/BipartiteGraph.h DESTINATION ${HUNGARIAN_SOURCE_DIR} )
	FILE ( COPY ${HUNGARIAN_SOURCE_DIR}/hungarian-v2.0/BipartiteGraph.cpp DESTINATION ${HUNGARIAN_SOURCE_DIR} )
	FILE ( COPY ${HUNGARIAN_SOURCE_DIR}/hungarian-v2.0/Matrix.h DESTINATION ${HUNGARIAN_SOURCE_DIR} )
	FILE ( RENAME ${HUNGARIAN_SOURCE_DIR}/Matrix.h ${HUNGARIAN_SOURCE_DIR}/HungarianMatrix.h )
	FILE ( COPY ${HUNGARIAN_SOURCE_DIR}/hungarian-v2.0/README DESTINATION ${HUNGARIAN_SOURCE_DIR} )

	# Remove the downloaded and extracted archive, we have extracted all the files we need.
    EXECUTE_PROCESS ( COMMAND ${CMAKE_COMMAND} -E remove_directory ${HUNGARIAN_SOURCE_DIR}/hungarian-v2.0 WORKING_DIRECTORY ${HUNGARIAN_SOURCE_DIR} )
    FILE ( REMOVE ${HUNGARIAN_ARCHIVE} )

    #Apply patch
    IF ( NOT WIN32 )
      EXECUTE_PROCESS ( COMMAND patch -p0 -i ${HUNGARIAN_SOURCE_DIR}/AssignmentH.patch ${HUNGARIAN_SOURCE_DIR}/Assignment.h )
      EXECUTE_PROCESS ( COMMAND patch -p0 -i ${HUNGARIAN_SOURCE_DIR}/AssignmentCPP.patch ${HUNGARIAN_SOURCE_DIR}/Assignment.cpp )
      EXECUTE_PROCESS ( COMMAND patch -p0 -i ${HUNGARIAN_SOURCE_DIR}/BipartiteGraphH.patch ${HUNGARIAN_SOURCE_DIR}/BipartiteGraph.h )
      EXECUTE_PROCESS ( COMMAND patch -p0 -i ${HUNGARIAN_SOURCE_DIR}/BipartiteGraphCPP.patch ${HUNGARIAN_SOURCE_DIR}/BipartiteGraph.cpp )
      EXECUTE_PROCESS ( COMMAND patch -p0 -i ${HUNGARIAN_SOURCE_DIR}/HungarianH.patch ${HUNGARIAN_SOURCE_DIR}/Hungarian.h )
      EXECUTE_PROCESS ( COMMAND patch -p0 -i ${HUNGARIAN_SOURCE_DIR}/HungarianCPP.patch ${HUNGARIAN_SOURCE_DIR}/Hungarian.cpp )
      EXECUTE_PROCESS ( COMMAND patch -p0 -i ${HUNGARIAN_SOURCE_DIR}/HungarianMatrixH.patch ${HUNGARIAN_SOURCE_DIR}/HungarianMatrix.h )
    ENDIF ( )
  ENDIF ( )

  ADD_SUBDIRECTORY( external/hungarian-source )
  INCLUDE_DIRECTORIES ( SYSTEM external/hungarian-source )
  LIST ( APPEND SYSTEM_LIBRARIES ${EXTERNAL_HUNGARIAN_LIBRARIES} )
  ADD_COMPILE_DEFINITIONS ( USE_LIB_HUNGARIAN )
ENDIF ( )

#! \cmakeoption{Use tiff library,OFF}
OPTION ( USE_TIFF "" OFF )
IF ( USE_TIFF )
  FIND_PACKAGE ( TIFF )
  IF ( TIFF_FOUND )
    INCLUDE_DIRECTORIES ( SYSTEM ${TIFF_INCLUDE_DIR} )
    LIST ( APPEND SYSTEM_LIBRARIES ${TIFF_LIBRARIES} )
    ADD_COMPILE_DEFINITIONS ( USE_LIB_TIFF )
  ENDIF ( )
ENDIF ( )

#! \cmakeoption{Use boost library,OFF}
OPTION ( USE_BOOST "" OFF )
IF ( USE_BOOST )
  SET ( BOOST_COMPONENTS "" CACHE STRING "Boost components that are not header only and should be linked (; seperated)" )
#  FIND_PACKAGE ( Boost COMPONENTS "${BOOST_COMPONENTS}" REQUIRED )
FIND_PACKAGE ( Boost COMPONENTS "${BOOST_COMPONENTS}" REQUIRED )
  IF ( Boost_FOUND )
    INCLUDE_DIRECTORIES ( SYSTEM ${Boost_INCLUDE_DIRS} )
    LINK_DIRECTORIES ( ${Boost_LIBRARY_DIR} )
    ADD_COMPILE_DEFINITIONS ( USE_BOOST )
    LIST ( APPEND SYSTEM_LIBRARIES ${Boost_LIBRARIES} )
    # Possibly this also needs to be activated for older GCC versions.
    IF ( CMAKE_COMPILER_IS_GNUCXX )
      IF ( NOT ( ${QUOC_GCC_VERSION} VERSION_LESS 6.3 ))
        SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fext-numeric-literals" )
      ENDIF ( )
    ENDIF ( )
  ENDIF ( )
ENDIF( USE_BOOST )

#! \cmakeoption{Use blas library,ON}
OPTION ( USE_BLAS "" OFF )
IF ( USE_BLAS )
  # Module of cmake 2.8.x does not really work, module of cmake 3.1.0 does well
  CMAKE_MINIMUM_REQUIRED ( VERSION 3.1 )
  FIND_PACKAGE ( BLAS )
  # Look for OpenBLASConfig.cmake config file, also try OpenBLASST
  FIND_PACKAGE( OpenBLAS QUIET NAMES OpenBLAS OpenBLASST )
  IF ( OpenBLAS_FOUND )
    SET ( OpenBLAS_INCLUDE_DIRS ${OpenBLAS_INCLUDE_DIRS} ${OpenBLASST_INCLUDE_DIRS} )
    SET ( OpenBLAS_LIBRARIES    ${OpenBLAS_LIBRARIES}    ${OpenBLASST_LIBRARIES} )
    FIND_PACKAGE_MESSAGE ( OpenBLAS "Found OpenBLAS: ${OpenBLAS_LIBRARIES}" "[${OpenBLAS_LIBRARIES}][${OpenBLAS_INCLUDE_DIRS}]" )
  ENDIF ( )
ENDIF( USE_BLAS )

#! \cmakeoption{Use lapack library,ON}
OPTION ( USE_LAPACK "" OFF )
IF ( USE_LAPACK )
  # Module of cmake 2.8.6 does not find local lib
  CMAKE_MINIMUM_REQUIRED ( VERSION 3.1 )
  FIND_PACKAGE ( LAPACK REQUIRED )
ENDIF ( USE_LAPACK )

#! \cmakeoption{Use fftw library,ON}
OPTION ( USE_FFTW "Look for libfftw (fastest Fourier tranform in the west) and use it if found" OFF )
IF ( USE_FFTW )
  FIND_PATH ( FFTW_INCLUDE_DIR fftw3.h )
  FIND_LIBRARY ( FFTW_LIBRARY NAMES fftw3 fftw3-3 libfftw3-3 )
  FIND_PACKAGE_HANDLE_STANDARD_ARGS( FFTW  DEFAULT_MSG  FFTW_LIBRARY  FFTW_INCLUDE_DIR )
  IF ( FFTW_FOUND )
    INCLUDE_DIRECTORIES ( ${FFTW_INCLUDE_DIR} )
    SET ( FFTW_LIBRARIES ${FFTW_LIBRARY} )
    # fftw needs additional libraries for float and long double precision. Also look for those.
    FIND_LIBRARY ( FFTW_FLOAT_LIBRARY NAMES fftw3f fftw3f-3 libfftw3f-3 )
    IF ( FFTW_FLOAT_LIBRARY )
      SET ( FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_FLOAT_LIBRARY} )
    ENDIF ( )
    FIND_LIBRARY ( FFTW_LONGDOUBLE_LIBRARY NAMES fftw3l fftw3l-3 libfftw3l-3 )
    IF ( FFTW_LONGDOUBLE_LIBRARY )
      SET ( FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_LONGDOUBLE_LIBRARY} )
      ADD_COMPILE_DEFINITIONS ( HAVE_LIB_FFTW_LONGDOUBLE )
    ENDIF ( )
    SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${FFTW_LIBRARIES} )
    ADD_COMPILE_DEFINITIONS ( USE_LIB_FFTW )
  ENDIF ( )
ENDIF ( USE_FFTW )

#! \cmakeoption{Build and use kissfft library,OFF}
OPTION ( BUILD_AND_USE_KISSFFT "Build and use KissFFT (Keep It Simple, Stupid. FFT). If the sources are not available try to download them." OFF )
IF ( USE_FFTW AND BUILD_AND_USE_KISSFFT )
  MESSAGE ( ERROR "USE_FFTW and BUILD_AND_USE_KISSFFT may not both be true!" )
ENDIF ( )
IF ( BUILD_AND_USE_KISSFFT )
  # If the KissFFT source is not in the external directory, try to download it.
  IF ( NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/kissfft-source/kiss_fft.c )
  	SET ( KISSFFT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/kissfft-source )
  	SET ( KISSFFT_ARCHIVE ${KISSFFT_SOURCE_DIR}/kiss_fft130.tar.gz )
  	FILE ( DOWNLOAD http://sourceforge.net/projects/kissfft/files/kissfft/v1_3_0/kiss_fft130.tar.gz/download ${KISSFFT_ARCHIVE} SHOW_PROGRESS )
  	EXECUTE_PROCESS ( COMMAND ${CMAKE_COMMAND} -E tar xzf ${KISSFFT_ARCHIVE} WORKING_DIRECTORY ${KISSFFT_SOURCE_DIR} )

  	# Find out in which subdirectory the source was extracted
    FILE ( GLOB KISSFFT_SUBDIRS RELATIVE ${KISSFFT_SOURCE_DIR} ${KISSFFT_SOURCE_DIR}/*)
    FOREACH ( KISSFFT_SUBDIR ${KISSFFT_SUBDIRS} )
      IF ( IS_DIRECTORY ${KISSFFT_SOURCE_DIR}/${KISSFFT_SUBDIR} )
        SET ( KISSFFT_EXTRACTED_SOURCE_DIR ${KISSFFT_SOURCE_DIR}/${KISSFFT_SUBDIR} )
        BREAK ( )
      ENDIF ( )
    ENDFOREACH ( )

  	# Copy all files into the directory we want to have the source in
	FILE ( GLOB KISSFFT_FILES ${KISSFFT_EXTRACTED_SOURCE_DIR}/*.c ${KISSFFT_EXTRACTED_SOURCE_DIR}/*.h ${KISSFFT_EXTRACTED_SOURCE_DIR}/tools/*.c ${KISSFFT_EXTRACTED_SOURCE_DIR}/tools/*.h ${KISSFFT_EXTRACTED_SOURCE_DIR}/*.hh )
	FOREACH ( KISSFFT_FILE ${KISSFFT_FILES} )
	  FILE ( COPY ${KISSFFT_FILE} DESTINATION ${KISSFFT_SOURCE_DIR} )
	ENDFOREACH ( )

	# Remove the extracted subdirectory and the downloaded archie, we have extracted all the files we need.
    EXECUTE_PROCESS ( COMMAND ${CMAKE_COMMAND} -E remove_directory ${KISSFFT_EXTRACTED_SOURCE_DIR} WORKING_DIRECTORY ${KISSFFT_SOURCE_DIR} )
    FILE ( REMOVE ${KISSFFT_ARCHIVE} )
  ENDIF ( )

  ADD_SUBDIRECTORY( external/kissfft-source )
  INCLUDE_DIRECTORIES ( SYSTEM external/kissfft-source ${CMAKE_CURRENT_BINARY_DIR}/external/kissfft-source )
  LIST ( APPEND SYSTEM_LIBRARIES ${EXTERNAL_KISSFFT_LIBRARIES} )
  ADD_COMPILE_DEFINITIONS ( USE_KISSFFT )
ENDIF ( BUILD_AND_USE_KISSFFT )

# We collect all quoc include directories (e.g. used module directories) in the global variable QUOCMESH_INCLUDE_DIRS.
# This macro adds the current source directory to these inlcude directories.
MACRO ( QUOC_APPEND_CURRENT_DIR_TO_INCLUDES )
  SET ( QUOCMESH_INCLUDE_DIRS ${QUOCMESH_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "" )
ENDMACRO ( QUOC_APPEND_CURRENT_DIR_TO_INCLUDES )

#! \cmakeoption{Use CImg library,ON}
OPTION ( USE_CIMG "" OFF )
IF ( USE_CIMG )
  FIND_PATH ( CIMG_INCLUDE_DIR CImg.h )
  FIND_PACKAGE_HANDLE_STANDARD_ARGS( CIMG  DEFAULT_MSG  CIMG_INCLUDE_DIR )
  IF ( CIMG_FOUND )
    INCLUDE_DIRECTORIES ( ${CIMG_INCLUDE_DIR} )
    ADD_SUBDIRECTORY ( external/cimg )
    ADD_COMPILE_DEFINITIONS ( USE_EXTERNAL_CIMG )
    IF ( TIFF_FOUND )
      ADD_COMPILE_DEFINITIONS ( cimg_use_tiff )
    ENDIF ( )
    IF ( UNIX )
      FIND_PACKAGE ( X11 )
      IF ( X11_FOUND )
        INCLUDE_DIRECTORIES ( ${X11_INCLUDE_DIR} )
        LIST ( APPEND SYSTEM_LIBRARIES ${X11_LIBRARIES} )
      ENDIF ( )
      FIND_LIBRARY ( PTHREAD_LIBRARY pthread )
      IF ( PTHREAD_LIBRARY )
        SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${PTHREAD_LIBRARY} )
      ENDIF ( )
    ENDIF ( )
  ELSE ( CIMG_FOUND )
    MESSAGE ( FATAL_ERROR "Could not find CImg.h." )
  ENDIF ( CIMG_FOUND )
ENDIF ( USE_CIMG )

#! \cmakeoption{Use NetCDF library,OFF}
OPTION ( USE_NETCDF "" OFF )
IF ( USE_NETCDF )
  FIND_PATH ( NETCDF_INCLUDE_DIR netcdf.h )
  FIND_LIBRARY ( NETCDF_LIBRARY NAMES netcdf )
  FIND_PACKAGE_HANDLE_STANDARD_ARGS( NETCDF DEFAULT_MSG NETCDF_LIBRARY NETCDF_INCLUDE_DIR )
  IF ( NETCDF_FOUND )
    INCLUDE_DIRECTORIES ( ${NETCDF_INCLUDE_DIR} )
    LIST ( APPEND SYSTEM_LIBRARIES ${NETCDF_LIBRARY} )
    ADD_COMPILE_DEFINITIONS ( USE_LIB_NETCDF )
  ELSE ( )
    MESSAGE ( FATAL_ERROR "Could not find NetCDF." )
  ENDIF ( )
ENDIF ( )

#! \cmakeoption{Use HDF5 library,OFF}
OPTION ( USE_HDF5 "" OFF )
IF ( USE_HDF5 )
  FIND_PACKAGE ( HDF5 COMPONENTS C CXX HL )
  IF ( HDF5_IS_PARALLEL AND NOT USE_MPI )
  	MESSAGE ( FATAL_ERROR "Parallel HDF5 found, which requires MPI." )
  ENDIF ( )
  IF ( HDF5_FOUND )
    INCLUDE_DIRECTORIES ( ${HDF5_INCLUDE_DIRS} )
    LIST ( APPEND SYSTEM_LIBRARIES ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} )
    ADD_COMPILE_OPTIONS ( ${HDF5_DEFINITIONS} )
    ADD_COMPILE_DEFINITIONS ( USE_LIB_HDF5 )
  ELSE ( )
    MESSAGE ( FATAL_ERROR "Could not find HDF5." )
  ENDIF ( )
ENDIF ( )

#! \cmakeoption{Use OpenCV library,ON}
OPTION ( USE_OPENCV "" OFF )
IF ( USE_OPENCV )
  FIND_PACKAGE ( OpenCV REQUIRED )
  IF ( OpenCV_FOUND )
    INCLUDE_DIRECTORIES ( ${OpenCV_INCLUDE_DIRS} )
    ADD_COMPILE_DEFINITIONS ( USE_OPENCV )
    LIST ( APPEND SYSTEM_LIBRARIES ${OpenCV_LIBS} )
  ENDIF ( )
ENDIF ( )

#! \cmakeoption{Use Gmm++ library,OFF}
OPTION ( USE_GMM "" OFF )
IF ( USE_GMM )
  FIND_PATH ( GMM_INCLUDE_DIR gmm/gmm.h )
  FIND_PACKAGE_HANDLE_STANDARD_ARGS( GMM DEFAULT_MSG GMM_INCLUDE_DIR )
  IF ( GMM_FOUND )
    INCLUDE_DIRECTORIES ( ${GMM_INCLUDE_DIR} )
    ADD_SUBDIRECTORY ( external/gmm )
    ADD_COMPILE_DEFINITIONS ( USE_EXTERNAL_GMM )
  ENDIF ( )
ENDIF ( )

#! \cmakeoption{Use Eigen library,OFF}
OPTION ( USE_EIGEN "" OFF )
IF ( USE_EIGEN )
  FIND_PATH ( EIGEN_INCLUDE_DIR Eigen/Core PATH_SUFFIXES eigen3 )
  FIND_PACKAGE_HANDLE_STANDARD_ARGS( EIGEN DEFAULT_MSG EIGEN_INCLUDE_DIR )
  IF ( EIGEN_FOUND )
    INCLUDE_DIRECTORIES ( ${EIGEN_INCLUDE_DIR} )
    ADD_SUBDIRECTORY ( external/eigen )
    ADD_COMPILE_DEFINITIONS ( USE_EXTERNAL_EIGEN )
  ENDIF ( )
ENDIF ( )

#! \cmakeoption{Use surface mesh,OFF}
OPTION ( USE_SURFACE_MESH "Build surface mesh generator from sources." OFF )
IF ( USE_SURFACE_MESH )
  ADD_SUBDIRECTORY ( external/surfaceMesh )
  INCLUDE_DIRECTORIES ( external/surfaceMesh )
  LIST ( APPEND SYSTEM_LIBRARIES ${EXTERNAL_SURFACE_MESH_LIBRARIES} )
  ADD_COMPILE_DEFINITIONS ( USE_EXTERNAL_SURFACE_MESH )
ENDIF ( USE_SURFACE_MESH )

#! \cmakeoption{Use OpenGL,ON}
OPTION ( USE_OPENGL "" OFF )
IF ( USE_OPENGL )
  FIND_PACKAGE ( OpenGL )
ENDIF ( USE_OPENGL )

IF ( NOT WIN32 )
  #! \cmakeoption{Use GRAPE,ON}
  OPTION ( USE_GRAPE "Look for grape and use it" OFF )
  IF ( USE_GRAPE AND OPENGL_FOUND )
    FIND_PATH ( GRAPE_INCLUDE_DIR grape.h )
    FIND_LIBRARY ( GRAPE_LIBRARY NAMES gr )
    FIND_PACKAGE_HANDLE_STANDARD_ARGS( GRAPE  DEFAULT_MSG  GRAPE_LIBRARY  GRAPE_INCLUDE_DIR )
    IF ( GRAPE_FOUND )
      INCLUDE_DIRECTORIES ( ${GRAPE_INCLUDE_DIR} )
      SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${GRAPE_LIBRARY} ${OPENGL_gl_LIBRARY} )
      ADD_COMPILE_DEFINITIONS ( USE_EXTERNAL_GRAPE )
    ENDIF ( )
  ENDIF ( )
ENDIF ( )

#! \cmakeoption{Use Ipopt,ON}
OPTION ( USE_IPOPT "Look for Ipopt (interior point optimizer) and use it if found. Relies on HSL solver being compiled as shared object" OFF )
IF ( USE_IPOPT )
  IF ( ( BLAS_FOUND AND LAPACK_FOUND ) OR OpenBLAS_FOUND )
    FIND_PATH ( IPOPT_INCLUDE_DIR IpTNLP.hpp PATH_SUFFIXES coin )
    FIND_LIBRARY ( IPOPT_LIBRARY NAMES ipopt PATH_SUFFIXES coin )
    FIND_LIBRARY ( IPOPT_LIBRARY_MUMPS NAMES coinmumps PATH_SUFFIXES coin )
    FIND_LIBRARY ( IPOPT_LIBRARY_METIS NAMES coinmetis PATH_SUFFIXES coin )
    FIND_LIBRARY ( IPOPT_LIBRARY_HSL NAMES coinhsl PATH_SUFFIXES coin/ThirdParty/ )
    FIND_PACKAGE_HANDLE_STANDARD_ARGS( IPOPT  DEFAULT_MSG  IPOPT_LIBRARY IPOPT_INCLUDE_DIR )
    IF ( IPOPT_FOUND )
      INCLUDE_DIRECTORIES ( ${IPOPT_INCLUDE_DIR} )
      SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${IPOPT_LIBRARY} )
      # OpenBLAS is sufficient
      IF ( OpenBLAS_FOUND )
        SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${OpenBLAS_LIBRARIES} )
      ELSE ( )
        SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
      ENDIF ( )
      INCLUDE_DIRECTORIES ( external/ipopt )
      ADD_COMPILE_DEFINITIONS ( USE_EXTERNAL_IPOPT )
      FIND_LIBRARY ( PTHREAD_LIBRARY pthread )
      IF ( PTHREAD_LIBRARY )
        SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${PTHREAD_LIBRARY} )
      ENDIF ( )
      # Link the HSL libraries and/or MUMPS, but print a warning if both are missing.
      IF ( IPOPT_LIBRARY_HSL )
        SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${IPOPT_LIBRARY_HSL} )
      ENDIF ( )
      IF ( IPOPT_LIBRARY_MUMPS )
        SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${IPOPT_LIBRARY_MUMPS} )
        IF ( IPOPT_LIBRARY_METIS )
          SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${IPOPT_LIBRARY_METIS} )
        ENDIF ( )
      ENDIF ( )
      IF ( ( NOT IPOPT_LIBRARY_MUMPS ) AND ( NOT IPOPT_LIBRARY_HSL ) )
        MESSAGE ( WARNING "Ipopt is activated, but no solver library was found. It is still possible to use the Linear Solver Loader." )
      ENDIF ( )
      # ldl not available under Windows?
      IF ( NOT WIN32 )
        SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} -ldl )
      ENDIF ( )
    ENDIF ( )
  ELSE ( )
    MESSAGE ( FATAL_ERROR "Ipopt requires Blas and Lapack or OpenBlas" )
  ENDIF ( )
ENDIF ( )

#! \cmakeoption{Use Metis,ON}
OPTION ( USE_METIS "" OFF )
IF ( USE_METIS )
  FIND_PATH ( METIS_INCLUDE_DIR metis.h )
  FIND_LIBRARY ( METIS_LIBRARY metis )
  FIND_PACKAGE_HANDLE_STANDARD_ARGS ( METIS DEFAULT_MSG METIS_LIBRARY METIS_INCLUDE_DIR )
  IF ( METIS_FOUND )
    INCLUDE_DIRECTORIES ( ${METIS_INCLUDE_DIR} )
    SET ( SYSTEN_LIBRARIES ${SYSTEM_LIBRARIES} ${METIS_LIBRARY} )
    ADD_COMPILE_DEFINITIONS ( ENABLE_METIS USE_EXTERNAL_METIS )
  ENDIF ( METIS_FOUND )
ENDIF ( )

#! \cmakeoption{Use AHMED,ON}
OPTION ( USE_AHMED "Look for AHMED (hierarchcical matrix library) and use it if found" OFF )

#! \cmakeoption{Ahmed interface version,0}
#! Set to 0 to enable the old interface for the AHMED library, referring to subversion revision 174. <br/>
#! Set to 1 to enable the new interface, referring to a more or less up-to-date git version
SET ( AHMED_INTERFACE_VERSION 0 CACHE STRING "Version of the QuocMesh AHMED interface" )
IF ( AHMED_INTERFACE_VERSION EQUAL 0 )
  IF ( USE_AHMED AND BLAS_FOUND AND LAPACK_FOUND )
    FIND_PATH ( AHMED_INCLUDE_DIR aca.h )
    FIND_LIBRARY ( AHMED_LIBRARY NAMES ahmed )
    FIND_PACKAGE_HANDLE_STANDARD_ARGS( AHMED  DEFAULT_MSG  AHMED_LIBRARY  AHMED_INCLUDE_DIR )
    IF ( AHMED_FOUND )
      INCLUDE_DIRECTORIES ( ${AHMED_INCLUDE_DIR} )
      SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${AHMED_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
      ADD_COMPILE_DEFINITIONS ( USE_EXTERNAL_AHMED )
      ADD_COMPILE_DEFINITIONS ( AHMED_INTERFACE_VERSION=0 )
    ENDIF ( AHMED_FOUND )
  ENDIF ( USE_AHMED AND BLAS_FOUND AND LAPACK_FOUND )
ELSEIF ( AHMED_INTERFACE_VERSION EQUAL 1 )
  IF ( USE_AHMED )
    IF ( ( BLAS_FOUND AND LAPACK_FOUND ) OR OpenBLAS_FOUND )
      #! \cmakeoption{INSTALL_AHMED_HEADERS,ON}
      OPTION ( INSTALL_AHMED_HEADERS "Install modified AHMED headers locally" ON )
      IF ( NOT INSTALL_AHMED_HEADERS )
        FIND_PATH ( AHMED_INCLUDE_DIR hmatrix/ACA.h )
      ENDIF ( )
      FIND_LIBRARY ( AHMED_LIBRARY NAMES ahmed )
      FIND_LIBRARY ( AHMED_BASMOD_LIBRARY NAMES basmod PATH_SUFFIXES Basmod )
      FIND_LIBRARY ( AHMED_MATRIX_LIBRARY NAMES matrix PATH_SUFFIXES Matrix )
      FIND_PACKAGE_HANDLE_STANDARD_ARGS( AHMED  DEFAULT_MSG  AHMED_LIBRARY )
      IF ( AHMED_FOUND )
        SET ( AHMED_EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/internal/external/ahmed/" )
        INCLUDE_DIRECTORIES ( ${AHMED_EXTERNAL_DIR} )
        IF ( NOT INSTALL_AHMED_HEADERS )
          INCLUDE_DIRECTORIES ( ${AHMED_INCLUDE_DIR} )
        ENDIF ( )
        SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${AHMED_LIBRARY} ${AHMED_BASMOD_LIBRARY} ${AHMED_MATRIX_LIBRARY} )
        # OpenBLAS is sufficient
        IF ( OpenBLAS_FOUND )
          SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${OpenBLAS_LIBRARIES} )
        ELSE ( )
          SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
        ENDIF ( )
        ADD_COMPILE_DEFINITIONS ( USE_EXTERNAL_AHMED )
        ADD_COMPILE_DEFINITIONS ( AHMED_INTERFACE_VERSION=1 )
        IF ( UNIX )
          GET_FILENAME_COMPONENT ( AHMED_BASE_DIR ${AHMED_LIBRARY} PATH )
          GET_FILENAME_COMPONENT ( AHMED_BASE_DIR ${AHMED_BASE_DIR} PATH )
          GET_FILENAME_COMPONENT ( AHMED_BASE_DIR ${AHMED_BASE_DIR} PATH )
          GET_FILENAME_COMPONENT ( AHMED_BASE_DIR ${AHMED_BASE_DIR} PATH )
          SET ( AHMED_INSTALL_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/internal/external/ahmed/install.sh" )
        ENDIF ( UNIX )
      ENDIF ( AHMED_FOUND )
    ELSE ( )
      MESSAGE ( FATAL_ERROR "Ahmed requires Blas and Lapack or OpenBlas" )
    ENDIF ( ( BLAS_FOUND AND LAPACK_FOUND ) OR OpenBLAS_FOUND )
  ENDIF ( USE_AHMED )
ELSEIF ( AHMED_INTERFACE_VERSION EQUAL 2 )
  IF ( USE_AHMED )
    IF ( ( BLAS_FOUND AND LAPACK_FOUND ) OR OpenBLAS_FOUND )
      #! \cmakeoption{AHMED_INSTALL_HEADERS,ON}
      OPTION ( AHMED_INSTALL_HEADERS "Install modified AHMED headers locally" ON )
      IF ( NOT AHMED_INSTALL_HEADERS )
        FIND_PATH ( AHMED_INCLUDE_DIR hmatrix/ACA.h )
      ENDIF ( )
      FIND_LIBRARY ( AHMED_LIBRARY NAMES ahmed PATH_SUFFIXES hmatrix )
      FIND_LIBRARY ( AHMED_BASMOD_LIBRARY NAMES basmod PATH_SUFFIXES basmod )
      FIND_LIBRARY ( AHMED_MATRIX_LIBRARY NAMES matrix PATH_SUFFIXES matrix )
      FIND_PACKAGE_HANDLE_STANDARD_ARGS( AHMED  DEFAULT_MSG  AHMED_LIBRARY )
      IF ( AHMED_FOUND )
        SET ( AHMED_EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/internal/external/ahmed/" )
        INCLUDE_DIRECTORIES ( ${AHMED_EXTERNAL_DIR} )
        IF ( NOT AHMED_INSTALL_HEADERS )
          INCLUDE_DIRECTORIES ( ${AHMED_INCLUDE_DIR} )
        ENDIF ( )
        SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${AHMED_LIBRARY} ${AHMED_BASMOD_LIBRARY} ${AHMED_MATRIX_LIBRARY} )
        # OpenBLAS is sufficient
        IF ( OpenBLAS_FOUND )
          SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${OpenBLAS_LIBRARIES} )
        ELSE ( )
          SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
        ENDIF ( )
        ADD_COMPILE_DEFINITIONS ( USE_EXTERNAL_AHMED )
        ADD_COMPILE_DEFINITIONS ( AHMED_INTERFACE_VERSION=2 )
        IF ( UNIX )
          GET_FILENAME_COMPONENT ( AHMED_BASE_DIR ${AHMED_LIBRARY} PATH )
          GET_FILENAME_COMPONENT ( AHMED_BASE_DIR ${AHMED_BASE_DIR} PATH )
          GET_FILENAME_COMPONENT ( AHMED_BASE_DIR ${AHMED_BASE_DIR} PATH )
          GET_FILENAME_COMPONENT ( AHMED_BASE_DIR ${AHMED_BASE_DIR} PATH )
          SET ( AHMED_INSTALL_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/internal/external/ahmed/install.sh" )
        ENDIF ( UNIX )
      ENDIF ( AHMED_FOUND )
    ELSE ( )
      MESSAGE ( FATAL_ERROR "Ahmed requires Blas and Lapack or OpenBlas" )
    ENDIF ( ( BLAS_FOUND AND LAPACK_FOUND ) OR OpenBLAS_FOUND )
  ENDIF ( USE_AHMED )
ELSE ( )
  MESSAGE ( FATAL_ERROR "Unknown AHMED_INTERFACE_VERSION" )
ENDIF ( )

#! \cmakeoption{Use SuiteSparse,ON}
OPTION ( USE_SUITESPARSE "" OFF )
IF ( USE_SUITESPARSE AND BLAS_FOUND AND LAPACK_FOUND )
  IF ( NOT WIN32 )
    FIND_PACKAGE ( SUITESPARSE )
    IF ( SUITESPARSE_FOUND OR SuiteSparse_FOUND )
      IF ( SUITESPARSE_FOUND )
        INCLUDE_DIRECTORIES ( ${SUITESPARSE_INCLUDE_DIRS} )
        LINK_DIRECTORIES ( ${SUITESPARSE_LIBRARY_DIR} )
        SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${SUITESPARSE_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
      ELSE ( )
        INCLUDE_DIRECTORIES ( ${SuiteSparse_INCLUDE_DIRS} )
        SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${SuiteSparse_LIBRARIES} )
      ENDIF ( )
      ADD_COMPILE_DEFINITIONS ( USE_EXTERNAL_SUITESPARSE )
    ELSE ( )
      MESSAGE ( FATAL_ERROR "Could NOT find SUITESPARSE" )
    ENDIF ( )
  ELSE ( )
    # This block links SuiteSparse on Windows
    # Download from: https://github.com/jlblancoc/suitesparse-metis-for-windows
    # Build with MinGW or MSVC and set SuiteSparse_DIR to the install directory
    SET ( SuiteSparse_USE_LAPACK_BLAS ON )
    FIND_PACKAGE ( SuiteSparse NO_MODULE )
    IF ( SuiteSparse_FOUND )
      INCLUDE ( ${USE_SuiteSparse} )
	  SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${SuiteSparse_LIBRARIES} ) # blas and lapack should be included with SuiteSparseWindows
	ADD_COMPILE_DEFINITIONS ( USE_EXTERNAL_SUITESPARSE )
      ELSE ( )
        MESSAGE ( STATUS "Could NOT find SUITESPARSE (hint: did you forget to specify SuiteSparse_DIR?)")
      ENDIF ( )
  ENDIF ( )
ELSE ( )
  IF ( USE_SUITESPARSE AND ( NOT BLAS_FOUND ) )
    MESSAGE ( FATAL_ERROR "Suitesparse requires blas" )
  ENDIF ( )
  IF ( USE_SUITESPARSE AND ( NOT LAPACK_FOUND ) )
    MESSAGE ( FATAL_ERROR "Suitesparse requires lapack" )
  ENDIF ( )
ENDIF ( )

#! \cmakeoption{Use FOX,ON}
OPTION ( USE_FOX "Look for fox-toolkit. If found, set the necessary include paths (including external/fox) and add the fox library." OFF )
IF ( USE_FOX )
  #! \cmakeoption{Build FOX,OFF}
OPTION ( BUILD_FOX "Build fox toolkit from sources." OFF )
  IF ( NOT BUILD_FOX )
    FIND_PATH ( FOX_INCLUDE_DIR fx.h PATH_SUFFIXES fox-1.6 )
    FIND_LIBRARY ( FOX_LIBRARY NAMES FOX-1.6 )
    FIND_PACKAGE_HANDLE_STANDARD_ARGS( FOX  DEFAULT_MSG  FOX_LIBRARY  FOX_INCLUDE_DIR )
    IF ( FOX_FOUND )
      INCLUDE_DIRECTORIES ( ${FOX_INCLUDE_DIR} )
      SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${FOX_LIBRARY} )
    ENDIF ( )
  ELSE ( )
    IF ( NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/fox-source/include/fx.h )
      SET ( FOX_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/fox-source )
      SET ( FOX_ARCHIVE ${FOX_SOURCE_DIR}/fox.tar.gz )
      FILE ( DOWNLOAD ftp://ftp.fox-toolkit.org/pub/fox-1.6.50.tar.gz ${FOX_ARCHIVE} SHOW_PROGRESS )
      EXECUTE_PROCESS ( COMMAND ${CMAKE_COMMAND} -E tar xzf ${FOX_ARCHIVE} WORKING_DIRECTORY ${FOX_SOURCE_DIR} )

      SET ( FOX_EXTRACTED_SOURCE_DIR ${FOX_SOURCE_DIR}/fox-1.6.50 )
      FILE ( MAKE_DIRECTORY ${FOX_SOURCE_DIR}/include )
      FILE ( MAKE_DIRECTORY ${FOX_SOURCE_DIR}/src )

      FILE ( GLOB FOX_HEADERS ${FOX_EXTRACTED_SOURCE_DIR}/include/*.h )
      FOREACH ( FOX_HEADER ${FOX_HEADERS} )
        FILE ( COPY ${FOX_HEADER} DESTINATION ${FOX_SOURCE_DIR}/include )
      ENDFOREACH ( )
      FILE ( GLOB FOX_SOURCES ${FOX_EXTRACTED_SOURCE_DIR}/src/*.cpp ${FOX_EXTRACTED_SOURCE_DIR}/src/*.h ${FOX_EXTRACTED_SOURCE_DIR}/src/*.xbm )
      FOREACH ( FOX_SOURCE ${FOX_SOURCES} )
        FILE ( COPY ${FOX_SOURCE} DESTINATION ${FOX_SOURCE_DIR}/src )
      ENDFOREACH ( )

      # Remove the extracted subdirectory and the downloaded archive, we have extracted all the files we need.
      EXECUTE_PROCESS ( COMMAND ${CMAKE_COMMAND} -E remove_directory ${FOX_EXTRACTED_SOURCE_DIR} WORKING_DIRECTORY ${FOX_SOURCE_DIR} )
      FILE ( REMOVE ${FOX_ARCHIVE} )
    ENDIF ( )

    ADD_SUBDIRECTORY ( external/fox-source )
    INCLUDE_DIRECTORIES ( ${EXTERNAL_FOX_INCLUDE_DIR} )
    SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${EXTERNAL_FOX_LIBRARIES} )
    SET ( FOX_FOUND 1 )
  ENDIF ( )
  IF ( FOX_FOUND OR BUILD_FOX )
    IF ( WIN32 )
      SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} opengl32 wsock32 )
    ENDIF ( )
    INCLUDE_DIRECTORIES ( external/fox )
    ADD_COMPILE_DEFINITIONS ( USE_EXTERNAL_FOX )
  ENDIF ( )
ENDIF ( USE_FOX )

# If desired, look for Qt4 and use it. Needs to be done before looking for VTK since VTK may use Qt.
#! \cmakeoption{Use Qt4,ON}
OPTION ( USE_QT "" OFF )
IF ( USE_QT )
  FIND_PACKAGE ( Qt4 REQUIRED QtCore QtGui QtXml )
  IF ( QT_FOUND )
    IF ( WIN32 )
      SET ( QT_USE_QTMAIN TRUE )
    ENDIF ( )
    INCLUDE ( ${QT_USE_FILE} )
    SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${QT_LIBRARIES} )
  ENDIF ( QT_FOUND )
ENDIF ( USE_QT )

#! \cmakeoption{Use Qt5,OFF}
OPTION ( USE_QT5 "If desired, look for Qt5 (but not for Qt5 and Qt4) and use it." OFF )
IF ( USE_QT AND USE_QT5 )
  MESSAGE ( ERROR "USE_QT and USE_QT5 may not both be true!" )
ENDIF ( )
IF ( USE_QT5 )
  FIND_PACKAGE ( Qt5 REQUIRED Core Gui Widgets PrintSupport )
ENDIF ( )

#! \cmakeoption{Use QCustomPlot,OFF}
OPTION ( USE_QCUSTOMPLOT "Build QCustomPlot from sources as lib." OFF )
IF ( USE_QCUSTOMPLOT )
  ADD_SUBDIRECTORY( external/qcustomplot-source )
  INCLUDE_DIRECTORIES ( ${EXTERNAL_QCUSTOMPLOT_INCLUDE_DIR} )
  LIST ( APPEND SYSTEM_LIBRARIES ${EXTERNAL_QCUSTOMPLOT_LIBRARIES} )
  ADD_COMPILE_DEFINITIONS ( USE_QCUSTOMPLOT )
ENDIF ( )

# Optional because vtkpng was not found for some vtk installations
#! \cmakeoption{Use VTK,ON}
OPTION ( USE_VTK "" OFF )
IF ( USE_VTK )
  FIND_PACKAGE ( VTK )
  IF ( VTK_FOUND )
    INCLUDE (${VTK_USE_FILE})
    INCLUDE_DIRECTORIES ( external/vtk )
    ADD_COMPILE_DEFINITIONS ( USE_EXTERNAL_VTK )
    FIND_LIBRARY ( QVTK_LIBRARY NAMES QVTK PATH_SUFFIXES vtk-5.10 )
    IF ( NOT VTK_LIBRARIES )
      SET ( VTK_LIBRARIES vtkRendering vtkGraphics vtkImaging vtkIO vtkFiltering vtkCommon vtksys vtkDICOMParser vtkpng vtktiff vtkzlib vtkjpeg vtkexpat vtkftgl vtkfreetype vtkVolumeRendering vtkRendering vtkHybrid vtkmetaio vtkverdict vtkNetCDF vtksqlite vtkexoIIc vtkViews vtkInfovis vtkWidgets vtklibxml2 )
    ENDIF ( )
    IF ( QT_FOUND AND QVTK_LIBRARY )
      SET ( VTK_LIBRARIES ${VTK_LIBRARIES} QVTK )
    ENDIF ( )

    # VC++ needs different static libraries for debug and release mode.
    IF ( MSVC )
      FOREACH ( VTK_LIB ${VTK_LIBRARIES} )
        SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} debug ${VTK_LIB}-debug optimized ${VTK_LIB} )
      ENDFOREACH ( )
    ELSE ( )
      SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${VTK_LIBRARIES} )
    ENDIF ( )
  ENDIF ( )
ENDIF ( USE_VTK )

# If fox-toolkit and VTK are found, use the vtkfox external
IF ( FOX_FOUND AND VTK_FOUND )
  ADD_SUBDIRECTORY ( external/vtkfox/src )
  INCLUDE_DIRECTORIES ( external/vtkfox/src )
  SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} vtkfox )
ENDIF ( )

#! \cmakeoption{Use CUDA,ON}
OPTION ( USE_CUDA "" OFF )
IF ( USE_CUDA )
  FIND_PACKAGE ( CUDA )
  IF ( CUDA_FOUND )
    # Buchner: Since we are using ADD_COMPILE_OPTIONS to add the -std=c++11 switch, we have to add
    # the flag manually to the nvcc invocation. (This could be fixed, in theory, by using
    # CMAKE_CXX_FLAGS instead, however doing this does not prompt g++ to parse the source code
    # as C++11. While the flag does get propagated this way, g++ seems to ignore it, probably
    # due to some position dependance of the -std flag.)
    # Note: The -std flag is only available since CUDA Toolkit 7.0!
    IF ( ( C++_STANDARD STREQUAL "C++11" ) AND ( NOT CUDA_VERSION_STRING VERSION_LESS 7.0 ) )
      LIST ( APPEND CUDA_NVCC_FLAGS -std=c++11 )
    ENDIF()
    IF ( CUDA_VERSION_STRING VERSION_LESS 7.0 )
      # ... then libcusparse needs to be searched manually
      IF ( CUDA_CUDART_LIBRARY AND NOT CUDA_cusparse_LIBRARY )
        # CUDA installs 32 and 64 bit libaries, make sure to use the correct one.
        # Here this is done by using the libary path an automatically detected CUDA library is in.
        GET_FILENAME_COMPONENT ( CUDA_LIBRARY_DIR ${CUDA_CUDART_LIBRARY} PATH )
        FIND_LIBRARY ( CUDA_cusparse_LIBRARY NAMES cusparse PATHS ${CUDA_LIBRARY_DIR} NO_DEFAULT_PATH )
      ENDIF ( )
      IF ( CUDA_cusparse_LIBRARY )
        LIST ( APPEND SYSTEM_LIBRARIES ${CUDA_cusparse_LIBRARY} )
      ENDIF ( )
    ELSE ( )
      IF ( CUDA_cusparse_LIBRARY )
        LIST ( APPEND SYSTEM_LIBRARIES ${CUDA_cusparse_LIBRARY} )
      ENDIF ( )
      IF ( CUDA_nppc_LIBRARY )
        LIST ( APPEND SYSTEM_LIBRARIES ${CUDA_nppc_LIBRARY} )
      ENDIF ( )
      IF ( CUDA_nppi_LIBRARY )
        LIST ( APPEND SYSTEM_LIBRARIES ${CUDA_nppi_LIBRARY} )
      ENDIF ( )
    ENDIF ( )

    # cutil.h needs to be located manually (Example code from an old FindCUDA.cmake version).
    # CUDA 5+ doesn't have cutil.h anymore, but replaces it with several headers like helper_cuda.h in the CUDA sample directory.
    FIND_PATH ( CUDA_CUT_INCLUDE_DIR
                NAMES cutil.h helper_cuda.h
                PATHS ${CUDA_SDK_SEARCH_PATH} ${CUDA_TOOLKIT_ROOT_DIR}/samples
                PATH_SUFFIXES common/inc
                DOC "Location of cutil.h" )
    IF ( NOT CUDA_CUT_INCLUDE_DIR )
      MESSAGE ( FATAL_ERROR "Could not find cutil.h or helper_cuda.h." )
    ENDIF ( )
    INCLUDE_DIRECTORIES ( ${CUDA_CUT_INCLUDE_DIR} )
    ADD_COMPILE_DEFINITIONS ( USE_EXTERNAL_CUDA )

    SET ( CUDA_CUSTOM_OPTIONS "" )
    IF ( CMAKE_COMPILER_IS_GNUCXX )
      # would like to check for old cuda version but version strings do not seem to work
      IF ( ( NOT ( ${QUOC_GCC_VERSION} VERSION_LESS 4.4 ) )  ) # AND ( ${CUDA_VERSION_MAJOR} LESS 2 )
        LIST ( APPEND CUDA_CUSTOM_OPTIONS --compiler-options -fno-inline )
      ENDIF()

      # The CUDA headers are full of old-style casts, don't let nvcc warn about them.
      LIST ( APPEND CUDA_CUSTOM_OPTIONS -Xcompiler -Wno-old-style-cast )
    ENDIF()
  ENDIF ( )
ENDIF ( )

#! \cmakeoption{Use OpenCL,OFF}
OPTION ( USE_OPENCL "" OFF )
IF ( USE_OPENCL )
  # cmake 3.1.0+ supports OpenCL.
  CMAKE_MINIMUM_REQUIRED ( VERSION 3.1 )
  FIND_PACKAGE ( OpenCL REQUIRED )
  INCLUDE_DIRECTORIES ( ${OpenCL_INCLUDE_DIRS} )
  LIST ( APPEND SYSTEM_LIBRARIES ${OpenCL_LIBRARY} )
  #! \cmakeoption{Use OpenCL BLAS,OFF}
  OPTION ( USE_CLBLAS "" ON )
  IF ( USE_CLBLAS )
    FIND_PATH ( CLBLAS_INCLUDE_DIR clBLAS.h )
    FIND_LIBRARY ( CLBLAS_LIBRARY NAMES clBLAS )
    FIND_PACKAGE_HANDLE_STANDARD_ARGS ( CLBLAS DEFAULT_MSG CLBLAS_LIBRARY CLBLAS_INCLUDE_DIR )
    IF ( CLBLAS_FOUND )
      INCLUDE_DIRECTORIES ( ${CLBLAS_INCLUDE_DIR} )
      LIST ( APPEND SYSTEM_LIBRARIES ${CLBLAS_LIBRARY} )
    ENDIF ( )
  ENDIF ( )
ENDIF ( )

#! \cmakeoption{Use OpenMesh,ON}
OPTION ( USE_OPENMESH "" OFF )
IF ( USE_OPENMESH )
  FIND_PATH ( OPENMESH_INCLUDE_DIR OpenMesh/Core/IO/MeshIO.hh )
  FIND_LIBRARY ( OPENMESH_CORE_LIBRARY OpenMesh_Core OpenMeshCore PATH_SUFFIXES lib/OpenMesh/ )
  FIND_PACKAGE_HANDLE_STANDARD_ARGS( OPENMESH  DEFAULT_MSG  OPENMESH_CORE_LIBRARY  OPENMESH_INCLUDE_DIR )
  IF ( OPENMESH_FOUND )
    INCLUDE_DIRECTORIES ( ${OPENMESH_INCLUDE_DIR} )
    SET ( SYSTEM_LIBRARIES ${SYSTEM_LIBRARIES} ${OPENMESH_CORE_LIBRARY} )
    IF ( MSVC )
      ADD_COMPILE_DEFINITIONS ( _USE_MATH_DEFINES )
    ENDIF ( )
    # Add the external which actually is just a header, also set the corresponing define.
    INCLUDE_DIRECTORIES ( external/openmesh )
    ADD_COMPILE_DEFINITIONS ( USE_EXTERNAL_OPENMESH )
  ENDIF ( )
ENDIF ( USE_OPENMESH )

SET ( QUOCMESH_VERSION "Version 1.6" )

#! \cmakeoption{Use mercurial,ON}
OPTION ( USE_MERCURIAL "" OFF )
IF ( USE_MERCURIAL )
  FIND_PACKAGE( Mercurial )
  IF( MERCURIAL_FOUND )
    MERCURIAL_HG_INFO( ${PROJECT_SOURCE_DIR} Quocmesh )
    SET ( QUOCMESH_VERSION "Latest HG changeset: ${Quocmesh_HG_CHANGESET}, ${Quocmesh_HG_DATE}" )
    ADD_COMPILE_DEFINITIONS ( HG_CHANGESET_ID="${Quocmesh_HG_CHANGESET}" )
  ENDIF ( MERCURIAL_FOUND )
ENDIF ( USE_MERCURIAL )

#! \cmakeoption{Use git,ON}
OPTION ( USE_GIT "" OFF )
IF ( USE_GIT )
  IF ( USE_MERCURIAL AND USE_GIT )
    MESSAGE ( ERROR "USE_MERCURIAL and USE_GIT may not both be true!" )
  ENDIF ( )
  FIND_PACKAGE( Git )
  IF( Git_FOUND )
    # Get the current working branch
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE GIT_BRANCH
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Get the latest abbreviated commit hash of the working branch
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE GIT_COMMIT_HASH
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Get the latest abbreviated commit hash of the working branch
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} log -1 --format=%cd
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE GIT_COMMIT_DATE
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    SET ( QUOCMESH_VERSION "Latest git working branch commit: ${GIT_BRANCH} ${GIT_COMMIT_HASH}, ${GIT_COMMIT_DATE}" )
    ADD_COMPILE_DEFINITIONS ( HG_CHANGESET_ID="${GIT_COMMIT_HASH}" )
  ENDIF ( Git_FOUND )
ENDIF ( USE_GIT )

#! \cmakeoption{Use subversion,OFF}
#! \warning use of subversion is deprecated
OPTION ( USE_SUBVERSION "" OFF )
IF ( USE_SUBVERSION )
  FIND_PACKAGE( Subversion )
  IF( Subversion_FOUND )
    EXECUTE_PROCESS( COMMAND ${Subversion_SVN_EXECUTABLE} info ${PROJECT_SOURCE_DIR} RESULT_VARIABLE SVN_INFO_RETURN OUTPUT_QUIET ERROR_QUIET )
    IF( ${SVN_INFO_RETURN} EQUAL 0 )
      Subversion_WC_INFO( ${PROJECT_SOURCE_DIR} Quocmesh )
      SET ( QUOCMESH_VERSION "Subversion revision ${Quocmesh_WC_REVISION}" )
    ENDIF ( )
  ENDIF ( )
ENDIF ( USE_SUBVERSION )

# Include additional cmake code only relevant for internal parts of the library
INCLUDE ( internal/CMakeLists.txt OPTIONAL )

# Find all files in the specified directory with the given suffix, while ignoring OS X
# resource forks. This is necessary under OS X when using a file system that treats
# the resource forks as "._" files.
MACRO ( QUOC_FIND_FILES_WITH_SUFFIX DirectoryName Suffix Files )
  FILE ( GLOB ${Files} ${DirectoryName}*.${Suffix} )
  FILE ( GLOB OSXResForks ${DirectoryName}._*.${Suffix} )
  FOREACH ( OSXResFork ${OSXResForks} )
    LIST ( REMOVE_ITEM ${Files} ${OSXResFork} )
  ENDFOREACH ( )
ENDMACRO ( QUOC_FIND_FILES_WITH_SUFFIX )

# Letting the executables under VC++, Xcode and Linux end up at the same place needs some special care.
# The custom macro QUOC_ADD_EXECUTABLE encapsulates this.
MACRO ( QUOC_ADD_EXECUTABLE ExecutableName ExecutableInputFiles )
  IF ( MSVC OR XCODE )
    SET ( CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR} )
    SET ( CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_BINARY_DIR} )
  ENDIF ( )

  # Handle that some of the input files may be cu files.
  IF ( CUDA_FOUND )
    CUDA_ADD_EXECUTABLE( ${ExecutableName} ${ExecutableInputFiles}  OPTIONS ${CUDA_CUSTOM_OPTIONS} )
  ELSE ( )
    ADD_EXECUTABLE( ${ExecutableName} ${ExecutableInputFiles} )
  ENDIF ( )

  IF ( MSVC OR XCODE )
    SET_TARGET_PROPERTIES( ${ExecutableName} PROPERTIES DEBUG_OUTPUT_NAME ${ExecutableName}d )
  ENDIF ( )
ENDMACRO ( QUOC_ADD_EXECUTABLE )

# Add a test that just calls a target.
MACRO ( QUOC_ADD_TEST TestTargetName )
  ADD_TEST ( NAME ${TestTargetName} COMMAND $<TARGET_FILE_DIR:${TestTargetName}>/$<TARGET_FILE_NAME:${TestTargetName}> )
  # Make sure that aol::callSystemPauseIfNecessaryOnPlatform doesn't stop the test from exiting automatically.
  SET_PROPERTY ( TEST ${TestTargetName} PROPERTY ENVIRONMENT "QUOC_NO_SYSTEM_PAUSE=ON" )
ENDMACRO ( QUOC_ADD_TEST )

# Writes the target name to the variable with the name specified by TargetName.
MACRO ( QUOC_GEN_TARGETNAME MainFile TargetName )
  STRING ( REPLACE ${CMAKE_SOURCE_DIR}/ "" ${TargetName} ${MainFile} )
  STRING ( REPLACE ${CMAKE_BINARY_DIR}/ "" ${TargetName} ${TARGETNAME} )
  STRING ( REPLACE ".cpp" "" ${TargetName} ${TARGETNAME} )
  STRING ( REPLACE "/" "_" ${TargetName} ${TARGETNAME} )
ENDMACRO ( QUOC_GEN_TARGETNAME )

# Add a test that just calls an automatically generated target from the name of the cpp main file with or without the .cpp suffix.
MACRO ( QUOC_ADD_TEST_FROM_CPP_BASENAME TargetCPPBaseName )
  QUOC_GEN_TARGETNAME ( ${CMAKE_CURRENT_SOURCE_DIR}/${TargetCPPBaseName} TARGETNAME )
  QUOC_ADD_TEST ( ${TARGETNAME} )
ENDMACRO ( QUOC_ADD_TEST_FROM_CPP_BASENAME )

# Add a test that depends on a target that can not be generated from the name of the main file (scripts etc.)
MACRO ( QUOC_ADD_TEST_WITH_DEPENDENCY TestCommand TargetCPPBaseName )
  QUOC_GEN_TARGETNAME ( ${CMAKE_CURRENT_SOURCE_DIR}/${TargetCPPBaseName} TARGETNAME )
  # TARGETNAME is just used as a name here. Could be something else, but this seems consistent.
  ADD_TEST ( NAME ${TARGETNAME} COMMAND ${TestCommand} ${ARGN} )
ENDMACRO ( QUOC_ADD_TEST_WITH_DEPENDENCY )

# Target to collect benchmarks.
ADD_CUSTOM_TARGET ( bench )
# Special target to initialize the benchmark-results file before any of the benchmarks are run.
ADD_CUSTOM_TARGET ( prebench )
# collect current flags and definitions
STRING ( TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}"        NAME_CXX_FLAGS_CONFIG )
STRING ( TOUPPER "CMAKE_EXE_LINKER_FLAGS_${CMAKE_BUILD_TYPE}" NAME_EXE_LINKER_FLAGS_CONFIG )
# this is the modern way to get compile definitions
GET_DIRECTORY_PROPERTY ( EFFECTIVE_DEFINITIONS        COMPILE_DEFINITIONS )
GET_DIRECTORY_PROPERTY ( EFFECTIVE_COMPILE_OPTIONS    COMPILE_OPTIONS )
ADD_CUSTOM_COMMAND ( TARGET prebench PRE_BUILD
                     COMMAND ${CMAKE_COMMAND} -DQUOC_HOSTNAME=${QUOC_HOSTNAME}
                                              -DQUOCMESH_VERSION=${QUOCMESH_VERSION}
                                              -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                                              -DEFFECTIVE_CXX_FLAGS=${CMAKE_CXX_FLAGS}\ "${EFFECTIVE_COMPILE_OPTIONS}"\ ${${NAME_CXX_FLAGS_CONFIG}}
                                              -DEFFECTIVE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}\ ${${NAME_EXE_LINKER_FLAGS_CONFIG}}
                                              -DEFFECTIVE_DEFINITIONS="${EFFECTIVE_DEFINITIONS}"
                                              -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                              -DMSYS=${MSYS}
                                              -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/InitBenchFile.cmake )

# Add a benchmark to the bench target.
MACRO ( QUOC_ADD_BENCH BenchName )
  SET ( BenchTargetName ${BenchName}-bench )
  SET ( QUOC_CUR_DIR_PREFIX "" )
  IF ( UNIX )
    SET ( QUOC_CUR_DIR_PREFIX "./" )
  ENDIF ( )

  ADD_CUSTOM_TARGET ( ${BenchTargetName} ${QUOC_CUR_DIR_PREFIX}${BenchName} bench file ${CMAKE_BINARY_DIR}/benchmark-results-${QUOC_HOSTNAME}.txt )
  ADD_DEPENDENCIES ( bench ${BenchTargetName} )
  ADD_DEPENDENCIES ( ${BenchTargetName} prebench )
ENDMACRO ( )

SET ( QUOCMESH_LIBRARIES ${SYSTEM_LIBRARIES} )

# Let the user configure which modules and projects shall be compiled with a local cmake.selection file.
INCLUDE ( ${CMAKE_CURRENT_BINARY_DIR}/cmake.selection )
IF ( GENERATE_INCLUDE_TEST )
  SET ( SELECTED_STDPROJECTS ${SELECTED_STDPROJECTS} ${CMAKE_BINARY_DIR}/selfTest/includeTest )
ENDIF ( )

IF ( WIN32 )
  #! \cmakeoption{Dynamic linking,OFF (Windows)\, ON (else)}
  OPTION ( DYNAMIC_LINKING "enable dynamic linking" OFF )
ELSE ( )
  OPTION ( DYNAMIC_LINKING "enable dynamic linking" ON )
ENDIF ( )
IF ( DYNAMIC_LINKING )
  # Since the modules depend on each other, we will just collect all sources of all modules to build one big shared library.
  SET ( LIB_SRCS )
ENDIF ( )

# Before building the modules we need to collect all module dirs and add them to the include dirs.
# This is necessary because modules include each other.
FOREACH ( SELECTED_MODULE ${SELECTED_MODULES} )
  # Add the module directory to the include dirs.
  SET ( QUOCMESH_INCLUDE_DIRS ${QUOCMESH_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/${SELECTED_MODULE} CACHE INTERNAL "" )
ENDFOREACH ( SELECTED_MODULE )
# Every time we run cmake, all selected module directories are appended to QUOCMESH_INCLUDE_DIRS,
# so we need to remove the duplicates.
LIST ( REMOVE_DUPLICATES QUOCMESH_INCLUDE_DIRS )
# Make sure QUOCMESH_INCLUDE_DIRS is cleared in the cache.
SET ( QUOCMESH_INCLUDE_DIRS ${QUOCMESH_INCLUDE_DIRS} CACHE INTERNAL "" )
INCLUDE_DIRECTORIES ( ${QUOCMESH_INCLUDE_DIRS} )

FOREACH ( SELECTED_MODULE ${SELECTED_MODULES} )
  # Find all headers and sources in the module directory and create a lib from those files.
  FILE ( GLOB HDRS ${SELECTED_MODULE}/*.h )
  FILE ( GLOB SRCS ${SELECTED_MODULE}/*.cpp )
  # If we are using CUDA, also compile all cu files.
  IF ( CUDA_FOUND )
    FILE ( GLOB CUDA_SRCS ${SELECTED_MODULE}/*.cu )
    LIST ( APPEND SRCS ${CUDA_SRCS} )
  ENDIF ( )
  IF ( QT_FOUND OR Qt5_FOUND )
    # Use uic on all .ui files.
    FILE ( GLOB QT_UIS ${SELECTED_MODULE}/*.ui ${SELECTED_MODULE}/ui/*.ui )
    SET ( QT_UI_HDRS "" )
    # QTX_WRAP_UI places generated header files into CMAKE_CURRENT_BINARY_DIR
    # thus we temporarily append SELECTED_MODULE to CMAKE_CURRENT_BINARY_DIR and include it (if required)
    SET ( CMAKE_CURRENT_BINARY_DIR_TMP ${CMAKE_CURRENT_BINARY_DIR} )
    SET ( CMAKE_CURRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR_TMP}/${SELECTED_MODULE} )
    IF ( QT_FOUND )
      QT4_WRAP_UI ( QT_UI_HDRS ${QT_UIS} )
    ELSE ( )
      QT5_WRAP_UI ( QT_UI_HDRS ${QT_UIS} )
    ENDIF ( )
    IF ( NOT QT_UI_HDRS STREQUAL "" )
	  INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_BINARY_DIR} )
	ENDIF ( )
    SET ( CMAKE_CURRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR_TMP} )
    # Use moc on all .hpp files.
    FILE ( GLOB MOC_HDRS ${SELECTED_MODULE}/*.hpp )
    SET ( MPC_SRCS "" )
    IF ( QT_FOUND )
      QT4_WRAP_CPP ( MPC_SRCS ${MOC_HDRS} )
    ELSE ( )
      QT5_WRAP_CPP ( MPC_SRCS ${MOC_HDRS} )
    ENDIF ( )
    LIST ( APPEND SRCS ${MPC_SRCS} ${QT_UI_HDRS} )
  ENDIF ( )
  GET_FILENAME_COMPONENT ( CURRENT_MODULE_BASE_NAME ${SELECTED_MODULE} NAME_WE )
  # Allows modules/xxx and internal/modules/xxx to coexist.
  # Other cases of identical module names in different directories are not supported.
  # If necessary, something similar to QUOC_GEN_TARGETNAME should be done.
  IF ( ${SELECTED_MODULE} MATCHES "internal" )
    SET ( CURRENT_MODULE_BASE_NAME ${CURRENT_MODULE_BASE_NAME}_int )
  ENDIF ( )
  IF ( MSVC )
    SET ( CURRENT_MODULE_LIB_NAME lib${CURRENT_MODULE_BASE_NAME} )
  ELSE ( )
    SET ( CURRENT_MODULE_LIB_NAME ${CURRENT_MODULE_BASE_NAME} )
  ENDIF ( )
  IF ( DYNAMIC_LINKING )
    SET ( LIB_SRCS ${LIB_SRCS} ${SRCS} ${HDRS} ${MOC_HDRS} )
    # IDE organization
    SOURCE_GROUP ( "${CURRENT_MODULE_LIB_NAME}\\Souce Files" FILES ${SRCS} )
    SOURCE_GROUP ( "${CURRENT_MODULE_LIB_NAME}\\Header Files" FILES ${HDRS} ${MOC_HDRS} )
  ELSE ( )
    # Don't try to compile a header-only module.
    IF ( NOT ( SRCS STREQUAL "" ) )
      # If we are using CUDA, there may be .cu files among the source.
      IF ( CUDA_FOUND )
        CUDA_ADD_LIBRARY ( ${CURRENT_MODULE_LIB_NAME} STATIC ${SRCS} ${HDRS} ${MOC_HDRS} OPTIONS ${CUDA_CUSTOM_OPTIONS} )
      ELSE ( )
        ADD_LIBRARY ( ${CURRENT_MODULE_LIB_NAME} STATIC ${SRCS} ${HDRS} ${MOC_HDRS} )
      ENDIF ( )

      # The module code possibly relies on some of the system libraries, e.g. bz2.
      TARGET_LINK_LIBRARIES( ${CURRENT_MODULE_LIB_NAME} ${SYSTEM_LIBRARIES} )
      IF ( Qt5_FOUND )
        TARGET_LINK_LIBRARIES( ${CURRENT_MODULE_LIB_NAME} Qt5::Widgets Qt5::PrintSupport )
      ENDIF ( )

      # Group all modules in a folder to improve the organization in an IDE.
      SET_TARGET_PROPERTIES( ${CURRENT_MODULE_LIB_NAME} PROPERTIES FOLDER modules )

      IF ( CMAKE_COMPILER_IS_GNUCXX )
        SET_TARGET_PROPERTIES( ${CURRENT_MODULE_LIB_NAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY lib )
      ENDIF ( )

      # Add the lib to the libs that are linked to our executables.
      SET ( QUOCMESH_LIBRARIES ${QUOCMESH_LIBRARIES} ${CURRENT_MODULE_LIB_NAME} )
    ELSE ( )
      # If we use an IDE and have a module with only headers, create a target with these headers so that we can access them in the IDE.
      IF ( ( MSVC OR XCODE ) AND ( HDRS OR MOC_HDRS ) )
        ADD_CUSTOM_TARGET ( ${CURRENT_MODULE_LIB_NAME} SOURCES ${HDRS} ${MOC_HDRS} )
        # IDE organization
        SET_TARGET_PROPERTIES( ${CURRENT_MODULE_LIB_NAME} PROPERTIES FOLDER modules )
      ENDIF ( )
    ENDIF ( )
  ENDIF ( )

  # Add a define for each selected module so that the source code can check whether a module is selected.
  STRING ( REPLACE "/" "_" SELECTED_MODULE_DEFINE_NAME ${SELECTED_MODULE} )
  STRING ( TOUPPER ${SELECTED_MODULE_DEFINE_NAME} SELECTED_MODULE_DEFINE_NAME  )
  ADD_COMPILE_DEFINITIONS ( USE_${SELECTED_MODULE_DEFINE_NAME} )

  IF ( GENERATE_INCLUDE_TEST )
    FOREACH ( HDR ${HDRS} )
      GET_FILENAME_COMPONENT ( HDR_BASENAME ${HDR} NAME_WE )
      SET ( HDR_INCLUDE_STRING "#include <${HDR_BASENAME}.h>" )
      SET ( ALL_HDRS_INCLUDE_STRING "${ALL_HDRS_INCLUDE_STRING}${HDR_INCLUDE_STRING}\n" )
      CONFIGURE_FILE ( ${CMAKE_SOURCE_DIR}/cmake.includeTest.cpp.in ${CMAKE_BINARY_DIR}/selfTest/includeTest/include_${HDR_BASENAME}_h.cpp @ONLY )
    ENDFOREACH ( )
  ENDIF ( )
ENDFOREACH ( SELECTED_MODULE )
IF ( GENERATE_INCLUDE_TEST )
  SET ( HDR_INCLUDE_STRING ${ALL_HDRS_INCLUDE_STRING} )
  CONFIGURE_FILE ( ${CMAKE_SOURCE_DIR}/cmake.includeTest.cpp.in ${CMAKE_BINARY_DIR}/selfTest/includeTest/include_all.cpp @ONLY )
ENDIF ( )

IF ( DYNAMIC_LINKING )
  # If we are using CUDA, there may be .cu files among the source.
  IF ( CUDA_FOUND )
    CUDA_ADD_LIBRARY ( quocmesh SHARED ${LIB_SRCS} OPTIONS ${CUDA_CUSTOM_OPTIONS} )
  ELSE ( )
    ADD_LIBRARY ( quocmesh SHARED ${LIB_SRCS} )
  ENDIF ( )
  INSTALL ( TARGETS quocmesh LIBRARY DESTINATION lib )
  TARGET_LINK_LIBRARIES( quocmesh ${SYSTEM_LIBRARIES} )
  IF ( Qt5_FOUND )
    TARGET_LINK_LIBRARIES( quocmesh Qt5::Widgets Qt5::PrintSupport )
  ENDIF ( )
  SET ( QUOCMESH_LIBRARIES ${QUOCMESH_LIBRARIES} quocmesh )
ELSE ( )
  # The quoc libraries depend on each other, specify them twice in order to avoid any dependence on the sorting of the libraries.
  # For some currently unknown reason tools/triangManipConv/LevelsetToTriangMesh.cpp needs this three times.
  SET ( QUOCMESH_LIBRARIES ${QUOCMESH_LIBRARIES} ${QUOCMESH_LIBRARIES} ${QUOCMESH_LIBRARIES} )
ENDIF ( )

# Hack to remove the red coloring of linking console output. cmake >= 3.3 has nice green color.
# ADD_CUSTOM_TARGET( color ALL COMMAND find ${CMAKE_BINARY_DIR} -name build.make -exec perl -pi -e "s/--red/--blue/" "{}" "\;")
# ADD_DEPENDENCIES( quocmesh color )

# Use CMakeLists.project.txt to handle all selected "standard" projects.
FOREACH ( SELECTED_STDPROJECT ${SELECTED_STDPROJECTS} )
  IF ( IS_ABSOLUTE ${SELECTED_STDPROJECT} )
    SET ( SELECTED_STDPROJECT_DIRECTORY ${SELECTED_STDPROJECT} )
  ELSE ( )
    SET ( SELECTED_STDPROJECT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${SELECTED_STDPROJECT} )
  ENDIF ( )
  # We need to create a temporary CMakeLists.txt. Refuse to do so if such a file already exists.
  IF ( EXISTS ${SELECTED_STDPROJECT_DIRECTORY}/CMakeLists.txt AND NOT ( IS_SYMLINK ${SELECTED_STDPROJECT_DIRECTORY}/CMakeLists.txt ) )
    MESSAGE ( FATAL_ERROR "${SELECTED_STDPROJECT} contains CMakeLists.txt!\nRemove this file if you want this directory to be used in SELECTED_STDPROJECTS." )
  ELSE ( )
    IF ( NOT EXISTS ${SELECTED_STDPROJECT_DIRECTORY}/CMakeLists.txt )
      # For some reason cmake won't generate Xcode projects, if we copy CMakeLists.project.txt and delete it directly after
      # adding the subdirectory. As workaround we create symlinks that are not deleted. For Unix based systems the symlink
      # approach could eventually replace the copying.
      IF ( XCODE )
        EXECUTE_PROCESS ( COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/CMakeLists.project.txt ${SELECTED_STDPROJECT_DIRECTORY}/CMakeLists.txt )
      ELSE ( )
        FILE ( COPY CMakeLists.project.txt DESTINATION ${SELECTED_STDPROJECT_DIRECTORY} )
        FILE ( RENAME ${SELECTED_STDPROJECT_DIRECTORY}/CMakeLists.project.txt ${SELECTED_STDPROJECT_DIRECTORY}/CMakeLists.txt )
      ENDIF ( )
    ENDIF ( )
    # If the path to the selected project is absolute, we have specify the binary dir manually.
    IF ( IS_ABSOLUTE ${SELECTED_STDPROJECT} )
      ADD_SUBDIRECTORY ( ${SELECTED_STDPROJECT} ${SELECTED_STDPROJECT}/build )
    ELSE ( )
      ADD_SUBDIRECTORY ( ${SELECTED_STDPROJECT} )
    ENDIF ( )
    # Leave the symlinks alone. Otherwise running something other than Xcode makes an existing
    # Xcode build dir unusable till cmake is manually run for Xcode again to restore the symlinks.
    IF ( NOT XCODE AND NOT ( IS_SYMLINK ${SELECTED_STDPROJECT_DIRECTORY}/CMakeLists.txt ) )
      FILE ( REMOVE ${SELECTED_STDPROJECT_DIRECTORY}/CMakeLists.txt )
    ENDIF ( )
  ENDIF ( )
ENDFOREACH ( SELECTED_STDPROJECT )

FOREACH ( SELECTED_PROJECT ${SELECTED_PROJECTS} )
  ADD_SUBDIRECTORY ( ${SELECTED_PROJECT} )
ENDFOREACH ( SELECTED_PROJECT )


IF ( USE_AHMED AND AHMED_FOUND AND NOT ( AHMED_INTERFACE_VERSION EQUAL 0 ) AND AHMED_INSTALL_HEADERS )
  MESSAGE ( STATUS "Using ${AHMED_BASE_DIR} as AHMED base directory" )
  ADD_CUSTOM_TARGET ( ahmed_install_headers ALL COMMAND "/bin/sh" ${AHMED_INSTALL_SCRIPT} ${AHMED_BASE_DIR} ${AHMED_LIBRARY} WORKING_DIRECTORY ${AHMED_EXTERNAL_DIR} )
  ADD_DEPENDENCIES ( quocmesh ahmed_install_headers )
ENDIF ( USE_AHMED AND AHMED_FOUND AND NOT ( AHMED_INTERFACE_VERSION EQUAL 0 ) AND AHMED_INSTALL_HEADERS )

#! \cmakeoption{Use doxygen,ON}
OPTION ( USE_DOXYGEN "" OFF )
IF ( USE_DOXYGEN )
  FIND_PACKAGE ( Doxygen )
  FIND_PACKAGE ( LATEX )
  IF ( DOXYGEN_FOUND AND PDFLATEX_COMPILER AND PYTHONINTERP_FOUND )
    #! \cmakeoption{Make Doxygen use MathJax to display Latex formulas instead of generating PNGs. Looks much better, in particular, on high-DPI displays.,OFF}
    OPTION ( USE_MATHJAX_FOR_DOXYGEN "" OFF )
      #! \cmakeoption{Relative (to doc/lib, doc/all) path of MathJax,""}
      SET ( DOXYGEN_MATHJAX_RELPATH "" CACHE STRING "relative (to doc/lib, doc/all) path of mathjax")
    IF ( USE_MATHJAX_FOR_DOXYGEN )
      IF ( DOXYGEN_MATHJAX_RELPATH STREQUAL "" )
	MESSAGE ( FATAL_ERROR "USE_MATHJAX_FOR_DOXYGEN is ON, but no DOXYGEN_MATHJAX_RELPATH has been set!" )
      ENDIF ( DOXYGEN_MATHJAX_RELPATH STREQUAL "" )
      SET( DOXYGEN_USE_MATHJAX "YES" )
    ELSE ( )
      SET( DOXYGEN_USE_MATHJAX "NO" )
    ENDIF ( )
    # doclib
    FILE ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/doc DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ )
    IF ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/internal/doc )
      FILE ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/internal/doc DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ )
    ENDIF ( )
    MAKE_DIRECTORY ( ${CMAKE_CURRENT_BINARY_DIR}/doc/lib )
    # input files to scan for documentation
    FOREACH ( SELECTED_MODULE ${SELECTED_MODULES} )
      # add module directory to the input file list
      SET ( DOXYGEN_INPUT ${DOXYGEN_INPUT} ${CMAKE_CURRENT_SOURCE_DIR}/${SELECTED_MODULE} )
    ENDFOREACH ( SELECTED_MODULE )
    # add doc folder to the input file list
    SET ( DOXYGEN_INPUT ${DOXYGEN_INPUT} ${CMAKE_CURRENT_SOURCE_DIR}/internal/doc )
    # if the cmake doxygen filter is available, use it for CMakeLists.txt files
    IF ( USE_CMAKE_DOXYGEN_FILTER )
      SET ( DOXYGEN_FILTERS "*.txt=${CMAKE_CURRENT_BINARY_DIR}/internal/external/cmakeDoxygenFilter/CMakeDoxygenFilter" )
      SET ( DOXYGEN_INPUT ${DOXYGEN_INPUT} ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt )
    ENDIF ( USE_CMAKE_DOXYGEN_FILTER )
    STRING( REGEX REPLACE ";" " " DOXYGEN_INPUT_LIST "${DOXYGEN_INPUT}" )
    # preprocessor defines
    GET_PROPERTY( DEFS DIRECTORY PROPERTY COMPILE_DEFINITIONS )
    STRING( REGEX REPLACE ";" " " DOXYGEN_PREDEFINED "${DEFS}" )
    # other stuff for doxygen
    SET( DOXYGEN_ENABLED_SECTIONS "DOC_LIB_ONLY" )
    SET( DOXYGEN_HTML_OUTPUT "lib" )
    SET( DOXYGEN_MAN_OUTPUT "lib/man" )
    # prepare doxygen configuration file
    CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen.conf.cmake ${CMAKE_CURRENT_BINARY_DIR}/doc/lib/doxygen.conf )
    # add target
    ADD_CUSTOM_TARGET ( doclib
                                ${PDFLATEX_COMPILER}  manual.tex
                        COMMAND ${PDFLATEX_COMPILER}  manual.tex
                        COMMAND ${PDFLATEX_COMPILER}  manual.tex
                        COMMAND ${PDFLATEX_COMPILER}  InputQuocConventions.tex
                        COMMAND ${PDFLATEX_COMPILER}  InputQuocConventions.tex
                        COMMAND ${PDFLATEX_COMPILER}  InputQuocConventions.tex
                        COMMAND ${DOXYGEN_EXECUTABLE} lib/doxygen.conf
                        WORKING_DIRECTORY doc
    )

    # docall
    MAKE_DIRECTORY ( ${CMAKE_CURRENT_BINARY_DIR}/doc/all )
    # additional input files to scan for documentation
    FOREACH ( SELECTED_PROJECT ${SELECTED_PROJECTS} ${SELECTED_STDPROJECTS} )
      # add project directory to the input file list
      SET ( DOXYGEN_INPUT ${DOXYGEN_INPUT} ${CMAKE_CURRENT_SOURCE_DIR}/${SELECTED_PROJECT} )
    ENDFOREACH ( SELECTED_PROJECT )
    STRING( REGEX REPLACE ";" " " DOXYGEN_INPUT_LIST "${DOXYGEN_INPUT}" )
    # other stuff for doxygen
    SET( DOXYGEN_ENABLED_SECTIONS "" )
    SET( DOXYGEN_HTML_OUTPUT "all" )
    SET( DOXYGEN_MAN_OUTPUT "all/man" )
    # prepare doxygen configuration file
    CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen.conf.cmake ${CMAKE_CURRENT_BINARY_DIR}/doc/all/doxygen.conf )
    # add target
    ADD_CUSTOM_TARGET ( docall
                                ${PDFLATEX_COMPILER}  manual.tex
                        COMMAND ${PDFLATEX_COMPILER}  manual.tex
                        COMMAND ${PDFLATEX_COMPILER}  manual.tex
                        COMMAND ${PDFLATEX_COMPILER}  InputQuocConventions.tex
                        COMMAND ${PDFLATEX_COMPILER}  InputQuocConventions.tex
                        COMMAND ${PDFLATEX_COMPILER}  InputQuocConventions.tex
                        COMMAND ${DOXYGEN_EXECUTABLE} all/doxygen.conf
                        WORKING_DIRECTORY doc
    )

    # if the cmake doxygen filter is enabled it has to build before the doc* targets are executed
    IF ( USE_CMAKE_DOXYGEN_FILTER )
      ADD_DEPENDENCIES( doclib CMakeDoxygenFilter )
      ADD_DEPENDENCIES( docall CMakeDoxygenFilter )
    ENDIF ()

    # Table of projects generation is only available for linux
    IF ( NOT WIN32 )
      ADD_CUSTOM_COMMAND( TARGET docall POST_BUILD COMMAND bash ARGS ${CMAKE_CURRENT_SOURCE_DIR}/util/makeTableOfProjects.sh ${SELECTED_PROJECTS} ${SELECTED_STDPROJECTS} )
    ENDIF ( )

    # cleanup folder doc on make clean
    SET_PROPERTY( DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES doc )
  ELSE ( )
    SET ( WARNING_TEXT "" )
    IF ( NOT DOXYGEN_FOUND  )
      SET ( WARNING_TEXT "${WARNING_TEXT}Doxygen not available, " )
    ENDIF ( )
    IF ( NOT PDFLATEX_COMPILER )
      SET ( WARNING_TEXT "${WARNING_TEXT}pdflatex compiler not available, " )
    ENDIF ( )
    IF ( NOT PYTHONINTERP_FOUND )
      SET ( WARNING_TEXT "${WARNING_TEXT}python interpreter not available, " )
    ENDIF ( )
    MESSAGE ( WARNING "${WARNING_TEXT}not generating doc targets." )
  ENDIF ( )
ENDIF ( )
