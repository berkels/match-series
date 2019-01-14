# Disable the VC++ warnings that occur in the default configuration when compiling zlib.
IF ( MSVC )
  ADD_DEFINITIONS ( /wd4996 )
ENDIF ( MSVC )

# Silence all GCC warnings triggered by zlib.
IF ( CMAKE_COMPILER_IS_GNUCXX OR ( CMAKE_CXX_COMPILER_ID STREQUAL "Clang" ) )
  ADD_DEFINITIONS ( -w )
ENDIF ( )

SET ( HDRS
  crc32.h
  deflate.h
  gzguts.h
  inffast.h
  inffixed.h
  inflate.h
  inftrees.h
  trees.h
  zlib.h
  zutil.h
)

SET ( SRCS
  adler32.c
  compress.c
  crc32.c
  deflate.c
  gzclose.c
  gzlib.c
  gzread.c
  gzwrite.c
  infback.c
  inffast.c
  inflate.c
  inftrees.c
  trees.c
  uncompr.c
  zutil.c
)

SET ( EXTERNAL_ZLIB_INCLUDE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}
  CACHE INTERNAL ""
)

ADD_LIBRARY ( zlib STATIC ${SRCS} ${HDRS} )
SET_TARGET_PROPERTIES( zlib PROPERTIES FOLDER external )

SET ( EXTERNAL_ZLIB_LIBRARIES zlib CACHE INTERNAL "" )