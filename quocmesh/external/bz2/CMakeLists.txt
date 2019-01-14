# Disable all VC++ warnings triggered by libbz2.
IF ( MSVC )
  ADD_DEFINITIONS ( /W0 )
ENDIF ( MSVC )

# Silence all GCC warnings triggered by libbz2.
IF ( CMAKE_COMPILER_IS_GNUCXX )
  ADD_DEFINITIONS ( -w )
  REMOVE_DEFINITIONS ( -pedantic )
ENDIF ( CMAKE_COMPILER_IS_GNUCXX )

SET ( HDRS
  bzlib.h
  bzlib_private.h
)

SET ( SRCS
  blocksort.c
  bzlib.c
  compress.c
  crctable.c
  decompress.c
  huffman.c
  randtable.c
)

SET ( EXTERNAL_BZIP2_INCLUDE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}
  CACHE INTERNAL ""
)

ADD_LIBRARY ( bzip2 STATIC ${SRCS} ${HDRS} )
SET_TARGET_PROPERTIES( bzip2 PROPERTIES FOLDER external )

SET ( EXTERNAL_BZIP2_LIBRARIES bzip2 CACHE INTERNAL "" )
