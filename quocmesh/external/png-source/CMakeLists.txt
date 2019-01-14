# Tested with libpng 1.2.56 and libpng-1.4.19

# Disable all VC++ warnings triggered by libpng.
IF ( MSVC )
  ADD_DEFINITIONS ( /W0 )
ENDIF ( MSVC )

# Silence all GCC warnings triggered by libpng.
IF ( CMAKE_COMPILER_IS_GNUCXX )
  ADD_DEFINITIONS ( -w )
  REMOVE_DEFINITIONS ( -pedantic )
ENDIF ( CMAKE_COMPILER_IS_GNUCXX )

SET ( HDRS
  png.h
  pngconf.h
)

SET ( SRCS
  png.c
  pngerror.c
  pngget.c
  pngmem.c
  pngpread.c
  pngread.c
  pngrio.c
  pngrtran.c
  pngrutil.c
  pngset.c
  pngtrans.c
  pngwio.c
  pngwrite.c
  pngwtran.c
  pngwutil.c
)

SET ( EXTERNAL_PNG_INCLUDE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}
  CACHE INTERNAL ""
)

ADD_LIBRARY ( png STATIC ${SRCS} ${HDRS} )
SET_TARGET_PROPERTIES( png PROPERTIES FOLDER external )

SET ( EXTERNAL_PNG_LIBRARIES png CACHE INTERNAL "" )