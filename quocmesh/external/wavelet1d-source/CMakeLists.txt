# You can download wavelet1d from:
# https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/wavelet1d/wavelet-03.rar
#
# TODO: some changes have to be made automatically to wavelet.cpp

# Silence all GCC warnings triggered by wavelet1d.
IF ( CMAKE_COMPILER_IS_GNUCXX )
  ADD_DEFINITIONS ( -w )
  REMOVE_DEFINITIONS ( -pedantic )
ENDIF ( CMAKE_COMPILER_IS_GNUCXX )

SET ( HDRS
  wavelet.h
)

SET ( SRCS
  wavelet.cpp
)

SET ( EXTERNAL_WAVELET_INCLUDE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}
  CACHE INTERNAL ""
)

ADD_LIBRARY ( wavelet1d STATIC ${SRCS} ${HDRS} )
SET_TARGET_PROPERTIES( wavelet1d PROPERTIES FOLDER external )

SET ( EXTERNAL_WAVELET_LIBRARIES wavelet1d CACHE INTERNAL "" )