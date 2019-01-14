# Silence all GCC warnings triggered by qcustomplot.
IF ( CMAKE_COMPILER_IS_GNUCXX )
  ADD_DEFINITIONS ( -w )
ENDIF ( CMAKE_COMPILER_IS_GNUCXX )

# If the QCustomPlot source is not in the external directory, try to download it.
IF ( NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/qcustomplot.h )
  SET ( QCUSTOMPLOT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
  SET ( QCUSTOMPLOT_ARCHIVE ${QCUSTOMPLOT_SOURCE_DIR}/wavelet1d.rar )
  FILE ( DOWNLOAD http://www.qcustomplot.com/release/1.3.2/QCustomPlot.tar.gz ${QCUSTOMPLOT_ARCHIVE} SHOW_PROGRESS )
  EXECUTE_PROCESS ( COMMAND ${CMAKE_COMMAND} -E tar xzf ${QCUSTOMPLOT_ARCHIVE} WORKING_DIRECTORY ${QCUSTOMPLOT_SOURCE_DIR} )
  	
  FILE ( COPY ${QCUSTOMPLOT_SOURCE_DIR}/qcustomplot/qcustomplot.cpp DESTINATION ${QCUSTOMPLOT_SOURCE_DIR} )
  FILE ( COPY ${QCUSTOMPLOT_SOURCE_DIR}/qcustomplot/qcustomplot.h DESTINATION ${QCUSTOMPLOT_SOURCE_DIR} )
	
  # Remove the downloaded and extracted archive, we have extracted all the files we need.
  EXECUTE_PROCESS ( COMMAND ${CMAKE_COMMAND} -E remove_directory ${QCUSTOMPLOT_SOURCE_DIR}/qcustomplot WORKING_DIRECTORY ${QCUSTOMPLOT_SOURCE_DIR} )
  FILE ( REMOVE ${QCUSTOMPLOT_ARCHIVE} )
ENDIF ( )

SET ( MOC_HDRS
  qcustomplot.h
)

SET ( SRCS
  qcustomplot.cpp
)

SET ( EXTERNAL_QCUSTOMPLOT_INCLUDE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}
  CACHE INTERNAL ""
)

IF ( QT_FOUND )
  QT4_WRAP_CPP ( MPC_SRCS ${MOC_HDRS} )
ELSE ( )
  QT5_WRAP_CPP ( MPC_SRCS ${MOC_HDRS} )
ENDIF ( )

ADD_LIBRARY ( qcustomplot STATIC ${SRCS} ${MPC_SRCS} )
SET_TARGET_PROPERTIES( qcustomplot PROPERTIES FOLDER external )

IF ( Qt5_FOUND )
  TARGET_LINK_LIBRARIES( qcustomplot Qt5::Widgets Qt5::PrintSupport )
ENDIF ( )

SET ( EXTERNAL_QCUSTOMPLOT_LIBRARIES qcustomplot CACHE INTERNAL "" )
