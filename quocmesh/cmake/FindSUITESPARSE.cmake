# /*===========================================================================*\
# *                                                                            *
# *                              OpenFlipper                                   *
# *      Copyright (C) 2001-2011 by Computer Graphics Group, RWTH Aachen       *
# *                           www.openflipper.org                              *
# *                                                                            *
# *--------------------------------------------------------------------------- *
# *  This file is part of OpenFlipper.                                         *
# *                                                                            *
# *  OpenFlipper is free software: you can redistribute it and/or modify       *
# *  it under the terms of the GNU Lesser General Public License as            *
# *  published by the Free Software Foundation, either version 3 of            *
# *  the License, or (at your option) any later version with the               *
# *  following exceptions:                                                     *
# *                                                                            *
# *  If other files instantiate templates or use macros                        *
# *  or inline functions from this file, or you compile this file and          *
# *  link it with other files to produce an executable, this file does         *
# *  not by itself cause the resulting executable to be covered by the         *
# *  GNU Lesser General Public License. This exception does not however        *
# *  invalidate any other reasons why the executable file might be             *
# *  covered by the GNU Lesser General Public License.                         *
# *                                                                            *
# *  OpenFlipper is distributed in the hope that it will be useful,            *
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
# *  GNU Lesser General Public License for more details.                       *
# *                                                                            *
# *  You should have received a copy of the GNU LesserGeneral Public           *
# *  License along with OpenFlipper. If not,                                   *
# *  see <http://www.gnu.org/licenses/>.                                       *
# *                                                                            *
# \*===========================================================================*/
#
# Modified for usage with QuocMesh


# - Try to find SUITESPARSE
# Once done this will define
#  
#  SUITESPARSE_FOUND            - system has SUITESPARSE
#  SUITESPARSE_INCLUDE_DIRS     - the SUITESPARSE include directory
#  SUITESPARSE_LIBRARIES        - Link these to use SUITESPARSE
#  SUITESPARSE_SPQR_LIBRARY     - name of spqr library (necessary due to error in debian package)
#  SUITESPARSE_SPQR_LIBRARY_DIR - name of spqr library (necessary due to error in debian package)
#  SUITESPARSE_LIBRARY_DIR      - Library main directory containing suitesparse libs
#  SUITESPARSE_LIBRARY_DIRS     - all Library directories containing suitesparse libs
#  SUITESPARSE_SPQR_VALID       - automatic identification whether or not spqr package is installed correctly

IF (SUITESPARSE_INCLUDE_DIRS)
  # Already in cache, be silent
  SET(SUITESPARSE_FIND_QUIETLY TRUE)
ENDIF (SUITESPARSE_INCLUDE_DIRS)

IF( WIN32 )
  # Find cholmod part of the suitesparse library collection

  FIND_PATH( CHOLMOD_INCLUDE_DIR cholmod.h
             PATHS "C:\\libs\\win32\\SuiteSparse\\Include"  )

  # Add cholmod include directory to collection include directories
  IF ( CHOLMOD_INCLUDE_DIR )
    LIST ( APPEND SUITESPARSE_INCLUDE_DIRS ${CHOLMOD_INCLUDE_DIR} )
  ENDIF( CHOLMOD_INCLUDE_DIR )


  # find path suitesparse library
  FIND_PATH( SUITESPARSE_LIBRARY_DIRS 
             amd.lib
             PATHS "C:\\libs\\win32\\SuiteSparse\\libs" )

  # if we found the library, add it to the defined libraries
  IF ( SUITESPARSE_LIBRARY_DIRS )
       LIST ( APPEND SUITESPARSE_LIBRARIES optimized;amd;optimized;camd;optimized;ccolamd;optimized;cholmod;optimized;colamd;optimized;metis;optimized;spqr;optimized;umfpack;debug;amdd;debug;camdd;debug;ccolamdd;debug;cholmodd;debug;spqrd;debug;umfpackd;debug;colamdd;debug;metisd;optimized;blas;optimized;libf2c;optimized;lapack;debug;blasd;debug;libf2cd;debug;lapackd )
  ENDIF( SUITESPARSE_LIBRARY_DIRS )  

ELSE( WIN32 )
  IF( APPLE)
    FIND_PATH( CHOLMOD_INCLUDE_DIR cholmod.h
               PATHS  /opt/local/include/ufsparse )

    FIND_PATH( SUITESPARSE_LIBRARY_DIR
               NAMES libcholmod.a 
               PATHS /opt/local/lib
               PATH_SUFFIXES lib )

  ELSE(APPLE)
    FIND_PATH( CHOLMOD_INCLUDE_DIR cholmod.h
               PATHS /usr/local/include 
                     /usr/include 
                     /usr/include/suitesparse/ 
                     ${CMAKE_SOURCE_DIR}/MacOS/Libs/cholmod
               PATH_SUFFIXES cholmod/ CHOLMOD/ )


    FIND_PATH ( SUITESPARSE_LIBRARY_DIR
                NAMES libcholmod.so libcholmod.a
                PATH_SUFFIXES lib lib64 )

ENDIF(APPLE)

 # Add cholmod include directory to collection include directories
IF ( CHOLMOD_INCLUDE_DIR )
  LIST ( APPEND SUITESPARSE_INCLUDE_DIRS ${CHOLMOD_INCLUDE_DIR} )
ENDIF( CHOLMOD_INCLUDE_DIR )


# if we found the library, add it to the defined libraries
IF ( SUITESPARSE_LIBRARY_DIR )
  LIST ( APPEND SUITESPARSE_LIBRARIES amd)
  LIST ( APPEND SUITESPARSE_LIBRARIES btf)
  LIST ( APPEND SUITESPARSE_LIBRARIES camd)
  LIST ( APPEND SUITESPARSE_LIBRARIES ccolamd)
  LIST ( APPEND SUITESPARSE_LIBRARIES cholmod)
  LIST ( APPEND SUITESPARSE_LIBRARIES colamd)
  #LIST ( APPEND SUITESPARSE_LIBRARIES csparse)
  LIST ( APPEND SUITESPARSE_LIBRARIES cxsparse)
  LIST ( APPEND SUITESPARSE_LIBRARIES klu)
  #LIST ( APPEND SUITESPARSE_LIBRARIES spqr)
  LIST ( APPEND SUITESPARSE_LIBRARIES umfpack)

  # from version 4 on an additional config library is needed
  FIND_LIBRARY( SUITESPARSE_CONFIG_LIB NAMES suitesparseconfig PATHS ${SUITESPARSE_LIBRARY_DIR})
  IF ( EXISTS ${SUITESPARSE_CONFIG_LIB} )
    LIST ( APPEND SUITESPARSE_LIBRARIES ${SUITESPARSE_CONFIG_LIB} )
  ENDIF ( EXISTS ${SUITESPARSE_CONFIG_LIB} )

  # Metis and spqr are optional
  FIND_LIBRARY( SUITESPARSE_METIS_LIBRARY
                NAMES metis
                PATHS ${SUITESPARSE_LIBRARY_DIR} )
  IF (SUITESPARSE_METIS_LIBRARY)                   
    LIST ( APPEND SUITESPARSE_LIBRARIES metis)
  ENDIF(SUITESPARSE_METIS_LIBRARY)

  IF(EXISTS  "${CHOLMOD_INCLUDE_DIR}/SuiteSparseQR.hpp")
    SET(SUITESPARSE_SPQR_VALID TRUE CACHE BOOL "SuiteSparseSPQR valid")
  ELSE()
    SET(SUITESPARSE_SPQR_VALID false CACHE BOOL "SuiteSparseSPQR valid")
  ENDIF()

  IF(SUITESPARSE_SPQR_VALID)
    FIND_LIBRARY( SUITESPARSE_SPQR_LIBRARY
                  NAMES spqr
                  PATHS ${SUITESPARSE_LIBRARY_DIR} )
    IF (SUITESPARSE_SPQR_LIBRARY)                 
      LIST ( APPEND SUITESPARSE_LIBRARIES spqr)
    ENDIF (SUITESPARSE_SPQR_LIBRARY)
  ENDIF()
       
ELSE( SUITESPARSE_LIBRARY_DIR )  
  MESSAGE ( WARNING "SUITESPARSE_LIBRARY_DIR not found!" )
ENDIF( SUITESPARSE_LIBRARY_DIR )  

ENDIF( WIN32 )


IF (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)
  IF(WIN32)
    LIST (APPEND SUITESPARSE_INCLUDE_DIRS ${CHOLMOD_INCLUDE_DIR}/../../UFconfig )
  ENDIF(WIN32)
    SET(SUITESPARSE_FOUND TRUE)
ELSE (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)
 SET( SUITESPARSE_FOUND FALSE )
ENDIF (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)

