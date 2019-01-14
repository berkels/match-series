#!/bin/bash
# Wrapper script for the exceptionTest. If the expected result changes, edit EXPECTED_RESULT below.

if [ $# -ne 2 ]; then
    printf 'Usage: runExceptionTest.sh PATH_TO_QUOCMESH_SOURCE_DIR PATH_TO_QUOCMESH_BUILD_DIR\n'
    exit 2
fi

QUOCMESH_SOURCE_DIR=$1
QUOCMESH_BUILD_DIR=$2

# If the expected result changes, edit this line!
EXPECTED_RESULT="UnimplementedCodeException in bar : file ${QUOCMESH_SOURCE_DIR}/selfTest/aol/exceptionTest.cpp, line 8"

result=$( ${QUOCMESH_BUILD_DIR}/selfTest/aol/exceptionTest 2>&1 )

if [[ "${result}" != "${EXPECTED_RESULT}" ]]
then
    echo "Unexpected result ${result}"
    echo "Should be ${EXPECTED_RESULT}"
    exit 1
fi

exit 0
