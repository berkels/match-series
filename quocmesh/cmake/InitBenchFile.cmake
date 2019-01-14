FILE ( APPEND ${CMAKE_BINARY_DIR}/benchmark-results-${QUOC_HOSTNAME}.txt
  "--------------------------------------------------------------------------------\n"
  "                               Benchmark results                                \n"
  "--------------------------------------------------------------------------------\n"
  "\n"
  "Computer name          : ${QUOC_HOSTNAME}\n"
)

# Platform dependent stuff ported from the old makefile.
IF ( UNIX )
  IF ( APPLE )
    EXECUTE_PROCESS ( COMMAND sysctl -n machdep.cpu.brand_string
                      OUTPUT_VARIABLE QUOC_CPU_MODEL_NAME )
  ELSE ( )
    EXECUTE_PROCESS ( COMMAND grep "model name" /proc/cpuinfo
                      COMMAND sed "s/model name\\s*:\\s*//"
                      COMMAND uniq
                      OUTPUT_VARIABLE QUOC_CPU_MODEL_NAME )
  ENDIF ( )
  EXECUTE_PROCESS ( COMMAND usleep 30 )
  IF ( APPLE )
    EXECUTE_PROCESS ( COMMAND sysctl hw.cpufrequency # get frequency in Hz
                      COMMAND sed "s/hw.cpufrequency: //"
                      COMMAND awk "{print $1/1000000}" # convert Hz to MHz
                      OUTPUT_VARIABLE QUOC_CPU_SPEED_MHZ )
    EXECUTE_PROCESS ( COMMAND sysctl -a
                      COMMAND grep machdep.cpu
                      COMMAND grep core_count
                      COMMAND sed "s/machdep.cpu.core_count: //"   
                      OUTPUT_VARIABLE QUOC_CPU_NUM_CORES )
    EXECUTE_PROCESS ( COMMAND sysctl -a
                      COMMAND grep hw.memsize:
                      COMMAND sed "s/hw.memsize: //"
                      COMMAND awk "{print $1/1024}" # convert B to KiB
                      OUTPUT_VARIABLE QUOC_MAIN_MEMORY_KB )
  ELSE ( )
    EXECUTE_PROCESS ( COMMAND grep "cpu MHz" /proc/cpuinfo
                      COMMAND sed "s/cpu MHz\\s*:\\s*//"
                      COMMAND sort
                      COMMAND uniq -c
                      COMMAND sed "s/\\([0-9]\\) \\([0-9]\\)/\\1 x \\2/"
                      COMMAND sed "s/$$/,/"
                      COMMAND xargs
                      COMMAND sed "s/,/, /"
                      COMMAND sed "s/, *$$//"
                      OUTPUT_VARIABLE QUOC_CPU_SPEED_MHZ )
    EXECUTE_PROCESS ( COMMAND grep "processor" /proc/cpuinfo
                      COMMAND wc -l
                      OUTPUT_VARIABLE QUOC_CPU_NUM_CORES )
    # Missing "Max. number of threads : "
    EXECUTE_PROCESS ( COMMAND grep "MemTotal" /proc/meminfo
                      COMMAND sed "s/MemTotal\\s*:\\s*//"
                      COMMAND sed "s/\\skB//"
                      OUTPUT_VARIABLE QUOC_MAIN_MEMORY_KB )
  ENDIF ( )
ENDIF ( )

IF ( UNIX OR MSYS )
  EXECUTE_PROCESS ( COMMAND ${CMAKE_CXX_COMPILER} --version
                    COMMAND head -n1
                    OUTPUT_VARIABLE QUOC_COMPILER_VERSION )
  EXECUTE_PROCESS ( COMMAND date
                    OUTPUT_VARIABLE QUOC_DATE )
ELSE ( )
  SET ( QUOC_COMPILER_VERSION "\n" )
  SET ( QUOC_DATE "\n" )
ENDIF ( )

IF ( WIN32 )
  MACRO ( QUOC_READ_WMIC_OUTPUT TypeName PropertyName WmicOutput )
    EXECUTE_PROCESS ( COMMAND wmic ${TypeName} get ${PropertyName}
                      OUTPUT_VARIABLE ${WmicOutput} )
    # Remove first line containing the property name
    STRING ( REGEX REPLACE "${PropertyName}[\r\n\t ]+" "" ${WmicOutput} ${${WmicOutput}} ) 
    # Remove trailing crap.
    STRING ( REGEX REPLACE "[\r\n\t ]+$" "" ${WmicOutput} ${${WmicOutput}} ) 
    # Properly end the string with a new line.
    SET ( ${WmicOutput} ${${WmicOutput}}\n ) 
  ENDMACRO ( )
  
  QUOC_READ_WMIC_OUTPUT ( cpu Name QUOC_CPU_MODEL_NAME )
  QUOC_READ_WMIC_OUTPUT ( cpu MaxClockSpeed QUOC_CPU_SPEED_MHZ )
  QUOC_READ_WMIC_OUTPUT ( cpu NumberOfCores QUOC_CPU_NUM_CORES )
  # ComputerSystem TotalPhysicalMemory is possibly more accurate, but returns the memory in byte
  # which may be too big to store as int32 and thus to convert to KB using MATH ( EXPR ... )
  QUOC_READ_WMIC_OUTPUT ( OS TotalVisibleMemorySize QUOC_MAIN_MEMORY_KB )
ENDIF ( )

STRING ( LENGTH "$ENV{OMP_NUM_THREADS}" QUOC_OMP_NUM_THREADS_STRING_LENGTH )
IF ( QUOC_OMP_NUM_THREADS_STRING_LENGTH GREATER 0 )
  SET ( QUOC_MAX_NUM_THREADS $ENV{OMP_NUM_THREADS} )
ELSE ( )
  SET ( QUOC_MAX_NUM_THREADS "unlimited" )
ENDIF ( )

FILE ( APPEND ${CMAKE_BINARY_DIR}/benchmark-results-${QUOC_HOSTNAME}.txt
  "CPU model              : ${QUOC_CPU_MODEL_NAME}"
  "Clock speed in MHz     : ${QUOC_CPU_SPEED_MHZ}"
  "Number of cores        : ${QUOC_CPU_NUM_CORES}"
  "Max. number of threads : ${QUOC_MAX_NUM_THREADS}\n"
  "Main memory in kB      : ${QUOC_MAIN_MEMORY_KB}"
  "Compiler version       : ${QUOC_COMPILER_VERSION}"
  "Date and time          : ${QUOC_DATE}"
)

# The compilation flags are still incomplete.
FILE ( APPEND ${CMAKE_BINARY_DIR}/benchmark-results-${QUOC_HOSTNAME}.txt
  "Quocmesh revision      : ${QUOCMESH_VERSION}\n"
  "Build type             : ${CMAKE_BUILD_TYPE}\n"
  "Compilation flags      : ${EFFECTIVE_CXX_FLAGS}\n"
  "Linker flags           : ${EFFECTIVE_LINKER_FLAGS}\n"
  "Defines                : ${EFFECTIVE_DEFINITIONS}\n"
  "\n"
  "benchmark program                                         |    nupsi |    wupsi\n"
  "--------------------------------------------------------------------------------\n"
)
