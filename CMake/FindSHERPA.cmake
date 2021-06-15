find_library(SHERPA_LIBRARY NAMES SherpaMain HINTS $ENV{SHERPA_ROOT_DIR}/lib/SHERPA-MC ${SHERPA_ROOT_DIR}/lib/SHERPA-MC)
IF(${SHERPA_LIBRARY} MATCHES "SHERPA_LIBRARY-NOTFOUND")
  FIND_PROGRAM(SHERPA_CONFIG_EXECUTABLE NAMES sherpa-config
    HINTS $ENV{SHERPA_ROOT_DIR}/bin ${SHERPA_ROOT_DIR}/bin)
  IF(${SHERPA_CONFIG_EXECUTABLE} MATCHES "SHERPA_CONFIG_EXECUTABLE-NOTFOUND")
    MESSAGE(STATUS "Looking for SHERPA... sherpa-config executable not found")
  ELSE(${SHERPA_CONFIG_EXECUTABLE} MATCHES "SHERPA_CONFIG_EXECUTABLE-NOTFOUND")
    MESSAGE(STATUS "Looking for SHERPA... using sherpa-config executable")
    EXEC_PROGRAM(${SHERPA_CONFIG_EXECUTABLE} ARGS "--prefix" OUTPUT_VARIABLE SHERPA_PREFIX)
    find_library(SHERPA_LIBRARY NAMES sherpa_v1 PATHS ${SHERPA_PREFIX}/lib)
  ENDIF(${SHERPA_CONFIG_EXECUTABLE} MATCHES "SHERPA_CONFIG_EXECUTABLE-NOTFOUND")
ENDIF(${SHERPA_LIBRARY} MATCHES "SHERPA_LIBRARY-NOTFOUND")

find_path(SHERPA_INCLUDE_DIR SHERPA/Main/Sherpa.H
  HINTS $ENV{SHERPA_ROOT_DIR}/include/SHERPA-MC ${SHERPA_ROOT_DIR}/include/SHERPA-MC ${SHERPA_PREFIX}/include/SHERPA-MC)

mark_as_advanced(SHERPA_LIBRARY SHERPA_INCLUDE_DIR)

# handle QUIETLY and REQUIRED arguments and set sherpa_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SHERPA DEFAULT_MSG SHERPA_INCLUDE_DIR SHERPA_LIBRARY)

set(SHERPA_LIBRARIES ${SHERPA_LIBRARY})
get_filename_component(SHERPA_LIBRARY_DIRS ${SHERPA_LIBRARY} PATH)

set(SHERPA_INCLUDE_DIRS ${SHERPA_INCLUDE_DIR})

mark_as_advanced(SHERPA_FOUND)