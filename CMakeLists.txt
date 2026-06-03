cmake_minimum_required(VERSION 3.17)

# Policy to address @foo@ variable expansion
if(POLICY CMP0053)
  cmake_policy(SET CMP0053 NEW)
endif()

# Set the project name and basic settings
project(ACHILLES LANGUAGES CXX VERSION 1.0.0)
include(CMake/StandardProjectSettings.cmake)

# RPATH handling
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Set RPATH for the install location
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
  if(APPLE)
    set(CMAKE_INSTALL_RPATH "@executable_path/../lib")
    list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
  elseif(UNIX)
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib;${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
  endif()
ENDIF()

include(GNUInstallDirs)

# Public usage requirements: everything a downstream consumer of Achilles
# needs (C++ standard, include path, public compile definitions). This
# target IS exported as part of the Achilles package.
add_library(achilles_public_headers INTERFACE)
add_library(achilles::PublicHeaders ALIAS achilles_public_headers)
target_compile_features(achilles_public_headers INTERFACE cxx_std_17)
target_include_directories(achilles_public_headers INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

# Internal-only build settings: sanitizers, coverage, PCH, TESTING define.
# Linked PRIVATE by Achilles libraries; never exported.
add_library(achilles_options INTERFACE)

enable_language(Fortran)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Internal-only compiler warnings; never exported.
add_library(achilles_warnings INTERFACE)
include(CMake/CompilerWarnings.cmake)
set_achilles_warnings(achilles_warnings)

# Sanitizer / coverage options if supported by compiler
include(CMake/Sanitizers.cmake)
enable_sanitizers(achilles_options)

# Allow for static analysis options
include(CMake/StaticAnalyzers.cmake)
include(CMake/ProgramOptions.cmake)

# Test if C++ compiler supports math special functions
include(CheckCXXSourceCompiles)

check_cxx_source_compiles("
#include <cmath>
int main() {
    std::sph_bessel(1, 0.0);
    return 0;
}" ACHILLES_HAS_MATH_SPECIAL_FUNCTIONS)
if(ACHILLES_HAS_MATH_SPECIAL_FUNCTIONS)
  # Public: this define gates code visible in installed headers.
  target_compile_definitions(achilles_public_headers INTERFACE ACHILLES_HAS_MATH_SPECIAL_FUNCTIONS)
endif()

option(ENABLE_AUTODIFF "Enable Autodiff" OFF)

# Allow build to be made even if not a git repo
SET(GIT_FAIL_IF_NONZERO_EXIT FALSE)

# Very basic PCH example
if(ACHILLES_ENABLE_PCH)
  # This sets a global PCH parameter, each project will build its own PCH,
  # which is a good idea if any #define's change
  # These should be headers included in many places
  target_precompile_headers(achilles_options INTERFACE <vector> <string> <map>)
endif()

# Check for python development environment
# find_package(Python REQUIRED COMPONENTS Development)

# Add CMake modules to CMAKE_MODULE_PATH
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMake)

# Ensure compiler supports filesystem and needed compiler flags
find_package(Filesystem REQUIRED)

# Find ZLIB to read gzip files
if(ACHILLES_ENABLE_GZIP)
find_package(ZLIB REQUIRED)
endif()

if(ACHILLES_ENABLE_SHERPA)
# Find Sherpa
find_package(SHERPA REQUIRED)
find_package(ZLIB REQUIRED)
endif()

# Find ROOT if needed
if(ACHILLES_ENABLE_ROOT)
  MESSAGE(STATUS "Building with ROOT flux files")
  find_package(ROOT 6 REQUIRED)
  string(REPLACE "." ";" VERSION_LIST ${ROOT_VERSION})
  list(GET VERSION_LIST 1 ROOT_MINOR_VERSION)
  if(${ROOT_MINOR_VERSION} LESS 14)
    # ROOT targets are missing includes and flags in ROOT 6.10 and 6.12
    set_property(TARGET ROOT::Core PROPERTY
                 INTERFACE_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIRS}")

    # Early ROOT does not include the flags required on targets
    add_library(ROOT::Flags_CXX IMPORTED INTERFACE)
  endif()

  if(${ROOT_MINOR_VERSION} LESS 16)
    # ROOT 6.14 and earlier have a spacing bug in the linker flags
    string(REPLACE "-L " "-L" ROOT_EXE_LINKER_FLAGS "${ROOT_EXE_LINKER_FLAGS}")

    # Fix for ROOT_CXX_FLAGS not actually being a CMake list
    separate_arguments(ROOT_CXX_FLAGS)
    set_property(TARGET ROOT::Flags_CXX APPEND PROPERTY
                 INTERFACE_COMPILE_OPTIONS ${ROOT_CXX_FLAGS})

    # Add definitions
    separate_arguments(ROOT_DEFINITIONS)
    foreach(flag ${ROOT_EXE_LINKER_FLAG_LIST})
      # Remove -D or /D if present
      string(REGEX REPLACE [=[^[-//]D]=] "" flag ${flag})
      set_property(TARGET ROOT::Flags
                   APPEND PROPERTY INTERFACE_LINK_LIBRARIES ${flag})
    endforeach()

    # This also fixes a bug in the linker flags
    separate_arguments(ROOT_EXE_LINKER_FLAGS)
    set_property(TARGET ROOT::Flags_CXX APPEND PROPERTY
                 INTERFACE_LINK_LIBRARIES ${ROOT_EXE_LINKER_FLAGS})
  endif()
endif()

# Add dependencies
include(CMake/CPM.cmake)
include(external/CMakeLists.txt)

# Configure system specific variables for loading shared libraries
set(LIBPREFIX ${CMAKE_SHARED_LIBRARY_PREFIX})
set(LIBSUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
set(BUILDPATH ${CMAKE_CURRENT_BINARY_DIR})
set(INSTALLPATH ${CMAKE_INSTALL_PREFIX})
configure_file(
    include/Achilles/System.hh.in
    ${CMAKE_CURRENT_BINARY_DIR}/generated/include/Achilles/System.hh
    @ONLY
)

# Copy over data files
file(COPY ${ACHILLES_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${ACHILLES_SOURCE_DIR}/flux DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${ACHILLES_SOURCE_DIR}/data/default/run.yml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${ACHILLES_SOURCE_DIR}/data/default/FormFactors.yml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${ACHILLES_SOURCE_DIR}/parameters.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${ACHILLES_SOURCE_DIR}/cascade.yml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy over example folder
file(COPY ${ACHILLES_SOURCE_DIR}/examples DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Unzip if GZIP not turned on
if(NOT ACHILLES_ENABLE_GZIP)
  file(GLOB zipfiles "${CMAKE_CURRENT_BINARY_DIR}/data/configurations/*.gz")
  MESSAGE(STATUS "Extracting configurations from zip files")
  foreach(file ${zipfiles})
    execute_process(COMMAND gunzip "${file}"
                    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data/configurations/")
  endforeach()
endif()

# Build Plugin Library
# add_subdirectory(src/Plugins)

# Testing
if(ACHILLES_ENABLE_TESTING)
  if(ACHILLES_COVERAGE)
    target_compile_options(achilles_options INTERFACE --coverage)
    target_link_libraries(achilles_options INTERFACE gcov)
  endif()
  target_compile_definitions(achilles_options INTERFACE TESTING)
  enable_testing()
  message(STATUS "Building Tests.")
  add_subdirectory(test)
endif()

# Main code
add_subdirectory(src)

# Install the public-headers INTERFACE target into the export set so that
# downstream consumers inherit the include path, C++ standard, and public
# compile definitions through find_package(Achilles).
install(TARGETS achilles_public_headers
        EXPORT AchillesTargets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Install the generated System.hh alongside the rest of the public headers.
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/generated/include/Achilles/System.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Achilles)

# Generate and install the package config files.
include(CMakePackageConfigHelpers)

# Capture build-time feature flags so the installed config file can issue
# the right find_dependency() calls (see AchillesConfig.cmake.in).
set(ACHILLES_ENABLE_SHERPA ${ACHILLES_ENABLE_SHERPA})
set(ACHILLES_ENABLE_ROOT   ${ACHILLES_ENABLE_ROOT})
set(ACHILLES_ENABLE_HEPMC3 ${ACHILLES_ENABLE_HEPMC3})
set(ACHILLES_ENABLE_GZIP   ${ACHILLES_ENABLE_GZIP})

configure_package_config_file(
    ${CMAKE_SOURCE_DIR}/CMake/AchillesConfig.cmake.in
    ${CMAKE_BINARY_DIR}/CMake/AchillesConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Achilles
)

write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/CMake/AchillesConfigVersion.cmake
    VERSION ${ACHILLES_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(
    FILES
        ${CMAKE_BINARY_DIR}/CMake/AchillesConfig.cmake
        ${CMAKE_BINARY_DIR}/CMake/AchillesConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Achilles
)

# Export the targets file. Library and executable targets register
# themselves into the AchillesTargets export set in src/Achilles/CMakeLists.txt.
install(EXPORT AchillesTargets
    FILE AchillesTargets.cmake
    NAMESPACE achilles::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Achilles
)

# Install default run card and form factors card
install(
    FILES
        ${CMAKE_BINARY_DIR}/data/default/run.yml
        ${CMAKE_BINARY_DIR}/data/default/FormFactors.yml
    DESTINATION ${CMAKE_INSTALL_PREFIX}
)

# Install the data files to be found
install(DIRECTORY
    ${CMAKE_BINARY_DIR}/data
    ${CMAKE_BINARY_DIR}/flux
    DESTINATION share/Achilles
)

# Install examples into the install directory
install(DIRECTORY
    ${CMAKE_BINARY_DIR}/examples
    DESTINATION ${CMAKE_INSTALL_PREFIX}
)
