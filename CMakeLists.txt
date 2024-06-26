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
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
ENDIF()

# Link this 'library' to set the c++ standard / compile-time options requested
# Additionally, link to get include and external dependencies
add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE cxx_std_17)
target_include_directories(project_options INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

enable_language(Fortran)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Link this 'library' to use the warnings specified in CompilerWarnings.cmake
add_library(project_warnings INTERFACE)

# Standard compiler warnings
include(CMake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

# Sanitizer options if supported by compiler
include(CMake/Sanitizers.cmake)
enable_sanitizers(project_options)

# Allow for static analysis options
include(CMake/StaticAnalyzers.cmake)
include(CMake/ProgramOptions.cmake)

option(ENABLE_AUTODIFF "Enable Autodiff" OFF)

# Allow build to be made even if not a git repo
SET(GIT_FAIL_IF_NONZERO_EXIT FALSE)

# Very basic PCH example
if(ACHILLES_ENABLE_PCH)
    # This sets a global PCH parameter, each project will build its own PCH,
    # which is a good idea if any #define's change
    # These should be headers included in many places
    target_precompile_headers(project_options INTERFACE <vector> <string> <map>)
endif()

# Check for python development environment
# find_package(Python REQUIRED COMPONENTS Development)

# Add CMake modules to CMAKE_MODULE_PATH
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMake)

# Ensure compiler supports filesystem and needed compiler flags
find_package(Filesystem REQUIRED)

# Find HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)

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
        foreach(_flag ${ROOT_EXE_LINKER_FLAG_LIST})
                # Remove -D or /D if present
                string(REGEX REPLACE [=[^[-//]D]=] "" _flag ${_flag})
                set_property(TARGET ROOT::Flags APPEND PROPERTY INTERFACE_LINK_LIBRARIES ${_flag})
        endforeach()
        
        # This also fixes a bug in the linker flags
        separate_arguments(ROOT_EXE_LINKER_FLAGS)
        set_property(TARGET ROOT::Flags_CXX APPEND PROPERTY
                     INTERFACE_LINK_LIBRARIES ${ROOT_EXE_LINKER_FLAGS})
    endif()
endif()

# Add dependencies
include(CMake/CPM.cmake)
add_subdirectory(external)

# Configure system specific variables for loading shared libraries
set(LIBPREFIX ${CMAKE_SHARED_LIBRARY_PREFIX})
set(LIBSUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
set(BUILDPATH ${CMAKE_CURRENT_BINARY_DIR})
set(INSTALLPATH ${CMAKE_INSTALL_PREFIX})
configure_file(include/Achilles/System.hh.in ${CMAKE_CURRENT_SOURCE_DIR}/include/Achilles/System.hh)

# Copy over data files
file(COPY ${ACHILLES_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${ACHILLES_SOURCE_DIR}/flux DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${ACHILLES_SOURCE_DIR}/data/default/run.yml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${ACHILLES_SOURCE_DIR}/FormFactors.yml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${ACHILLES_SOURCE_DIR}/parameters.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${ACHILLES_SOURCE_DIR}/cascade.yml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Setup RPATH to ensure the right libraries are used
file(RELATIVE_PATH relative_rpath
    ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
    ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
)
set(CMAKE_INSTALL_RPATH $ORIGIN $ORIGIN/${relative_rpath})

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
        target_compile_options(project_options INTERFACE --coverage) 
        target_link_libraries(project_options INTERFACE gcov)
    endif()
    target_compile_definitions(project_options INTERFACE TESTING)
    enable_testing()
    message(STATUS "Building Tests.")
    add_subdirectory(test)
endif()

# Main code
add_subdirectory(src)

# Install CMake find_package files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_SOURCE_DIR}/CMake/achilles-config.cmake.in
    ${CMAKE_BINARY_DIR}/CMake/achilles-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/achilles
)

write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/CMake/achilles-config-version.cmake
    VERSION ${ACHILLES_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(
    FILES
        ${CMAKE_BINARY_DIR}/CMake/achilles-config.cmake
        ${CMAKE_BINARY_DIR}/CMake/achilles-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/achilles
)

# Install the data files to be found
install(DIRECTORY
    ${CMAKE_BINARY_DIR}/data
    ${CMAKE_BINARY_DIR}/flux
    DESTINATION share/Achilles
)
