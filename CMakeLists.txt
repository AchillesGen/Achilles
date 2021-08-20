cmake_minimum_required(VERSION 3.12)

# Policy to address @foo@ variable expansion
if(POLICY CMP0053)
    cmake_policy(SET CMP0053 NEW)
endif()

# Set the project name and basic settings
project(NUCHIC LANGUAGES CXX VERSION 1.0.0)
include(CMake/StandardProjectSettings.cmake)

# Link this 'library' to set the c++ standard / compile-time options requested
# Additionally, link to get include and external dependencies
add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE cxx_std_17)
target_include_directories(project_options INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

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

option(BUILD_SHARED_LIBS "Enable compilation of shared libraries" ON)
option(ENABLE_TESTING "Enable Test Builds" OFF)
option(ENABLE_AUTODIFF "Enable Autodiff" OFF)
option(ENABLE_GZIP "Enable compression of event files" ON)
option(ENABLE_CASCADE_TEST "Enable executables for testing the cascade" OFF)

# Very basic PCH example
option(ENABLE_PCH "Enable Precompiled Headers" OFF)
if(ENABLE_PCH)
    # This sets a global PCH parameter, each project will build its own PCH,
    # which is a good idea if any #define's change
    # These should be headers included in many places
    target_precompile_headers(project_options INTERFACE <vector> <string> <map>)
endif()

# Check for python development environment
# find_package(Python REQUIRED COMPONENTS Development)

# Find HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)

# Find ZLIB to read gzip files
if(ENABLE_GZIP)
find_package(ZLIB REQUIRED)
endif()

# Add dependencies
include(CMake/CPM.cmake)
add_subdirectory(external)

# Configure system specific variables for loading shared libraries
set(LIBPREFIX ${CMAKE_SHARED_LIBRARY_PREFIX})
set(LIBSUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
set(BUILDPATH ${CMAKE_CURRENT_BINARY_DIR})
set(INSTALLPATH ${CMAKE_INSTALL_PREFIX})
configure_file(include/nuchic/System.hh.in ${CMAKE_CURRENT_SOURCE_DIR}/include/nuchic/System.hh)

# Copy over data files
file(COPY ${NUCHIC_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${NUCHIC_SOURCE_DIR}/run.yml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Unzip if GZIP not turned on
if(NOT ENABLE_GZIP)
    file(GLOB zipfiles "${CMAKE_CURRENT_BINARY_DIR}/data/configurations/*.gz")
    MESSAGE(STATUS "Extracting configurations from zip files")
    foreach(file ${zipfiles})
        execute_process(COMMAND gunzip "${file}"
                        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data/configurations/")
    endforeach()
endif()

# Testing
if(ENABLE_TESTING)
    if(COVERAGE) 
        target_compile_options(project_options INTERFACE --coverage) 
        target_link_libraries(project_options INTERFACE gcov)
    endif()
    target_compile_definitions(project_options INTERFACE TESTING)
    enable_testing()
    message(STATUS "Building Tests.")
    add_subdirectory(test)
endif()

if(ENABLE_FUZZING)
    message(STATUS 
        "Building Fuzz Tests, using fuzzing sanitizer https://www.llvm.org/docs/LibFuzzer.html"
    )
    add_subdirectory(fuzz_test)
endif()

# Main code
add_subdirectory(src)

# Install CMake find_package files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_SOURCE_DIR}/CMake/nuchic-config.cmake.in
    ${CMAKE_BINARY_DIR}/CMake/nuchic-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nuchic
)

write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/CMake/nuchic-config-version.cmake
    VERSION ${NUCHIC_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(
    FILES
        ${CMAKE_BINARY_DIR}/CMake/nuchic-config.cmake
        ${CMAKE_BINARY_DIR}/CMake/nuchic-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nuchic
)
