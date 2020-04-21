cmake_minimum_required(VERSION 3.12)

# Policy to address @foo@ variable expansion
if(POLICY CMP0053)
    cmake_policy(SET CMP0053 NEW)
endif()

# Set the project name and basic settings
project(NUCHIC CXX)
include(CMake/StandardProjectSettings.cmake)

# Link this 'library' to set the c++ standard / compile-time options requested
# Additionally, link to get include and external dependencies
add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE cxx_std_11)
target_include_directories(project_options INTERFACE include external)


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

# Very basic PCH example
option(ENABLE_PCH "Enable Precompiled Headers" OFF)
if(ENABLE_PCH)
    # This sets a global PCH parameter, each project will build its own PCH,
    # which is a good idea if any #define's change
    # These should be headers included in many places
    target_precompile_headers(project_options INTERFACE <vector> <string> <map>)
endif()

# Check for python development environment
find_package(Python REQUIRED COMPONENTS Development)

# Find HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)

# Add dependencies
include(CMake/CPM.cmake)
add_subdirectory(external)

# Main code
add_subdirectory(src/nuchic)

# Testing
if(ENABLE_TESTING)
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
