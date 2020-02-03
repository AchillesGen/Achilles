cmake_minimum_required(VERSION 3.9)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS TRUE)

project(NUCHIC)

# Policy to address @foo@ variable expansion
if(POLICY CMP0053)
    cmake_policy(SET CMP0053 NEW)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Setup CMake Directory
set(NUCHIC_CMAKE_DIR "${NUCHIC_SOURCE_DIR}/CMake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake")
include(PyLink)

# Check for python development environment
find_package(Python REQUIRED COMPONENTS Development)

# Find HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)

set(NUCHIC_INCLUDES ${NUCHIC_BUILD_INSTALL_PREFIX}/include ${NUCHIC_SOURCE_DIR}/include ${HDF5_INCLUDE_DIRS})
link_directories(${NUCHIC_BUILD_INSTALL_PREFIX}/lib)

# Pull down pybind11 submodule
execute_process(COMMAND git submodule update --init -- external/pybind11
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(external/pybind11)
include_directories(${PYBIND11_INCLUDE_DIR})

add_subdirectory(src/nuchic)
