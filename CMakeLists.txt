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

# Check for python development environment
find_package(Python REQUIRED COMPONENTS Development)

# Pull down pybind11 submodule
execute_process(COMMAND git submodule update --init -- external/pybind11
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(external/pybind11)


# Find HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)

include_directories(${NUCHIC_BUILD_INSTALL_PREFIX}/include)
link_directories(${NUCHIC_BUILD_INSTALL_PREFIX}/lib)
include_directories(${NUCHIC_SOURCE_DIR}/include)
include_directories(${HDF5_INCLUDE_DIRS})

add_library(utilities SHARED
    src/nuchic/FourVector.cc
    src/nuchic/Interpolation.cc
    src/nuchic/Particle.cc
    src/nuchic/ThreeVector.cc
    src/nuchic/Utilities.cc
)
target_include_directories(utilities PUBLIC ${Python_INCLUDE_DIRS})

add_library(physics SHARED
    src/nuchic/Cascade.cc
    src/nuchic/Interactions.cc
    src/nuchic/InteractionsFactory.cc
    src/nuchic/Nucleus.cc
)

target_link_libraries(physics PRIVATE utilities ${HDF5_LIBRARIES})

pybind11_add_module(vectors MODULE
    src/nuchic/VectorModule.cc
)
target_link_libraries(vectors PRIVATE utilities physics ${HDF5_LIBRARIES})

pybind11_add_module(interpolation MODULE
    src/nuchic/InterpolationModule.cc
)
target_link_libraries(interpolation PRIVATE utilities physics ${HDF5_LIBRARIES})

pybind11_add_module(particle MODULE
    src/nuchic/ParticleModule.cc
)
target_link_libraries(particle PRIVATE utilities physics ${HDF5_LIBRARIES})

pybind11_add_module(nucleus MODULE
    src/nuchic/NucleusModule.cc
)
target_link_libraries(nucleus PRIVATE utilities physics ${HDF5_LIBRARIES})

pybind11_add_module(interactions MODULE
    src/nuchic/InteractionsModule.cc
)
target_link_libraries(interactions PRIVATE utilities physics ${HDF5_LIBRARIES})

pybind11_add_module(cascade MODULE
    src/nuchic/CascadeModule.cc
)
target_link_libraries(cascade PRIVATE utilities physics ${HDF5_LIBRARIES})
