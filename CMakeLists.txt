cmake_minimum_required(VERSION 3.9)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(NUCHIC)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(pybind11 REQUIRED)
find_package(Python REQUIRED COMPONENTS Development)
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
message(${Python_INCLUDE_DIRS})
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
