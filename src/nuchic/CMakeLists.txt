add_library(utilities SHARED
    FourVector.cc
    Interpolation.cc
    Particle.cc
    ThreeVector.cc
    Utilities.cc
    ParticleInfo.cc
)
target_include_directories(utilities PUBLIC ${yaml-cpp_INCLUDE_DIRS})
target_link_libraries(utilities PRIVATE project_options project_warnings
                                PUBLIC spdlog::spdlog pybind11::pybind11 yaml-cpp)

add_library(physics SHARED
    Cascade.cc
    Interactions.cc
    InteractionsFactory.cc
    Nucleus.cc
)
target_include_directories(physics SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(physics PRIVATE project_options project_warnings
                              PUBLIC utilities ${HDF5_LIBRARIES})

pybind11_add_module(_nuchic MODULE
    # Main file
    PyBindings.cc

    # Utilities
    LoggingModule.cc
    ConstantsModule.cc
    InterpolationModule.cc

    # Physics objects
    VectorModule.cc
    ParticleModule.cc
    ParticleInfoModule.cc
    NucleusModule.cc

    # Calculation modules
    InteractionsModule.cc
    CascadeModule.cc
)
target_link_libraries(_nuchic PRIVATE project_options project_warnings
                              PUBLIC spdlog::spdlog utilities physics)
