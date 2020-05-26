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

pybind11_add_module(logger MODULE
    LoggingModule.cc
)
target_link_libraries(logger PUBLIC spdlog::spdlog)

pybind11_add_module(constants MODULE
    ConstantsModule.cc
)
target_link_libraries(constants PRIVATE project_options)

pybind11_add_module(vectors MODULE
    VectorModule.cc
)
target_link_libraries(vectors PRIVATE project_options project_warnings
                              PUBLIC utilities)

pybind11_add_module(interpolation MODULE
    InterpolationModule.cc
)
target_link_libraries(interpolation PRIVATE project_options project_warnings
                                    PUBLIC utilities)

pybind11_add_module(particle MODULE
    ParticleModule.cc
)
target_link_libraries(particle PRIVATE project_options project_warnings
                               PUBLIC utilities)

pybind11_add_module(particle_info MODULE
    ParticleInfoModule.cc
)
target_link_libraries(particle_info PRIVATE project_options project_warnings
                                    PUBLIC utilities)

pybind11_add_module(nucleus MODULE
    NucleusModule.cc
)
target_link_libraries(nucleus PRIVATE project_options project_warnings 
                              PUBLIC physics utilities)

pybind11_add_module(interactions MODULE
    InteractionsModule.cc
)
target_link_libraries(interactions PRIVATE project_options project_warnings 
                                   PUBLIC physics utilities)

pybind11_add_module(cascade MODULE
    CascadeModule.cc
)
target_link_libraries(cascade PRIVATE project_options project_warnings 
                              PUBLIC physics utilities)
