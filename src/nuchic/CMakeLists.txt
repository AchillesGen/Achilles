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
set_target_properties(utilities PROPERTIES LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})

add_library(interaction_plugin SHARED
    InteractionLoader.cc
    Interactions.cc
    InteractionsFactory.cc
)
target_link_libraries(interaction_plugin PRIVATE project_options project_warnings
                                         PUBLIC spdlog::spdlog stdc++fs dl
)

add_library(physics SHARED
    Cascade.cc
    Nucleus.cc
)
target_include_directories(physics SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(physics PRIVATE project_options project_warnings
                              PUBLIC utilities interaction_plugin ${HDF5_LIBRARIES})

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
