add_library(utilities SHARED
    Constants.cc
    FourVector.cc
    Interpolation.cc
    Particle.cc
    ThreeVector.cc
    Utilities.cc
    ParticleInfo.cc
    Vegas.cc
    AdaptiveMap.cc
    Multichannel.cc
    Histogram.cc
)
target_include_directories(utilities PUBLIC $<BUILD_INTERFACE:${yaml-cpp_INCLUDE_DIRS}>)
target_link_libraries(utilities PRIVATE project_options project_warnings
                                PUBLIC spdlog::spdlog yaml::cpp) #pybind11::pybind11 

add_library(cuts SHARED
    Cuts.cc
    OneParticleCuts.cc
    TwoParticleCuts.cc
    CombinedCuts.cc
)
target_link_libraries(cuts PRIVATE project_options project_warnings
                           PUBLIC utilities) #pybind11::pybind11 

add_library(physics SHARED
    Cascade.cc
    Interactions.cc
    InteractionsFactory.cc
    Nucleus.cc
    FormFactor.cc
    Beams.cc
    HardScattering.cc
    HardScatteringFactory.cc
    Configuration.cc
    Event.cc
)
target_include_directories(physics SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(physics PRIVATE project_options project_warnings
                              PUBLIC cuts ${HDF5_CXX_LIBRARIES} gzstream sherpa)

add_library(mappers SHARED
    PhaseSpaceMapper.cc
    FinalStateMapper.cc
    BeamMapper.cc
    HadronicMapper.cc
    QuasielasticTestMapper.cc
)
target_link_libraries(mappers PRIVATE project_options project_warnings
                              PUBLIC physics)

add_subdirectory(fortran)

add_library(event_gen SHARED
    EventGen.cc
    EventWriter.cc
)
target_link_libraries(event_gen PRIVATE project_options project_warnings
                                PUBLIC physics mappers fortran_interface)

                            # pybind11_add_module(_nuchic MODULE
                            #     # Main file
                            #     PyBindings.cc
                            # 
                            #     # Utilities
                            #     LoggingModule.cc
                            #     ConstantsModule.cc
                            #     InterpolationModule.cc
                            # 
                            #     # Physics objects
                            #     VectorModule.cc
                            #     ParticleModule.cc
                            #     ParticleInfoModule.cc
                            #     NucleusModule.cc
                            # 
                            #     # Calculation modules
                            #     InteractionsModule.cc
                            #     CascadeModule.cc
                            # )
                            # target_link_libraries(_nuchic PRIVATE project_options project_warnings
                            #                               PUBLIC spdlog::spdlog utilities physics)

add_executable(nuchic main.cc)
target_link_libraries(nuchic PRIVATE project_options project_warnings
                             PUBLIC event_gen sherpa docopt::docopt)

set_target_properties(nuchic event_gen fortran_interface physics cuts utilities gzstream sherpa mappers PROPERTIES # fortran_interface_f
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

include(GNUInstallDirs)

install(TARGETS nuchic event_gen fortran_interface physics cuts utilities gzstream spdlog yaml-cpp fmt sherpa mappers #pybind11 fortran_interface_f 
        EXPORT nuchic-targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${HDF5_INCLUDE_DIRS}
)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/nuchic
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT nuchic-targets
    FILE nuchic-targets.cmake
    NAMESPACE nuchic::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nuchic)
