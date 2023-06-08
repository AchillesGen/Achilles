set(achilles_targets "")

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
    MomSolver.cc
    Autodiff.cc
    SymplecticIntegrator.cc
    Potential.cc
    Spinor.cc
    ProcessInfo.cc
    Poincare.cc
    Unweighter.cc
    Settings.cc
)
target_include_directories(utilities PUBLIC $<BUILD_INTERFACE:${yaml-cpp_INCLUDE_DIRS}>)
target_link_libraries(utilities PRIVATE project_options project_warnings
                                PUBLIC spdlog::spdlog yaml::cpp) #pybind11::pybind11 
list(APPEND achilles_targets utilities)
if(ENABLE_AUTODIFF)
target_compile_definitions(utilities PUBLIC AUTODIFF)
target_link_libraries(utilities PUBLIC autodiff::autodiff)
endif()

# add_library(interaction_plugin SHARED
#     InteractionLoader.cc
# )
# target_link_libraries(interaction_plugin PRIVATE project_options project_warnings
#                                          PUBLIC spdlog::spdlog stdc++fs dl yaml::cpp
# )

add_library(cuts SHARED
    OneParticleCuts.cc
    TwoParticleCuts.cc
    CombinedCuts.cc
)
target_link_libraries(cuts PRIVATE project_options project_warnings
                           PUBLIC utilities) #pybind11::pybind11 
list(APPEND achilles_targets cuts)

add_library(physics SHARED
    Cascade.cc
    Nucleus.cc
    FormFactor.cc
    FormFactorBuilder.cc
    Beams.cc
    HardScattering.cc
    Configuration.cc
    Event.cc
    Interactions.cc
    InteractionsFactory.cc
    SpectralFunction.cc
    Channels.cc
    ElectronPDF.cc
    EventHistory.cc
)
# target_include_directories(physics SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS})
set(physics_libs "")
if(ENABLE_GZIP)
    list(APPEND achilles_targets gzstream)
    target_compile_definitions(physics PUBLIC GZIP)
    target_link_libraries(physics PUBLIC gzstream)
endif()
if(ENABLE_BSM)
    list(APPEND achilles_targets sherpa)
    target_compile_definitions(physics PUBLIC ENABLE_BSM)
    target_link_libraries(physics PUBLIC sherpa)
endif()
if(USE_ROOT)
    target_link_libraries(physics PUBLIC ROOT::Core ROOT::Hist)
    target_compile_definitions(physics PUBLIC -DUSE_ROOT)
    if(${ROOT_MINOR_VERSION} LESS 16)
        target_link_libraries(physics PUBLIC ROOT::Flags_CXX)
    endif()
endif()
if(ACHILLES_LOW_MEMORY)
    target_compile_definitions(physics PUBLIC -DACHILLES_LOW_MEMORY)
endif()
target_link_libraries(physics PRIVATE project_options project_warnings
                      PUBLIC utilities cuts HighFive)
list(APPEND achilles_targets physics)

add_library(mappers SHARED
    PhaseSpaceMapper.cc
    FinalStateMapper.cc
    BeamMapper.cc
    HadronicMapper.cc
    QuasielasticTestMapper.cc
    PhaseSpaceBuilder.cc
)
target_link_libraries(mappers PRIVATE project_options project_warnings
                              PUBLIC physics)
list(APPEND achilles_targets mappers)

# add_subdirectory(fortran)

add_library(event_gen SHARED
    # TODO: Move to its own library
    NuclearModel.cc
    EventGen.cc
    EventWriter.cc
)
if(ENABLE_HEPMC3)
list(APPEND achilles_targets hepmc3)
target_compile_definitions(physics PUBLIC ENABLE_HEPMC3)
target_link_libraries(event_gen PUBLIC hepmc3)
endif()
target_link_libraries(event_gen PRIVATE project_options project_warnings
                                PUBLIC physics mappers) #fortran_interface
list(APPEND achilles_targets event_gen)

                            # pybind11_add_module(_achilles MODULE
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
                            # target_link_libraries(_achilles PRIVATE project_options project_warnings
                            #                               PUBLIC spdlog::spdlog utilities physics)

add_executable(achilles main.cc)
target_link_libraries(achilles PRIVATE project_options project_warnings
                               PUBLIC event_gen docopt::docopt dl)
list(APPEND achilles_targets achilles)

if(ENABLE_CASCADE_TEST)
    add_executable(achilles-cascade CascadeMain.cc RunCascade.cc)
    target_link_libraries(achilles-cascade PRIVATE project_options project_warnings
                                         PUBLIC physics docopt::docopt)
    list(APPEND achilles_targets achilles-cascade)
endif()

if(ENABLE_POTENTIAL_TEST)
    add_executable(achilles-potential PotentialMain.cc RunPotential.cc)
    target_link_libraries(achilles-potential PRIVATE project_options project_warnings
                                           PUBLIC physics docopt::docopt)
    list(APPEND achilles_targets achilles-potential)
endif()

set_target_properties(${achilles_targets} PROPERTIES # fortran_interface_f
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

include(GNUInstallDirs)

install(TARGETS ${achilles_targets} spdlog yaml-cpp fmt #pybind11 fortran_interface_f 
        # EXPORT achilles-targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/achilles
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
# install(EXPORT achilles-targets
#     FILE achilles-targets.cmake
#     NAMESPACE achilles::
#     DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/achilles)
