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
    System.cc
    Legendre.cc
    ClebschGordan.cc
    ReferenceHandler.cc
    XSecSpline.cc
)
target_link_libraries(utilities PRIVATE achilles_options achilles_warnings
                                PUBLIC achilles_public_headers spdlog::spdlog yaml-cpp::yaml-cpp std::filesystem) #pybind11::pybind11
if(ENABLE_AUTODIFF)
target_compile_definitions(utilities PUBLIC AUTODIFF)
target_link_libraries(utilities PUBLIC autodiff::autodiff)
endif()
if(ACHILLES_EVENT_DETAILS)
target_compile_definitions(utilities PUBLIC ACHILLES_EVENT_DETAILS)
endif()
list(APPEND achilles_targets utilities)

add_subdirectory(Integrators)

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
target_link_libraries(cuts PRIVATE achilles_options achilles_warnings
                           PUBLIC utilities) #pybind11::pybind11
list(APPEND achilles_targets cuts)

add_library(mappers SHARED
    PhaseSpaceMapper.cc
    FinalStateMapper.cc
    BeamMapper.cc
    HadronicMapper.cc
    QuasielasticTestMapper.cc
    PhaseSpaceBuilder.cc
    Channels.cc
)
target_link_libraries(mappers PRIVATE achilles_options achilles_warnings
                              PUBLIC utilities)
if(ACHILLES_ENABLE_SHERPA)
    list(APPEND achilles_targets sherpa)
    target_compile_definitions(mappers PUBLIC ACHILLES_SHERPA_INTERFACE)
    target_link_libraries(mappers PUBLIC sherpa)
endif()
list(APPEND achilles_targets mappers)

add_library(physics SHARED
    Cascade.cc
    Nucleus.cc
    FormFactor.cc
    FormFactorBuilder.cc
    Beams.cc
    FluxLoaders.cc
    LeptonicCurrent.cc
    Configuration.cc
    Event.cc
    Interactions.cc
    InteractionHandler.cc
    DecayHandler.cc
    SpectralFunction.cc
    ElectronPDF.cc
    EventHistory.cc
    MesonBaryonAmplitudes.cc
    OsetCrossSections.cc
    ResonanceHelper.cc
)
set(physics_libs "")
if(ACHILLES_ENABLE_GZIP)
    list(APPEND achilles_targets gzstream)
    target_compile_definitions(physics PUBLIC GZIP)
    target_link_libraries(physics PUBLIC gzstream)
endif()
if(ACHILLES_ENABLE_ROOT)
    target_link_libraries(physics PUBLIC ROOT::Core ROOT::Hist)
    target_compile_definitions(physics PUBLIC -DUSE_ROOT)
    if(${ROOT_MINOR_VERSION} LESS 16)
        target_link_libraries(physics PUBLIC ROOT::Flags_CXX)
    endif()
endif()
target_link_libraries(physics PRIVATE achilles_options achilles_warnings
                      PUBLIC utilities cuts AchillesIntegrator)
list(APPEND achilles_targets physics)

add_subdirectory(CascadeInteractions)

add_library(nuclear_models SHARED
    NuclearModel.cc
    XSecBackend.cc
    Process.cc
)
target_link_libraries(nuclear_models PRIVATE achilles_options achilles_warnings
                                     PUBLIC physics mappers)
list(APPEND achilles_targets nuclear_models)

add_subdirectory(fortran)
list(APPEND achilles_targets fortran_interface)

add_subdirectory(Writers)

add_library(event_gen SHARED
    EventGen.cc
    Debug.cc
)
target_link_libraries(event_gen PRIVATE achilles_options achilles_warnings
                                PUBLIC physics mappers fortran_interface nuclear_models
                                       AchillesCascadeInteractions AchillesWriters)
list(APPEND achilles_targets event_gen)

# Aggregate interface target for downstream consumers: find_package(Achilles)
# exposes achilles::Achilles. event_gen links the physics / mappers / writers /
# cascade stack PUBLIC, so this pulls the whole event-generation API plus the
# public headers.
add_library(Achilles INTERFACE)
add_library(achilles::Achilles ALIAS Achilles)
target_link_libraries(Achilles INTERFACE achilles_public_headers event_gen)
list(APPEND achilles_targets Achilles)

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
                            # target_link_libraries(_achilles PRIVATE project_warnings
                            #                               PUBLIC spdlog::spdlog utilities physics)

add_executable(achilles main.cc)
target_link_libraries(achilles PRIVATE achilles_options achilles_warnings
                                       event_gen docopt::docopt dl cmake_git_version_tracking plugin_manager)
list(APPEND achilles_targets achilles plugin_manager)

# TODO: Change the name of the option to be more clear
if(ACHILLES_ENABLE_CASCADE_TEST)
    add_executable(achilles-cascade CascadeMain.cc RunCascade.cc)
    # TODO: Refactor the event_gen library to pull pieces unused by the cascade executable
    target_link_libraries(achilles-cascade PRIVATE achilles_options achilles_warnings
                                                   event_gen docopt::docopt cmake_git_version_tracking)
    list(APPEND achilles_targets achilles-cascade)
endif()

if(ACHILLES_ENABLE_POTENTIAL_TEST)
    add_executable(achilles-potential PotentialMain.cc RunPotential.cc)
    target_link_libraries(achilles-potential PRIVATE achilles_options achilles_warnings
                                                     physics docopt::docopt)
    list(APPEND achilles_targets achilles-potential)
endif()

if(ACHILLES_ENABLE_PRECOMPUTED)
    add_executable(achilles-precomputed Precomputed.cc PrecomputedMain.cc)
    target_link_libraries(achilles-precomputed PRIVATE achilles_options achilles_warnings
                                                       physics event_gen docopt::docopt)
    list(APPEND achilles_targets achilles-precomputed)
endif()

set_target_properties(${achilles_targets} PROPERTIES # fortran_interface_f
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

include(GNUInstallDirs)

# Install Achilles targets and register them in the AchillesTargets export
# set. The export(EXPORT ...) and install(EXPORT ...) calls live in the
# top-level CMakeLists.txt so all subdirectories can contribute.
install(TARGETS ${achilles_targets}
        EXPORT AchillesTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Achilles
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
