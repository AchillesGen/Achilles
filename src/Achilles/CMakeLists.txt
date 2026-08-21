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
    Integrator.cc
    FlowIntegrator.cc
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
    Profiler.cc
    Legendre.cc
    ClebschGordan.cc
    ReferenceHandler.cc
)
target_link_libraries(utilities PRIVATE project_warnings
                                PUBLIC project_options spdlog::spdlog yaml::cpp std::filesystem) #pybind11::pybind11
if(ENABLE_AUTODIFF)
target_compile_definitions(utilities PUBLIC AUTODIFF)
target_link_libraries(utilities PUBLIC autodiff::autodiff)
endif()
if(ACHILLES_EVENT_DETAILS)
target_compile_definitions(utilities PUBLIC ACHILLES_EVENT_DETAILS)
endif()
if(ACHILLES_ENABLE_PYTHON)
# The Flow integrator backend drives a Python flow object via pybind11. Compiling
# PyFlowSampler here keeps the MakeFlowIntegrator factory in the same library as the
# Integrator factory that calls it.
target_sources(utilities PRIVATE PyFlowSampler.cc)
target_compile_definitions(utilities PUBLIC ACHILLES_ENABLE_PYTHON)
# The standalone achilles executable, which truly embeds an interpreter, links
# pybind11::embed itself (see below).
target_link_libraries(utilities PUBLIC pybind11::module)
# pybind11 types carry hidden visibility; holding a py::object in a default-visibility
# class trips -Wattributes, which the project promotes to an error. Silence it just for
# this pybind glue file.
set_source_files_properties(PyFlowSampler.cc PROPERTIES COMPILE_OPTIONS "-Wno-attributes")
# Because utilities links pybind11::module (not ::embed), libutilities.so references
# libpython symbols it does not itself link -- they are resolved at runtime by the
# interpreter embedded in the standalone `achilles` executable, or by the host
# interpreter when loaded as the `_achilles` extension module. When GNU ld links an
# *executable* against such a shared library it reports those deferred symbols as
# undefined by default, breaking every executable that pulls in utilities transitively
# but never drives the Python flow path (the Fortran/unit tests, achilles-cascade, ...).
# Propagate --allow-shlib-undefined to all consumers so ld tolerates them; the achilles
# executable additionally links pybind11::embed and resolves them for real.
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_link_options(utilities INTERFACE "LINKER:--allow-shlib-undefined")
endif()
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
target_link_libraries(cuts PRIVATE project_warnings
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
target_link_libraries(mappers PRIVATE project_warnings
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
target_link_libraries(physics PRIVATE project_warnings
                      PUBLIC utilities cuts HighFive AchillesIntegrator)
list(APPEND achilles_targets physics)

add_subdirectory(CascadeInteractions)

add_library(nuclear_models SHARED
    NuclearModel.cc
    XSecBackend.cc
    Process.cc
)
target_link_libraries(nuclear_models PRIVATE project_options project_warnings
                                     PUBLIC physics mappers)
list(APPEND achilles_targets nuclear_models)

add_subdirectory(fortran)
list(APPEND achilles_targets fortran_interface)

add_library(event_gen SHARED
    EventGen.cc
    Debug.cc
    EventWriter.cc
)
if(ACHILLES_ENABLE_HEPMC3)
list(APPEND achilles_targets AchillesHepMC3 nuhepmc)
target_compile_definitions(physics PUBLIC ACHILLES_ENABLE_HEPMC3)
target_link_libraries(event_gen PUBLIC AchillesHepMC3 nuhepmc)
endif()
target_link_libraries(event_gen PRIVATE project_warnings
                                PUBLIC physics mappers fortran_interface nuclear_models
                                       AchillesCascadeInteractions) #fortran_interface
list(APPEND achilles_targets event_gen)

if(ACHILLES_ENABLE_PYTHON)
    pybind11_add_module(_achilles MODULE
        # Main file
        PyBindings.cc
        AchillesInit.cc
    )
target_link_libraries(_achilles PRIVATE project_warnings
                                PUBLIC event_gen cmake_git_version_tracking plugin_manager dl)
set_target_properties(_achilles PROPERTIES
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "$<IF:$<PLATFORM_ID:Darwin>,@loader_path/../lib,$ORIGIN/../lib>")
install(TARGETS _achilles DESTINATION Achilles)
endif(ACHILLES_ENABLE_PYTHON)

add_executable(achilles main.cc AchillesInit.cc)
target_link_libraries(achilles PRIVATE project_options project_warnings
                               PUBLIC event_gen docopt::docopt dl cmake_git_version_tracking plugin_manager)
if(ACHILLES_ENABLE_PYTHON)
    # The standalone executable embeds its own interpreter (scoped_interpreter in
    # PyFlowSampler), so it -- and only it -- needs libpython linked.
    target_link_libraries(achilles PUBLIC pybind11::embed)
endif()
list(APPEND achilles_targets achilles plugin_manager)

# TODO: Change the name of the option to be more clear
if(ACHILLES_ENABLE_CASCADE_TEST)
    add_executable(achilles-cascade CascadeMain.cc RunCascade.cc)
    # TODO: Refactor the event_gen library to pull pieces unused by the cascade executable
    target_link_libraries(achilles-cascade PRIVATE project_warnings
                                           PUBLIC event_gen docopt::docopt cmake_git_version_tracking)
    if(ACHILLES_ENABLE_PYTHON)
        # Pulls in libutilities.so, which references libpython without linking it.
        target_link_libraries(achilles-cascade PUBLIC pybind11::embed)
    endif()
    list(APPEND achilles_targets achilles-cascade)
endif()

if(ACHILLES_ENABLE_POTENTIAL_TEST)
    add_executable(achilles-potential PotentialMain.cc RunPotential.cc)
    target_link_libraries(achilles-potential PRIVATE project_warnings
                                             PUBLIC physics docopt::docopt)
    if(ACHILLES_ENABLE_PYTHON)
        # Pulls in libutilities.so, which references libpython without linking it.
        target_link_libraries(achilles-potential PUBLIC pybind11::embed)
    endif()
    list(APPEND achilles_targets achilles-potential)
endif()

if(ACHILLES_ENABLE_PRECOMPUTED)
    add_executable(achilles-precomputed Precomputed.cc PrecomputedMain.cc)
    target_link_libraries(achilles-precomputed PRIVATE project_warnings
                                               PUBLIC physics event_gen docopt::docopt)
    if(ACHILLES_ENABLE_PYTHON)
        # Pulls in libutilities.so, which references libpython without linking it.
        target_link_libraries(achilles-precomputed PUBLIC pybind11::embed)
    endif()
    list(APPEND achilles_targets achilles-precomputed)
endif()

set_target_properties(${achilles_targets} PROPERTIES # fortran_interface_f
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

include(GNUInstallDirs)

# TODO: Ensure that all packages installed by CPM are included when installing
install(TARGETS ${achilles_targets} spdlog yaml-cpp fmt docopt #pybind11 fortran_interface_f
        # EXPORT achilles-targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION include
)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Achilles
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
# install(EXPORT achilles-targets
#     FILE achilles-targets.cmake
#     NAMESPACE achilles::
#     DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/achilles)
