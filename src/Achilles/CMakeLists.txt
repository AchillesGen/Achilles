# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

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

add_library(event_gen SHARED
    EventGen.cc
    Debug.cc
    EventWriter.cc
)
if(ACHILLES_ENABLE_HEPMC3)
target_compile_definitions(physics PUBLIC ACHILLES_ENABLE_HEPMC3)
# EventGen.cc instantiates the writer classes directly, so event_gen link-depends
# on them. They pull in HepMC3's own (non-exported) build targets, which is why
# event_gen itself is kept out of the AchillesTargets export set below.
target_link_libraries(event_gen PRIVATE AchillesHepMC3 nuhepmc)
endif()
target_link_libraries(event_gen PRIVATE achilles_options achilles_warnings
                                PUBLIC physics mappers fortran_interface nuclear_models
                                       ${AchillesCascadeInteractionsLink}) #fortran_interface
list(APPEND achilles_targets event_gen)

if(ACHILLES_ENABLE_PYTHON)
    pybind11_add_module(_achilles MODULE
        # Main file
        PyBindings.cc
        AchillesInit.cc
    )
target_link_libraries(_achilles PRIVATE achilles_options achilles_warnings
                                PUBLIC event_gen cmake_git_version_tracking plugin_manager dl)
# _achilles lives in the achilles package dir; the shared libraries are in the
# lib/ subdir (achilles/lib), so point the rpath there ($ORIGIN/lib), not at a
# top-level ../lib.
set_target_properties(_achilles PROPERTIES
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "$<IF:$<PLATFORM_ID:Darwin>,@loader_path/lib,$ORIGIN/lib>")
install(TARGETS _achilles DESTINATION achilles)
endif(ACHILLES_ENABLE_PYTHON)

add_executable(achilles main.cc AchillesInit.cc)
target_link_libraries(achilles PRIVATE achilles_options achilles_warnings
                                       event_gen docopt::docopt dl cmake_git_version_tracking plugin_manager)
list(APPEND achilles_targets achilles plugin_manager)
# Like _achilles, the executable lives in the package dir and needs achilles/lib.
# Bake the rpath at build time: macOS install-time rpath rewriting of executables
# is unreliable and otherwise leaves the binary with no rpath (delocate then
# cannot resolve its @rpath dependencies).
if(SKBUILD)
    set_target_properties(achilles PROPERTIES
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_RPATH "$<IF:$<PLATFORM_ID:Darwin>,@loader_path/lib,$ORIGIN/lib>")
endif()

# TODO: Change the name of the option to be more clear
if(ACHILLES_ENABLE_CASCADE_TEST)
    add_executable(achilles-cascade CascadeMain.cc RunCascade.cc)
    # TODO: Refactor the event_gen library to pull pieces unused by the cascade executable
    target_link_libraries(achilles-cascade PRIVATE achilles_options achilles_warnings
                                           PUBLIC event_gen docopt::docopt cmake_git_version_tracking)
    list(APPEND achilles_targets achilles-cascade)
endif()

if(ACHILLES_ENABLE_POTENTIAL_TEST)
    add_executable(achilles-potential PotentialMain.cc RunPotential.cc)
    target_link_libraries(achilles-potential PRIVATE achilles_options achilles_warnings
                                             PUBLIC physics docopt::docopt)
    list(APPEND achilles_targets achilles-potential)
endif()

if(ACHILLES_ENABLE_PRECOMPUTED)
    add_executable(achilles-precomputed Precomputed.cc PrecomputedMain.cc)
    target_link_libraries(achilles-precomputed PRIVATE achilles_options achilles_warnings
                                               PUBLIC physics event_gen docopt::docopt)
    list(APPEND achilles_targets achilles-precomputed)
endif()

set_target_properties(${achilles_targets} PROPERTIES # fortran_interface_f
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

include(GNUInstallDirs)

# In a wheel, install everything inside the achilles package directory:
# This keeps the libraries inside the package so auditwheel does not orphan them
if(SKBUILD)
    set(achilles_runtime_dest achilles)
    set(achilles_library_dest achilles/lib)
else()
    set(achilles_runtime_dest ${CMAKE_INSTALL_BINDIR})
    set(achilles_library_dest ${CMAKE_INSTALL_LIBDIR})
endif()

# The HepMC/NuHepMC writer plugins link HepMC3's own (non-exported) build
# targets, so they must not go into the AchillesTargets export set. Install
# them on their own; the install(EXPORT) below only sees the core libraries.
if(ACHILLES_ENABLE_HEPMC3)
    set_target_properties(AchillesHepMC3 nuhepmc PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    install(TARGETS AchillesHepMC3 nuhepmc
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${achilles_library_dest}
        RUNTIME DESTINATION ${achilles_runtime_dest})
endif()

# Split the targets into the exported library API and the rest. The generator
# driver (event_gen), the executables and plugin_manager link the HepMC writers
# -- which drag in HepMC3's own non-exported build targets -- so they cannot go
# into the export set. They are still installed, just without EXPORT. What
# remains (utilities, cuts, mappers, physics, nuclear_models, integrators,
# cascade interactions, fortran_interface, gzstream, sherpa) is exactly the
# library surface a downstream consumer or plugin links against.
set(achilles_export_targets ${achilles_targets})
list(REMOVE_ITEM achilles_export_targets
     event_gen plugin_manager achilles
     achilles-cascade achilles-potential achilles-precomputed)
set(achilles_noexport_targets ${achilles_targets})
list(REMOVE_ITEM achilles_noexport_targets ${achilles_export_targets})

# Register the library API in the AchillesTargets export set so a consumer's
# find_package(Achilles) resolves the imported targets (achilles::utilities,
# achilles::cuts, ...). The install(EXPORT) that writes AchillesTargets.cmake
# lives in the top-level CMakeLists.txt after every subdirectory has been added.
# spdlog/fmt/docopt are NOT listed here: they already register themselves in
# their own CPM export sets (a target may live in only one) and are re-resolved
# by find_dependency() in AchillesConfig.cmake.in. yaml-cpp IS listed because,
# unlike them, it does not register itself in an export set yet utilities links
# it publicly.
install(TARGETS ${achilles_export_targets} yaml-cpp #pybind11 fortran_interface_f
        EXPORT AchillesTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${achilles_library_dest}
        RUNTIME DESTINATION ${achilles_runtime_dest}
        INCLUDES DESTINATION include
)
# Install the remaining (non-exported) Achilles targets and the CPM dependency
# runtimes so a standalone install tree is self-contained.
install(TARGETS ${achilles_noexport_targets} spdlog fmt docopt
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${achilles_library_dest}
        RUNTIME DESTINATION ${achilles_runtime_dest})
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Achilles
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Some CPM dependencies hardcode a top-level lib/ in their own install rules and
# ignore CMAKE_INSTALL_LIBDIR, so they land outside achilles/lib where the
# $ORIGIN rpath cannot reach them. Install those targets into achilles/lib too.
if(SKBUILD)
    foreach(_straggler nuhepmc_cpputils yaml-fortran)
        if(TARGET ${_straggler})
            install(TARGETS ${_straggler} LIBRARY DESTINATION achilles/lib)
        endif()
    endforeach()
endif()
