# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

set(achilles_targets "")

add_library(AchillesUtilities SHARED
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
# achilles_options carries build-only flags (sanitizers, coverage), so keep it out of the export.
target_link_libraries(AchillesUtilities PRIVATE achilles_warnings
                                        PUBLIC $<BUILD_INTERFACE:achilles_options> spdlog::spdlog yaml-cpp::yaml-cpp std::filesystem) #pybind11::pybind11
target_compile_features(AchillesUtilities PUBLIC cxx_std_17)
if(ENABLE_AUTODIFF)
target_compile_definitions(AchillesUtilities PUBLIC AUTODIFF)
target_link_libraries(AchillesUtilities PUBLIC autodiff::autodiff)
endif()
if(ACHILLES_EVENT_DETAILS)
target_compile_definitions(AchillesUtilities PUBLIC ACHILLES_EVENT_DETAILS)
endif()
list(APPEND achilles_targets AchillesUtilities)

add_subdirectory(Integrators)

# add_library(interaction_plugin SHARED
#     InteractionLoader.cc
# )
# target_link_libraries(interaction_plugin PRIVATE achilles_options achilles_warnings
#                                          PUBLIC spdlog::spdlog stdc++fs dl yaml-cpp::yaml-cpp
# )

add_library(AchillesCuts SHARED
    OneParticleCuts.cc
    TwoParticleCuts.cc
    CombinedCuts.cc
)
target_link_libraries(AchillesCuts PRIVATE achilles_warnings
                                   PUBLIC AchillesUtilities) #pybind11::pybind11
list(APPEND achilles_targets AchillesCuts)

add_library(AchillesMappers SHARED
    PhaseSpaceMapper.cc
    FinalStateMapper.cc
    BeamMapper.cc
    HadronicMapper.cc
    QuasielasticTestMapper.cc
    PhaseSpaceBuilder.cc
    Channels.cc
)
target_link_libraries(AchillesMappers PRIVATE achilles_warnings
                                      PUBLIC AchillesUtilities)
if(ACHILLES_ENABLE_SHERPA)
    list(APPEND achilles_targets AchillesSherpa)
    target_compile_definitions(AchillesMappers PUBLIC ACHILLES_SHERPA_INTERFACE)
    target_link_libraries(AchillesMappers PUBLIC AchillesSherpa)
endif()
list(APPEND achilles_targets AchillesMappers)

add_library(AchillesPhysics SHARED
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
    ResonanceHelper.cc
)
set(physics_libs "")
if(ACHILLES_ENABLE_GZIP)
    target_compile_definitions(AchillesPhysics PUBLIC $<BUILD_INTERFACE:GZIP>)
    target_link_libraries(AchillesPhysics PRIVATE $<BUILD_INTERFACE:gzstream>)
endif()
if(ACHILLES_ENABLE_ROOT)
    target_link_libraries(AchillesPhysics PUBLIC ROOT::Core ROOT::Hist)
    target_compile_definitions(AchillesPhysics PUBLIC -DUSE_ROOT)
    if(${ROOT_MINOR_VERSION} LESS 16)
        target_link_libraries(AchillesPhysics PUBLIC ROOT::Flags_CXX)
    endif()
endif()
target_link_libraries(AchillesPhysics PRIVATE achilles_warnings
                                      PUBLIC AchillesUtilities AchillesCuts HighFive AchillesIntegrator)
list(APPEND achilles_targets AchillesPhysics)

add_subdirectory(CascadeInteractions)

add_library(AchillesNuclearModels SHARED
    NuclearModel.cc
    XSecBackend.cc
    Process.cc
)
target_link_libraries(AchillesNuclearModels PRIVATE achilles_options achilles_warnings
                                            PUBLIC AchillesPhysics AchillesMappers)
list(APPEND achilles_targets AchillesNuclearModels)

add_subdirectory(fortran)
list(APPEND achilles_targets AchillesFortranInterface)

add_library(AchillesEventGen SHARED
    EventGen.cc
    Debug.cc
    EventWriter.cc
)
if(ACHILLES_ENABLE_HEPMC3)
list(APPEND achilles_targets AchillesHepMC3 AchillesNuHepMC)
target_compile_definitions(AchillesPhysics PUBLIC ACHILLES_ENABLE_HEPMC3)
target_link_libraries(AchillesEventGen PUBLIC AchillesHepMC3 AchillesNuHepMC)
endif()
target_link_libraries(AchillesEventGen PRIVATE achilles_warnings
                                       PUBLIC AchillesPhysics AchillesMappers AchillesFortranInterface AchillesNuclearModels
                                              ${AchillesCascadeInteractionsLink}) #AchillesFortranInterface
list(APPEND achilles_targets AchillesEventGen)
if(ACHILLES_ENABLE_GZIP)
    target_link_libraries(AchillesEventGen PRIVATE $<BUILD_INTERFACE:gzstream>)
endif()

if(ACHILLES_ENABLE_PYTHON)
    pybind11_add_module(_achilles MODULE
        # Main file
        PyBindings.cc
        AchillesInit.cc
    )
target_link_libraries(_achilles PRIVATE achilles_warnings
                                PUBLIC AchillesEventGen cmake_git_version_tracking AchillesManager dl)
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
                               PUBLIC AchillesEventGen docopt::docopt dl cmake_git_version_tracking AchillesManager)
list(APPEND achilles_targets achilles AchillesManager)
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
    # TODO: Refactor the AchillesEventGen library to pull pieces unused by the cascade executable
    target_link_libraries(achilles-cascade PRIVATE achilles_warnings
                                           PUBLIC AchillesEventGen docopt::docopt cmake_git_version_tracking)
    list(APPEND achilles_targets achilles-cascade)
endif()

if(ACHILLES_ENABLE_POTENTIAL_TEST)
    add_executable(achilles-potential PotentialMain.cc RunPotential.cc)
    target_link_libraries(achilles-potential PRIVATE achilles_warnings
                                             PUBLIC AchillesPhysics docopt::docopt)
    list(APPEND achilles_targets achilles-potential)
endif()

if(ACHILLES_ENABLE_PRECOMPUTED)
    add_executable(achilles-precomputed Precomputed.cc PrecomputedMain.cc)
    target_link_libraries(achilles-precomputed PRIVATE achilles_warnings
                                               PUBLIC AchillesPhysics AchillesEventGen docopt::docopt)
    list(APPEND achilles_targets achilles-precomputed)
endif()

set_target_properties(${achilles_targets} PROPERTIES # fortran_interface_f
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Single entry point for find_package(Achilles) consumers
add_library(Achilles INTERFACE)
add_library(achilles::Achilles ALIAS Achilles)
target_link_libraries(Achilles INTERFACE AchillesEventGen)
list(APPEND achilles_targets Achilles)

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

install(TARGETS ${achilles_targets}
        EXPORT AchillesTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${achilles_library_dest}
        RUNTIME DESTINATION ${achilles_runtime_dest}
        INCLUDES DESTINATION include
)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Achilles
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Some CPM dependencies hardcode a top-level lib/ in their own install rules and
# ignore CMAKE_INSTALL_LIBDIR, so they land outside achilles/lib where the
# $ORIGIN rpath cannot reach them. Install those targets into achilles/lib too.
if(SKBUILD)
    foreach(_straggler nuhepmc_cpputils)
        if(TARGET ${_straggler})
            install(TARGETS ${_straggler} LIBRARY DESTINATION achilles/lib)
        endif()
    endforeach()
endif()
