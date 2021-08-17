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
    Histogram.cc
    Cuts.cc
    MomSolver.cc
    Autodiff.cc
    SymplecticIntegrator.cc
    Potential.cc
)
target_include_directories(utilities PUBLIC $<BUILD_INTERFACE:${yaml-cpp_INCLUDE_DIRS}>)
target_link_libraries(utilities PRIVATE project_options project_warnings
                                PUBLIC spdlog::spdlog yaml::cpp) #pybind11::pybind11 
if(AUTODIFF)
target_link_libraries(utilities PUBLIC autodiff::autodiff)
endif()

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
if(ENABLE_GZIP)
target_compile_definitions(physics PUBLIC GZIP) 
target_link_libraries(physics PRIVATE project_options project_warnings
                              PUBLIC utilities ${HDF5_CXX_LIBRARIES} gzstream)
else()
target_link_libraries(physics PRIVATE project_options project_warnings
                              PUBLIC utilities ${HDF5_CXX_LIBRARIES})
endif()

add_subdirectory(fortran)

add_library(event_gen SHARED
    EventGen.cc
    EventWriter.cc
)
target_link_libraries(event_gen PRIVATE project_options project_warnings
                                PUBLIC physics fortran_interface)

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
                             PUBLIC event_gen docopt::docopt)

set(nuchic_targets nuchic event_gen fortran_interface physics utilities)
if(ENABLE_GZIP)
list(APPEND nuchic_targets gzstream)
endif()
set_target_properties(${nuchic_targets} PROPERTIES # fortran_interface_f
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

include(GNUInstallDirs)

install(TARGETS ${nuchic_targets} spdlog fmt yaml-cpp #pybind11 fortran_interface_f 
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
