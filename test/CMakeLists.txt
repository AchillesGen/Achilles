# Load Catch2 information into variables
include(FetchContent)
FetchContent_GetProperties(catch2)
include(${catch2_SOURCE_DIR}/contrib/Catch.cmake)

add_library(catch_main STATIC catch_main.cpp)
target_link_libraries(catch_main PUBLIC Catch2::Catch2 trompeloeil)
target_compile_definitions(catch_main PUBLIC CATCH_CONFIG_ENABLE_BENCHMARKING)

add_executable(achilles-testsuite 
    # Files with tests 
    test_interp.cc
    test_vectors.cc
    test_utils.cc
    test_particle_info.cc
    test_particle.cc
    test_nucleus.cc
    test_form_factor.cc
    test_cascade.cc
    test_beams.cc
    test_event.cc
    test_cuts.cc
    test_mom_solver.cc
    test_autodiff.cc
    test_symplectic_integrator.cc
    test_potential.cc
    test_adaptive_map.cc
    test_stats.cc
    test_vegas.cc
    test_multichannel.cc
    # test_integrand.cc
    test_spectral.cc
    test_spinor.cc
    test_histogram.cc
    test_configuration.cc
    test_nuclear_model.cc
    test_event_writer.cc
    test_process_info.cc
    test_hadronic_mapper.cc
    test_final_state_mapper.cc
    test_quasielastic_mapper.cc
    test_beam_mapper.cc
    test_ps_mapper.cc
    test_event_history.cc
    test_poincare.cc
    test_process.cc
    test_process_group.cc
    # NOTE: NuHepMC validator is broken, so the tests are broken
    # test_nuhepmc.cc
    test_cache.cc
)
target_link_libraries(achilles-testsuite PRIVATE project_options project_warnings catch_main 
    PUBLIC physics mappers event_gen)
if(ENABLE_BSM)
    target_link_libraries(achilles-testsuite PUBLIC -L${SHERPA_ROOT_DIR}/lib/SHERPA-MC -lToolsOrg -lMEToolsMain -lToolsPhys)
endif()

set_target_properties(NuHepMCReferenceWriter NuHepMCReferenceValidator PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)


add_executable(achilles-fortran-yaml test_yaml_fortran.f90)
target_link_libraries(achilles-fortran-yaml PUBLIC yaml::fortran)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -ffree-line-length-0")
add_executable(achilles-fortran-testsuite
    test_fortran.f90
    test_spectral_interface.f90
    test_fortran_map.f90
    test_fortran_interference.f90
)
message(STATUS ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(achilles-fortran-testsuite PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_link_libraries(achilles-fortran-testsuite PRIVATE project_options
    PUBLIC physics fortran_interface FUT::FUT)

include(CTest)
# add_test(NAME fortran-testsuite
#          COMMAND achilles-fortran-testsuite)
catch_discover_tests(achilles-testsuite
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
