# Load Catch2 information into variables
include(FetchContent)
FetchContent_GetProperties(catch2)
include(${catch2_SOURCE_DIR}/contrib/Catch.cmake)

add_library(catch_main STATIC catch_main.cpp)
target_link_libraries(catch_main PUBLIC Catch2::Catch2 trompeloeil)

add_executable(nuchic-testsuite 
    # Files with tests 
    test_interp.cc
    test_vectors.cc
    # test_utils.cc
    # test_nucleus.cc
    # test_form_factor.cc
    test_beams.cc
    test_event.cc
)
target_link_libraries(nuchic-testsuite PRIVATE project_options project_warnings catch_main 
                                       PUBLIC physics)

enable_language(Fortran)
add_executable(nuchic-fortran-yaml test_yaml_fortran.f90)
target_link_libraries(nuchic-fortran-yaml PUBLIC yaml::fortran)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -ffree-line-length-0")
add_executable(nuchic-fortran-testsuite
    test_fortran.f90
)
target_link_libraries(nuchic-fortran-testsuite PRIVATE project_options
    PUBLIC physics fortran_interface FUT::FUT)

include(CTest)
add_test(NAME fortran-testsuite
         COMMAND nuchic-fortran-testsuite)
catch_discover_tests(nuchic-testsuite)
