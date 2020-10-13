# Load Catch2 information into variables
include(FetchContent)
FetchContent_GetProperties(catch2)
include(${catch2_SOURCE_DIR}/contrib/Catch.cmake)

add_library(catch_main STATIC catch_main.cpp)
target_link_libraries(catch_main PUBLIC Catch2::Catch2 trompeloeil)

add_executable(nuchic-testsuite 
    # Files with tests 
    test_vectors.cc
    # test_utils.cc
    # test_nucleus.cc
    # test_form_factor.cc
    test_beams.cc
)
target_link_libraries(nuchic-testsuite PRIVATE project_options project_warnings catch_main 
                                       PUBLIC physics)

add_executable(test_delta delta_test.cc)
target_link_libraries(test_delta PRIVATE project_options project_warnings
                                 PUBLIC physics)

add_executable(test_scattering test_hard_scattering.cc)
target_link_libraries(test_scattering PRIVATE project_options project_warnings
                                      PUBLIC physics)

include(CTest)
catch_discover_tests(nuchic-testsuite)
