#pybind11
# CPMAddPackage(
#     NAME pybind11
#     GITHUB_REPOSITORY pybind/pybind11
#     GIT_TAG v2.4.3
#     PATCH_COMMAND git apply "${CMAKE_SOURCE_DIR}/patches/pybind11.patch"
# )
# add_library(pybind11::pybind11 ALIAS pybind11)

CPMFindPackage(
    NAME fmt
    GIT_TAG 7.1.3
    GITHUB_REPOSITORY fmtlib/fmt
)
add_library(fmt::fmt ALIAS fmt)


CPMFindPackage(
    NAME spdlog
    GIT_TAG "v1.8.5"
    GITHUB_REPOSITORY gabime/spdlog
    OPTIONS
        "SPDLOG_COMPILED_LIB ON"
        "SPDLOG_FMT_EXTERNAL ON"
        "SPDLOG_BUILD_SHARED ON"
)
add_library(spdlog::spdlog ALIAS spdlog)

CPMAddPackage(
    NAME docopt
    GIT_TAG v0.6.3
    GITHUB_REPOSITORY docopt/docopt.cpp
)
add_library(docopt::docopt ALIAS docopt)

CPMAddPackage(
    NAME yaml-cpp
    GITHUB_REPOSITORY jbeder/yaml-cpp
    GIT_TAG 4edff1fa5dbfca16fc72d89870841bee89f8ef89
    OPTIONS
        "YAML_CPP_BUILD_TESTS OFF"
        "YAML_CPP_BUILD_CONTRIB OFF"
        "YAML_CPP_BUILD_TOOLS OFF"
)
add_library(yaml::cpp ALIAS yaml-cpp)

CPMAddPackage(
    NAME yaml-fortran
    GITHUB_REPOSITORY BoldingBruggeman/fortran-yaml
    GIT_TAG e465992dab9bcc147abd4c758fa05ee419f50e12
    PATCH_COMMAND git apply "${CMAKE_SOURCE_DIR}/patches/yaml-fortran.patch"
)
add_library(yaml::fortran ALIAS yaml-fortran)

# Install testing framework
if(ENABLE_TESTING OR ENABLE_FUZZING)
    # Catch2
    CPMAddPackage(
        NAME catch2
        GITHUB_REPOSITORY catchorg/Catch2
        GIT_TAG cfb695669830ceb73b0875a20af4241944eb48e0
    )
    add_library(Catch2::Catch2 ALIAS Catch2)

    # Trompeloeil
    CPMAddPackage(
        NAME trompeloeil
        GITHUB_REPOSITORY rollbear/trompeloeil
        VERSION 37 
    )

    # Fortran Unit Tests
    CPMAddPackage(
        NAME fortran_unit_test
        GITHUB_REPOSITORY dongli/fortran-unit-test
        GIT_TAG 3e291757c7ef028a93dccaac8393cc0371641dbe
        PATCH_COMMAND git apply "${CMAKE_SOURCE_DIR}/patches/fortran_unit_test.patch"
    )
    add_library(FUT::FUT ALIAS fortran_unit_test)
endif()
