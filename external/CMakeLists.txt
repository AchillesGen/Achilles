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

if(AUTODIFF)
CPMFindPackage(
    NAME eigen3
    GIT_TAG "3.4"
    GIT_REPOSITORY "https://gitlab.com/libeigen/eigen"
    OPTIONS
        "BUILD_TESTING OFF"
)

CPMFindPackage(
    NAME autodiff
    GIT_TAG "v0.6.1"
    GITHUB_REPOSITORY autodiff/autodiff
    OPTIONS
        "AUTODIFF_BUILD_TESTS OFF"
        "AUTODIFF_BUILD_PYTHON OFF"
        "AUTODIFF_BUILD_EXAMPLES OFF"
        "AUTODIFF_BUILD_DOCS OFF"
)
add_library(autodiff::autodiff ALIAS autodiff)
endif()

CPMAddPackage(
    NAME docopt
    GIT_TAG v0.6.3
    GITHUB_REPOSITORY docopt/docopt.cpp
)
add_library(docopt::docopt ALIAS docopt)

CPMAddPackage(
    NAME yaml-cpp
    GITHUB_REPOSITORY jxi24/yaml-cpp
    GIT_TAG master
    OPTIONS
        "YAML_CPP_BUILD_TESTS OFF"
        "YAML_CPP_BUILD_CONTRIB OFF"
        "YAML_CPP_BUILD_TOOLS OFF"
)
add_library(yaml::cpp ALIAS yaml-cpp)

CPMAddPackage(
    NAME yaml-fortran
    GITHUB_REPOSITORY jxi24/fortran-yaml
    GIT_TAG master
)
add_library(yaml::fortran ALIAS yaml-fortran)

CPMAddPackage(
    NAME HepMC3
    VERSION 3.2.6
    GIT_REPOSITORY "https://gitlab.cern.ch/hepmc/HepMC3.git"
    GIT_TAG 3.2.6
    OPTIONS
      "HEPMC3_CXX_STANDARD ${CMAKE_CXX_STANDARD}"
      "HEPMC3_ENABLE_SEARCH OFF"
      "HEPMC3_ENABLE_ROOTIO ${USE_ROOT}"
      "HEPMC3_ENABLE_PROTOBUFIO OFF"
      "HEPMC3_ENABLE_PYTHON OFF"
      "HEPMC3_BUILD_DOCS OFF"
      "HEPMC3_BUILD_EXAMPLES OFF"
      "HEPMC3_INSTALL_EXAMPLES OFF"
      "HEPMC3_ENABLE_TEST OFF"
      "HEPMC3_INSTALL_INTERFACES OFF"
      "HEPMC3_BUILD_STATIC_LIBS OFF"
)

CPMAddPackage(
    NAME HighFive
    VERSION 2.6.2
    GITHUB_REPOSITORY BlueBrain/HighFive
    OPTIONS
        "HIGHFIVE_USE_BOOST OFF"
        "HIGHFIVE_EXAMPLES OFF"
        "HIGHFIVE_BUILD_DOCS OFF"
        "HIGHFIVE_UNIT_TESTS OFF"
)

# Install testing framework
if(ENABLE_TESTING OR ENABLE_FUZZING)
    # Catch2
    CPMAddPackage(
        NAME catch2
        GITHUB_REPOSITORY catchorg/Catch2
        GIT_TAG v2.13.6
        OPTIONS
            "CATCH_ENABLE_COVERAGE"
    )
    add_library(Catch2::Catch2 ALIAS Catch2)

    # Trompeloeil
    CPMAddPackage(
        NAME trompeloeil
        GITHUB_REPOSITORY rollbear/trompeloeil
        VERSION 44 
    )

    # Fortran Unit Tests
    CPMAddPackage(
        NAME fortran_unit_test
        GITHUB_REPOSITORY jxi24/fortran-unit-test
        GIT_TAG master
    )
    add_library(FUT::FUT ALIAS fortran_unit_test)
endif()
