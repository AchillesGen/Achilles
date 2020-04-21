#pybind11
CPMAddPackage(
    NAME pybind11
    GITHUB_REPOSITORY pybind/pybind11
    GIT_TAG v2.4.0
    PATCH_COMMAND git apply "${CMAKE_SOURCE_DIR}/patches/pybind11.patch"
)
add_library(pybind11::pybind11 ALIAS pybind11)

CPMFindPackage(
    NAME fmt
    # Use this commit which eliminates variable shadowing but is a stable point
    GIT_TAG 63479c8519c7fa3571707ec7f40b87438b9b9c51
    GITHUB_REPOSITORY fmtlib/fmt
)
add_library(fmt::fmt ALIAS fmt)


CPMFindPackage(
    NAME spdlog
    GIT_TAG "v1.5.0"
    GITHUB_REPOSITORY gabime/spdlog
    OPTIONS 
        "SPDLOG_FMT_EXTERNAL ON"
        "SPDLOG_BUILD_SHARED ON"
)
add_library(spdlog::spdlog ALIAS spdlog)

CPMFindPackage(
    NAME yaml-cpp
    GITHUB_REPOSITORY jbeder/yaml-cpp
    GIT_TAG 4edff1fa5dbfca16fc72d89870841bee89f8ef89
    OPTIONS
        "YAML_CPP_BUILD_TESTS OFF"
        "YAML_CPP_BUILD_CONTRIB OFF"
        "YAML_CPP_BUILD_TOOLS OFF"
)
