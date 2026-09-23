# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

CPMAddPackage(
    NAME fmt
    GIT_TAG 11.0.0
    GITHUB_REPOSITORY fmtlib/fmt
)
add_library(fmt::fmt ALIAS fmt)

CPMAddPackage(
    NAME spdlog
    GIT_TAG "v1.14.1"
    GITHUB_REPOSITORY gabime/spdlog
    OPTIONS
        "SPDLOG_COMPILED_LIB ON"
        "SPDLOG_FMT_EXTERNAL ON"
        "SPDLOG_BUILD_SHARED ON"
        "SPDLOG_INSTALL ON"
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
  NAME NuHepMC_CPPUtils
  GIT_TAG v1-RC18
  GIT_REPOSITORY "https://github.com/NuHepMC/cpputils.git"
  OPTIONS "NuHepMC_CPPUtils_BUILTIN_HEPMC3 ON"
          "NuHepMC_CPPUtils_PYTHON_ENABLED OFF"
          "NuHepMC_CPPUtils_PROTOBUF_INTERFACE OFF"
)
# cpputils hardcodes -Werror, and the HepMC3 headers it bundles trigger
# warnings on newer compilers; demote those back to warnings.
if(TARGET nuhepmc_private_compile_options)
    target_compile_options(nuhepmc_private_compile_options INTERFACE -Wno-error)
endif()

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

CPMAddPackage(
    Name CmakeGitVersion
    GITHUB_REPOSITORY andrew-hardin/cmake-git-version-tracking
    GIT_TAG master
    OPTIONS
        "GIT_FAIL_IF_NONZERO_EXIT OFF"
)

# This static library gets linked into the _achilles Python extension (a shared
# object), so its objects must be position-independent.
if(TARGET cmake_git_version_tracking)
    set_target_properties(cmake_git_version_tracking PROPERTIES
        POSITION_INDEPENDENT_CODE ON)
endif()

if(ACHILLES_ENABLE_PYTHON)
find_package(python COMPONENTS Interpreter Development)
CPMFindPackage(
    NAME pybind11
    VERSION 2.13.6
    GITHUB_REPOSITORY pybind/pybind11
    GIT_TAG "v2.13.6"
)
endif(ACHILLES_ENABLE_PYTHON)

# Install testing framework
if(ACHILLES_ENABLE_TESTING)
    # Catch2
    # SYSTEM YES treats Catch2's bundled headers as system headers so that
    # the project's -Wconversion -Werror does not fail on warnings triggered
    # inside Catch2's own templates (e.g. the random/floating-point generators).
    CPMAddPackage(
        NAME catch2
        GITHUB_REPOSITORY catchorg/Catch2
        GIT_TAG v3.15.1
        SYSTEM YES
        OPTIONS
            "CATCH_ENABLE_COVERAGE"
    )
    add_library(Catch2::Catch2 ALIAS Catch2)

    # Trompeloeil
    CPMAddPackage(
        NAME trompeloeil
        GITHUB_REPOSITORY rollbear/trompeloeil
        VERSION 48
    )

    # Fortran Unit Tests
    CPMAddPackage(
        NAME fortran_unit_test
        GITHUB_REPOSITORY jxi24/fortran-unit-test
        GIT_TAG master
    )
    add_library(FUT::FUT ALIAS fortran_unit_test)

    # NuHepMC Reference Implementation
    CPMAddPackage(
        NAME NuHepMC
        GITHUB_REPOSITORY NuHepMC/ReferenceImplementation
        GIT_TAG main
    )
endif()
