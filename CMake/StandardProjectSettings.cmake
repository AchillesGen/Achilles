# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(
        STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
    set(CMAKE_BUILD_TYPE
        RelWithDebInfo
        CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui, ccmake
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                                 "MinSizeRel" "RelWithDebInfo")
endif()

# Debug info in optimized builds. Full `-g` costs ~25% compile time and produces
# ~4x larger objects than `-g1`; `-g1` still gives line tables, so backtraces,
# perf and addr2line keep working. Debug builds are left at full `-g`.
set(ACHILLES_DEBUG_INFO_LEVEL "1" CACHE STRING
    "C++ debug info level for RelWithDebInfo: 0 (none), 1 (line tables, fast), 2 (full)")
set_property(CACHE ACHILLES_DEBUG_INFO_LEVEL PROPERTY STRINGS 0 1 2)
string(REPLACE "-g " "-g${ACHILLES_DEBUG_INFO_LEVEL} "
       CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ")
string(STRIP "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CMAKE_CXX_FLAGS_RELWITHDEBINFO)

find_program(CCACHE ccache)
if(CCACHE)
    message(STATUS "Using ccache")
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE})
else()
    message(STATUS "ccache not found, cannot use")
endif()

# Generate compile_commands.json to make it easier to work with clang based tools
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ENABLE_IPO
       "Enable Interprocedural Optimization, aka Link Time Optimization (LTO)"
       OFF)

if(ENABLE_IPO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT result OUTPUT output)
    if(result)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
        message(SEND_ERROR "IPO is not supported: ${output}")
    endif()
endif()
