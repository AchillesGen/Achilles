include(FetchContent)

# pybind11
FetchContent_Declare(
    extern_pybind11

    GIT_REPOSITORY "https://github.com/pybind/pybind11.git"
    GIT_TAG "v2.4"
    PATCH_COMMAND git apply "${CMAKE_SOURCE_DIR}/patches/pybind11.patch" 
    )

FetchContent_GetProperties(extern_pybind11)
if(NOT extern_pybind11_POPULATED)
    FetchContent_Populate(extern_pybind11)
    add_subdirectory(
        ${extern_pybind11_SOURCE_DIR}
        ${extern_pybind11_BINARY_DIR}
        EXCLUDE_FROM_ALL)
endif()

# spdlog
FetchContent_Declare(
    extern_spdlog

    GIT_REPOSITORY "https://github.com/gabime/spdlog.git"
    GIT_TAG "v1.x"
    PATCH_COMMAND git apply "${CMAKE_SOURCE_DIR}/patches/spdlog.patch"
    )

FetchContent_GetProperties(extern_spdlog)
if(NOT extern_spdlog_POPULATED)
    FetchContent_Populate(extern_spdlog)
    set(SPDLOG_BUILD_SHARED ON CACHE INTERNAL "")
    add_subdirectory(
        ${extern_spdlog_SOURCE_DIR}
        ${extern_spdlog_BINARY_DIR}
        EXCLUDE_FROM_ALL)
endif()

# Install testing framework
if(ENABLE_TESTING OR ENABLE_FUZZING)
    # Catch2
    FetchContent_Declare(
        extern_catch2

        GIT_REPOSITORY "https://github.com/catchorg/Catch2.git"
    )

    FetchContent_GetProperties(extern_catch2)
    if(NOT extern_catch2_POPULATED)
        FetchContent_Populate(extern_catch2)
        add_subdirectory(
            ${extern_catch2_SOURCE_DIR}
            ${extern_catch2_BINARY_DIR})
    endif()

    # Trompeloeil
    FetchContent_Declare(
        extern_trompeloeil

        GIT_REPOSITORY "https://github.com/rollbear/trompeloeil.git"
    )

    FetchContent_GetProperties(extern_trompeloeil)
    if(NOT extern_trompeloeil_POPULATED)
        FetchContent_Populate(extern_trompeloeil)
        add_subdirectory(
            ${extern_trompeloeil_SOURCE_DIR}
            ${extern_trompeloeil_BINARY_DIR})
    endif()
endif()
