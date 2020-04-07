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

# yaml-cpp
if(yaml-cpp_FOUND)
    MESSAGE(STATUS "Found yaml-cpp at ${YAML_CPP_INCLUDE_DIR}")
else()
    MESSAGE(STATUS "yaml-cpp not found. Cloning from github as external dependency")
    FetchContent_Declare(
        extern_yamlcpp
    
        GIT_REPOSITORY "https://github.com/jbeder/yaml-cpp.git"
        GIT_TAG "yaml-cpp-0.6.3"
        )
    
    FetchContent_GetProperties(extern_yamlcpp)
    if(NOT extern_yamlcpp_POPULATED)
        FetchContent_Populate(extern_yamlcpp)
        set(YAML_CPP_BUILD_TESTS OFF CACHE INTERNAL "")
        set(YAML_BUILD_SHARED_LIBS ON CACHE INTERNAL "")
        add_subdirectory(
            ${extern_yamlcpp_SOURCE_DIR}
            ${extern_yamlcpp_BINARY_DIR}
            )
    endif()
endif()
