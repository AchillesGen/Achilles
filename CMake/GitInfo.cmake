# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

# achilles_git_info: a small static library with the state of the git checkout this
# copy of Achilles was built from (see include/Achilles/GitInfo.hh). The information
# is refreshed on every build by CMake/WriteGitInfo.cmake.

find_package(Git QUIET)

set(ACHILLES_GIT_INFO_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/generated/GitInfo.cc)

add_custom_target(achilles_git_info_refresh
    COMMAND ${CMAKE_COMMAND}
        -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
        -DOUTPUT=${ACHILLES_GIT_INFO_SOURCE}
        -DGIT_EXECUTABLE=${GIT_EXECUTABLE}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/CMake/WriteGitInfo.cmake
    BYPRODUCTS ${ACHILLES_GIT_INFO_SOURCE}
    COMMENT "Recording git state"
    VERBATIM
)

add_library(achilles_git_info STATIC ${ACHILLES_GIT_INFO_SOURCE})
add_dependencies(achilles_git_info achilles_git_info_refresh)
target_link_libraries(achilles_git_info PRIVATE project_options)
# Linked into the _achilles Python extension (a shared object) as well
set_target_properties(achilles_git_info PROPERTIES POSITION_INDEPENDENT_CODE ON)
