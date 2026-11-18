# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

add_library(AchillesPlugin_cuts TestCut.cc)
target_link_libraries(AchillesPlugin_cuts PRIVATE achilles_options achilles_warnings
                                          PUBLIC cuts)

set_target_properties(AchillesPlugin_cuts PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
