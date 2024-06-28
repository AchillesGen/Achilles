add_library(AchillesPlugin_cuts TestCut.cc)
target_link_libraries(AchillesPlugin_cuts PRIVATE project_options project_warnings
                                          PUBLIC cuts)

set_target_properties(AchillesPlugin_cuts PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
