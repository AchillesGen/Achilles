add_library(ConstantInteractions MODULE
    ConstantInteractions.cc
)
target_link_libraries(ConstantInteractions PRIVATE project_options project_warnings
                                           PUBLIC spdlog::spdlog utilities interaction_plugin
)
set_target_properties(ConstantInteractions PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/../Plugins/lib)

add_library(GeantInteractions MODULE
    GeantInteractions.cc
)
target_include_directories(GeantInteractions SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(GeantInteractions PRIVATE project_options project_warnings
                                        PUBLIC spdlog::spdlog utilities ${HDF5_LIBRARIES}
                                               interaction_plugin
)
set_target_properties(GeantInteractions PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/../Plugins/lib)

add_library(NasaInteractions MODULE
    NasaInteractions.cc
)
target_link_libraries(NasaInteractions PRIVATE project_options project_warnings
                                       PUBLIC spdlog::spdlog utilities interaction_plugin
)
set_target_properties(NasaInteractions PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/../Plugins/lib)

add_executable(test_loading test.cc)
target_link_libraries(test_loading PRIVATE project_options project_warnings
                                   PUBLIC utilities interaction_plugin stdc++fs dl)
