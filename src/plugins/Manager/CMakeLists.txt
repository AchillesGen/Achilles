add_library(plugin_manager PluginManager.cc)
target_link_libraries(plugin_manager PRIVATE project_options project_warnings
                                     PUBLIC event_gen)

install(TARGETS plugin_manager
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION include
)
