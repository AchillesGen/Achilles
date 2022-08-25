add_library(plugin_manager PluginManager.cc)
target_link_libraries(plugin_manager PRIVATE project_options project_warnings
                                     PUBLIC event_gen)
