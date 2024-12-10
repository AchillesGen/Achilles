add_library(nuhepmc SHARED NuHepMCWriter.cc)
target_link_libraries(nuhepmc PRIVATE project_options
                               PUBLIC NuHepMC::CPPUtils utilities physics nuclear_models)
