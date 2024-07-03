add_library(nuhepmc SHARED NuHepMCWriter.cc)
target_link_libraries(nuhepmc PRIVATE project_options
                               PUBLIC NuHepMC::CPPUtils fmt::fmt spdlog::spdlog yaml::cpp utilities physics nuclear_models)
