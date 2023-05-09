add_library(nuhepmc SHARED NuHepMCWriter.cc)
target_link_libraries(nuhepmc PRIVATE project_options
                               PUBLIC HepMC3 fmt::fmt spdlog::spdlog yaml::cpp)
