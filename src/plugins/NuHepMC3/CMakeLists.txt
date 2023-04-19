add_library(nuhepmc3 SHARED NuHepMC3Writer.cc)
target_link_libraries(nuhepmc3 PRIVATE project_options
                               PUBLIC HepMC3 fmt::fmt spdlog::spdlog yaml::cpp)
