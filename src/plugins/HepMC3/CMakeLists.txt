add_library(hepmc3 SHARED HepMC3EventWriter.cc)
target_link_libraries(hepmc3 PRIVATE project_options
                             PUBLIC HepMC3::All fmt::fmt spdlog::spdlog yaml::cpp)
