add_library(AchillesHepMC3 SHARED HepMC3EventWriter.cc)
target_link_libraries(AchillesHepMC3 PRIVATE project_options
                             PUBLIC HepMC3::All fmt::fmt spdlog::spdlog yaml::cpp utilities physics nuclear_models)
