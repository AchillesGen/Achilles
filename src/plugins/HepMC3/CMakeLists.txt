#had to change this lib name as it was colliding on a linuxFS on hfs+ with libHepMC3.so
add_library(ach_hepmc3 SHARED HepMC3EventWriter.cc)
target_link_libraries(ach_hepmc3 PRIVATE project_options
                             PUBLIC HepMC3::All fmt::fmt spdlog::spdlog yaml::cpp utilities physics)
