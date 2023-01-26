add_library(hepmc3 SHARED HepMC3EventWriter.cc)
target_include_directories(hepmc3 SYSTEM PUBLIC ${HEPMC3_INCLUDE_DIR})
target_link_libraries(hepmc3 PRIVATE project_options
                             PUBLIC physics ${HEPMC3_LIB} fmt::fmt spdlog::spdlog yaml::cpp)
