add_library(hepmc3 SHARED HepMC3EventWriter.cc)
# target_include_directories(hepmc3 PUBLIC ${HEPMC3_INCLUDE_DIR})
target_link_libraries(hepmc3 PRIVATE project_options
                             PUBLIC HepMC3 fmt::fmt spdlog::spdlog yaml::cpp)
if(ENABLE_GZIP)
    target_compile_definitions(hepmc3 PUBLIC GZIP)
    target_link_libraries(hepmc3 PUBLIC gzstream)
endif()
