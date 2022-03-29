add_executable(convert convert.cc)
target_include_directories(convert SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(convert PRIVATE project_options project_warnings
                              PUBLIC ${HDF5_LIBRARIES} gzstream)

add_executable(evt_read ReadHepMC3.cc)
target_include_directories(evt_read PUBLIC ${HEPMC3_INCLUDE_DIR})
target_link_libraries(evt_read PRIVATE project_options
                               PUBLIC ${HEPMC3_LIB})
