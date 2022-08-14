add_subdirectory(Manager)
add_subdirectory(Cuts)
# add_subdirectory(NucleonInteractions)

if(ENABLE_BSM)
add_subdirectory(Sherpa)
endif()

if(ENABLE_HEPMC3)
add_subdirectory(HepMC3)
endif()
