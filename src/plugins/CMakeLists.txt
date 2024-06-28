add_subdirectory(Manager)
add_subdirectory(Cuts)
# add_subdirectory(NucleonInteractions)

if(ACHILLES_ENABLE_SHERPA)
add_subdirectory(Sherpa)
endif()

if(ACHILLES_ENABLE_HEPMC3)
add_subdirectory(HepMC3)
add_subdirectory(NuHepMC)
endif()
