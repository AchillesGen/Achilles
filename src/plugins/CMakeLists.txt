add_subdirectory(Manager)
add_subdirectory(Cuts)
# add_subdirectory(NucleonInteractions)

if(ACHILLES_ENABLE_SHERPA)
add_subdirectory(Sherpa)
endif()

