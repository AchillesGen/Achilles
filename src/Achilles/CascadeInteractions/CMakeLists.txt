add_library(AchillesCascadeInteractions SHARED
    GeantInteractions.cc
    NasaInteractions.cc
    ConstantInteractions.cc
    MesonBaryonInteractions.cc
    OsetMesonBaryonInteractions.cc
    PionAbsorption.cc
    PionInteractions.cc
    DeltaInteractions.cc
)
target_link_libraries(AchillesCascadeInteractions PRIVATE project_options project_warnings
                                                  PUBLIC physics AchillesIntegrator)
get_directory_property(SUBDIR_TARGETS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DEFINITON achilles_targets)
set(achilles_targets ${achilles_targets} AchillesCascadeInteractions PARENT_SCOPE)
